#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GWAS Summary Statistics Standardizer & Allele Aligner
=====================================================
Author : Kesar-Chi
Date   : 2026-03-09
Description:
    交互式 GWAS 摘要统计数据标准化工具。
    - 交互式列映射（questionary + rich）
    - 等位基因对齐（含链翻转、回文位点处理）
    - 质量控制（OR→BETA, -log10P→P, 过滤, 去重）
    - 参考面板强制引入

Usage:
    python gwas_standardizer.py
"""

import argparse
import os
import sys
import math
import gzip
import logging
import shutil
import subprocess
import tempfile
from enum import Enum
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple
from pathlib import Path
from datetime import datetime

try:
    import polars as pl
except ImportError:
    sys.exit("[ERROR] 请安装 polars: pip install polars")

try:
    import questionary
    from questionary import Style
except ImportError:
    sys.exit("[ERROR] 请安装 questionary: pip install questionary")

try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    from rich.logging import RichHandler
    from rich import print as rprint
except ImportError:
    sys.exit("[ERROR] 请安装 rich: pip install rich")

# ─────────────────────────────────────────────────────────────
# 全局配置
# ─────────────────────────────────────────────────────────────
console = Console()

CUSTOM_STYLE = Style([
    ("qmark", "fg:cyan bold"),
    ("question", "fg:white bold"),
    ("answer", "fg:green bold"),
    ("pointer", "fg:cyan bold"),
    ("highlighted", "fg:cyan bold"),
    ("selected", "fg:green"),
])

# 互补碱基映射
COMPLEMENT: Dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}

# 回文碱基对
PALINDROMIC_PAIRS = {frozenset({"A", "T"}), frozenset({"C", "G"})}

# EAF 阈值
EAF_THRESHOLD = 0.15
EAF_MID_THRESHOLD = 0.10  # |EAF - 0.5| < 0.10 视为中间频率


def _env_float(name: str, default: float) -> float:
    value = os.environ.get(name)
    if value is None or value == "":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def _env_int(name: str, default: int) -> int:
    value = os.environ.get(name)
    if value is None or value == "":
        return default
    try:
        return int(value)
    except ValueError:
        return default


DEFAULT_R_LIB_PATH = os.environ.get("MR_PIPELINE_R_LIB_PATH", "/home/ding/R/4.4.1_MR")
DEFAULT_REFERENCE_PANEL = os.environ.get("MR_PIPELINE_REFERENCE_PANEL", "")
DEFAULT_ORG_DIR = os.environ.get("MR_PIPELINE_ORG_DIR", "")
DEFAULT_STANDARDIZED_OUTPUT_DIR = os.environ.get("MR_PIPELINE_STANDARDIZED_OUTPUT_DIR", "")
DEFAULT_EXPOSURE_OUTPUT_DIR = os.environ.get("MR_PIPELINE_EXP_DIR", os.environ.get("MR_PIPELINE_EXPOSURE_DIR", ""))
DEFAULT_OUTCOME_OUTPUT_DIR = os.environ.get("MR_PIPELINE_OUT_DIR", os.environ.get("MR_PIPELINE_OUTCOME_DIR", ""))
DEFAULT_CLUMP_PLINK = os.environ.get("MR_PIPELINE_CLUMP_PLINK", "/home/ding/miniconda3/envs/GWAS/bin/plink")
DEFAULT_CLUMP_BFILE = os.environ.get("MR_PIPELINE_CLUMP_BFILE", "/home/ding/MR_LPA/Ref/g1000_eur/g1000_eur")
DEFAULT_CLUMP_R2 = _env_float("MR_PIPELINE_CLUMP_R2", 0.1)
DEFAULT_CLUMP_KB = _env_int("MR_PIPELINE_CLUMP_KB", 500)
DEFAULT_CLUMP_P1 = _env_float("MR_PIPELINE_CLUMP_P1", 1e-4)
DEFAULT_CLUMP_POP = os.environ.get("MR_PIPELINE_CLUMP_POP", "EUR")


# ─────────────────────────────────────────────────────────────
# 数据模型
# ─────────────────────────────────────────────────────────────
class AlleleMode(Enum):
    """等位基因记录模式"""
    REF_ALT = "A"       # Mode A: REF/ALT (VCF 规范)
    EFFECT_OTHER = "B"  # Mode B: Effect Allele / Other Allele
    A1_A2 = "C"         # Mode C: Allele1 / Allele2 (方向二义性)


class StatType(Enum):
    """效应统计量类型"""
    BETA = "BETA"
    OR = "OR"


class PvalFormat(Enum):
    """P 值格式"""
    RAW = "raw"         # 原始 P 值
    NEGLOG10 = "neglog10"  # -log10(P)


class FreqType(Enum):
    """频率类型"""
    EAF = "EAF"  # 效应等位基因频率
    MAF = "MAF"  # 次等位基因频率


@dataclass
class ColumnMapping:
    """列映射配置"""
    # 基础坐标
    snp_col: Optional[str] = None
    chr_col: str = ""
    pos_col: str = ""

    # 等位基因列
    allele_mode: AlleleMode = AlleleMode.REF_ALT
    allele1_col: str = ""   # REF / EA / A1
    allele2_col: str = ""   # ALT / OA / A2

    # Mode C 特有：效应方向
    effect_on: Optional[str] = None  # "A1", "A2", "Unknown"

    # 统计量
    stat_type: StatType = StatType.BETA
    stat_col: str = ""      # BETA 或 OR 列名
    se_col: str = ""
    p_col: str = ""
    pval_format: PvalFormat = PvalFormat.RAW

    # 频率
    freq_type: FreqType = FreqType.EAF
    freq_col: Optional[str] = None

    # 文件属性
    separator: str = "\t"


@dataclass
class AlignmentStats:
    """对齐统计"""
    total_input: int = 0
    matched: int = 0
    swapped: int = 0
    flipped: int = 0
    flipped_swapped: int = 0
    dropped_mismatch: int = 0
    dropped_ambiguous_palindromic: int = 0
    dropped_unresolved_mode_c: int = 0
    qc_failed: int = 0
    duplicates_removed: int = 0
    not_in_ref: int = 0
    total_output: int = 0


# ─────────────────────────────────────────────────────────────
# 工具函数
# ─────────────────────────────────────────────────────────────
def complement_base(base: str) -> str:
    """返回互补碱基"""
    return COMPLEMENT.get(base.upper(), "N")


def complement_alleles(a1: str, a2: str) -> Tuple[str, str]:
    """返回互补碱基对"""
    return complement_base(a1), complement_base(a2)


def is_palindromic(a1: str, a2: str) -> bool:
    """判断是否为回文位点 (A/T 或 C/G)"""
    return frozenset({a1.upper(), a2.upper()}) in PALINDROMIC_PAIRS


def make_sort_snpid(chrom: str, pos: str, a1: str, a2: str) -> str:
    """生成标准化 Sort_SNPID (碱基字母排序)"""
    a1, a2 = a1.upper(), a2.upper()
    if a1 <= a2:
        return f"{chrom}:{pos}:{a1}:{a2}"
    else:
        return f"{chrom}:{pos}:{a2}:{a1}"


def make_bim_id(chrom: str, pos: str, a1: str, a2: str) -> str:
    """生成与 BIM 一致的字典序 ID: CHR_POS_A1_A2"""
    return make_sort_snpid(chrom, pos, a1, a2).replace(":", "_")


def detect_separator(filepath: str) -> str:
    """自动检测文件分隔符"""
    open_fn = gzip.open if filepath.endswith(".gz") else open
    with open_fn(filepath, "rt") as f:
        first_line = f.readline().strip()

    if "\t" in first_line:
        return "\t"
    elif "," in first_line:
        return ","
    elif " " in first_line:
        return " "
    return "\t"


def read_file_header(filepath: str, n_rows: int = 100) -> pl.DataFrame:
    """读取文件的前 N 行"""
    sep = detect_separator(filepath)
    try:
        df = pl.read_csv(
            filepath,
            separator=sep,
            n_rows=n_rows,
            infer_schema_length=n_rows,
            truncate_ragged_lines=True,
            ignore_errors=True,
        )
        return df
    except Exception as e:
        console.print(f"[red]读取文件失败: {e}[/red]")
        sys.exit(1)


def get_input_base(input_path: str) -> str:
    """去掉常见扩展名后的输入基名"""
    input_base = input_path
    if input_base.endswith(".gz"):
        input_base = input_base.rsplit(".", 2)[0]
    else:
        input_base = input_base.rsplit(".", 1)[0]
    return input_base


def get_configured_output_dir(output_format: str, mr_role: Optional[str] = None) -> Optional[str]:
    """根据输出类型选择 config 中的默认输出目录"""
    if output_format == "mr":
        role = normalize_mr_role(mr_role)
        if role == "outcome" and DEFAULT_OUTCOME_OUTPUT_DIR:
            return DEFAULT_OUTCOME_OUTPUT_DIR
        if role == "exposure" and DEFAULT_EXPOSURE_OUTPUT_DIR:
            return DEFAULT_EXPOSURE_OUTPUT_DIR
        return None

    return DEFAULT_STANDARDIZED_OUTPUT_DIR or None


def derive_default_output_path(
    input_path: str,
    output_format: str = "standardized",
    mr_role: Optional[str] = None,
    output_dir: Optional[str] = None,
) -> str:
    """根据输入路径和输出模式生成默认输出路径"""
    input_base = get_input_base(input_path)
    if output_format == "mr":
        role_suffix_map = {
            "out": "_outcome.csv",
            "outcome": "_outcome.csv",
            "exp": "_exposure.csv",
            "exposure": "_exposure.csv",
        }
        suffix = role_suffix_map.get((mr_role or "").lower(), "_mr.csv")
    else:
        suffix = "_standardized.tsv.gz"

    target_dir = output_dir or get_configured_output_dir(output_format, mr_role)
    if target_dir:
        return str(Path(target_dir) / f"{Path(input_base).name}{suffix}")
    return input_base + suffix


def parse_args() -> argparse.Namespace:
    """解析命令行参数，支持交互式和非交互式两种模式"""
    parser = argparse.ArgumentParser(
        description="GWAS 摘要统计数据标准化与等位基因对齐工具",
    )
    parser.add_argument("--input", help="输入 GWAS 摘要统计文件")
    parser.add_argument("--reference", help="参考面板文件")
    parser.add_argument("--output", help="输出文件路径")
    parser.add_argument(
        "--output-format",
        choices=["standardized", "mr"],
        default="standardized",
        help="输出格式: standardized 或 mr (TwoSampleMR 输入友好格式)",
    )
    parser.add_argument("--mode", choices=[m.value for m in AlleleMode], help="等位基因模式: A/B/C")
    parser.add_argument("--snp-col", help="SNP/variant ID 列名")
    parser.add_argument("--chr-col", help="染色体列名")
    parser.add_argument("--pos-col", help="位置列名")
    parser.add_argument("--allele1-col", help="等位基因 1 列名")
    parser.add_argument("--allele2-col", help="等位基因 2 列名")
    parser.add_argument("--effect-on", choices=["A1", "A2", "Unknown"], help="Mode C 时效应作用在哪个等位基因上")
    parser.add_argument("--stat-type", choices=[s.value for s in StatType], help="统计量类型: BETA/OR")
    parser.add_argument("--stat-col", help="效应统计量列名")
    parser.add_argument("--se-col", help="标准误列名")
    parser.add_argument("--p-col", help="P 值列名")
    parser.add_argument("--pval-format", choices=[p.value for p in PvalFormat], help="P 值格式: raw/neglog10")
    parser.add_argument("--freq-type", choices=[f.value for f in FreqType], help="频率类型: EAF/MAF")
    parser.add_argument("--freq-col", help="频率列名")
    parser.add_argument("--phenotype", help="MR 输出时附加的 phenotype 列值")
    parser.add_argument("--sample-size", type=int, help="MR 输出时附加的样本量 N")
    parser.add_argument(
        "--mr-role",
        choices=["out", "outcome", "exp", "exposure"],
        help="MR 输出角色: out/outcome 或 exp/exposure；exp 会自动执行 clump",
    )
    parser.add_argument(
        "--r-lib-path",
        default=DEFAULT_R_LIB_PATH,
        help="R 包库路径，用于加载 data.table 和 TwoSampleMR",
    )
    parser.add_argument("--clump-r2", type=float, default=DEFAULT_CLUMP_R2, help="Exposure clump 的 r2 阈值")
    parser.add_argument("--clump-kb", type=int, default=DEFAULT_CLUMP_KB, help="Exposure clump 的窗口大小（kb）")
    parser.add_argument("--clump-p1", type=float, default=DEFAULT_CLUMP_P1, help="Exposure clump 的 p 值阈值")
    parser.add_argument("--clump-pop", default=DEFAULT_CLUMP_POP, help="Exposure clump 的参考人群")
    parser.add_argument(
        "--clump-plink",
        default=DEFAULT_CLUMP_PLINK,
        help="Exposure clump 使用的 plink 可执行文件",
    )
    parser.add_argument(
        "--clump-bfile",
        default=DEFAULT_CLUMP_BFILE,
        help="Exposure clump 使用的 PLINK bfile 前缀",
    )
    parser.add_argument(
        "--non-interactive",
        action="store_true",
        help="禁止进入交互式提问；缺少必要参数时直接报错",
    )
    return parser.parse_args()


def is_cli_mode(args: argparse.Namespace) -> bool:
    """判断当前是否使用命令行参数模式"""
    cli_fields = [
        args.input,
        args.reference,
        args.output,
        args.mode,
        args.snp_col,
        args.chr_col,
        args.pos_col,
        args.allele1_col,
        args.allele2_col,
        args.effect_on,
        args.stat_type,
        args.stat_col,
        args.se_col,
        args.p_col,
        args.pval_format,
        args.freq_type,
        args.freq_col,
        args.mr_role,
    ]
    return args.non_interactive or any(value is not None for value in cli_fields)


def build_mapping_from_args(args: argparse.Namespace) -> ColumnMapping:
    """根据命令行参数构建列映射"""
    required = [
        ("--input", args.input),
        ("--reference", args.reference),
        ("--mode", args.mode),
        ("--chr-col", args.chr_col),
        ("--pos-col", args.pos_col),
        ("--allele1-col", args.allele1_col),
        ("--allele2-col", args.allele2_col),
        ("--stat-type", args.stat_type),
        ("--stat-col", args.stat_col),
        ("--se-col", args.se_col),
        ("--p-col", args.p_col),
        ("--pval-format", args.pval_format),
        ("--freq-type", args.freq_type),
    ]
    missing = [flag for flag, value in required if value is None]
    if missing:
        console.print(f"[red]非交互模式缺少必要参数: {', '.join(missing)}[/red]")
        sys.exit(1)

    if args.mode == AlleleMode.A1_A2.value and args.effect_on is None:
        console.print("[red]Mode C 需要提供 --effect-on (A1/A2/Unknown)[/red]")
        sys.exit(1)

    if args.mode != AlleleMode.A1_A2.value and args.effect_on is not None:
        console.print("[yellow]--effect-on 仅在 Mode C 下生效，当前将忽略。[/yellow]")

    mapping = ColumnMapping(
        snp_col=args.snp_col,
        chr_col=args.chr_col or "",
        pos_col=args.pos_col or "",
        allele_mode=AlleleMode(args.mode),
        allele1_col=args.allele1_col or "",
        allele2_col=args.allele2_col or "",
        effect_on=args.effect_on if args.mode == AlleleMode.A1_A2.value else None,
        stat_type=StatType(args.stat_type),
        stat_col=args.stat_col or "",
        se_col=args.se_col or "",
        p_col=args.p_col or "",
        pval_format=PvalFormat(args.pval_format),
        freq_type=FreqType(args.freq_type),
        freq_col=args.freq_col,
    )
    return mapping


# ─────────────────────────────────────────────────────────────
# 交互式列映射模块
# ─────────────────────────────────────────────────────────────
def display_preview(df: pl.DataFrame):
    """使用 rich 展示数据预览"""
    table = Table(title="📊 输入文件预览 (前 5 行)", show_lines=True)
    for col in df.columns:
        table.add_column(col, style="cyan", no_wrap=True)
    for row in df.head(5).iter_rows():
        table.add_row(*[str(v) for v in row])
    console.print(table)


def ask_file_path(prompt: str, must_exist: bool = True, default: str = "") -> str:
    """交互式询问文件路径"""
    while True:
        path = questionary.path(
            prompt,
            default=default,
            style=CUSTOM_STYLE,
        ).ask()
        if path is None:
            console.print("[red]操作已取消[/red]")
            sys.exit(0)
        path = os.path.expanduser(path.strip())
        if must_exist and not os.path.isfile(path):
            console.print(f"[red]文件不存在: {path}[/red]")
            continue
        return path


def ask_output_path(default_path: str) -> str:
    """交互式询问输出路径"""
    path = questionary.text(
        "📁 请输入输出文件路径:",
        default=default_path,
        style=CUSTOM_STYLE,
    ).ask()
    if path is None:
        console.print("[red]操作已取消[/red]")
        sys.exit(0)
    return path.strip()


def ask_text(prompt: str, default: str = "") -> str:
    """通用文本输入"""
    answer = questionary.text(
        prompt,
        default=default,
        style=CUSTOM_STYLE,
    ).ask()
    if answer is None:
        console.print("[red]操作已取消[/red]")
        sys.exit(0)
    return answer.strip()


def ask_reference_path() -> str:
    """交互式获取参考面板路径，优先使用 config 默认值"""
    if DEFAULT_REFERENCE_PANEL and os.path.isfile(DEFAULT_REFERENCE_PANEL):
        use_default = questionary.confirm(
            f"📄 使用 config 中的参考面板?\n{DEFAULT_REFERENCE_PANEL}",
            default=True,
            style=CUSTOM_STYLE,
        ).ask()
        if use_default:
            return DEFAULT_REFERENCE_PANEL

    default_path = DEFAULT_REFERENCE_PANEL if DEFAULT_REFERENCE_PANEL else ""
    return ask_file_path("📄 请输入参考面板 (Reference Panel) 文件路径:", default=default_path)


def ask_input_path() -> str:
    """交互式获取输入 GWAS 路径，优先提示 Org 目录"""
    default_path = str(Path(DEFAULT_ORG_DIR)) if DEFAULT_ORG_DIR else ""
    return ask_file_path("📄 请输入 GWAS 摘要统计文件路径:", default=default_path)


def interactive_runtime_options(input_path: str) -> Dict[str, object]:
    """交互式询问输出格式与 clump 配置"""
    input_name = Path(get_input_base(input_path)).name
    output_format = questionary.select(
        "请选择输出格式:",
        choices=[
            questionary.Choice("standardized (标准化 TSV)", value="standardized"),
            questionary.Choice("mr (TwoSampleMR 友好格式)", value="mr"),
        ],
        style=CUSTOM_STYLE,
    ).ask()
    if output_format is None:
        console.print("[red]操作已取消[/red]")
        sys.exit(0)

    runtime: Dict[str, object] = {
        "output_format": output_format,
        "mr_role": None,
        "phenotype": None,
        "sample_size": None,
        "clump_r2": DEFAULT_CLUMP_R2,
        "clump_kb": DEFAULT_CLUMP_KB,
        "clump_p1": DEFAULT_CLUMP_P1,
        "clump_pop": DEFAULT_CLUMP_POP,
        "clump_plink": DEFAULT_CLUMP_PLINK,
        "clump_bfile": DEFAULT_CLUMP_BFILE,
        "r_lib_path": DEFAULT_R_LIB_PATH,
    }

    if output_format != "mr":
        return runtime

    mr_role = questionary.select(
        "请选择 MR 数据角色:",
        choices=[
            questionary.Choice("Outcome", value="outcome"),
            questionary.Choice("Exposure (自动 clump)", value="exposure"),
        ],
        style=CUSTOM_STYLE,
    ).ask()
    if mr_role is None:
        console.print("[red]操作已取消[/red]")
        sys.exit(0)
    runtime["mr_role"] = mr_role

    phenotype = ask_text("📛 请输入 phenotype 名称 (可留空):", default=input_name)
    runtime["phenotype"] = phenotype or None

    sample_size_text = ask_text("👥 请输入样本量 N (可留空):", default="")
    runtime["sample_size"] = int(sample_size_text) if sample_size_text else None

    if mr_role != "exposure":
        return runtime

    console.print("[bold yellow]━━━ Clump 参数配置 ━━━[/bold yellow]")
    runtime["clump_r2"] = float(ask_text("clump_r2:", default=str(DEFAULT_CLUMP_R2)))
    runtime["clump_kb"] = int(ask_text("clump_kb:", default=str(DEFAULT_CLUMP_KB)))
    runtime["clump_p1"] = float(ask_text("clump_p1:", default=str(DEFAULT_CLUMP_P1)))
    runtime["clump_pop"] = ask_text("pop:", default=DEFAULT_CLUMP_POP)
    runtime["clump_plink"] = ask_text("plink 路径:", default=DEFAULT_CLUMP_PLINK)
    runtime["clump_bfile"] = ask_text("clump 参考 bfile 前缀:", default=DEFAULT_CLUMP_BFILE)

    return runtime


def select_column(columns: List[str], prompt: str, allow_none: bool = False) -> Optional[str]:
    """从列名列表中选择一个列"""
    choices = list(columns)
    if allow_none:
        choices.append("⏭  跳过 (无此列)")

    answer = questionary.select(
        prompt,
        choices=choices,
        style=CUSTOM_STYLE,
    ).ask()

    if answer is None:
        console.print("[red]操作已取消[/red]")
        sys.exit(0)
    if answer == "⏭  跳过 (无此列)":
        return None
    return answer


def interactive_mapping(df: pl.DataFrame) -> ColumnMapping:
    """交互式列映射主流程"""
    columns = df.columns
    mapping = ColumnMapping()

    console.print(Panel.fit(
        "[bold cyan]🧬 GWAS 摘要统计数据标准化工具[/bold cyan]\n"
        "[dim]交互式列名映射 & 等位基因对齐[/dim]",
        border_style="cyan",
    ))

    # 展示预览
    display_preview(df)
    console.print()

    # ── Step 1: 选择等位基因模式 ──
    console.print("[bold yellow]━━━ Step 1: 选择等位基因记录模式 ━━━[/bold yellow]")
    mode_answer = questionary.select(
        "请选择输入文件的等位基因记录规范:",
        choices=[
            questionary.Choice("Mode A: REF/ALT (VCF 规范, BETA 针对 ALT)", value="A"),
            questionary.Choice("Mode B: Effect/Other (明确效应等位基因 EA/OA)", value="B"),
            questionary.Choice("Mode C: Allele1/Allele2 (方向可能有二义性)", value="C"),
        ],
        style=CUSTOM_STYLE,
    ).ask()
    if mode_answer is None:
        sys.exit(0)

    mapping.allele_mode = AlleleMode(mode_answer)
    console.print(f"  ✅ 已选择: Mode {mode_answer}\n")

    # ── Step 2: 基础坐标映射 ──
    console.print("[bold yellow]━━━ Step 2: 映射基础坐标列 ━━━[/bold yellow]")
    mapping.snp_col = select_column(columns, "🔹 SNP ID 列 (rsID / variant_id):", allow_none=True)
    mapping.chr_col = select_column(columns, "🔹 染色体 (CHR) 列:")
    mapping.pos_col = select_column(columns, "🔹 位置 (POS) 列:")
    console.print()

    # ── Step 3: 等位基因映射 ──
    console.print("[bold yellow]━━━ Step 3: 映射等位基因列 ━━━[/bold yellow]")
    if mapping.allele_mode == AlleleMode.REF_ALT:
        mapping.allele1_col = select_column(columns, "🔹 REF (参考等位基因) 列:")
        mapping.allele2_col = select_column(columns, "🔹 ALT (替代等位基因) 列:")
    elif mapping.allele_mode == AlleleMode.EFFECT_OTHER:
        mapping.allele1_col = select_column(columns, "🔹 Effect Allele (EA, 效应等位基因) 列:")
        mapping.allele2_col = select_column(columns, "🔹 Other Allele (OA, 非效应等位基因) 列:")
    else:  # Mode C
        mapping.allele1_col = select_column(columns, "🔹 Allele1 (A1) 列:")
        mapping.allele2_col = select_column(columns, "🔹 Allele2 (A2) 列:")
        # 追加提问
        effect_answer = questionary.select(
            "效应值是针对哪个等位基因计算的?",
            choices=[
                questionary.Choice("A1 (Allele1)", value="A1"),
                questionary.Choice("A2 (Allele2)", value="A2"),
                questionary.Choice("Unknown (将通过频率反推)", value="Unknown"),
            ],
            style=CUSTOM_STYLE,
        ).ask()
        if effect_answer is None:
            sys.exit(0)
        mapping.effect_on = effect_answer
    console.print()

    # ── Step 4: 统计量映射 ──
    console.print("[bold yellow]━━━ Step 4: 映射统计量列 ━━━[/bold yellow]")

    # BETA 或 OR
    stat_answer = questionary.select(
        "效应统计量类型:",
        choices=[
            questionary.Choice("BETA (效应值)", value="BETA"),
            questionary.Choice("OR (比值比, 将自动转换为 ln(OR))", value="OR"),
        ],
        style=CUSTOM_STYLE,
    ).ask()
    if stat_answer is None:
        sys.exit(0)
    mapping.stat_type = StatType(stat_answer)
    mapping.stat_col = select_column(columns, f"🔹 {stat_answer} 列:")

    # SE
    mapping.se_col = select_column(columns, "🔹 SE (标准误) 列:")

    # P
    mapping.p_col = select_column(columns, "🔹 P (P 值) 列:")
    pval_answer = questionary.select(
        "P 值格式:",
        choices=[
            questionary.Choice("原始 P 值 (如 0.05)", value="raw"),
            questionary.Choice("-log10(P) 格式 (如 1.3)", value="neglog10"),
        ],
        style=CUSTOM_STYLE,
    ).ask()
    if pval_answer is None:
        sys.exit(0)
    mapping.pval_format = PvalFormat(pval_answer)

    # EAF / MAF
    freq_answer = questionary.select(
        "频率列类型:",
        choices=[
            questionary.Choice("EAF (效应等位基因频率)", value="EAF"),
            questionary.Choice("MAF (次等位基因频率)", value="MAF"),
        ],
        style=CUSTOM_STYLE,
    ).ask()
    if freq_answer is None:
        sys.exit(0)
    mapping.freq_type = FreqType(freq_answer)
    mapping.freq_col = select_column(columns, f"🔹 {freq_answer} 列:", allow_none=True)
    console.print()

    # ── Step 5: 确认 ──
    display_mapping_summary(mapping)

    confirm = questionary.confirm(
        "以上映射配置是否正确?",
        default=True,
        style=CUSTOM_STYLE,
    ).ask()
    if not confirm:
        console.print("[yellow]请重新运行程序配置映射。[/yellow]")
        sys.exit(0)

    return mapping


def display_mapping_summary(mapping: ColumnMapping):
    """展示映射配置摘要"""
    table = Table(title="📋 列映射配置摘要", show_lines=True)
    table.add_column("配置项", style="cyan bold", width=20)
    table.add_column("值", style="green")

    mode_labels = {
        AlleleMode.REF_ALT: "Mode A (REF/ALT)",
        AlleleMode.EFFECT_OTHER: "Mode B (Effect/Other)",
        AlleleMode.A1_A2: "Mode C (A1/A2)",
    }
    table.add_row("等位基因模式", mode_labels[mapping.allele_mode])
    table.add_row("SNP 列", mapping.snp_col or "(无)")
    table.add_row("CHR 列", mapping.chr_col)
    table.add_row("POS 列", mapping.pos_col)
    table.add_row("等位基因 1", mapping.allele1_col)
    table.add_row("等位基因 2", mapping.allele2_col)

    if mapping.allele_mode == AlleleMode.A1_A2:
        table.add_row("效应方向", mapping.effect_on or "Unknown")

    table.add_row("统计量类型", mapping.stat_type.value)
    table.add_row("统计量列", mapping.stat_col)
    table.add_row("SE 列", mapping.se_col)
    table.add_row("P 值列", mapping.p_col)
    table.add_row("P 值格式", mapping.pval_format.value)
    table.add_row("频率类型", mapping.freq_type.value)
    table.add_row("频率列", mapping.freq_col or "(无)")

    console.print(table)


def add_frequency_columns(
    df: pl.DataFrame,
    mapping: ColumnMapping,
    *,
    flip_eaf: bool = False,
) -> pl.DataFrame:
    """根据频率类型统一生成 _EAF / _MAF 列"""
    null_float = pl.lit(None).cast(pl.Float64)

    if "_FREQ" not in df.columns:
        return df.with_columns([
            null_float.alias("_EAF"),
            null_float.alias("_MAF"),
        ])

    if mapping.freq_type == FreqType.EAF:
        eaf_expr = pl.col("_FREQ")
        if flip_eaf:
            eaf_expr = 1.0 - eaf_expr
        return df.with_columns([
            eaf_expr.alias("_EAF"),
            null_float.alias("_MAF"),
        ])

    return df.with_columns([
        null_float.alias("_EAF"),
        pl.col("_FREQ").alias("_MAF"),
    ])


def normalize_mr_role(mr_role: Optional[str]) -> Optional[str]:
    """规范化 MR 角色命名"""
    if mr_role is None:
        return None

    role = mr_role.lower()
    if role in ("out", "outcome"):
        return "outcome"
    if role in ("exp", "exposure"):
        return "exposure"
    raise ValueError(f"Unsupported MR role: {mr_role}")


# ─────────────────────────────────────────────────────────────
# QC & 数据转换
# ─────────────────────────────────────────────────────────────
def transform_and_qc(df: pl.DataFrame, mapping: ColumnMapping) -> Tuple[pl.DataFrame, int]:
    """
    执行数据转换和 QC 过滤。
    返回 (清洗后的 DataFrame, QC 失败的行数)
    """
    initial_count = df.height
    missing_tokens = ["", ".", "NA", "NAN", "NULL"]

    # ── 重命名列 ──
    rename_map = {
        mapping.chr_col: "_CHR",
        mapping.pos_col: "_POS",
        mapping.allele1_col: "_A1",
        mapping.allele2_col: "_A2",
        mapping.stat_col: "_STAT",
        mapping.se_col: "_SE",
        mapping.p_col: "_P",
    }
    if mapping.snp_col:
        rename_map[mapping.snp_col] = "_SNP"
    if mapping.freq_col:
        rename_map[mapping.freq_col] = "_FREQ"

    df = df.rename(rename_map)

    # ── 类型转换 ──
    df = df.with_columns([
        pl.col("_CHR").cast(pl.Utf8).str.replace("chr", "").str.strip_chars().alias("_CHR"),
        pl.col("_POS").cast(pl.Utf8).str.strip_chars().alias("_POS"),
        pl.col("_A1").cast(pl.Utf8).str.to_uppercase().str.strip_chars().alias("_A1"),
        pl.col("_A2").cast(pl.Utf8).str.to_uppercase().str.strip_chars().alias("_A2"),
        pl.col("_STAT").cast(pl.Float64, strict=False).alias("_STAT"),
        pl.col("_SE").cast(pl.Float64, strict=False).alias("_SE"),
        pl.col("_P").cast(pl.Float64, strict=False).alias("_P"),
    ])

    if "_FREQ" in df.columns:
        df = df.with_columns(
            pl.col("_FREQ").cast(pl.Float64, strict=False).alias("_FREQ")
        )

    # ── OR → BETA ──
    if mapping.stat_type == StatType.OR:
        df = df.with_columns(
            pl.col("_STAT").log().alias("_STAT")  # ln(OR) → BETA
        )

    # ── -log10(P) → P ──
    if mapping.pval_format == PvalFormat.NEGLOG10:
        df = df.with_columns(
            (pl.lit(10.0).pow(-pl.col("_P"))).alias("_P")
        )

    # ── QC 过滤 ──
    # 剔除核心列的缺失值
    qc_filter = (
        pl.col("_CHR").is_not_null()
        & pl.col("_POS").is_not_null()
        & pl.col("_A1").is_not_null()
        & pl.col("_A2").is_not_null()
        & (~pl.col("_CHR").str.to_uppercase().is_in(missing_tokens))
        & (~pl.col("_POS").str.to_uppercase().is_in(missing_tokens))
        & (~pl.col("_A1").str.to_uppercase().is_in(missing_tokens))
        & (~pl.col("_A2").str.to_uppercase().is_in(missing_tokens))
        & pl.col("_STAT").is_not_null()
        & pl.col("_SE").is_not_null()
        & pl.col("_P").is_not_null()
        # P ∈ (0, 1]
        & (pl.col("_P") > 0)
        & (pl.col("_P") <= 1)
        # SE > 0
        & (pl.col("_SE") > 0)
    )

    if "_FREQ" in df.columns:
        qc_filter = qc_filter & (
            pl.col("_FREQ").is_null()
            | ((pl.col("_FREQ") >= 0) & (pl.col("_FREQ") <= 1))
        )

    df_clean = df.filter(qc_filter)
    qc_failed = initial_count - df_clean.height

    return df_clean, qc_failed


# ─────────────────────────────────────────────────────────────
# 等位基因对齐算法
# ─────────────────────────────────────────────────────────────
def load_reference_panel(filepath: str) -> pl.DataFrame:
    """加载参考面板并生成 Sort_SNPID"""
    sep = detect_separator(filepath)
    console.print(f"[cyan]🔧 正在加载参考面板: {filepath}[/cyan]")

    ref = pl.read_csv(
        filepath,
        separator=sep,
        infer_schema_length=10000,
        truncate_ragged_lines=True,
        ignore_errors=True,
    )

    # 标准化列名（不区分大小写匹配）
    col_map = {}
    for c in ref.columns:
        cl = c.upper().strip()
        if cl in ("CHR", "#CHR", "CHROM", "#CHROM"):
            col_map[c] = "REF_CHR"
        elif cl in ("POS", "BP", "POSITION"):
            col_map[c] = "REF_POS"
        elif cl in ("REF", "A1", "ALLELE1"):
            col_map[c] = "REF_REF"
        elif cl in ("ALT", "A2", "ALLELE2"):
            col_map[c] = "REF_ALT"
        elif cl in ("AF", "MAF", "A1F", "FRQ", "FREQ", "EAF"):
            col_map[c] = "REF_AF"
        elif cl in ("RSID", "ID", "SNP_ID"):
            col_map[c] = "REF_RSID"

    ref = ref.rename(col_map)

    # 确保必要列存在
    required = ["REF_CHR", "REF_POS", "REF_REF", "REF_ALT", "REF_AF"]
    missing = [c for c in required if c not in ref.columns]
    if missing:
        console.print(f"[red]参考面板缺少必要列: {missing}[/red]")
        console.print(f"[dim]当前列: {ref.columns}[/dim]")
        sys.exit(1)

    # 类型标准化
    ref = ref.with_columns([
        pl.col("REF_CHR").cast(pl.Utf8).str.replace("chr", "").str.strip_chars(),
        pl.col("REF_POS").cast(pl.Utf8).str.strip_chars(),
        pl.col("REF_REF").cast(pl.Utf8).str.to_uppercase().str.strip_chars(),
        pl.col("REF_ALT").cast(pl.Utf8).str.to_uppercase().str.strip_chars(),
        pl.col("REF_AF").cast(pl.Float64, strict=False),
    ])

    # 生成 Sort_SNPID
    ref = ref.with_columns([
        pl.struct(["REF_CHR", "REF_POS", "REF_REF", "REF_ALT"])
        .map_elements(
            lambda row: make_sort_snpid(
                str(row["REF_CHR"]), str(row["REF_POS"]),
                str(row["REF_REF"]), str(row["REF_ALT"])
            ),
            return_dtype=pl.Utf8,
        )
        .alias("REF_SORT_SNPID"),
        pl.struct(["REF_CHR", "REF_POS", "REF_REF", "REF_ALT"])
        .map_elements(
            lambda row: make_bim_id(
                str(row["REF_CHR"]), str(row["REF_POS"]),
                str(row["REF_REF"]), str(row["REF_ALT"])
            ),
            return_dtype=pl.Utf8,
        )
        .alias("REF_BIM_ID"),
    ])

    # 去重（保留第一条）
    ref = ref.unique(subset=["REF_SORT_SNPID"], keep="first")

    select_cols = ["REF_SORT_SNPID", "REF_BIM_ID", "REF_CHR", "REF_POS", "REF_REF", "REF_ALT", "REF_AF"]
    if "REF_RSID" in ref.columns:
        select_cols.append("REF_RSID")

    console.print(f"  ✅ 参考面板加载完成: {ref.height:,} 个变异位点\n")
    return ref.select(select_cols)


def resolve_effect_allele(df: pl.DataFrame, mapping: ColumnMapping) -> pl.DataFrame:
    """
    根据等位基因模式，统一列名为 _Aeff (效应等位基因) 和 _Aref (非效应等位基因)。
    同时生成 _EAF 列。
    """
    if mapping.allele_mode == AlleleMode.REF_ALT:
        # Mode A: ALT 是效应等位基因
        df = df.with_columns([
            pl.col("_A2").alias("_Aeff"),  # ALT = effect
            pl.col("_A1").alias("_Aref"),  # REF = reference
        ])
        df = add_frequency_columns(df, mapping)

    elif mapping.allele_mode == AlleleMode.EFFECT_OTHER:
        # Mode B: A1 = EA (effect), A2 = OA (other)
        df = df.with_columns([
            pl.col("_A1").alias("_Aeff"),
            pl.col("_A2").alias("_Aref"),
        ])
        df = add_frequency_columns(df, mapping)

    elif mapping.allele_mode == AlleleMode.A1_A2:
        # Mode C: 需要根据 effect_on 决定
        if mapping.effect_on == "A1":
            df = df.with_columns([
                pl.col("_A1").alias("_Aeff"),
                pl.col("_A2").alias("_Aref"),
            ])
            df = add_frequency_columns(df, mapping)
        elif mapping.effect_on == "A2":
            df = df.with_columns([
                pl.col("_A2").alias("_Aeff"),
                pl.col("_A1").alias("_Aref"),
            ])
            df = add_frequency_columns(df, mapping, flip_eaf=True)
        else:  # Unknown - 先保留原始，后续通过参考面板反推
            df = df.with_columns([
                pl.col("_A1").alias("_Aeff"),
                pl.col("_A2").alias("_Aref"),
            ])
            df = add_frequency_columns(df, mapping)

    return df


def align_alleles(
    df: pl.DataFrame,
    ref: pl.DataFrame,
    mapping: ColumnMapping,
) -> Tuple[pl.DataFrame, AlignmentStats]:
    """
    核心等位基因对齐逻辑。
    返回 (对齐后的 DataFrame, 对齐统计)
    """
    stats = AlignmentStats()
    stats.total_input = df.height

    # 生成 Sort_SNPID
    df = df.with_columns(
        pl.struct(["_CHR", "_POS", "_Aeff", "_Aref"])
        .map_elements(
            lambda row: make_sort_snpid(
                str(row["_CHR"]), str(row["_POS"]),
                str(row["_Aeff"]), str(row["_Aref"])
            ),
            return_dtype=pl.Utf8,
        )
        .alias("_SORT_SNPID")
    )

    # ── 去重：保留 P 值最小的 ──
    before_dedup = df.height
    df = df.sort("_P").unique(subset=["_SORT_SNPID"], keep="first")
    stats.duplicates_removed = before_dedup - df.height

    # ── 与参考面板连接 ──
    df = df.join(ref, left_on="_SORT_SNPID", right_on="REF_SORT_SNPID", how="left")

    # 未在参考面板中找到的位点
    not_in_ref = df.filter(pl.col("REF_REF").is_null())
    stats.not_in_ref = not_in_ref.height
    df = df.filter(pl.col("REF_REF").is_not_null())

    # ── Mode C Unknown 反推 ──
    if mapping.allele_mode == AlleleMode.A1_A2 and mapping.effect_on == "Unknown":
        df, unresolved_count = _resolve_mode_c(df)
        stats.dropped_unresolved_mode_c = unresolved_count

    # ── 逐行对齐判定 ──
    # 使用 polars map_elements 进行判定
    result = df.with_columns(
        pl.struct([
            "_Aeff", "_Aref", "REF_REF", "REF_ALT", "REF_AF",
            "_EAF", "_STAT",
        ]).map_elements(
            _classify_alignment,
            return_dtype=pl.Utf8,
        ).alias("_ALIGN_ACTION")
    )

    # 统计各种结果
    action_counts = result.group_by("_ALIGN_ACTION").len().to_dict(as_series=False)
    action_map = dict(zip(action_counts["_ALIGN_ACTION"], action_counts["len"]))
    stats.matched = action_map.get("match", 0)
    stats.swapped = action_map.get("swap", 0)
    stats.flipped = action_map.get("flip", 0)
    stats.flipped_swapped = action_map.get("flip_swap", 0)
    stats.dropped_mismatch = action_map.get("drop_mismatch", 0)
    stats.dropped_ambiguous_palindromic = action_map.get("drop_ambiguous", 0)

    # ── 应用对齐操作 ──
    # 过滤掉需要剔除的行
    result = result.filter(
        pl.col("_ALIGN_ACTION").is_in(["match", "swap", "flip", "flip_swap"])
    )

    # 对 swap 和 flip_swap 的行执行 BETA 取反和 EAF 翻转
    result = result.with_columns([
        pl.when(pl.col("_ALIGN_ACTION").is_in(["swap", "flip_swap"]))
          .then(-pl.col("_STAT"))
          .otherwise(pl.col("_STAT"))
          .alias("_STAT"),
        pl.when(pl.col("_ALIGN_ACTION").is_in(["swap", "flip_swap"]))
          .then(
              pl.when(pl.col("_EAF").is_not_null())
                .then(1.0 - pl.col("_EAF"))
                .otherwise(pl.col("_EAF"))
          )
          .otherwise(pl.col("_EAF"))
          .alias("_EAF"),
    ])

    if "_MAF" in result.columns:
        result = result.with_columns(
            pl.when(pl.col("_EAF").is_not_null())
              .then(pl.col("_EAF"))
              .otherwise(
                  pl.when(pl.col("_MAF").is_not_null() & pl.col("REF_AF").is_not_null())
                    .then(
                        pl.when(pl.col("REF_AF") <= 0.5)
                          .then(pl.col("_MAF"))
                          .otherwise(1.0 - pl.col("_MAF"))
                    )
                    .otherwise(pl.lit(None).cast(pl.Float64))
              )
              .alias("_EAF")
        )

    stats.total_output = result.height
    return result, stats


def _classify_alignment(row: dict) -> str:
    """
    对单行数据进行对齐分类。
    返回: match, swap, flip, flip_swap, drop_mismatch, drop_ambiguous
    """
    a_eff = str(row["_Aeff"]).upper()
    a_ref = str(row["_Aref"]).upper()
    p_ref = str(row["REF_REF"]).upper()
    p_alt = str(row["REF_ALT"]).upper()
    ref_af = row["REF_AF"]
    eaf = row["_EAF"]

    # 检查是否为回文位点
    palindromic = is_palindromic(a_eff, a_ref)

    if palindromic:
        return _handle_palindromic(a_eff, a_ref, p_ref, p_alt, ref_af, eaf)

    # ── 直接匹配 ──
    if a_ref == p_ref and a_eff == p_alt:
        return "match"

    # ── 交换 ──
    if a_ref == p_alt and a_eff == p_ref:
        return "swap"

    # ── 互补链匹配 ──
    c_eff = complement_base(a_eff)
    c_ref = complement_base(a_ref)

    if c_ref == p_ref and c_eff == p_alt:
        return "flip"

    if c_ref == p_alt and c_eff == p_ref:
        return "flip_swap"

    return "drop_mismatch"


def _handle_palindromic(
    a_eff: str, a_ref: str,
    p_ref: str, p_alt: str,
    ref_af: float, eaf: float,
) -> str:
    """处理回文位点 (A/T, C/G)"""
    if eaf is None or ref_af is None:
        return "drop_ambiguous"

    try:
        eaf = float(eaf)
        ref_af = float(ref_af)
    except (ValueError, TypeError):
        return "drop_ambiguous"

    # 中间频率，无法区分
    if abs(eaf - 0.5) < EAF_MID_THRESHOLD:
        return "drop_ambiguous"

    # 方向一致
    if abs(eaf - ref_af) < EAF_THRESHOLD:
        if a_ref == p_ref and a_eff == p_alt:
            return "match"
        elif a_ref == p_alt and a_eff == p_ref:
            return "swap"

    # 方向相反
    if abs((1 - eaf) - ref_af) < EAF_THRESHOLD:
        if a_ref == p_ref and a_eff == p_alt:
            return "swap"
        elif a_ref == p_alt and a_eff == p_ref:
            return "match"

    return "drop_ambiguous"


def _resolve_mode_c(df: pl.DataFrame) -> Tuple[pl.DataFrame, int]:
    """
    Mode C Unknown 情况下，通过参考面板频率反推效应等位基因方向。
    """
    df = df.with_columns(
        pl.struct(["_Aeff", "_Aref", "_EAF", "REF_AF", "_STAT"])
        .map_elements(
            lambda row: _mode_c_resolve_row(row),
            return_dtype=pl.Struct({
                "_Aeff": pl.Utf8,
                "_Aref": pl.Utf8,
                "_EAF": pl.Float64,
                "_STAT": pl.Float64,
                "_resolved": pl.Boolean,
            }),
        )
        .alias("_resolved_struct")
    )

    unresolved_mask = ~pl.col("_resolved_struct").struct.field("_resolved")
    unresolved_count = df.filter(unresolved_mask).height

    df = df.filter(~unresolved_mask).with_columns([
        pl.col("_resolved_struct").struct.field("_Aeff").alias("_Aeff"),
        pl.col("_resolved_struct").struct.field("_Aref").alias("_Aref"),
        pl.col("_resolved_struct").struct.field("_EAF").alias("_EAF"),
        pl.col("_resolved_struct").struct.field("_STAT").alias("_STAT"),
    ]).drop("_resolved_struct")

    return df, unresolved_count


def _mode_c_resolve_row(row: dict) -> dict:
    """单行 Mode C 反推"""
    a_eff = str(row["_Aeff"])
    a_ref = str(row["_Aref"])
    eaf = row["_EAF"]
    ref_af = row["REF_AF"]
    stat = row["_STAT"]

    if eaf is None or ref_af is None:
        return {"_Aeff": a_eff, "_Aref": a_ref, "_EAF": eaf, "_STAT": stat, "_resolved": False}

    try:
        eaf = float(eaf)
        ref_af = float(ref_af)
    except (ValueError, TypeError):
        return {"_Aeff": a_eff, "_Aref": a_ref, "_EAF": eaf, "_STAT": stat, "_resolved": False}

    # 比较频率，判断 A1 vs A2 谁是效应碱基
    diff_a1 = abs(eaf - ref_af)
    diff_a2 = abs((1 - eaf) - ref_af)

    if diff_a1 < EAF_THRESHOLD and diff_a1 < diff_a2:
        # A1 的频率更接近 REF_ALT 的频率 → A1 就是效应等位基因，保持不变
        return {"_Aeff": a_eff, "_Aref": a_ref, "_EAF": eaf, "_STAT": stat, "_resolved": True}
    elif diff_a2 < EAF_THRESHOLD and diff_a2 < diff_a1:
        # A2 的频率更接近 → 翻转
        return {"_Aeff": a_ref, "_Aref": a_eff, "_EAF": 1 - eaf, "_STAT": -stat, "_resolved": True}
    else:
        # 无法区分
        return {"_Aeff": a_eff, "_Aref": a_ref, "_EAF": eaf, "_STAT": stat, "_resolved": False}


# ─────────────────────────────────────────────────────────────
# 输出模块
# ─────────────────────────────────────────────────────────────
def write_output(
    df: pl.DataFrame,
    output_path: str,
    mapping: ColumnMapping,
    output_format: str = "standardized",
    phenotype: Optional[str] = None,
    sample_size: Optional[int] = None,
):
    """写出标准化结果文件"""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    bim_id = pl.col("REF_BIM_ID")
    standardized_id = (
        pl.col("_CHR")
        + pl.lit(":")
        + pl.col("_POS")
        + pl.lit(":")
        + pl.col("REF_REF")
        + pl.lit(":")
        + pl.col("REF_ALT")
    )
    rsid_expr: Optional[pl.Expr] = None
    if "_SNP" in df.columns:
        rsid_expr = pl.col("_SNP")
    elif "REF_RSID" in df.columns:
        rsid_expr = pl.col("REF_RSID")

    # 构建最终输出列
    if output_format == "mr_raw":
        out_cols = [
            bim_id.alias("bim_id"),
            standardized_id.alias("variant_id"),
            pl.col("_CHR").alias("chromosome"),
            pl.col("_POS").alias("base_pair_location"),
            pl.col("REF_ALT").alias("effect_allele"),
            pl.col("REF_REF").alias("other_allele"),
            pl.col("_EAF").alias("effect_allele_frequency"),
            pl.col("_STAT").alias("beta"),
            pl.col("_SE").alias("standard_error"),
            pl.col("_P").alias("p_value"),
        ]
        if rsid_expr is not None:
            out_cols.append(rsid_expr.alias("rsid"))
        if phenotype is not None:
            out_cols.append(pl.lit(phenotype).alias("phenotype"))
        if sample_size is not None:
            out_cols.append(pl.lit(sample_size).alias("N"))
    else:
        out_cols = [
            bim_id.alias("BIM_ID"),
            pl.col("_CHR").alias("CHR"),
            pl.col("_POS").alias("POS"),
            standardized_id.alias("VARIANT_ID"),
        ]
        if rsid_expr is not None:
            out_cols.append(rsid_expr.alias("RSID"))
        out_cols.extend([
            pl.col("REF_REF").alias("REF"),
            pl.col("REF_ALT").alias("ALT"),
            pl.col("_STAT").alias("BETA"),
            pl.col("_SE").alias("SE"),
            pl.col("_P").alias("P"),
            pl.col("_EAF").alias("EAF"),
        ])

    output = df.select(out_cols)

    # 写出
    if output_path.endswith(".gz"):
        tmp_output_path = output_path[:-3]
        output.write_csv(tmp_output_path, separator="\t")
        with open(tmp_output_path, "rb") as src, gzip.open(output_path, "wb") as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
        os.remove(tmp_output_path)
    else:
        output.write_csv(output_path, separator="\t")

    console.print(f"\n[green]✅ 标准化文件已输出: {output_path}[/green]")
    console.print(f"   共 {output.height:,} 个变异位点\n")


def format_mr_output(
    raw_output_path: str,
    final_output_path: str,
    mr_role: str,
    r_lib_path: str,
    phenotype: Optional[str],
    sample_size: Optional[int],
    clump_r2: float,
    clump_kb: int,
    clump_p1: float,
    clump_pop: str,
    clump_plink: str,
    clump_bfile: str,
):
    """调用 R/TwoSampleMR 将中间表转换为最终 MR 输入，并在 exposure 模式下执行 clump"""
    role = normalize_mr_role(mr_role)
    script_path = Path(__file__).with_name("mr_format_and_clump.R")

    if not script_path.exists():
        console.print(f"[red]缺少 R 后处理脚本: {script_path}[/red]")
        sys.exit(1)

    Path(final_output_path).parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "Rscript",
        str(script_path),
        "--r-lib-path", r_lib_path,
        "--input", raw_output_path,
        "--output", final_output_path,
        "--mr-role", role,
        "--clump-r2", str(clump_r2),
        "--clump-kb", str(clump_kb),
        "--clump-p1", str(clump_p1),
        "--clump-pop", clump_pop,
        "--plink", clump_plink,
        "--bfile", clump_bfile,
    ]

    if phenotype is not None:
        cmd.extend(["--phenotype", phenotype])
    if sample_size is not None:
        cmd.extend(["--sample-size", str(sample_size)])

    console.print(
        f"[cyan]🔧 正在调用 TwoSampleMR 生成 {role} 数据"
        + (" 并执行 clump..." if role == "exposure" else "...") 
        + "[/cyan]"
    )

    try:
        completed = subprocess.run(
            cmd,
            check=True,
            text=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as exc:
        if exc.stdout:
            console.print(exc.stdout)
        if exc.stderr:
            console.print(f"[red]{exc.stderr}[/red]")
        sys.exit(1)

    if completed.stdout.strip():
        console.print(completed.stdout.strip())
    if completed.stderr.strip():
        console.print(f"[yellow]{completed.stderr.strip()}[/yellow]")


def generate_report(
    stats: AlignmentStats,
    qc_failed: int,
    output_path: str,
    input_path: str,
    ref_path: str,
):
    """生成审计报告"""
    report_path = (
        output_path
        .replace(".tsv.gz", ".log")
        .replace(".csv.gz", ".log")
        .replace(".tsv", ".log")
        .replace(".csv", ".log")
    )
    if report_path == output_path:
        report_path = output_path + ".log"

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines = [
        "=" * 60,
        "GWAS Standardizer - 对齐审计报告",
        "=" * 60,
        f"时间: {timestamp}",
        f"输入文件: {input_path}",
        f"参考面板: {ref_path}",
        f"输出文件: {output_path}",
        "-" * 60,
        f"初始输入变异总数:           {stats.total_input:>10,}",
        f"QC 不达标剔除:              {qc_failed:>10,}",
        f"去重剔除:                   {stats.duplicates_removed:>10,}",
        f"参考面板未找到:             {stats.not_in_ref:>10,}",
        f"Mode C 未解析剔除:          {stats.dropped_unresolved_mode_c:>10,}",
        "-" * 60,
        f"对齐成功 (方向一致):        {stats.matched:>10,}",
        f"效应反转 (Swapped):         {stats.swapped:>10,}",
        f"链翻转 (Strand Flipped):    {stats.flipped:>10,}",
        f"链翻转+效应反转:            {stats.flipped_swapped:>10,}",
        f"碱基不匹配剔除:             {stats.dropped_mismatch:>10,}",
        f"回文位点歧义剔除:           {stats.dropped_ambiguous_palindromic:>10,}",
        "-" * 60,
        f"最终输出变异总数:           {stats.total_output:>10,}",
        "=" * 60,
    ]

    report_text = "\n".join(lines)

    # 打印到终端
    console.print(Panel(report_text, title="📊 对齐审计报告", border_style="green"))

    # 写入日志文件
    with open(report_path, "w") as f:
        f.write(report_text + "\n")

    console.print(f"[dim]日志已保存至: {report_path}[/dim]")


# ─────────────────────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────────────────────
def main():
    """主程序入口"""
    args = parse_args()

    console.print(Panel.fit(
        "[bold magenta]🧬 GWAS Summary Statistics Standardizer[/bold magenta]\n"
        "[dim]GWAS 摘要统计数据标准化与等位基因对齐工具[/dim]\n"
        "[dim]Version 1.0 | hg19[/dim]",
        border_style="magenta",
    ))
    console.print()

    # ── Step 1: 获取文件路径 ──
    console.print("[bold yellow]━━━ 文件路径配置 ━━━[/bold yellow]")
    cli_mode = is_cli_mode(args)
    output_format = args.output_format
    mr_role = args.mr_role
    phenotype = args.phenotype
    sample_size = args.sample_size
    r_lib_path = args.r_lib_path
    clump_r2 = args.clump_r2
    clump_kb = args.clump_kb
    clump_p1 = args.clump_p1
    clump_pop = args.clump_pop
    clump_plink = args.clump_plink
    clump_bfile = args.clump_bfile

    if cli_mode:
        mapping = build_mapping_from_args(args)
        input_path = args.input or ""
        ref_path = args.reference or ""
        if output_format == "mr" and mr_role is None:
            console.print("[red]--output-format mr 时必须提供 --mr-role out/exp[/red]")
            sys.exit(1)

        output_path = args.output or derive_default_output_path(
            input_path,
            output_format=output_format,
            mr_role=mr_role,
        )
        if not os.path.isfile(input_path):
            console.print(f"[red]输入文件不存在: {input_path}[/red]")
            sys.exit(1)
        if not os.path.isfile(ref_path):
            console.print(f"[red]参考面板不存在: {ref_path}[/red]")
            sys.exit(1)
    else:
        input_path = ask_input_path()
        ref_path = ask_reference_path()
        runtime = interactive_runtime_options(input_path)
        output_format = str(runtime["output_format"])
        mr_role = runtime["mr_role"]
        phenotype = runtime["phenotype"]
        sample_size = runtime["sample_size"]
        r_lib_path = str(runtime["r_lib_path"])
        clump_r2 = float(runtime["clump_r2"])
        clump_kb = int(runtime["clump_kb"])
        clump_p1 = float(runtime["clump_p1"])
        clump_pop = str(runtime["clump_pop"])
        clump_plink = str(runtime["clump_plink"])
        clump_bfile = str(runtime["clump_bfile"])
        output_default = derive_default_output_path(input_path, output_format=output_format, mr_role=mr_role)
        output_input = ask_output_path(output_default)
        output_path = output_input if output_input else output_default

    console.print(f"  📂 输出文件: [green]{output_path}[/green]\n")

    # ── Step 2: 读取预览 & 交互映射 ──
    console.print("[cyan]🔧 正在读取输入文件...[/cyan]")
    df_preview = read_file_header(input_path)

    if cli_mode:
        display_mapping_summary(mapping)
    else:
        mapping = interactive_mapping(df_preview)

    mapping.separator = detect_separator(input_path)

    # ── Step 3: 加载完整数据 ──
    console.print("[cyan]🔧 正在加载完整数据...[/cyan]")
    df_full = pl.read_csv(
        input_path,
        separator=mapping.separator,
        infer_schema_length=10000,
        truncate_ragged_lines=True,
        ignore_errors=True,
    )
    console.print(f"  ✅ 已加载 {df_full.height:,} 行 × {df_full.width} 列\n")

    # ── Step 4: QC & 转换 ──
    console.print("[cyan]🔧 正在执行数据转换 & QC...[/cyan]")
    df_clean, qc_failed = transform_and_qc(df_full, mapping)
    console.print(f"  ✅ QC 通过 {df_clean.height:,} 行 (剔除 {qc_failed:,} 行)\n")

    # ── Step 5: 解析效应等位基因 ──
    console.print("[cyan]🔧 正在解析效应等位基因方向...[/cyan]")
    df_clean = resolve_effect_allele(df_clean, mapping)

    # ── Step 6: 加载参考面板 ──
    ref = load_reference_panel(ref_path)

    # ── Step 7: 等位基因对齐 ──
    console.print("[cyan]🔧 正在执行等位基因对齐...[/cyan]")
    df_aligned, stats = align_alleles(df_clean, ref, mapping)
    stats.qc_failed = qc_failed

    # ── Step 8: 输出 ──
    if output_format == "mr":
        console.print("[cyan]🔧 正在写出 MR 中间结果...[/cyan]")
        tmp_file = tempfile.NamedTemporaryFile(
            suffix=".mr_raw.tsv.gz",
            prefix="standardlizer_",
            dir=str(Path(output_path).parent),
            delete=False,
        )
        tmp_output_path = tmp_file.name
        tmp_file.close()

        try:
            write_output(
                df_aligned,
                tmp_output_path,
                mapping,
                output_format="mr_raw",
                phenotype=phenotype,
                sample_size=sample_size,
            )
            format_mr_output(
                tmp_output_path,
                output_path,
                mr_role=mr_role or "outcome",
                r_lib_path=r_lib_path,
                phenotype=phenotype,
                sample_size=sample_size,
                clump_r2=clump_r2,
                clump_kb=clump_kb,
                clump_p1=clump_p1,
                clump_pop=clump_pop,
                clump_plink=clump_plink,
                clump_bfile=clump_bfile,
            )
        finally:
            if os.path.exists(tmp_output_path):
                os.remove(tmp_output_path)
    else:
        console.print("[cyan]🔧 正在写出标准化结果...[/cyan]")
        write_output(
            df_aligned,
            output_path,
            mapping,
            output_format="standardized",
            phenotype=phenotype,
            sample_size=sample_size,
        )

    # ── Step 9: 报告 ──
    generate_report(stats, qc_failed, output_path, input_path, ref_path)

    console.print("[bold green]🎉 标准化流程完成！[/bold green]\n")


if __name__ == "__main__":
    main()
