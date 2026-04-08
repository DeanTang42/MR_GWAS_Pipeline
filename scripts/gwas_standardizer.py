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
    - canonical ID 生成
    - 质量控制（OR→BETA, -log10P→P, 过滤, 去重）

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
import tempfile
from enum import Enum
from dataclasses import dataclass, field
from typing import Optional, List, Tuple
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

QUESTIONARY_FALLBACK_WARNED = False


def _cancel_interaction() -> None:
    console.print("[red]操作已取消[/red]")
    sys.exit(0)


def _warn_questionary_fallback(exc: Exception) -> None:
    global QUESTIONARY_FALLBACK_WARNED
    if QUESTIONARY_FALLBACK_WARNED:
        return
    console.print(
        "[yellow]检测到 questionary 菜单不可用，已自动回退为纯文本输入模式。[/yellow]"
    )
    console.print(f"[dim]{exc.__class__.__name__}: {exc}[/dim]")
    QUESTIONARY_FALLBACK_WARNED = True


def _plain_text_input(prompt: str, default: str = "") -> str:
    prompt_text = prompt
    if default:
        prompt_text += f" [默认: {default}]"
    prompt_text += " "
    try:
        answer = input(prompt_text)
    except (EOFError, KeyboardInterrupt):
        _cancel_interaction()
    answer = answer.strip()
    return answer or default


def _normalize_choice(choice) -> Tuple[str, object]:
    if hasattr(choice, "value"):
        label = getattr(choice, "title", None)
        value = choice.value
        return str(label if label is not None else value), value
    return str(choice), choice


def ask_select(prompt: str, choices: List[object]):
    try:
        answer = questionary.select(
            prompt,
            choices=choices,
            style=CUSTOM_STYLE,
        ).ask()
        if answer is None:
            _cancel_interaction()
        return answer
    except Exception as exc:
        _warn_questionary_fallback(exc)

    normalized_choices = [_normalize_choice(choice) for choice in choices]
    console.print(f"[bold cyan]{prompt}[/bold cyan]")
    for idx, (label, _) in enumerate(normalized_choices, start=1):
        console.print(f"  {idx}. {label}")

    while True:
        raw = _plain_text_input("请输入编号:", default="1")
        try:
            selected_idx = int(raw)
        except ValueError:
            console.print(f"[red]请输入 1 到 {len(normalized_choices)} 之间的数字。[/red]")
            continue
        if 1 <= selected_idx <= len(normalized_choices):
            return normalized_choices[selected_idx - 1][1]
        console.print(f"[red]请输入 1 到 {len(normalized_choices)} 之间的数字。[/red]")


def ask_confirm(prompt: str, default: bool = True) -> bool:
    try:
        answer = questionary.confirm(
            prompt,
            default=default,
            style=CUSTOM_STYLE,
        ).ask()
        if answer is None:
            _cancel_interaction()
        return bool(answer)
    except Exception as exc:
        _warn_questionary_fallback(exc)

    default_value = "y" if default else "n"
    while True:
        answer = _plain_text_input(f"{prompt} [y/n]:", default=default_value).lower()
        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        console.print("[red]请输入 y 或 n。[/red]")


DEFAULT_ORG_DIR = os.environ.get("MR_PIPELINE_ORG_DIR", "")
DEFAULT_STANDARDIZED_OUTPUT_DIR = os.environ.get("MR_PIPELINE_STANDARDIZED_OUTPUT_DIR", "")


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
    effect_on: Optional[str] = None  # "A1", "A2"

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
class ProcessingStats:
    """标准化统计"""
    total_input: int = 0
    qc_failed: int = 0
    duplicates_removed: int = 0
    maf_without_eaf: int = 0
    eaf_missing: int = 0
    total_output: int = 0


# ─────────────────────────────────────────────────────────────
# 工具函数
# ─────────────────────────────────────────────────────────────
def make_sort_snpid(chrom: str, pos: str, a1: str, a2: str) -> str:
    """生成标准化 Sort_SNPID (碱基字母排序)"""
    a1, a2 = a1.upper(), a2.upper()
    if a1 <= a2:
        return f"{chrom}:{pos}:{a1}:{a2}"
    else:
        return f"{chrom}:{pos}:{a2}:{a1}"


def make_bim_id(chrom: str, pos: str, a1: str, a2: str) -> str:
    """生成与 BIM 一致的字典序 ID: CHR:POS:A1:A2"""
    return make_sort_snpid(chrom, pos, a1, a2)


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


def needs_whitespace_normalization(filepath: str, separator: str) -> bool:
    """识别表头和数据行空白分隔不一致的文件。"""
    if separator not in {"\t", " "}:
        return False

    open_fn = gzip.open if filepath.endswith(".gz") else open
    with open_fn(filepath, "rt") as f:
        header = f.readline().strip()
        data_line = ""
        for line in f:
            if line.strip():
                data_line = line.strip()
                break

    if not header or not data_line:
        return False

    header_fields = header.split()
    data_fields = data_line.split()
    if len(header_fields) <= 1 or len(data_fields) <= 1:
        return False

    if separator == "\t":
        header_tab_fields = header.split("\t")
        data_tab_fields = data_line.split("\t")
        return (
            len(header_tab_fields) > 1
            and len(data_tab_fields) == 1
            and len(data_fields) == len(header_fields)
        )

    header_space_fields = header.split(" ")
    data_space_fields = data_line.split(" ")
    return (
        len(header_space_fields) > 1
        and len(data_space_fields) == 1
        and len(data_fields) == len(header_fields)
    )


def normalize_whitespace_delimited_file(filepath: str, output_dir: str) -> str:
    """把任意空白分隔的 GWAS 文件规范化为 tab 分隔临时文件。"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    tmp_file = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        suffix=".normalized.tsv",
        prefix="standardlizer_input_",
        dir=output_dir,
        delete=False,
    )
    tmp_path = tmp_file.name
    open_fn = gzip.open if filepath.endswith(".gz") else open

    try:
        with tmp_file as out, open_fn(filepath, "rt") as src:
            for line in src:
                fields = line.strip().split()
                if not fields:
                    continue
                out.write("\t".join(fields) + "\n")
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise

    return tmp_path


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


def derive_default_output_path(
    input_path: str,
    output_dir: Optional[str] = None,
) -> str:
    """根据输入路径生成默认标准化输出路径"""
    input_base = get_input_base(input_path)
    suffix = "_standardized.tsv.gz"
    target_dir = output_dir or DEFAULT_STANDARDIZED_OUTPUT_DIR or None
    if target_dir:
        return str(Path(target_dir) / f"{Path(input_base).name}{suffix}")
    return input_base + suffix


def resolve_output_path(
    input_path: str,
    output_path: Optional[str],
) -> str:
    """解析输出路径；如果用户传入目录，则在目录内生成默认文件名。"""
    if output_path:
        output_path = os.path.expanduser(output_path.strip())
        if output_path.endswith(os.sep) or os.path.isdir(output_path):
            return derive_default_output_path(
                input_path,
                output_dir=output_path,
            )
        return output_path

    return derive_default_output_path(input_path)


def parse_args() -> argparse.Namespace:
    """解析命令行参数，支持交互式和非交互式两种模式"""
    parser = argparse.ArgumentParser(
        description="GWAS 摘要统计数据标准化工具",
    )
    parser.add_argument("--input", help="输入 GWAS 摘要统计文件")
    parser.add_argument("--output", help="输出文件路径")
    parser.add_argument(
        "--output-format",
        choices=["standardized"],
        default="standardized",
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--mode", choices=[m.value for m in AlleleMode], help="等位基因模式: A/B/C")
    parser.add_argument("--snp-col", help="SNP/variant ID 列名")
    parser.add_argument("--chr-col", help="染色体列名")
    parser.add_argument("--pos-col", help="位置列名")
    parser.add_argument("--allele1-col", help="等位基因 1 列名")
    parser.add_argument("--allele2-col", help="等位基因 2 列名")
    parser.add_argument("--effect-on", choices=["A1", "A2", "Unknown"], help="Mode C 时效应作用在哪个等位基因上；Unknown 当前不支持")
    parser.add_argument("--stat-type", choices=[s.value for s in StatType], help="统计量类型: BETA/OR")
    parser.add_argument("--stat-col", help="效应统计量列名")
    parser.add_argument("--se-col", help="标准误列名")
    parser.add_argument("--p-col", help="P 值列名")
    parser.add_argument("--pval-format", choices=[p.value for p in PvalFormat], help="P 值格式: raw/neglog10")
    parser.add_argument("--freq-type", choices=[f.value for f in FreqType], help="频率类型: EAF/MAF")
    parser.add_argument("--freq-col", help="频率列名")
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
    ]
    return args.non_interactive or any(value is not None for value in cli_fields)


def build_mapping_from_args(args: argparse.Namespace) -> ColumnMapping:
    """根据命令行参数构建列映射"""
    required = [
        ("--input", args.input),
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
        console.print("[red]Mode C 需要提供 --effect-on (A1/A2)[/red]")
        sys.exit(1)

    if args.mode == AlleleMode.A1_A2.value and args.effect_on == "Unknown":
        console.print("[red]Mode C 的 Unknown 方向当前无法自动解析，请显式指定 A1 或 A2。[/red]")
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
        try:
            path = questionary.path(
                prompt,
                default=default,
                style=CUSTOM_STYLE,
            ).ask()
            if path is None:
                _cancel_interaction()
        except Exception as exc:
            _warn_questionary_fallback(exc)
            path = _plain_text_input(prompt, default=default)
        path = os.path.expanduser(path.strip())
        if must_exist and not os.path.isfile(path):
            console.print(f"[red]文件不存在: {path}[/red]")
            continue
        return path


def ask_output_path(default_path: str, prompt: str = "📁 请输入输出文件路径:") -> str:
    """交互式询问输出路径"""
    path = ask_text(prompt, default=default_path)
    return path.strip()


def ask_text(prompt: str, default: str = "") -> str:
    """通用文本输入"""
    try:
        answer = questionary.text(
            prompt,
            default=default,
            style=CUSTOM_STYLE,
        ).ask()
        if answer is None:
            _cancel_interaction()
    except Exception as exc:
        _warn_questionary_fallback(exc)
        answer = _plain_text_input(prompt, default=default)
    return answer.strip()


def ask_input_path() -> str:
    """交互式获取输入 GWAS 路径，优先提示 Org 目录"""
    default_path = str(Path(DEFAULT_ORG_DIR)) if DEFAULT_ORG_DIR else ""
    return ask_file_path("📄 请输入 GWAS 摘要统计文件路径:", default=default_path)


def select_column(columns: List[str], prompt: str, allow_none: bool = False) -> Optional[str]:
    """从列名列表中选择一个列"""
    choices = list(columns)
    if allow_none:
        choices.append("⏭  跳过 (无此列)")

    answer = ask_select(
        prompt,
        choices,
    )
    if answer == "⏭  跳过 (无此列)":
        return None
    return answer


def interactive_mapping(df: pl.DataFrame) -> ColumnMapping:
    """交互式列映射主流程"""
    columns = df.columns
    mapping = ColumnMapping()

    console.print(Panel.fit(
        "[bold cyan]🧬 GWAS 摘要统计数据标准化工具[/bold cyan]\n"
        "[dim]交互式列名映射与标准化[/dim]",
        border_style="cyan",
    ))

    # 展示预览
    display_preview(df)
    console.print()

    # ── Step 1: 选择等位基因模式 ──
    console.print("[bold yellow]━━━ Step 1: 选择等位基因记录模式 ━━━[/bold yellow]")
    mode_answer = ask_select(
        "请选择输入文件的等位基因记录规范:",
        [
            questionary.Choice("Mode A: REF/ALT (VCF 规范, BETA 针对 ALT)", value="A"),
            questionary.Choice("Mode B: Effect/Other (明确效应等位基因 EA/OA)", value="B"),
            questionary.Choice("Mode C: Allele1/Allele2 (方向可能有二义性)", value="C"),
        ],
    )

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
        effect_answer = ask_select(
            "效应值是针对哪个等位基因计算的?",
            [
                questionary.Choice("A1 (Allele1)", value="A1"),
                questionary.Choice("A2 (Allele2)", value="A2"),
            ],
        )
        mapping.effect_on = effect_answer
    console.print()

    # ── Step 4: 统计量映射 ──
    console.print("[bold yellow]━━━ Step 4: 映射统计量列 ━━━[/bold yellow]")

    # BETA 或 OR
    stat_answer = ask_select(
        "效应统计量类型:",
        [
            questionary.Choice("BETA (效应值)", value="BETA"),
            questionary.Choice("OR (比值比, 将自动转换为 ln(OR))", value="OR"),
        ],
    )
    mapping.stat_type = StatType(stat_answer)
    mapping.stat_col = select_column(columns, f"🔹 {stat_answer} 列:")

    # SE
    mapping.se_col = select_column(columns, "🔹 SE (标准误) 列:")

    # P
    mapping.p_col = select_column(columns, "🔹 P (P 值) 列:")
    pval_answer = ask_select(
        "P 值格式:",
        [
            questionary.Choice("原始 P 值 (如 0.05)", value="raw"),
            questionary.Choice("-log10(P) 格式 (如 1.3)", value="neglog10"),
        ],
    )
    mapping.pval_format = PvalFormat(pval_answer)

    # EAF / MAF
    freq_answer = ask_select(
        "频率列类型:",
        [
            questionary.Choice("EAF (效应等位基因频率)", value="EAF"),
            questionary.Choice("MAF (次等位基因频率)", value="MAF"),
        ],
    )
    mapping.freq_type = FreqType(freq_answer)
    mapping.freq_col = select_column(columns, f"🔹 {freq_answer} 列:", allow_none=True)
    console.print()

    # ── Step 5: 确认 ──
    display_mapping_summary(mapping)

    confirm = ask_confirm("以上映射配置是否正确?", default=True)
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
        table.add_row("效应方向", mapping.effect_on or "(未设置)")

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
def resolve_effect_allele(df: pl.DataFrame, mapping: ColumnMapping) -> pl.DataFrame:
    """
    根据等位基因模式，统一列名为 _Aeff (效应等位基因) 和 _Aref (非效应等位基因)。
    同时生成 _EAF 列。
    """
    if mapping.allele_mode == AlleleMode.REF_ALT:
        # Mode A: ALT 是效应等位基因
        df = df.with_columns([
            pl.col("_A2").alias("_Aeff"),  # ALT = effect
            pl.col("_A1").alias("_Aref"),  # REF = non-effect
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
        else:
            console.print("[red]Mode C 的 Unknown 方向当前不支持。[/red]")
            sys.exit(1)

    return df


def canonicalize_variants(df: pl.DataFrame) -> Tuple[pl.DataFrame, ProcessingStats]:
    """生成 canonical ID，并按位点去重。"""
    stats = ProcessingStats(total_input=df.height)

    sort_a1 = (
        pl.when(pl.col("_Aeff") <= pl.col("_Aref"))
        .then(pl.col("_Aeff"))
        .otherwise(pl.col("_Aref"))
    )
    sort_a2 = (
        pl.when(pl.col("_Aeff") <= pl.col("_Aref"))
        .then(pl.col("_Aref"))
        .otherwise(pl.col("_Aeff"))
    )

    df = df.with_columns([
        sort_a1.alias("_SORT_A1"),
        sort_a2.alias("_SORT_A2"),
    ]).with_columns([
        (
            pl.col("_CHR")
            + pl.lit(":")
            + pl.col("_POS")
            + pl.lit(":")
            + pl.col("_SORT_A1")
            + pl.lit(":")
            + pl.col("_SORT_A2")
        ).alias("_SORT_SNPID")
    ]).with_columns(
        pl.col("_SORT_SNPID").alias("_BIM_ID")
    )

    before_dedup = df.height
    df = df.sort("_P").unique(subset=["_SORT_SNPID"], keep="first")
    stats.duplicates_removed = before_dedup - df.height
    stats.maf_without_eaf = df.filter(pl.col("_MAF").is_not_null() & pl.col("_EAF").is_null()).height
    stats.eaf_missing = df.filter(pl.col("_EAF").is_null()).height
    stats.total_output = df.height
    return df, stats


# ─────────────────────────────────────────────────────────────
# 输出模块
# ─────────────────────────────────────────────────────────────
def write_output(
    df: pl.DataFrame,
    output_path: str,
    mapping: ColumnMapping,
):
    """写出标准化结果文件"""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    bim_id = pl.col("_BIM_ID")
    standardized_id = pl.col("_SORT_SNPID")
    eaf_out = (
        pl.when(pl.col("_EAF").is_null() | pl.col("_EAF").is_nan())
        .then(pl.lit("NA"))
        .otherwise(pl.col("_EAF").cast(pl.Utf8))
        .alias("EAF")
    )
    maf_out = (
        pl.when(pl.col("_MAF").is_null() | pl.col("_MAF").is_nan())
        .then(pl.lit("NA"))
        .otherwise(pl.col("_MAF").cast(pl.Utf8))
        .alias("MAF")
    )
    rsid_expr: Optional[pl.Expr] = None
    if "_SNP" in df.columns:
        rsid_expr = pl.col("_SNP")

    out_cols = [
        bim_id.alias("BIM_ID"),
        pl.col("_CHR").alias("CHR"),
        pl.col("_POS").alias("POS"),
        standardized_id.alias("VARIANT_ID"),
    ]
    if rsid_expr is not None:
        out_cols.append(rsid_expr.alias("RSID"))
    out_cols.extend([
        pl.col("_Aeff").alias("EFFECT_ALLELE"),
        pl.col("_Aref").alias("OTHER_ALLELE"),
        pl.col("_STAT").alias("BETA"),
        pl.col("_SE").alias("SE"),
        pl.col("_P").alias("P"),
        eaf_out,
        maf_out,
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


def generate_report(
    stats: ProcessingStats,
    output_path: str,
    input_path: str,
    output_paths: Optional[List[str]] = None,
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

    if output_paths and len(output_paths) > 1:
        output_lines = ["输出文件:"]
        output_lines.extend([f"  - {path}" for path in output_paths])
    else:
        output_lines = [f"输出文件: {output_path}"]

    lines = [
        "=" * 60,
        "GWAS Standardizer 审计报告",
        "=" * 60,
        f"时间: {timestamp}",
        f"输入文件: {input_path}",
        *output_lines,
        "-" * 60,
        f"初始输入变异总数:           {stats.total_input:>10,}",
        f"QC 不达标剔除:              {stats.qc_failed:>10,}",
        f"去重剔除:                   {stats.duplicates_removed:>10,}",
        f"仅有 MAF 待后续注释:        {stats.maf_without_eaf:>10,}",
        f"EAF 缺失位点数:             {stats.eaf_missing:>10,}",
        "-" * 60,
        f"最终输出变异总数:           {stats.total_output:>10,}",
        "=" * 60,
    ]

    report_text = "\n".join(lines)

    # 打印到终端
    console.print(Panel(report_text, title="📊 标准化审计报告", border_style="green"))

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
        "[dim]GWAS 摘要统计数据标准化与 canonical 位点 ID 工具[/dim]\n"
        "[dim]Version 1.0 | hg19[/dim]",
        border_style="magenta",
    ))
    console.print()

    # ── Step 1: 获取文件路径 ──
    console.print("[bold yellow]━━━ 文件路径配置 ━━━[/bold yellow]")
    cli_mode = is_cli_mode(args)

    if cli_mode:
        mapping = build_mapping_from_args(args)
        input_path = args.input or ""
        output_path = resolve_output_path(input_path, args.output)
        if not os.path.isfile(input_path):
            console.print(f"[red]输入文件不存在: {input_path}[/red]")
            sys.exit(1)
    else:
        input_path = ask_input_path()
        output_default = derive_default_output_path(input_path)
        output_input = ask_output_path(output_default)
        output_path = resolve_output_path(input_path, output_input)

    console.print(f"  📂 输出文件: [green]{output_path}[/green]\n")

    read_input_path = input_path
    input_separator = detect_separator(input_path)
    normalized_input_path: Optional[str] = None
    if needs_whitespace_normalization(input_path, input_separator):
        console.print(
            "[yellow]⚠ 检测到表头和数据行使用了不同的空白分隔符，"
            "正在生成临时 tab 分隔输入文件...[/yellow]"
        )
        normalized_input_path = normalize_whitespace_delimited_file(
            input_path,
            str(Path(output_path).parent),
        )
        read_input_path = normalized_input_path
        input_separator = "\t"
        console.print(f"  ✅ 临时规范化文件: [dim]{normalized_input_path}[/dim]\n")

    # ── Step 2: 读取预览 & 交互映射 ──
    try:
        console.print("[cyan]🔧 正在读取输入文件...[/cyan]")
        df_preview = read_file_header(read_input_path)

        if cli_mode:
            display_mapping_summary(mapping)
        else:
            mapping = interactive_mapping(df_preview)

        mapping.separator = input_separator

        # ── Step 3: 加载完整数据 ──
        console.print("[cyan]🔧 正在加载完整数据...[/cyan]")
        df_full = pl.read_csv(
            read_input_path,
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

        # ── Step 6: 生成 canonical ID & 去重 ──
        console.print("[cyan]🔧 正在生成 canonical 位点 ID...[/cyan]")
        df_aligned, stats = canonicalize_variants(df_clean)
        stats.qc_failed = qc_failed

        if stats.maf_without_eaf > 0:
            console.print(
                f"[yellow]⚠ 检测到 {stats.maf_without_eaf:,} 个位点只有 MAF，"
                "标准化阶段不会直接把 MAF 转成 EAF；如需补充 EAF，请在后续运行 annotate_eaf.sh。[/yellow]\n"
            )

        # ── Step 7: 输出 ──
        console.print("[cyan]🔧 正在写出标准化结果...[/cyan]")
        write_output(
            df_aligned,
            output_path,
            mapping,
        )

        # ── Step 8: 报告 ──
        generate_report(
            stats,
            output_path,
            input_path,
        )
    finally:
        if normalized_input_path and os.path.exists(normalized_input_path):
            os.remove(normalized_input_path)

    console.print("[bold green]🎉 标准化流程完成！[/bold green]\n")


if __name__ == "__main__":
    main()
