#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把标准化 GWAS TSV 转换为 TwoSampleMR exposure/outcome 输入。"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

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
    from rich.panel import Panel
except ImportError:
    sys.exit("[ERROR] 请安装 rich: pip install rich")


console = Console()
QUESTIONARY_FALLBACK_WARNED = False

CUSTOM_STYLE = Style([
    ("qmark", "fg:cyan bold"),
    ("question", "fg:white bold"),
    ("answer", "fg:green bold"),
    ("pointer", "fg:cyan bold"),
    ("highlighted", "fg:cyan bold"),
    ("selected", "fg:green"),
])

DEFAULT_R_LIB_PATH = os.environ.get("MR_PIPELINE_R_LIB_PATH", "/home/ding/R/4.4.1_MR")
DEFAULT_STANDARDIZED_DIR = os.environ.get("MR_PIPELINE_STANDARDIZED_OUTPUT_DIR", "")
DEFAULT_EXP_DIR = os.environ.get("MR_PIPELINE_EXP_DIR", os.environ.get("MR_PIPELINE_EXPOSURE_DIR", ""))
DEFAULT_OUT_DIR = os.environ.get("MR_PIPELINE_OUT_DIR", os.environ.get("MR_PIPELINE_OUTCOME_DIR", ""))
DEFAULT_CLUMP_PLINK = os.environ.get("MR_PIPELINE_CLUMP_PLINK", "/home/ding/miniconda3/envs/GWAS/bin/plink")
DEFAULT_CLUMP_BFILE = os.environ.get("MR_PIPELINE_CLUMP_BFILE", "/home/ding/MR_LPA/Ref/g1000_eur/g1000_eur_colon")
DEFAULT_CLUMP_R2 = float(os.environ.get("MR_PIPELINE_CLUMP_R2", "0.1"))
DEFAULT_CLUMP_KB = int(os.environ.get("MR_PIPELINE_CLUMP_KB", "500"))
DEFAULT_CLUMP_P1 = float(os.environ.get("MR_PIPELINE_CLUMP_P1", "1e-4"))
DEFAULT_CLUMP_POP = os.environ.get("MR_PIPELINE_CLUMP_POP", "EUR")
DEFAULT_EAF_THRESHOLD_RAW = os.environ.get("MR_PIPELINE_FORMAT_EAF_THRESHOLD", "").strip()
DEFAULT_EAF_THRESHOLD = DEFAULT_EAF_THRESHOLD_RAW if DEFAULT_EAF_THRESHOLD_RAW else None


def _cancel_interaction() -> None:
    console.print("[red]操作已取消[/red]")
    sys.exit(0)


def _warn_questionary_fallback(exc: Exception) -> None:
    global QUESTIONARY_FALLBACK_WARNED
    if QUESTIONARY_FALLBACK_WARNED:
        return
    console.print("[yellow]检测到 questionary 菜单不可用，已自动回退为纯文本输入模式。[/yellow]")
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


def _is_separator(choice: object) -> bool:
    return choice.__class__.__name__.lower() == "separator"


def _normalize_choice(choice: object) -> Tuple[str, object]:
    if _is_separator(choice):
        return "", None
    if hasattr(choice, "value"):
        label = getattr(choice, "title", None)
        value = choice.value
        return str(label if label is not None else value), value
    return str(choice), choice


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="标准化 GWAS -> TwoSampleMR exposure/outcome 转换工具")
    parser.add_argument("--input", help="data/standardized 下的标准化 TSV/TSV.GZ 文件")
    parser.add_argument("--output", help="输出 CSV 路径；如果是目录则自动生成文件名")
    parser.add_argument("--role", "--mr-role", choices=["exp", "exposure", "out", "outcome"], help="输出角色")
    parser.add_argument("--phenotype", help="TwoSampleMR phenotype 名称")
    parser.add_argument("--sample-size", type=int, help="样本量 N，可选")
    parser.add_argument("--r-lib-path", default=DEFAULT_R_LIB_PATH, help="R 包库路径")
    parser.add_argument("--clump-r2", type=float, default=DEFAULT_CLUMP_R2)
    parser.add_argument("--clump-kb", type=int, default=DEFAULT_CLUMP_KB)
    parser.add_argument("--clump-p1", type=float, default=DEFAULT_CLUMP_P1)
    parser.add_argument("--clump-pop", default=DEFAULT_CLUMP_POP)
    parser.add_argument("--clump-plink", default=DEFAULT_CLUMP_PLINK)
    parser.add_argument("--clump-bfile", default=DEFAULT_CLUMP_BFILE)
    parser.add_argument("--region", help="可选区域筛选，仅 exposure 建议使用；格式 CHR:START-END")
    parser.add_argument("--region-chr", help="可选区域染色体")
    parser.add_argument("--region-start", type=int, help="可选区域起点")
    parser.add_argument("--region-end", type=int, help="可选区域终点")
    parser.add_argument("--use-maf-as-eaf", action="store_true", help="当 EAF 缺失时用 MAF 填充 EAF")
    parser.add_argument("--eaf-threshold", default=DEFAULT_EAF_THRESHOLD, help="仅 exposure 使用；按对称阈值过滤 EAF，保留 threshold <= EAF <= 1-threshold")
    parser.add_argument("--non-interactive", action="store_true", help="缺参数时直接报错，不进入交互")
    return parser.parse_args()


def normalize_role(role: str) -> str:
    role = role.lower()
    if role in {"exp", "exposure"}:
        return "exposure"
    if role in {"out", "outcome"}:
        return "outcome"
    raise ValueError(f"Unsupported role: {role}")


def output_suffix(role: str) -> str:
    return "_exposure.csv" if role == "exposure" else "_outcome.csv"


def input_base(input_path: str) -> str:
    name = Path(input_path).name
    for suffix in (".tsv.gz", ".csv.gz", ".tsv", ".csv", ".gz"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.endswith("_standardized"):
        name = name[: -len("_standardized")]
    return name


def resolve_output(input_path: str, output: Optional[str], role: str) -> str:
    default_dir = DEFAULT_EXP_DIR if role == "exposure" else DEFAULT_OUT_DIR
    default_dir = default_dir or str(Path(input_path).parent)
    default_name = f"{input_base(input_path)}{output_suffix(role)}"

    if not output:
        return str(Path(default_dir) / default_name)

    output = os.path.expanduser(output.strip())
    if output.endswith(os.sep) or os.path.isdir(output):
        return str(Path(output) / default_name)
    return output


def list_standardized_files() -> List[str]:
    if not DEFAULT_STANDARDIZED_DIR or not os.path.isdir(DEFAULT_STANDARDIZED_DIR):
        return []
    suffixes = (".tsv.gz", ".tsv", ".csv.gz", ".csv")
    return sorted(
        str(path)
        for path in Path(DEFAULT_STANDARDIZED_DIR).iterdir()
        if path.is_file() and path.name.endswith(suffixes)
    )


def ask_text(prompt: str, default: str = "") -> str:
    try:
        answer = questionary.text(prompt, default=default, style=CUSTOM_STYLE).ask()
        if answer is None:
            _cancel_interaction()
        return answer.strip()
    except Exception as exc:
        _warn_questionary_fallback(exc)
        return _plain_text_input(prompt, default=default)


def ask_confirm(prompt: str, default: bool = False) -> bool:
    try:
        answer = questionary.confirm(prompt, default=default, style=CUSTOM_STYLE).ask()
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


def ask_select(prompt: str, choices: List[object]):
    try:
        answer = questionary.select(prompt, choices=choices, style=CUSTOM_STYLE).ask()
        if answer is None:
            _cancel_interaction()
        return answer
    except Exception as exc:
        _warn_questionary_fallback(exc)

    normalized_choices = []
    for choice in choices:
        label, value = _normalize_choice(choice)
        if value is None:
            continue
        normalized_choices.append((label, value))

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


def parse_region(region: Optional[str], region_chr: Optional[str], region_start: Optional[int], region_end: Optional[int]) -> Optional[Tuple[str, int, int]]:
    if region:
        match = re.match(r"^(?:chr)?([^:]+):([0-9]+)-([0-9]+)$", region.strip(), re.IGNORECASE)
        if not match:
            raise ValueError("--region 格式应为 CHR:START-END，例如 1:55000000-56000000")
        region_chr, start, end = match.groups()
        region_start, region_end = int(start), int(end)

    provided = [region_chr is not None, region_start is not None, region_end is not None]
    if any(provided) and not all(provided):
        raise ValueError("区域参数必须同时提供 --region-chr/--region-start/--region-end，或直接提供 --region")
    if not any(provided):
        return None
    if int(region_start) > int(region_end):
        raise ValueError("region-start 不能大于 region-end")
    return str(region_chr).replace("chr", "").replace("CHR", ""), int(region_start), int(region_end)


def validate_eaf_threshold(eaf_threshold: Optional[float]) -> Optional[float]:
    if eaf_threshold is None:
        return None
    eaf_threshold = float(eaf_threshold)
    if not 0 <= eaf_threshold < 0.5:
        raise ValueError("EAF 阈值必须满足 0 <= threshold < 0.5")
    return eaf_threshold


def interactive_fill(args: argparse.Namespace) -> argparse.Namespace:
    if not args.input:
        candidates = list_standardized_files()
        if candidates:
            choices = candidates + [questionary.Separator(), "手动输入路径"]
            selected = ask_select("请选择标准化文件:", choices)
            args.input = ask_text("请输入标准化文件路径:") if selected == "手动输入路径" else selected
        else:
            args.input = ask_text("请输入标准化文件路径:", default=DEFAULT_STANDARDIZED_DIR)

    if not args.role:
        args.role = ask_select(
            "请选择输出角色:",
            [
                questionary.Choice("Exposure (可 clump)", value="exposure"),
                questionary.Choice("Outcome", value="outcome"),
            ],
        )

    role = normalize_role(args.role)
    default_pheno = input_base(args.input)
    if not args.phenotype:
        args.phenotype = ask_text("请输入 phenotype 名称:", default=default_pheno)

    if args.sample_size is None:
        sample_size = ask_text("请输入样本量 N (可留空):", default="")
        args.sample_size = int(sample_size) if sample_size else None

    if role == "exposure":
        eaf_default = "" if args.eaf_threshold is None else str(args.eaf_threshold)
        eaf_threshold = ask_text("请输入 EAF 过滤阈值 (可留空，不过滤):", default=eaf_default)
        args.eaf_threshold = float(eaf_threshold) if eaf_threshold else None
        if not any([args.region, args.region_chr, args.region_start, args.region_end]):
            if ask_confirm("是否先按基因/染色体区域筛选后 clump?", default=False):
                args.region = ask_text("请输入区域 (CHR:START-END):")
        args.clump_r2 = float(ask_text("clump_r2:", default=str(args.clump_r2)))
        args.clump_kb = int(ask_text("clump_kb:", default=str(args.clump_kb)))
        args.clump_p1 = float(ask_text("clump_p1:", default=str(args.clump_p1)))
        args.clump_pop = ask_text("pop:", default=args.clump_pop)
        args.clump_plink = ask_text("plink 路径:", default=args.clump_plink)
        args.clump_bfile = ask_text("clump 参考 bfile 前缀:", default=args.clump_bfile)

    if args.output is None:
        default_output = resolve_output(args.input, None, role)
        args.output = ask_text("请输入输出 CSV 路径:", default=default_output)

    return args


def build_command(args: argparse.Namespace, role: str, output: str, region: Optional[Tuple[str, int, int]]) -> List[str]:
    script_path = Path(__file__).with_name("mr_format_and_clump.R")
    cmd = [
        "Rscript", str(script_path),
        "--input-format", "standardized",
        "--r-lib-path", args.r_lib_path,
        "--input", args.input,
        "--output", output,
        "--mr-role", role,
        "--clump-r2", str(args.clump_r2),
        "--clump-kb", str(args.clump_kb),
        "--clump-p1", str(args.clump_p1),
        "--clump-pop", args.clump_pop,
        "--plink", args.clump_plink,
        "--bfile", args.clump_bfile,
        "--use-maf-as-eaf", "true" if args.use_maf_as_eaf else "false",
    ]
    if args.phenotype:
        cmd.extend(["--phenotype", args.phenotype])
    if args.sample_size is not None:
        cmd.extend(["--sample-size", str(args.sample_size)])
    if region:
        cmd.extend(["--region-chr", region[0], "--region-start", str(region[1]), "--region-end", str(region[2])])
    if args.eaf_threshold is not None:
        cmd.extend(["--eaf-threshold", str(args.eaf_threshold)])
    return cmd


def standardized_schema(input_path: str) -> List[str]:
    return pl.scan_csv(input_path, separator="\t", infer_schema_length=1000).collect_schema().names()


def build_tsmr_lazyframe(
    input_path: str,
    role: str,
    phenotype: str,
    sample_size: Optional[int],
    use_maf_as_eaf: bool,
    region: Optional[Tuple[str, int, int]] = None,
    p_threshold: Optional[float] = None,
    eaf_threshold: Optional[float] = None,
) -> pl.LazyFrame:
    schema = standardized_schema(input_path)
    required = {"BIM_ID", "VARIANT_ID", "CHR", "POS", "EFFECT_ALLELE", "OTHER_ALLELE", "BETA", "SE", "P"}
    missing = sorted(required - set(schema))
    if missing:
        raise ValueError(f"标准化文件缺少必要列: {', '.join(missing)}")

    lf = pl.scan_csv(
        input_path,
        separator="\t",
        infer_schema_length=10000,
        ignore_errors=True,
    )

    if region:
        region_chr, region_start, region_end = region
        lf = lf.filter(
            (pl.col("CHR").cast(pl.Utf8).str.replace("chr", "").str.replace("CHR", "") == str(region_chr))
            & (pl.col("POS").cast(pl.Int64, strict=False) >= int(region_start))
            & (pl.col("POS").cast(pl.Int64, strict=False) <= int(region_end))
        )

    if p_threshold is not None:
        lf = lf.filter(pl.col("P").cast(pl.Float64, strict=False) <= float(p_threshold))

    eaf_expr = pl.col("EAF").cast(pl.Float64, strict=False) if "EAF" in schema else pl.lit(None).cast(pl.Float64)
    if use_maf_as_eaf and "MAF" in schema:
        eaf_expr = pl.coalesce([eaf_expr, pl.col("MAF").cast(pl.Float64, strict=False)])
    lf = lf.with_columns(eaf_expr.alias("__EAF_VALUE"))

    if eaf_threshold is not None:
        threshold = float(eaf_threshold)
        lf = lf.filter(
            pl.col("__EAF_VALUE").is_not_null()
            & (pl.col("__EAF_VALUE") >= threshold)
            & (pl.col("__EAF_VALUE") <= (1.0 - threshold))
        )

    role_label = "exposure" if role == "exposure" else "outcome"
    out_cols = [
        pl.col("BIM_ID").cast(pl.Utf8).str.to_uppercase().alias("SNP"),
        pl.col("EFFECT_ALLELE").cast(pl.Utf8).str.to_uppercase().alias(f"effect_allele.{role_label}"),
        pl.col("OTHER_ALLELE").cast(pl.Utf8).str.to_uppercase().alias(f"other_allele.{role_label}"),
        pl.col("__EAF_VALUE").alias(f"eaf.{role_label}"),
        pl.col("BETA").cast(pl.Float64, strict=False).alias(f"beta.{role_label}"),
        pl.col("SE").cast(pl.Float64, strict=False).alias(f"se.{role_label}"),
        pl.col("P").cast(pl.Float64, strict=False).alias(f"pval.{role_label}"),
    ]
    if sample_size is not None:
        out_cols.append(pl.lit(sample_size).alias(f"samplesize.{role_label}"))
    out_cols.extend([
        pl.lit(phenotype).alias(role_label),
        pl.lit(True).alias(f"mr_keep.{role_label}"),
        pl.lit("reported").alias(f"pval_origin.{role_label}"),
        pl.lit(phenotype).alias(f"id.{role_label}"),
        pl.col("VARIANT_ID").cast(pl.Utf8).str.to_uppercase().alias("variant_id"),
    ])
    if "RSID" in schema:
        out_cols.append(pl.col("RSID").cast(pl.Utf8).alias("rsid"))
    else:
        out_cols.append(pl.lit(None).cast(pl.Utf8).alias("rsid"))

    return lf.select(out_cols)


def sink_csv(lf: pl.LazyFrame, output: str) -> None:
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    lf.sink_csv(output)


def write_empty_like(lf: pl.LazyFrame, output: str) -> None:
    lf.head(0).collect().write_csv(output)


def parse_clumped_snps(clumped_path: str) -> List[str]:
    if not os.path.exists(clumped_path) or os.path.getsize(clumped_path) == 0:
        return []
    with open(clumped_path, "r", encoding="utf-8", errors="replace") as handle:
        header: Optional[List[str]] = None
        snp_idx: Optional[int] = None
        snps: List[str] = []
        for line in handle:
            fields = line.strip().split()
            if not fields:
                continue
            if header is None:
                header = fields
                if "SNP" not in header:
                    return []
                snp_idx = header.index("SNP")
                continue
            if snp_idx is not None and len(fields) > snp_idx:
                snps.append(fields[snp_idx].upper())
    return snps


def run_plink_clump(args: argparse.Namespace, candidates_path: str) -> List[str]:
    with tempfile.TemporaryDirectory(prefix="format_tsmr_clump_") as tmp_dir:
        clump_input = str(Path(tmp_dir) / "clump_input.tsv")
        clump_out = str(Path(tmp_dir) / "clump")
        (
            pl.scan_csv(candidates_path)
            .select([
                pl.col("SNP"),
                pl.col("pval.exposure").alias("P"),
            ])
            .sink_csv(clump_input, separator="\t")
        )
        cmd = [
            args.clump_plink,
            "--bfile", args.clump_bfile,
            "--clump", clump_input,
            "--clump-snp-field", "SNP",
            "--clump-field", "P",
            "--clump-kb", str(args.clump_kb),
            "--clump-r2", str(args.clump_r2),
            "--clump-p1", str(args.clump_p1),
            "--out", clump_out,
        ]
        completed = subprocess.run(cmd, text=True, capture_output=True)
        if completed.stdout:
            console.print(completed.stdout.strip())
        if completed.stderr:
            console.print(completed.stderr.strip(), style="yellow" if completed.returncode == 0 else "red")
        if completed.returncode != 0:
            raise RuntimeError(f"PLINK clump failed with exit code {completed.returncode}")
        return parse_clumped_snps(clump_out + ".clumped")


def filter_candidates_by_snps(candidates_path: str, output: str, snps: Sequence[str]) -> int:
    candidates = pl.read_csv(candidates_path)
    if not snps:
        candidates.head(0).write_csv(output)
        return 0
    result = candidates.filter(pl.col("SNP").is_in(list(snps)))
    result.write_csv(output)
    return result.height


def convert_standardized(args: argparse.Namespace, role: str, output: str, region: Optional[Tuple[str, int, int]]) -> None:
    phenotype = args.phenotype or input_base(args.input)
    if role == "outcome":
        lf = build_tsmr_lazyframe(
            args.input,
            role=role,
            phenotype=phenotype,
            sample_size=args.sample_size,
            use_maf_as_eaf=args.use_maf_as_eaf,
            region=None,
            eaf_threshold=args.eaf_threshold,
        )
        sink_csv(lf, output)
        return

    console.print(f"按 p <= {args.clump_p1:g} 预筛 exposure 候选位点")
    candidate_lf = build_tsmr_lazyframe(
        args.input,
        role=role,
        phenotype=phenotype,
        sample_size=args.sample_size,
        use_maf_as_eaf=args.use_maf_as_eaf,
        region=region,
        p_threshold=args.clump_p1,
        eaf_threshold=args.eaf_threshold,
    )
    with tempfile.TemporaryDirectory(prefix="format_tsmr_candidates_") as tmp_dir:
        candidates_path = str(Path(tmp_dir) / "candidates.csv")
        sink_csv(candidate_lf, candidates_path)
        candidate_count = pl.scan_csv(candidates_path).select(pl.len()).collect().item()
        console.print(f"进入 clump 的候选位点: [green]{candidate_count:,}[/green]")
        if candidate_count == 0:
            write_empty_like(candidate_lf, output)
            return

        clumped_snps = run_plink_clump(args, candidates_path)
        final_count = filter_candidates_by_snps(candidates_path, output, clumped_snps)
        console.print(f"clump 后保留位点: [green]{final_count:,}[/green]")


def main() -> None:
    args = parse_args()
    cli_mode = bool(args.input and args.role)
    if not cli_mode:
        if args.non_interactive:
            console.print("[red]--non-interactive 模式必须提供 --input 和 --role[/red]")
            sys.exit(1)
        console.print(Panel.fit("标准化 GWAS -> TwoSampleMR 格式转换", border_style="cyan"))
        args = interactive_fill(args)

    if not args.input or not os.path.isfile(args.input):
        console.print(f"[red]输入文件不存在: {args.input}[/red]")
        sys.exit(1)

    role = normalize_role(args.role)
    region = parse_region(args.region, args.region_chr, args.region_start, args.region_end)
    try:
        args.eaf_threshold = validate_eaf_threshold(args.eaf_threshold)
    except ValueError as exc:
        console.print(f"[red]{exc}[/red]")
        sys.exit(1)
    if role == "outcome" and args.eaf_threshold is not None:
        console.print("[yellow]已忽略 --eaf-threshold：EAF 过滤仅在 exposure 转换时生效。[/yellow]")
        args.eaf_threshold = None
    if region and role != "exposure":
        console.print("[red]区域筛选/clump 只建议用于 exposure；outcome 不支持区域参数。[/red]")
        sys.exit(1)

    output = resolve_output(args.input, args.output, role)
    Path(output).parent.mkdir(parents=True, exist_ok=True)

    console.print(f"输入文件: [green]{args.input}[/green]")
    console.print(f"输出角色: [green]{role}[/green]")
    console.print(f"输出文件: [green]{output}[/green]")
    if region:
        console.print(f"区域筛选: [green]{region[0]}:{region[1]}-{region[2]}[/green]")
    if args.eaf_threshold is not None:
        console.print(f"EAF 过滤: [green]{args.eaf_threshold:g} <= EAF <= {1 - args.eaf_threshold:g}[/green]")

    try:
        convert_standardized(args, role, output, region)
    except Exception as exc:
        console.print(f"[red]{exc}[/red]")
        sys.exit(1)

    console.print(f"转换完成: [green]{output}[/green]")


if __name__ == "__main__":
    main()
