#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用 1KG PLINK .frq 和原始 MAF 为标准化 GWAS 补充 EAF。"""

import argparse
import gzip
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

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

CUSTOM_STYLE = Style([
    ("qmark", "fg:cyan bold"),
    ("question", "fg:white bold"),
    ("answer", "fg:green bold"),
    ("pointer", "fg:cyan bold"),
    ("highlighted", "fg:cyan bold"),
    ("selected", "fg:green"),
])

DEFAULT_STANDARDIZED_DIR = os.environ.get("MR_PIPELINE_STANDARDIZED_OUTPUT_DIR", "")
DEFAULT_REFERENCE_DIR = os.environ.get("MR_PIPELINE_REFERENCE_DIR", "")
DEFAULT_FRQ_REFERENCE = os.environ.get("MR_PIPELINE_FRQ_REFERENCE", "")
LEGACY_AUDIT_COLUMNS = {"EAF_1KG", "EAF_FROM_MAF", "EAF_SOURCE", "EAF_ABS_DIFF"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="为标准化 GWAS 中缺失的 EAF 添加 MAF 定向注释")
    parser.add_argument("--input", help="标准化 GWAS TSV/TSV.GZ")
    parser.add_argument("--frq", default=DEFAULT_FRQ_REFERENCE, help="PLINK .frq 或规范化 .frq.tsv.gz")
    parser.add_argument("--output", help="完整输出文件；默认替换原 input")
    parser.add_argument("--matched-output", help="只包含 BIM_ID 匹配到参考的行的输出文件")
    parser.add_argument("--diff-output", help="记录 EAF_1KG 和 EAF_FROM_MAF 差异超过阈值的位点")
    parser.add_argument("--log", help="注释日志路径")
    parser.add_argument("--diff-threshold", type=float, default=0.01, help="频率差异记录阈值，默认 0.01")
    parser.add_argument("--replace-input", dest="replace_input", action="store_true", default=True, help="替换原输入文件，默认启用")
    parser.add_argument("--no-replace-input", dest="replace_input", action="store_false", help="不替换原输入，写到 --output 或默认派生文件")
    parser.add_argument("--non-interactive", action="store_true", help="缺少参数时直接报错")
    return parser.parse_args()


def strip_table_suffix(path: str) -> str:
    for suffix in (".tsv.gz", ".csv.gz", ".tsv", ".csv", ".gz"):
        if path.endswith(suffix):
            return path[: -len(suffix)]
    return path


def default_outputs(input_path: str, replace_input: bool) -> Tuple[str, str, str, str]:
    base = strip_table_suffix(input_path)
    output = input_path if replace_input else f"{base}.eaf_annotated.tsv.gz"
    return (
        output,
        f"{base}.eaf_matched.tsv.gz",
        f"{base}.eaf_frequency_diff.tsv.gz",
        f"{base}.eaf_annotation.log",
    )


def list_files(directory: str, suffixes: Tuple[str, ...]) -> List[str]:
    if not directory or not os.path.isdir(directory):
        return []
    return sorted(
        str(path)
        for path in Path(directory).iterdir()
        if path.is_file() and path.name.endswith(suffixes)
    )


def ask_text(prompt: str, default: str = "") -> str:
    answer = questionary.text(prompt, default=default, style=CUSTOM_STYLE).ask()
    if answer is None:
        sys.exit(0)
    return answer.strip()


def ask_select(prompt: str, choices: List[object]):
    answer = questionary.select(prompt, choices=choices, style=CUSTOM_STYLE).ask()
    if answer is None:
        sys.exit(0)
    return answer


def ask_confirm(prompt: str, default: bool = True) -> bool:
    answer = questionary.confirm(prompt, default=default, style=CUSTOM_STYLE).ask()
    if answer is None:
        sys.exit(0)
    return bool(answer)


def interactive_fill(args: argparse.Namespace) -> argparse.Namespace:
    if not args.input:
        candidates = list_files(DEFAULT_STANDARDIZED_DIR, (".tsv.gz", ".tsv"))
        if candidates:
            selected = ask_select("请选择需要注释 EAF 的标准化文件:", candidates + ["手动输入路径"])
            args.input = ask_text("请输入标准化文件路径:") if selected == "手动输入路径" else selected
        else:
            args.input = ask_text("请输入标准化文件路径:", DEFAULT_STANDARDIZED_DIR)

    if not args.frq:
        candidates = list_files(DEFAULT_REFERENCE_DIR, (".frq.tsv.gz", ".frq.tsv", ".frq.gz", ".frq"))
        if candidates:
            selected = ask_select("请选择 1KG .frq 参考文件:", candidates + ["手动输入路径"])
            args.frq = ask_text("请输入 .frq 路径:") if selected == "手动输入路径" else selected
        else:
            args.frq = ask_text("请输入 .frq 路径:", DEFAULT_REFERENCE_DIR)

    if args.output is None:
        args.replace_input = ask_confirm("是否直接替换原标准化文件?", default=True)
        if not args.replace_input:
            args.output = ask_text("请输入完整输出文件路径:", default=default_outputs(args.input, False)[0])

    if args.matched_output is None:
        args.matched_output = ask_text("请输入 BIM_ID 匹配行输出路径:", default=default_outputs(args.input, args.replace_input)[1])
    if args.diff_output is None:
        args.diff_output = ask_text("请输入频率差异位点输出路径:", default=default_outputs(args.input, args.replace_input)[2])
    if args.log is None:
        args.log = ask_text("请输入日志输出路径:", default=default_outputs(args.input, args.replace_input)[3])

    return args


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8", errors="replace")


def normalize_frq_to_tsv(frq_path: str, tmp_dir: str) -> str:
    """把 PLINK 空白分隔 .frq 规范化为 SNP/A1/A2/MAF TSV。"""
    normalized = str(Path(tmp_dir) / "reference.frq.tsv")
    with open_text(frq_path) as src, open(normalized, "w", encoding="utf-8") as dst:
        header = None
        indices = {}
        for line in src:
            fields = line.strip().split()
            if not fields:
                continue
            header = fields
            indices = {name: i for i, name in enumerate(header)}
            break
        if header is None:
            raise ValueError(f"空的 .frq 文件: {frq_path}")
        missing = [name for name in ("SNP", "A1", "A2", "MAF") if name not in indices]
        if missing:
            raise ValueError(f".frq 缺少必要列: {', '.join(missing)}")
        dst.write("SNP\tA1\tA2\tMAF\n")
        for line in src:
            fields = line.strip().split()
            if not fields:
                continue
            try:
                dst.write(
                    f"{fields[indices['SNP']]}\t{fields[indices['A1']]}\t"
                    f"{fields[indices['A2']]}\t{fields[indices['MAF']]}\n"
                )
            except IndexError:
                continue
    return normalized


def resolve_reference(frq_path: str, tmp_dir: str) -> str:
    """优先使用 build 脚本生成的 .frq.tsv.gz；否则临时规范化 .frq。"""
    if frq_path.endswith((".tsv", ".tsv.gz", ".csv", ".csv.gz")):
        return frq_path

    sibling = f"{frq_path}.tsv.gz"
    if os.path.exists(sibling):
        return sibling

    return normalize_frq_to_tsv(frq_path, tmp_dir)


def scan_reference(reference_path: str) -> pl.LazyFrame:
    schema = pl.scan_csv(
        reference_path,
        separator="\t",
        null_values=["NA", ""],
        infer_schema_length=1000,
    ).collect_schema().names()
    missing = [name for name in ("SNP", "A1", "A2", "MAF") if name not in schema]
    if missing:
        raise ValueError(f"参考频率文件缺少必要列: {', '.join(missing)}")

    return (
        pl.scan_csv(
            reference_path,
            separator="\t",
            null_values=["NA", ""],
            infer_schema_length=1000,
        )
        .select([
            pl.col("SNP").cast(pl.Utf8).str.to_uppercase().alias("BIM_ID"),
            pl.col("A1").cast(pl.Utf8).str.to_uppercase().alias("REF_A1"),
            pl.col("A2").cast(pl.Utf8).str.to_uppercase().alias("REF_A2"),
            pl.col("MAF").cast(pl.Float64, strict=False).alias("REF_MAF"),
        ])
        .unique(subset=["BIM_ID"], keep="first")
    )


def required_standardized_columns(schema: List[str]) -> None:
    required = {"BIM_ID", "EFFECT_ALLELE", "OTHER_ALLELE", "EAF", "MAF"}
    missing = sorted(required - set(schema))
    if missing:
        raise ValueError(f"标准化 GWAS 缺少必要列: {', '.join(missing)}")


def build_annotated_lf(input_path: str, reference_path: str, diff_threshold: float) -> Tuple[pl.LazyFrame, List[str]]:
    schema = pl.scan_csv(
        input_path,
        separator="\t",
        null_values=["NA", ""],
        infer_schema_length=10000,
    ).collect_schema().names()
    required_standardized_columns(schema)

    gwas = pl.scan_csv(
        input_path,
        separator="\t",
        null_values=["NA", ""],
        infer_schema_length=10000,
        ignore_errors=True,
    ).with_row_index("__ROW_NR")
    ref = scan_reference(reference_path)
    ref_keys = ref.select("BIM_ID")

    eaf_num = pl.col("EAF").cast(pl.Float64, strict=False)
    maf_num = pl.col("MAF").cast(pl.Float64, strict=False)
    ea = pl.col("EFFECT_ALLELE").cast(pl.Utf8).str.to_uppercase()
    oa = pl.col("OTHER_ALLELE").cast(pl.Utf8).str.to_uppercase()
    effect_is_ref_a1 = (ea == pl.col("REF_A1")) & (oa == pl.col("REF_A2"))
    effect_is_ref_a2 = (ea == pl.col("REF_A2")) & (oa == pl.col("REF_A1"))
    allele_match = effect_is_ref_a1 | effect_is_ref_a2

    matched = gwas.join(ref, on="BIM_ID", how="inner")
    unmatched = gwas.join(ref_keys, on="BIM_ID", how="anti")

    matched_annotated = matched.with_columns([
        eaf_num.alias("__EAF_NUM"),
        maf_num.alias("__MAF_NUM"),
        eaf_num.is_null().alias("__NEEDS_EAF"),
        pl.lit(True).alias("__BIM_ID_MATCHED"),
        pl.lit(True).alias("__HAS_REFERENCE"),
        effect_is_ref_a1.alias("__EFFECT_IS_REF_A1"),
        effect_is_ref_a2.alias("__EFFECT_IS_REF_A2"),
        allele_match.alias("__ALLELE_MATCH"),
    ]).with_columns([
        pl.when(pl.col("__EFFECT_IS_REF_A1"))
        .then(pl.col("REF_MAF"))
        .when(pl.col("__EFFECT_IS_REF_A2"))
        .then(1.0 - pl.col("REF_MAF"))
        .otherwise(None)
        .alias("__EAF_1KG"),
        pl.when(pl.col("__EFFECT_IS_REF_A1"))
        .then(pl.col("__MAF_NUM"))
        .when(pl.col("__EFFECT_IS_REF_A2"))
        .then(1.0 - pl.col("__MAF_NUM"))
        .otherwise(None)
        .alias("__EAF_FROM_MAF"),
    ]).with_columns([
        (
            pl.col("__NEEDS_EAF")
            & pl.col("__EAF_FROM_MAF").is_not_null()
            & pl.col("__ALLELE_MATCH")
        ).alias("__ANNOTATED_FROM_MAF"),
        (pl.col("__EAF_1KG") - pl.col("__EAF_FROM_MAF")).abs().alias("__EAF_ABS_DIFF"),
    ]).with_columns([
        pl.when(pl.col("__EAF_NUM").is_not_null())
        .then(pl.col("__EAF_NUM"))
        .when(pl.col("__ANNOTATED_FROM_MAF"))
        .then(pl.col("__EAF_FROM_MAF"))
        .otherwise(None)
        .alias("__EAF_FILLED"),
        (pl.col("__EAF_ABS_DIFF") > float(diff_threshold)).fill_null(False).alias("__EAF_DIFF_GT_THRESHOLD"),
        (pl.col("__EAF_ABS_DIFF") > 0.0).fill_null(False).alias("__EAF_DIFF_NONZERO"),
    ])

    unmatched_annotated = unmatched.with_columns([
        eaf_num.alias("__EAF_NUM"),
        maf_num.alias("__MAF_NUM"),
        eaf_num.is_null().alias("__NEEDS_EAF"),
        pl.lit(False).alias("__BIM_ID_MATCHED"),
        pl.lit(False).alias("__HAS_REFERENCE"),
        pl.lit(False).alias("__EFFECT_IS_REF_A1"),
        pl.lit(False).alias("__EFFECT_IS_REF_A2"),
        pl.lit(False).alias("__ALLELE_MATCH"),
        pl.lit(None, dtype=pl.Float64).alias("__EAF_1KG"),
        pl.lit(None, dtype=pl.Float64).alias("__EAF_FROM_MAF"),
        pl.lit(False).alias("__ANNOTATED_FROM_MAF"),
        pl.lit(None, dtype=pl.Float64).alias("__EAF_ABS_DIFF"),
        pl.lit(False).alias("__EAF_DIFF_GT_THRESHOLD"),
        pl.lit(False).alias("__EAF_DIFF_NONZERO"),
    ]).with_columns([
        pl.col("__EAF_NUM").alias("__EAF_FILLED"),
    ])

    annotated = pl.concat([matched_annotated, unmatched_annotated], how="diagonal_relaxed").sort("__ROW_NR")
    return annotated, schema


def final_columns(schema: List[str]) -> List[pl.Expr]:
    columns: List[pl.Expr] = []
    for name in schema:
        if name in LEGACY_AUDIT_COLUMNS:
            continue
        if name == "EAF":
            columns.append(pl.col("__EAF_FILLED").alias("EAF"))
        elif name == "MAF":
            columns.append(pl.col("__MAF_NUM").alias("MAF"))
        else:
            columns.append(pl.col(name))
    columns.extend([
        pl.col("__EAF_1KG").alias("EAF_1KG"),
        pl.col("__EAF_FROM_MAF").alias("EAF_FROM_MAF"),
        pl.col("__EAF_ABS_DIFF").alias("EAF_ABS_DIFF"),
    ])
    return columns


def stats_frame(lf: pl.LazyFrame) -> pl.DataFrame:
    return lf.select([
        pl.len().alias("total_rows"),
        pl.col("__BIM_ID_MATCHED").sum().alias("bimid_matched_rows"),
        (~pl.col("__BIM_ID_MATCHED")).sum().alias("bimid_unmatched_rows"),
        pl.col("__NEEDS_EAF").sum().alias("eaf_missing_before"),
        (pl.col("__NEEDS_EAF") & pl.col("__MAF_NUM").is_not_null()).sum().alias("eaf_missing_with_maf"),
        (pl.col("__NEEDS_EAF") & pl.col("__HAS_REFERENCE")).sum().alias("reference_snp_matched"),
        (pl.col("__NEEDS_EAF") & pl.col("__ALLELE_MATCH")).sum().alias("reference_allele_matched"),
        pl.col("__ANNOTATED_FROM_MAF").sum().alias("eaf_annotated_from_input_maf"),
        (pl.col("__NEEDS_EAF") & ~pl.col("__HAS_REFERENCE")).sum().alias("missing_reference"),
        (
            pl.col("__NEEDS_EAF")
            & pl.col("__HAS_REFERENCE")
            & ~pl.col("__ALLELE_MATCH")
        ).sum().alias("allele_mismatch"),
        (pl.col("__NEEDS_EAF") & pl.col("__MAF_NUM").is_null()).sum().alias("missing_input_maf"),
        (
            pl.col("__ANNOTATED_FROM_MAF")
            & pl.col("__EAF_DIFF_NONZERO")
        ).sum().alias("eaf_1kg_vs_input_maf_diff_nonzero"),
        (
            pl.col("__ANNOTATED_FROM_MAF")
            & pl.col("__EAF_DIFF_GT_THRESHOLD")
        ).sum().alias("eaf_1kg_vs_input_maf_diff_gt_threshold"),
    ]).collect()


def sink_tsv_atomic(lf: pl.LazyFrame, output: str) -> None:
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_plain = tempfile.mkstemp(prefix=f".{output_path.name}.", suffix=".tsv.tmp", dir=str(output_path.parent))
    os.close(fd)
    tmp_final = tmp_plain
    try:
        lf.sink_csv(
            tmp_plain,
            separator="\t",
            null_value="NA",
        )
        if output.endswith(".gz"):
            tmp_final = f"{tmp_plain}.gz"
            with open(tmp_plain, "rb") as src, gzip.open(tmp_final, "wb") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)
            with gzip.open(tmp_final, "rb") as handle:
                while handle.read(1024 * 1024):
                    pass
        os.replace(tmp_final, output)
    except Exception:
        if os.path.exists(tmp_final):
            os.remove(tmp_final)
        raise
    finally:
        if os.path.exists(tmp_plain):
            os.remove(tmp_plain)


def write_log(log_path: str, input_path: str, reference_path: str, output_path: str, matched_path: str, diff_path: str, diff_threshold: float, stats: dict) -> None:
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "EAF annotation report",
        "=" * 60,
        f"input: {input_path}",
        f"reference: {reference_path}",
        f"output: {output_path}",
        f"matched_output: {matched_path}",
        f"diff_output: {diff_path}",
        f"diff_threshold: {diff_threshold}",
        "-" * 60,
    ]
    lines.extend(f"{key}: {value}" for key, value in stats.items())
    with open(log_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    if not args.non_interactive:
        args = interactive_fill(args)

    if not args.input:
        console.print("[red]缺少 --input[/red]")
        sys.exit(1)
    if not args.frq:
        console.print("[red]缺少 --frq；请先运行 build_1kg_frq_reference.sh 或提供参考频率文件[/red]")
        sys.exit(1)
    if not os.path.exists(args.input):
        console.print(f"[red]输入文件不存在: {args.input}[/red]")
        sys.exit(1)
    if not os.path.exists(args.frq):
        console.print(f"[red]参考频率文件不存在: {args.frq}[/red]")
        sys.exit(1)

    default_output, default_matched, default_diff, default_log = default_outputs(args.input, args.replace_input)
    output = args.output or default_output
    matched_output = args.matched_output or default_matched
    diff_output = args.diff_output or default_diff
    log_path = args.log or default_log

    console.print(Panel.fit(
        "[bold cyan]EAF 注释[/bold cyan]\n"
        "[dim]用 1KG .frq 为原文件 MAF 定向，并保留 1KG 参考频率用于比较[/dim]",
        border_style="cyan",
    ))
    console.print(f"输入文件: [green]{args.input}[/green]")
    console.print(f"参考频率: [green]{args.frq}[/green]")
    console.print(f"完整输出: [green]{output}[/green]")
    console.print(f"BIM_ID 匹配输出: [green]{matched_output}[/green]")
    console.print(f"差异位点: [green]{diff_output}[/green]")
    console.print(f"日志: [green]{log_path}[/green]\n")

    with tempfile.TemporaryDirectory(prefix="eaf_frq_reference_") as tmp_dir:
        reference_path = resolve_reference(args.frq, tmp_dir)
        if reference_path != args.frq:
            console.print(f"使用规范化参考文件: [dim]{reference_path}[/dim]")

        annotated_lf, schema = build_annotated_lf(args.input, reference_path, args.diff_threshold)
        stats = stats_frame(annotated_lf).to_dicts()[0]

        full_lf = annotated_lf.select(final_columns(schema))
        matched_lf = annotated_lf.filter(pl.col("__BIM_ID_MATCHED")).select(final_columns(schema))
        diff_lf = (
            annotated_lf
            .filter(pl.col("__BIM_ID_MATCHED") & pl.col("__ALLELE_MATCH") & pl.col("__EAF_DIFF_GT_THRESHOLD"))
            .select(final_columns(schema))
        )

        sink_tsv_atomic(matched_lf, matched_output)
        sink_tsv_atomic(diff_lf, diff_output)
        sink_tsv_atomic(full_lf, output)
        write_log(
            log_path,
            input_path=args.input,
            reference_path=reference_path,
            output_path=output,
            matched_path=matched_output,
            diff_path=diff_output,
            diff_threshold=args.diff_threshold,
            stats=stats,
        )

    console.print("[green]EAF 注释完成[/green]")
    console.print(f"  成功注释 EAF: {stats['eaf_annotated_from_input_maf']:,}")
    console.print(f"  1KG 与原文件 MAF 推导 EAF 差异 > {args.diff_threshold}: {stats['eaf_1kg_vs_input_maf_diff_gt_threshold']:,}")


if __name__ == "__main__":
    main()
