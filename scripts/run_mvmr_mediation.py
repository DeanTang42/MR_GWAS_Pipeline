#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""交互式构建并运行 X-M-Y 中介 MVMR。"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple

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
DEFAULT_EXP_DIR = os.environ.get("MR_PIPELINE_EXP_DIR", os.environ.get("MR_PIPELINE_EXPOSURE_DIR", ""))
DEFAULT_STANDARDIZED_DIR = os.environ.get("MR_PIPELINE_STANDARDIZED_OUTPUT_DIR", "")
DEFAULT_RESULTS_DIR = os.environ.get("MR_PIPELINE_RESULTS_DIR", "")


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="交互式运行 X-M-Y 中介 MVMR 分析")
    parser.add_argument("--x-exp", help="X 对应的 exposure 文件名（不含后缀），固定从 data/exp 选择")
    parser.add_argument("--m-exp", help="M 对应的 exposure 文件名（不含后缀），固定从 data/exp 选择")
    parser.add_argument("--x-std", help="X 对应的 standardized .gz 文件")
    parser.add_argument("--m-std", help="M 对应的 standardized .gz 文件")
    parser.add_argument("--y-std", help="Y 对应的 standardized .gz 文件")
    parser.add_argument("--x-label", help="X 的展示名称，默认使用 exposure 文件名")
    parser.add_argument("--m-label", help="M 的展示名称，默认使用 exposure 文件名")
    parser.add_argument("--y-label", help="Y 的展示名称，默认使用 standardized 文件名")
    parser.add_argument("--output-dir", help="结果输出目录")
    parser.add_argument("--r-lib-path", default=DEFAULT_R_LIB_PATH, help="R 包库路径")
    parser.add_argument("--non-interactive", action="store_true", help="缺参数时直接报错")
    return parser.parse_args()


def list_exp_datasets() -> List[str]:
    if not DEFAULT_EXP_DIR or not os.path.isdir(DEFAULT_EXP_DIR):
        return []
    items = set()
    for path in Path(DEFAULT_EXP_DIR).iterdir():
        if not path.is_file():
            continue
        name = path.name
        for suffix in (".csv.gz", ".csv"):
            if name.endswith(suffix):
                items.add(name[: -len(suffix)])
                break
    return sorted(items)


def list_standardized_gz_files() -> List[str]:
    if not DEFAULT_STANDARDIZED_DIR or not os.path.isdir(DEFAULT_STANDARDIZED_DIR):
        return []
    files = []
    for path in Path(DEFAULT_STANDARDIZED_DIR).iterdir():
        if path.is_file() and path.name.endswith(".gz"):
            files.append(str(path))
    return sorted(files)


def std_label(path: str) -> str:
    name = Path(path).name
    for suffix in (".tsv.gz", ".csv.gz", ".gz"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    for suffix in ("_standardized", ".eaf_matched", ".eaf_annotated"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def default_output_dir(x_label: str, m_label: str, y_label: str) -> str:
    base = DEFAULT_RESULTS_DIR or str(Path.cwd() / "results")
    return str(Path(base) / f"{x_label}_{m_label}_{y_label}_mvmr_mediation")


def choose_exp_dataset(prompt: str, exclude: Optional[str] = None) -> str:
    candidates = list_exp_datasets()
    if exclude:
        candidates = [item for item in candidates if item != exclude]
    if candidates:
        selected = ask_select(prompt, candidates + [questionary.Separator(), "手动输入文件名"])
        return ask_text("请输入 exposure 文件名（不含后缀）:") if selected == "手动输入文件名" else selected
    return ask_text("请输入 exposure 文件名（不含后缀）:", default=DEFAULT_EXP_DIR)


def choose_standardized_file(prompt: str) -> str:
    candidates = list_standardized_gz_files()
    if candidates:
        selected = ask_select(prompt, candidates + [questionary.Separator(), "手动输入路径"])
        return ask_text("请输入 standardized .gz 文件路径:") if selected == "手动输入路径" else selected
    return ask_text("请输入 standardized .gz 文件路径:", default=DEFAULT_STANDARDIZED_DIR)


def interactive_fill(args: argparse.Namespace) -> argparse.Namespace:
    if not args.x_exp:
        args.x_exp = choose_exp_dataset("请选择 X 对应的 exposure 文件:")
    if not args.m_exp:
        args.m_exp = choose_exp_dataset("请选择 M 对应的 exposure 文件:", exclude=args.x_exp)
    if not args.x_std:
        args.x_std = choose_standardized_file("请选择 X 对应的 standardized .gz 文件:")
    if not args.m_std:
        args.m_std = choose_standardized_file("请选择 M 对应的 standardized .gz 文件:")
    if not args.y_std:
        args.y_std = choose_standardized_file("请选择 Y 对应的 standardized .gz 文件:")

    args.x_label = args.x_label or args.x_exp
    args.m_label = args.m_label or args.m_exp
    args.y_label = args.y_label or std_label(args.y_std)

    if not args.output_dir:
        args.output_dir = ask_text(
            "请输入结果输出目录:",
            default=default_output_dir(args.x_label, args.m_label, args.y_label),
        )
    return args


def validate_args(args: argparse.Namespace) -> None:
    required = {
        "--x-exp": args.x_exp,
        "--m-exp": args.m_exp,
        "--x-std": args.x_std,
        "--m-std": args.m_std,
        "--y-std": args.y_std,
    }
    missing = [name for name, value in required.items() if not value]
    if missing:
        raise ValueError(f"缺少必要参数: {', '.join(missing)}")
    for path in (args.x_std, args.m_std, args.y_std):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"文件不存在: {path}")


def build_command(args: argparse.Namespace) -> List[str]:
    project_dir = Path(__file__).resolve().parent.parent
    output_dir = os.path.abspath(os.path.expanduser(args.output_dir))
    return [
        "Rscript",
        str(project_dir / "scripts" / "MVMR_mediation_pipeline.R"),
        "--project-dir", str(project_dir),
        "--r-lib-path", args.r_lib_path,
        "--x-exp-name", args.x_exp,
        "--m-exp-name", args.m_exp,
        "--x-std", os.path.abspath(args.x_std),
        "--m-std", os.path.abspath(args.m_std),
        "--y-std", os.path.abspath(args.y_std),
        "--x-label", args.x_label,
        "--m-label", args.m_label,
        "--y-label", args.y_label,
        "--output-dir", output_dir,
    ]


def print_summary(args: argparse.Namespace) -> None:
    lines = [
        "MVMR 中介分析配置",
        "",
        f"X exposure: {args.x_exp}",
        f"M exposure: {args.m_exp}",
        f"X standardized: {args.x_std}",
        f"M standardized: {args.m_std}",
        f"Y standardized: {args.y_std}",
        f"X label: {args.x_label}",
        f"M label: {args.m_label}",
        f"Y label: {args.y_label}",
        f"R_LIB_PATH: {args.r_lib_path}",
        f"输出目录: {os.path.abspath(os.path.expanduser(args.output_dir))}",
    ]
    console.print(Panel.fit("\n".join(lines), border_style="cyan"))


def main() -> None:
    args = parse_args()
    cli_mode = all([args.x_exp, args.m_exp, args.x_std, args.m_std, args.y_std])

    if not cli_mode:
        if args.non_interactive:
            console.print("[red]--non-interactive 模式必须提供 --x-exp/--m-exp/--x-std/--m-std/--y-std[/red]")
            sys.exit(1)
        console.print(Panel.fit("X-M-Y 中介 MVMR 分析", border_style="cyan"))
        args = interactive_fill(args)
    else:
        args.x_label = args.x_label or args.x_exp
        args.m_label = args.m_label or args.m_exp
        args.y_label = args.y_label or std_label(args.y_std)
        args.output_dir = args.output_dir or default_output_dir(args.x_label, args.m_label, args.y_label)

    try:
        validate_args(args)
    except Exception as exc:
        console.print(f"[red]{exc}[/red]")
        sys.exit(1)

    print_summary(args)
    if not ask_confirm("以上配置是否正确并开始分析?", default=True):
        _cancel_interaction()

    cmd = build_command(args)
    completed = subprocess.run(cmd)
    sys.exit(completed.returncode)


if __name__ == "__main__":
    main()
