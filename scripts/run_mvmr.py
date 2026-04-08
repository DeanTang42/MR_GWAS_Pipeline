#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""交互式构建并运行可配置 exposure 数量的 MVMR。"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

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


def ask_select(prompt: str, choices: Sequence[object]):
    try:
        answer = questionary.select(prompt, choices=list(choices), style=CUSTOM_STYLE).ask()
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


def ask_int(prompt: str, default: int, minimum: int = 1) -> int:
    while True:
        raw = ask_text(prompt, default=str(default))
        try:
            value = int(raw)
        except ValueError:
            console.print("[red]请输入整数。[/red]")
            continue
        if value < minimum:
            console.print(f"[red]请输入大于等于 {minimum} 的整数。[/red]")
            continue
        return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="交互式运行可配置 exposure 数量的 MVMR")
    parser.add_argument("--exp-count", type=int, help="纳入分析的 exposure 数量，至少为 2")
    parser.add_argument("--exp", dest="exp_names", action="append", default=[], help="exposure 文件名（不含后缀），可重复传入")
    parser.add_argument("--exp-std", dest="exp_stds", action="append", default=[], help="对应 exposure 的 standardized .gz 文件，可重复传入")
    parser.add_argument("--exp-label", dest="exp_labels", action="append", default=[], help="对应 exposure 的展示名称，可重复传入")
    parser.add_argument("--y-std", help="outcome 对应的 standardized .gz 文件")
    parser.add_argument("--y-label", help="outcome 的展示名称")
    parser.add_argument("--output-dir", help="结果输出目录")
    parser.add_argument("--r-lib-path", default=DEFAULT_R_LIB_PATH, help="R 包库路径")
    parser.add_argument("--non-interactive", action="store_true", help="缺参数时直接报错，且不再进行确认")
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
                items.add(name[:-len(suffix)])
                break
    return sorted(items)


def list_standardized_gz_files() -> List[str]:
    if not DEFAULT_STANDARDIZED_DIR or not os.path.isdir(DEFAULT_STANDARDIZED_DIR):
        return []
    return sorted(str(path) for path in Path(DEFAULT_STANDARDIZED_DIR).iterdir() if path.is_file() and path.name.endswith(".gz"))


def std_label(path: str) -> str:
    name = Path(path).name
    for suffix in (".tsv.gz", ".csv.gz", ".gz"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    for suffix in ("_standardized", ".eaf_matched", ".eaf_annotated"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name


def default_output_dir(exp_labels: Sequence[str], y_label: str) -> str:
    base = DEFAULT_RESULTS_DIR or str(Path.cwd() / "results")
    exp_part = "_".join(exp_labels)
    return str(Path(base) / f"{exp_part}_{y_label}_mvmr")


def choose_exp_dataset(prompt: str, exclude: Sequence[str]) -> str:
    candidates = [item for item in list_exp_datasets() if item not in exclude]
    if candidates:
        selected = ask_select(prompt, list(candidates) + [questionary.Separator(), "手动输入文件名"])
        return ask_text("请输入 exposure 文件名（不含后缀）:") if selected == "手动输入文件名" else selected
    return ask_text("请输入 exposure 文件名（不含后缀）:", default=DEFAULT_EXP_DIR)


def choose_standardized_file(prompt: str) -> str:
    candidates = list_standardized_gz_files()
    if candidates:
        selected = ask_select(prompt, candidates + [questionary.Separator(), "手动输入路径"])
        return ask_text("请输入 standardized .gz 文件路径:") if selected == "手动输入路径" else selected
    return ask_text("请输入 standardized .gz 文件路径:", default=DEFAULT_STANDARDIZED_DIR)


def interactive_fill(args: argparse.Namespace) -> argparse.Namespace:
    if not args.exp_count:
        args.exp_count = ask_int("请输入纳入分析的 exposure 数量:", default=2, minimum=2)

    while len(args.exp_names) < args.exp_count:
        idx = len(args.exp_names) + 1
        args.exp_names.append(choose_exp_dataset(f"请选择第 {idx} 个 exposure 文件:", exclude=args.exp_names))

    while len(args.exp_stds) < args.exp_count:
        idx = len(args.exp_stds) + 1
        args.exp_stds.append(choose_standardized_file(f"请选择第 {idx} 个 exposure 对应的 standardized .gz 文件:"))

    while len(args.exp_labels) < args.exp_count:
        idx = len(args.exp_labels)
        default_label = args.exp_names[idx]
        args.exp_labels.append(ask_text(f"请输入第 {idx + 1} 个 exposure 的展示名称:", default=default_label))

    if not args.y_std:
        args.y_std = choose_standardized_file("请选择 outcome 对应的 standardized .gz 文件:")
    args.y_label = args.y_label or ask_text("请输入 outcome 的展示名称:", default=std_label(args.y_std))

    if not args.output_dir:
        args.output_dir = ask_text(
            "请输入结果输出目录:",
            default=default_output_dir(args.exp_labels, args.y_label),
        )
    return args


def validate_args(args: argparse.Namespace) -> None:
    if args.exp_count is None:
        raise ValueError("缺少必要参数: --exp-count")
    if args.exp_count < 2:
        raise ValueError("--exp-count 至少为 2")
    if len(args.exp_names) != args.exp_count:
        raise ValueError(f"--exp 需要提供 {args.exp_count} 次，当前为 {len(args.exp_names)} 次")
    if len(args.exp_stds) != args.exp_count:
        raise ValueError(f"--exp-std 需要提供 {args.exp_count} 次，当前为 {len(args.exp_stds)} 次")
    if args.exp_labels and len(args.exp_labels) != args.exp_count:
        raise ValueError(f"--exp-label 需要提供 {args.exp_count} 次，当前为 {len(args.exp_labels)} 次")
    if len(set(args.exp_names)) != len(args.exp_names):
        raise ValueError("exposure 文件名存在重复，请为每个 exposure 选择不同文件")
    if not args.y_std:
        raise ValueError("缺少必要参数: --y-std")
    if not args.y_label:
        raise ValueError("缺少必要参数: --y-label")
    if not args.output_dir:
        raise ValueError("缺少必要参数: --output-dir")
    for path in list(args.exp_stds) + [args.y_std]:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"文件不存在: {path}")


def build_command(args: argparse.Namespace) -> List[str]:
    project_dir = Path(__file__).resolve().parent.parent
    output_dir = os.path.abspath(os.path.expanduser(args.output_dir))
    cmd = [
        "Rscript",
        str(project_dir / "scripts" / "MVMR_pipeline.R"),
        "--project-dir", str(project_dir),
        "--r-lib-path", args.r_lib_path,
        "--y-std", os.path.abspath(args.y_std),
        "--y-label", args.y_label,
        "--output-dir", output_dir,
    ]
    for item in args.exp_names:
        cmd.extend(["--exp-name", item])
    for item in args.exp_stds:
        cmd.extend(["--exp-std", os.path.abspath(item)])
    for item in args.exp_labels:
        cmd.extend(["--exp-label", item])
    return cmd


def print_summary(args: argparse.Namespace) -> None:
    lines = [
        "通用 MVMR 配置",
        "",
        f"Exposure 数量: {args.exp_count}",
    ]
    for idx, (name, std_path, label) in enumerate(zip(args.exp_names, args.exp_stds, args.exp_labels), start=1):
        lines.extend([
            f"EXP {idx}: {name}",
            f"  standardized: {std_path}",
            f"  label: {label}",
        ])
    lines.extend([
        f"Outcome standardized: {args.y_std}",
        f"Outcome label: {args.y_label}",
        f"R_LIB_PATH: {args.r_lib_path}",
        f"输出目录: {os.path.abspath(os.path.expanduser(args.output_dir))}",
    ])
    console.print(Panel.fit("\n".join(lines), border_style="cyan"))


def main() -> None:
    args = parse_args()
    cli_complete = (
        args.exp_count is not None and
        len(args.exp_names) == args.exp_count and
        len(args.exp_stds) == args.exp_count and
        (len(args.exp_labels) in (0, args.exp_count)) and
        bool(args.y_std)
    )

    if not cli_complete:
        if args.non_interactive:
            console.print("[red]--non-interactive 模式下必须完整提供 exposure 数量、所有 exposure 文件、所有 standardized 文件和 outcome standardized 文件[/red]")
            sys.exit(1)
        console.print(Panel.fit("可配置 exposure 数量的 MVMR 分析", border_style="cyan"))
        args = interactive_fill(args)
    else:
        if not args.exp_labels:
            args.exp_labels = list(args.exp_names)
        args.y_label = args.y_label or std_label(args.y_std)
        args.output_dir = args.output_dir or default_output_dir(args.exp_labels, args.y_label)

    try:
        validate_args(args)
    except Exception as exc:
        console.print(f"[red]{exc}[/red]")
        sys.exit(1)

    print_summary(args)
    if not args.non_interactive and not ask_confirm("以上配置是否正确并开始分析?", default=True):
        _cancel_interaction()

    cmd = build_command(args)
    completed = subprocess.run(cmd)
    sys.exit(completed.returncode)


if __name__ == "__main__":
    main()
