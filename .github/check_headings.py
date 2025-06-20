#!/usr/bin/env python
# /// script
# dependencies = ["typer","rich"]
# ///

import os
import re
import sys
from pathlib import Path
from typing import List, Tuple

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

app = typer.Typer()
console = Console(force_terminal=True if os.getenv("GITHUB_ACTIONS") else False)


def is_heading(line: str) -> bool:
    """
    Check if a line is a numeric markdown heading.
    """
    return bool(re.match(r"^#+\s+\d+(\.\d+)*\.?\s+", line))


def has_trailing_period(heading_text: str) -> bool:
    """Check if heading ends with a period."""
    return bool(re.match(r"^#+\s+\d+(\.\d+)*\.\s+", heading_text))


def extract_heading_number(heading_text: str) -> str:
    """Extract the number part from a heading."""
    match = re.match(r"^#+\s+(\d+(\.\d+)*)", heading_text)
    return match.group(1) if match else ""


def extract_heading_level(heading_text: str) -> int:
    """Extract the number part from a heading."""
    match = re.match(r"^(#+)\s+\d+(\.\d+)*", heading_text)
    return len(match.group(1)) if match else 0


def check_heading_numbering(content: str, fix: bool = False) -> Tuple[List[tuple], str]:
    """
    Check if headings follow proper numbering conventions and optionally fix issues.

    Args:
        content: The markdown content
        fix: Whether to fix issues

    Returns:
        A tuple of (errors, fixed_content)
    """
    errors = []
    lines = content.split("\n")
    fixed_lines = []
    last_depth_nums = {}
    last_level_depth = 0

    for i, line in enumerate(lines):
        if not is_heading(line):
            fixed_lines.append(line)
            continue

        number_part = extract_heading_number(line)
        heading_level = extract_heading_level(line)
        number_sections = number_part.split(".")
        level_depth = len(number_sections)

        # Check for trailing period
        if not has_trailing_period(line):
            errors.append(
                (
                    line,
                    "Missing trailing period",
                    "Heading number should end with a period",
                )
            )
            if fix:
                line = line.replace(number_part, f"{number_part}.")

        # Check if heading level matches numbering level
        # Top level headings (1.) should be h2 (##), so expected level is depth + 1
        expected_level = level_depth + 1
        if heading_level != expected_level:
            errors.append(
                (
                    line,
                    "Heading level mismatch",
                    f"Numbering suggests h{expected_level} but found {heading_level}",
                )
            )
            if fix:
                line = ("#" * expected_level) + line.lstrip("#")

        # Check sections start at 0 or 1
        if (
            level_depth not in last_depth_nums
            and int(number_sections[-1]) != 0
            and int(number_sections[-1]) != 1
        ):
            errors.append(
                (
                    line,
                    "Numbering must start at 0 or 1",
                    f"Expected 0 or 1, but got {number_sections[-1]}.",
                )
            )
            if fix:
                fixed_numbers = number_sections[:-1] + ["1"]
                line = line.replace(number_part, ".".join(fixed_numbers))

        # Going back down a level, reset
        if int(last_level_depth) > int(level_depth):
            for i in range(level_depth + 1, 5):
                last_depth_nums[i] = 0

        # Check sequential numbering
        if (
            level_depth in last_depth_nums
            and int(number_sections[-1]) != int(last_depth_nums[level_depth]) + 1
        ):
            errors.append(
                (
                    line,
                    "Non-sequential heading",
                    f"Expected {'.'.join(number_sections[:-1])}.{last_depth_nums[level_depth] + 1}. but got {'.'.join(number_sections[:-1])}.{number_sections[-1]}.",
                )
            )
            if fix:
                fixed_numbers = number_sections[:-1] + [
                    str(last_depth_nums[level_depth] + 1)
                ]
                line = line.replace(number_part, ".".join(fixed_numbers))

        last_level_depth = level_depth
        last_depth_nums[level_depth] = int(number_sections[-1])

        fixed_lines.append(line)

    fixed_content = "\n".join(fixed_lines)
    return errors, fixed_content


@app.command()
def check(
    markdown_files: List[Path],
    fix: bool = typer.Option(False, "--fix", help="Automatically fix detected issues"),
):
    """
    Check markdown files for proper heading numbering.

    Validates:
    - Sequential numbering at each level
    - Trailing period after numbers
    - Heading level matches number level (## for 1., ### for 1.1., etc.)

    With --fix flag, automatically corrects issues in the files.
    """
    # Filter out files in transcripts directories
    markdown_files = [f for f in markdown_files if "transcripts" not in f.parts]

    has_errors = False
    all_error_messages = []
    total_errors = 0
    fixed_files = 0

    for file_path in markdown_files:
        content = file_path.read_text(encoding="utf-8")
        errors, fixed_content = check_heading_numbering(content, fix)

        if errors:
            has_errors = True
            total_errors += len(errors)

            # Create table for errors
            table = Table(show_header=False, box=None)
            table.add_column("Heading", style="green")
            table.add_column("Reason", style="bold")
            table.add_column("Description")

            # Add each error directly to the table
            for heading_text, error_type, error_details in errors:
                table.add_row(f"'{heading_text}'", error_type, error_details)

            panel = Panel(
                table,
                title=f"[bold red]{file_path}[/bold red]",
                border_style="red",
                title_align="left",
            )
            all_error_messages.append(panel)

            # Write fixed content if needed
            if fix and fixed_content != content:
                file_path.write_text(fixed_content, encoding="utf-8")
                fixed_files += 1

    if has_errors:
        for panel in all_error_messages:
            console.print(panel)

        if fix:
            msg = f":wrench: [green bold]Checked {len(markdown_files)} files, found {total_errors} errors, "
            msg += f"fixed {fixed_files} files"
            console.print(msg)
        else:
            console.print(
                f":x: [red]Checked {len(markdown_files)} files and found {total_errors} errors"
            )

        if not fix:
            sys.exit(1)
        elif fixed_files < len(all_error_messages):
            # Some files couldn't be fully fixed
            sys.exit(1)
    else:
        console.print(
            f":white_check_mark: [green]Checked {len(markdown_files)} files and found {total_errors} errors"
        )


if __name__ == "__main__":
    app()
