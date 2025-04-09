#!/usr/bin/env python

import re
import sys
from pathlib import Path
from typing import List

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

app = typer.Typer()
console = Console()

def extract_headings(content: str) -> List[tuple]:
    """Extract all markdown headings with their level and text."""
    heading_pattern = r'^(#{1,6})\s+(.*?)$'
    return [(len(match.group(1)), match.group(2).strip()) for match in re.finditer(heading_pattern, content, re.MULTILINE)]

def is_numbered_heading(heading_text: str) -> bool:
    """Check if heading starts with a number pattern like 1., 1.1., etc."""
    return bool(re.match(r'^\d+(\.\d+)*', heading_text))

def has_trailing_period(heading_text: str) -> bool:
    """Check if heading ends with a period."""
    return bool(re.match(r'^\d+(\.\d+)*\. ', heading_text))

def extract_heading_number(heading_text: str) -> str:
    """Extract the number part from a heading."""
    match = re.match(r'^(\d+(\.\d+)*)', heading_text)
    return match.group(1) if match else ""

def check_heading_numbering(headings: List[tuple]) -> List[tuple]:
    """
    Check if headings follow proper numbering conventions.
    Returns a list of tuples (heading_text, error_type, error_details).
    """
    errors = []
    current_numbers = {}  # Track the current number at each level

    for heading_level, heading_text in headings:
        if not is_numbered_heading(heading_text):
            # Check if it's a numbered heading but missing the trailing period
            if re.match(r'^\d+(\.\d+)*\s', heading_text):
                errors.append((
                    heading_text,
                    "Missing trailing period",
                    "Number part should end with a period"
                ))
            continue

        if not has_trailing_period(heading_text):
            errors.append((
                heading_text,
                "Missing trailing period",
                "Heading number should end with a period"
            ))

        number_part = extract_heading_number(heading_text)
        number_sections = number_part.split('.')
        level_depth = len(number_sections)

        # Check if heading level matches numbering level
        # Top level headings (1.) should be h2 (##), so expected level is depth + 1
        expected_level = level_depth + 1
        if heading_level != expected_level:
            errors.append((
                heading_text,
                "Heading level mismatch",
                f"Numbering suggests h{expected_level} but found {heading_level}"
            ))

        # Check sequential numbering
        parent_key = '.'.join(number_sections[:-1])
        current_level = int(number_sections[-1])

        # For sub-levels, we don't need to verify parent level exists
        # as it might come later in the document
        if level_depth > 1:
            # Initialize parent level if not seen yet
            if parent_key not in current_numbers:
                current_numbers[parent_key] = 0

        # Check if numbering is sequential
        expected = current_numbers.get(parent_key, 0) + 1
        if current_level != expected:
            # Special case: allow starting with 0
            if level_depth == 1 and current_level == 0 and expected == 1:
                # This is the first heading and it starts with 0, which is valid
                pass
            else:
                if level_depth == 1:
                    errors.append((
                        heading_text,
                        "Non-sequential top-level heading",
                        f"Expected {expected}. but got {current_level}."
                    ))
                else:
                    errors.append((
                        heading_text,
                        "Non-sequential heading",
                        f"Expected {parent_key}.{expected}. but got {number_part}"
                    ))

        # Update the current number for this level
        current_numbers[parent_key] = current_level

        # Reset all deeper levels when encountering a new parent
        deeper_levels = [
            k for k in current_numbers
            if k.startswith(f"{parent_key}.") or (parent_key == "" and '.' in k)
        ]
        for k in deeper_levels:
            del current_numbers[k]

    return errors

@app.command()
def check(markdown_files: List[Path]):
    """
    Check markdown files for proper heading numbering.

    Validates:
    - Sequential numbering at each level
    - Trailing period after numbers
    - Heading level matches number level (## for 1., ### for 1.1., etc.)
    """
    has_errors = False
    all_error_messages = []
    total_errors = 0

    for file_path in markdown_files:
        content = file_path.read_text(encoding="utf-8")
        headings = extract_headings(content)
        errors = check_heading_numbering(headings)

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
                title_align="left"
            )
            all_error_messages.append(panel)

    if has_errors:
        for panel in all_error_messages:
            console.print(panel)
        console.print(f":x: Checked {len(markdown_files)} files and found {total_errors} errors")
        sys.exit(1)
    else:
        console.print(f":white_check_mark: Checked {len(markdown_files)} files and found {total_errors} errors")

if __name__ == "__main__":
    app()
