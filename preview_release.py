#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rich",
#     "rich-click",
# ]
# ///
"""
Preview Release Script

Serves the training docs locally at https://training.nextflow.io/ with the
current branch appearing as a specified version release.

Note: This script should NOT be run with sudo. It will prompt for sudo
privileges only when needed (for /etc/hosts modification).

See CONTRIBUTING.md for full documentation.
"""

import atexit
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

try:
    import rich_click as click
    from rich.console import Console
    from rich.panel import Panel
    from rich.progress import Progress, SpinnerColumn, TextColumn
    from rich.table import Table
except ImportError:
    print("Error: Required packages not found.")
    print("Run this script with uv: uv run ./preview_release.py --version 3.0")
    sys.exit(1)

# Configuration
DOMAIN = "training.nextflow.io"
WORK_DIR = Path(".preview-release")
CERTS_DIR = WORK_DIR / "certs"
SITE_DIR = WORK_DIR / "site"
CADDYFILE = WORK_DIR / "Caddyfile"

# Rich console
console = Console()

# Global state for cleanup
_cleanup_done = False
_hosts_modified = False
_caddy_process = None
_keep_files = True

# Configure rich-click
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True
click.rich_click.STYLE_ERRORS_SUGGESTION = "dim"


def run_cmd(cmd, check=True, capture=False, **kwargs):
    """Run a shell command."""
    if capture:
        result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
        if check and result.returncode != 0:
            console.print(f"[red]Error running {' '.join(cmd)}[/red]")
            console.print(result.stderr)
            sys.exit(1)
        return result
    else:
        return subprocess.run(cmd, check=check, **kwargs)


def status_msg(phase: str, message: str, success: bool | None = None):
    """Print a status message."""
    phase_style = "cyan"
    if success is True:
        icon = "[green]✓[/green]"
    elif success is False:
        icon = "[red]✗[/red]"
    else:
        icon = "[dim]•[/dim]"
    console.print(f"[{phase_style}][{phase}][/{phase_style}] {icon} {message}")


def check_sudo():
    """Check if running with sudo privileges."""
    return os.geteuid() == 0


def run_with_sudo(cmd, check=True, capture=False, **kwargs):
    """Run a command with sudo, prompting for password if needed."""
    return run_cmd(["sudo"] + cmd, check=check, capture=capture, **kwargs)


def cleanup():
    """Clean up everything on exit."""
    global _cleanup_done, _hosts_modified, _caddy_process, _keep_files

    if _cleanup_done:
        return
    _cleanup_done = True

    console.print()
    status_msg("cleanup", "Shutting down...")

    # Stop Caddy
    if _caddy_process and _caddy_process.poll() is None:
        status_msg("cleanup", "Stopping server...")
        _caddy_process.terminate()
        try:
            _caddy_process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            _caddy_process.kill()
        status_msg("cleanup", "Server stopped", success=True)

    # Remove hosts entry
    if _hosts_modified:
        status_msg("cleanup", "Removing hosts entry...")
        run_with_sudo(["sed", "-i", "", f"/{DOMAIN}/d", "/etc/hosts"], check=False)
        status_msg("cleanup", "Hosts entry removed", success=True)

    # Remove work directory
    if not _keep_files and WORK_DIR.exists():
        status_msg("cleanup", "Removing work directory...")
        shutil.rmtree(WORK_DIR)
        status_msg("cleanup", "Work directory removed", success=True)
    elif _keep_files:
        status_msg("cleanup", "Keeping work directory for faster restart", success=True)

    console.print()
    console.print("[green bold]Done![/green bold]")


def check_dependencies():
    """Check and install required dependencies."""
    status_msg("setup", "Checking dependencies...")

    # Check for Homebrew
    if not shutil.which("brew"):
        status_msg(
            "setup", "Homebrew not found. Please install it first.", success=False
        )
        sys.exit(1)

    # Check for Docker
    result = run_cmd(["docker", "info"], check=False, capture=True)
    if result.returncode != 0:
        status_msg(
            "setup",
            "Docker is not running. Please start Docker Desktop.",
            success=False,
        )
        sys.exit(1)
    status_msg("setup", "Docker running", success=True)

    # Check/install mkcert
    if not shutil.which("mkcert"):
        status_msg("setup", "Installing mkcert...")
        run_cmd(["brew", "install", "mkcert"])
    status_msg("setup", "mkcert installed", success=True)

    # Check/install caddy
    if not shutil.which("caddy"):
        status_msg("setup", "Installing caddy...")
        run_cmd(["brew", "install", "caddy"])
    status_msg("setup", "caddy installed", success=True)

    # Check if mkcert CA is installed in the system keychain
    result = run_cmd(
        [
            "security",
            "find-certificate",
            "-c",
            "mkcert",
            "/Library/Keychains/System.keychain",
        ],
        check=False,
        capture=True,
    )
    ca_in_keychain = result.returncode == 0 and "mkcert" in result.stdout

    if not ca_in_keychain:
        status_msg("setup", "mkcert CA not found in system keychain.", success=False)
        console.print()
        console.print(
            Panel(
                "[yellow]The mkcert root CA must be installed for browsers to trust the certificate.[/yellow]\n\n"
                "This is a one-time setup. Run this command as your regular user (not sudo):\n\n"
                "    [bold cyan]mkcert -install[/bold cyan]\n\n"
                "Then restart your browser and re-run this script.",
                title="Action Required",
                border_style="yellow",
            )
        )
        console.print()
        sys.exit(1)

    status_msg("setup", "mkcert CA trusted by system", success=True)


def generate_certificates():
    """Generate TLS certificates for the domain."""
    status_msg("setup", "Generating certificates...")

    CERTS_DIR.mkdir(parents=True, exist_ok=True)

    cert_file = (CERTS_DIR / f"{DOMAIN}.pem").resolve()
    key_file = (CERTS_DIR / f"{DOMAIN}-key.pem").resolve()

    if cert_file.exists() and key_file.exists():
        status_msg("setup", "Certificates already exist", success=True)
        return

    run_cmd(
        ["mkcert", "-cert-file", str(cert_file), "-key-file", str(key_file), DOMAIN]
    )
    status_msg("setup", "Certificates created", success=True)


def fetch_gh_pages():
    """Fetch gh-pages content from upstream."""
    # Skip if already fetched
    versions_file = SITE_DIR / "versions.json"
    if versions_file.exists():
        versions = json.loads(versions_file.read_text())
        status_msg(
            "setup", f"Using cached gh-pages ({len(versions)} versions)", success=True
        )
        return

    SITE_DIR.mkdir(parents=True, exist_ok=True)

    # Check if upstream remote exists
    result = run_cmd(
        ["git", "remote", "get-url", "upstream"], check=False, capture=True
    )
    if result.returncode != 0:
        status_msg("setup", "upstream remote not found.", success=False)
        console.print(
            "  Please add it: [cyan]git remote add upstream https://github.com/nextflow-io/training.git[/cyan]"
        )
        sys.exit(1)

    with Progress(
        TextColumn(""),
        SpinnerColumn(),
        TextColumn("[cyan][setup][/cyan] {task.description}"),
        console=console,
        transient=True,
    ) as progress:
        progress.add_task("Fetching gh-pages content...", total=None)

        # Fetch latest from upstream
        run_cmd(["git", "fetch", "upstream", "gh-pages"], capture=True)

        # Extract gh-pages content using git archive
        archive_proc = subprocess.Popen(
            ["git", "archive", "upstream/gh-pages"], stdout=subprocess.PIPE
        )
        subprocess.run(
            ["tar", "-x", "-C", str(SITE_DIR)], stdin=archive_proc.stdout, check=True
        )
        archive_proc.wait()

    # Count versions fetched
    if versions_file.exists():
        versions = json.loads(versions_file.read_text())
        status_msg("setup", f"Fetched {len(versions)} existing versions", success=True)
    else:
        status_msg(
            "setup", "Warning: versions.json not found in gh-pages", success=False
        )


def compute_source_hash():
    """Compute a hash of source files to detect changes."""
    hasher = hashlib.md5()

    # Hash key source directories and files
    source_paths = [
        Path("docs"),
        Path("mkdocs.yml"),
    ]

    for source_path in source_paths:
        if not source_path.exists():
            continue
        if source_path.is_file():
            hasher.update(source_path.read_bytes())
        else:
            # For directories, hash all markdown and config files
            for file in sorted(source_path.rglob("*")):
                if file.is_file() and file.suffix in (
                    ".md",
                    ".yml",
                    ".yaml",
                    ".css",
                    ".js",
                    ".html",
                ):
                    try:
                        hasher.update(str(file).encode())
                        hasher.update(file.read_bytes())
                    except (IOError, OSError):
                        pass

    return hasher.hexdigest()[:12]


def build_docs(version: str):
    """Build docs for the specified version using Docker."""
    version_dir = SITE_DIR / version
    hash_file = WORK_DIR / f"{version}.hash"

    # Compute current source hash
    current_hash = compute_source_hash()

    # Check if we can skip the build
    if (
        version_dir.exists()
        and (version_dir / "index.html").exists()
        and hash_file.exists()
    ):
        stored_hash = hash_file.read_text().strip()
        if stored_hash == current_hash:
            status_msg(
                "setup",
                f"v{version} already built [dim](source unchanged)[/dim]",
                success=True,
            )
            return
        else:
            status_msg("setup", "Source files changed, rebuilding...")
            shutil.rmtree(version_dir)

    # Ensure clean build directory (use rm -rf to handle root-owned files from Docker)
    if version_dir.exists():
        run_cmd(["rm", "-rf", str(version_dir)], check=False)
    version_dir.mkdir(parents=True, exist_ok=True)
    # Make directory world-writable so Docker container can write to it
    run_cmd(["chmod", "-R", "777", str(version_dir)], check=False)

    # Get absolute paths
    repo_root = Path.cwd().resolve()

    with Progress(
        TextColumn(""),
        SpinnerColumn(),
        TextColumn("[cyan][setup][/cyan] {task.description}"),
        console=console,
        transient=True,
    ) as progress:
        progress.add_task(
            f"Building v{version} docs (this may take a few minutes)...", total=None
        )

        result = run_cmd(
            [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{repo_root}:/docs",
                "-e",
                "CARDS=false",
                "ghcr.io/nextflow-io/training-mkdocs:latest",
                "build",
                "-d",
                f"/docs/{WORK_DIR}/site/{version}",
            ],
            check=False,
            capture=True,
        )

    if result.returncode != 0:
        status_msg("setup", "Docker build failed", success=False)
        console.print(result.stderr)
        console.print(result.stdout)
        sys.exit(1)

    # Save hash for future cache checks
    hash_file.write_text(current_hash)

    status_msg("setup", f"v{version} built successfully", success=True)


def update_versions_json(version: str):
    """Update versions.json to make the specified version the latest."""
    status_msg("setup", "Updating versions.json...")

    versions_file = SITE_DIR / "versions.json"

    if not versions_file.exists():
        versions = []
    else:
        versions = json.loads(versions_file.read_text())

    # Remove 'latest' alias from all existing versions
    for v in versions:
        if "latest" in v.get("aliases", []):
            v["aliases"].remove("latest")

    # Check if our version already exists
    v_exists = any(v["version"] == version for v in versions)

    if v_exists:
        for v in versions:
            if v["version"] == version:
                v["aliases"] = ["latest"]
    else:
        versions.insert(
            0, {"version": version, "title": version, "aliases": ["latest"]}
        )

    versions_file.write_text(json.dumps(versions, indent=2))

    # Create/update 'latest' symlink
    latest_link = SITE_DIR / "latest"
    if latest_link.is_dir() and not latest_link.is_symlink():
        shutil.rmtree(latest_link)
    elif latest_link.exists() or latest_link.is_symlink():
        latest_link.unlink()
    latest_link.symlink_to(version)

    status_msg("setup", f"v{version} set as latest", success=True)


def create_caddyfile():
    """Create Caddyfile for serving the site."""
    status_msg("setup", "Creating Caddyfile...")

    cert_path = (CERTS_DIR / f"{DOMAIN}.pem").resolve()
    key_path = (CERTS_DIR / f"{DOMAIN}-key.pem").resolve()
    site_path = SITE_DIR.resolve()

    caddyfile_content = f"""{{
    auto_https off
    admin off
}}

{DOMAIN}:443 {{
    tls {cert_path} {key_path}
    root * {site_path}
    file_server
    try_files {{path}} {{path}}/ {{path}}.html {{path}}/index.html
}}
"""

    CADDYFILE.write_text(caddyfile_content)
    status_msg("setup", "Caddyfile created", success=True)


def modify_hosts():
    """Add entry to /etc/hosts."""
    global _hosts_modified

    status_msg("setup", "Configuring /etc/hosts...")

    with open("/etc/hosts") as f:
        hosts_content = f.read()

    if DOMAIN in hosts_content:
        status_msg("setup", "Hosts entry already exists", success=True)
        _hosts_modified = True
        return

    entry = f"127.0.0.1 {DOMAIN}"
    run_with_sudo(["sh", "-c", f'echo "{entry}" >> /etc/hosts'])

    _hosts_modified = True
    status_msg("setup", "Hosts entry added", success=True)


def run_server(version: str):
    """Run Caddy server in foreground."""
    global _caddy_process

    status_msg("run", "Starting server...")

    caddyfile_path = CADDYFILE.resolve()
    _caddy_process = subprocess.Popen(
        [
            "caddy",
            "run",
            "--config",
            str(caddyfile_path),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    # Give it a moment to start
    time.sleep(1)

    if _caddy_process.poll() is not None:
        status_msg("run", "Failed to start server", success=False)
        console.print(_caddy_process.stderr.read().decode())
        sys.exit(1)

    status_msg("run", "Server running", success=True)
    console.print()
    console.print(
        Panel(
            f"[bold green]Site available at:[/bold green] [link=https://{DOMAIN}/]https://{DOMAIN}/[/link]\n\n"
            f"Serving current branch as [cyan]v{version}[/cyan] (latest)\n\n"
            "[dim]Press Ctrl+C to stop and clean up[/dim]",
            title="Preview Server Running",
            border_style="green",
        )
    )
    console.print()

    # Wait for the process (blocks until Ctrl+C or process exits)
    try:
        _caddy_process.wait()
    except KeyboardInterrupt:
        pass


@click.group(invoke_without_command=True)
@click.option(
    "--version",
    "-v",
    "version",
    type=str,
    help="Version number to publish current branch as (e.g., 3.0)",
)
@click.option(
    "--clean",
    "-c",
    is_flag=True,
    help="Delete work directory on exit (default: keep for faster restart)",
)
@click.pass_context
def cli(ctx, version: str | None, clean: bool):
    """
    **Preview Release** - Serve training docs locally at production URL.

    Builds the current branch and serves it at https://training.nextflow.io/
    as if it were a released version. Useful for recording videos or previewing
    how a release will look.

    **First-time setup:** Run `mkcert -install` once as your regular user.

    **Note:** Do not run with sudo. The script will prompt for sudo when needed.
    """
    global _keep_files
    _keep_files = not clean

    # If a subcommand was invoked, let it handle things
    if ctx.invoked_subcommand is not None:
        return

    # Ensure we're in the repo root
    if not Path("mkdocs.yml").exists():
        console.print("[red]Error:[/red] Must be run from the training repository root")
        sys.exit(1)

    if version:
        if check_sudo():
            console.print("[red]Error:[/red] Do not run this script with sudo.")
            console.print("  The script will request sudo privileges only when needed.")
            console.print(
                f"  Please re-run: [cyan]./preview_release.py --version {version}[/cyan]"
            )
            sys.exit(1)

        # Register cleanup handlers
        atexit.register(cleanup)
        signal.signal(signal.SIGINT, lambda s, f: sys.exit(0))
        signal.signal(signal.SIGTERM, lambda s, f: sys.exit(0))

        # Setup
        check_dependencies()
        generate_certificates()
        fetch_gh_pages()
        build_docs(version)
        update_versions_json(version)
        create_caddyfile()
        modify_hosts()

        # Run (blocks until Ctrl+C)
        run_server(version)
    else:
        # No version specified, show help
        click.echo(ctx.get_help())


@cli.command()
def status():
    """Show current preview release status."""
    # Ensure we're in the repo root
    if not Path("mkdocs.yml").exists():
        console.print("[red]Error:[/red] Must be run from the training repository root")
        sys.exit(1)

    table = Table(title="Preview Release Status", show_header=False, box=None)
    table.add_column("Key", style="cyan")
    table.add_column("Value")

    # Check hosts entry
    with open("/etc/hosts") as f:
        hosts_has_entry = DOMAIN in f.read()
    table.add_row(
        "Hosts entry", "[green]yes[/green]" if hosts_has_entry else "[dim]no[/dim]"
    )

    # Check work directory
    if WORK_DIR.exists():
        table.add_row("Work directory", "[green]exists[/green]")
        versions_file = SITE_DIR / "versions.json"
        if versions_file.exists():
            versions = json.loads(versions_file.read_text())
            latest_version = next(
                (v["version"] for v in versions if "latest" in v.get("aliases", [])),
                None,
            )
            table.add_row(
                "Versions cached", ", ".join(v["version"] for v in versions[:5]) + "..."
            )
            if latest_version:
                table.add_row("Latest points to", f"[cyan]{latest_version}[/cyan]")
    else:
        table.add_row("Work directory", "[dim]not present[/dim]")

    # Check if server is running
    result = run_cmd(["pgrep", "-f", "caddy"], check=False, capture=True)
    server_running = result.returncode == 0
    table.add_row(
        "Server running", "[green]yes[/green]" if server_running else "[dim]no[/dim]"
    )

    console.print(table)

    # Check if site is accessible
    if server_running and hosts_has_entry:
        result = run_cmd(
            [
                "curl",
                "-s",
                "-o",
                "/dev/null",
                "-w",
                "%{http_code}",
                f"https://{DOMAIN}/",
            ],
            check=False,
            capture=True,
        )
        if result.stdout.strip() == "200":
            console.print(
                f"\n[green]✓[/green] Site accessible at [link=https://{DOMAIN}/]https://{DOMAIN}/[/link]"
            )
        else:
            console.print(
                f"\n[red]✗[/red] Site not responding (HTTP {result.stdout.strip()})"
            )


if __name__ == "__main__":
    cli()
