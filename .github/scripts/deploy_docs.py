#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.11"
# dependencies = ["gitpython"]
# ///

"""Deploy documentation to gh-pages branch."""

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import git


def discover_languages(docs_dir: Path = Path("docs")) -> list[str]:
    """Find all language directories containing mkdocs.yml."""
    return sorted(p.parent.name for p in docs_dir.glob("*/mkdocs.yml"))


def get_version_info(
    event_name: str,
    release_tag: str | None,
    input_version: str | None,
    input_latest: bool,
) -> dict:
    """Determine version, alias, and release status from trigger context."""
    if event_name == "release":
        return {"version": release_tag, "alias": "latest", "is_release": True}
    elif input_version:
        return {
            "version": input_version,
            "alias": "latest" if input_latest else "",
            "is_release": True,
        }
    else:
        return {"version": "0.dev", "alias": "development", "is_release": False}


def combine_artifacts(artifacts_dir: Path, output_dir: Path) -> None:
    """Combine per-language build artifacts into single directory."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # English goes to root
    en_dir = artifacts_dir / "site-en"
    if en_dir.exists():
        shutil.copytree(en_dir, output_dir, dirs_exist_ok=True)

    # Other languages go to subdirectories
    for site_dir in artifacts_dir.glob("site-*"):
        lang = site_dir.name.removeprefix("site-")
        if lang != "en":
            shutil.copytree(site_dir, output_dir / lang, dirs_exist_ok=True)


def update_versions_json(version: str, alias: str, is_release: bool) -> None:
    """Update versions.json for mike compatibility."""
    versions_file = Path("versions.json")
    versions = json.loads(versions_file.read_text()) if versions_file.exists() else []

    # Remove existing entry for this version
    versions = [v for v in versions if v.get("version") != version]

    # Remove 'latest' alias from other versions if setting new latest
    if alias == "latest":
        for v in versions:
            if "aliases" in v:
                v["aliases"] = [a for a in v["aliases"] if a != "latest"]

    # Create new entry
    new_entry = {
        "version": version,
        "title": version,
        "aliases": [alias] if alias else [],
    }

    # Insert at correct position (releases after dev versions, dev at start)
    if is_release:
        insert_idx = sum(1 for v in versions if v.get("version", "").endswith(".dev"))
        versions.insert(insert_idx, new_entry)
    else:
        versions.insert(0, new_entry)

    versions_file.write_text(json.dumps(versions, indent=2) + "\n")


def pin_codespaces_links(version: str, directories: list[str]) -> None:
    """Replace ref=master with ref=VERSION in Codespaces links."""
    for directory in directories:
        if Path(directory).exists():
            subprocess.run(
                [
                    "find",
                    directory,
                    "-type",
                    "f",
                    "-exec",
                    "sed",
                    "-i",
                    f"s|ref=master|ref={version}|g",
                    "{}",
                    "+",
                ],
                check=True,
            )


def deploy(
    version: str, alias: str, is_release: bool, artifacts_dir: Path, force: bool = False
) -> None:
    """Deploy combined artifacts to gh-pages branch."""
    repo = git.Repo(".")

    # Combine artifacts to temp location
    combined = Path("/tmp/combined")
    shutil.rmtree(combined, ignore_errors=True)
    combine_artifacts(artifacts_dir, combined)

    # Configure git
    with repo.config_writer() as cw:
        cw.set_value("user", "name", "github-actions[bot]")
        cw.set_value(
            "user", "email", "41898282+github-actions[bot]@users.noreply.github.com"
        )

    # Fetch and checkout gh-pages
    try:
        repo.remotes.origin.fetch("gh-pages:gh-pages")
        repo.git.checkout("gh-pages")
    except git.GitCommandError:
        # Branch doesn't exist, create orphan
        repo.git.checkout("--orphan", "gh-pages")
        repo.git.rm("-rf", ".", "--ignore-unmatch")

    # Remove old version directories
    for d in [version, alias] if alias else [version]:
        shutil.rmtree(d, ignore_errors=True)

    # Copy combined site
    shutil.copytree(combined, version)
    if alias:
        shutil.copytree(combined, alias)

    # Update versions.json
    update_versions_json(version, alias, is_release)

    # Commit and push
    paths_to_add = [version, "versions.json"]
    if alias:
        paths_to_add.append(alias)
    repo.index.add(paths_to_add)
    repo.index.commit(f"Deploy {version} docs")

    repo.remotes.origin.push("gh-pages", force=force)

    # Pin codespaces links for releases
    if is_release:
        dirs_to_pin = [version, alias] if alias else [version]
        pin_codespaces_links(version, dirs_to_pin)
        repo.index.add(dirs_to_pin)
        try:
            repo.index.commit("[automated] Pin GitHub Codespaces link versions")
            repo.remotes.origin.push("gh-pages")
        except git.GitCommandError:
            pass  # Nothing to commit


def main():
    parser = argparse.ArgumentParser(description="Deploy documentation to gh-pages")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # discover-languages
    subparsers.add_parser("discover-languages", help="Find language directories")

    # get-version-info
    info_parser = subparsers.add_parser(
        "get-version-info", help="Get version/alias from event context"
    )
    info_parser.add_argument("--event-name", required=True)
    info_parser.add_argument("--release-tag", default="")
    info_parser.add_argument("--input-version", default="")
    info_parser.add_argument(
        "--input-latest", type=lambda x: x.lower() == "true", default=True
    )
    info_parser.add_argument("--github-output", help="Path to GITHUB_OUTPUT file")

    # deploy
    deploy_parser = subparsers.add_parser("deploy", help="Deploy docs to gh-pages")
    deploy_parser.add_argument("--version", required=True)
    deploy_parser.add_argument("--alias", default="")
    deploy_parser.add_argument(
        "--is-release", type=lambda x: x.lower() == "true", default=False
    )
    deploy_parser.add_argument("--artifacts-dir", type=Path, default=Path("artifacts"))

    args = parser.parse_args()

    if args.command == "discover-languages":
        print(json.dumps(discover_languages()))

    elif args.command == "get-version-info":
        info = get_version_info(
            args.event_name,
            args.release_tag or None,
            args.input_version or None,
            args.input_latest,
        )
        if args.github_output:
            with open(args.github_output, "a") as f:
                f.write(f"version={info['version']}\n")
                f.write(f"alias={info['alias']}\n")
                f.write(f"is_release={str(info['is_release']).lower()}\n")
        else:
            print(json.dumps(info))

    elif args.command == "deploy":
        force = not args.is_release  # Dev builds use force push
        deploy(args.version, args.alias, args.is_release, args.artifacts_dir, force)


if __name__ == "__main__":
    main()
