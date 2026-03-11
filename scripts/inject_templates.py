#!/usr/bin/env python
import argparse
import re
from pathlib import Path


TAG_RE = re.compile(r"<template[\\/](?P<path>[^>]+)>")


def resolve_repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_wiki_dir(repo_root: Path) -> Path:
    return repo_root.with_name(repo_root.name + ".wiki")


def load_template(repo_root: Path, rel_path: str) -> str:
    # Support both forward and backslashes in tag paths.
    rel_path = rel_path.replace("\\", "/")
    path = (repo_root / rel_path).resolve()
    if not path.exists():
        raise FileNotFoundError(f"Template not found: {path}")
    code = path.read_text(encoding="utf-8").rstrip()
    # Strip verbose header block if present.
    if code.startswith("% TEMPLATE-HEADER-BEGIN"):
        end = code.find("% TEMPLATE-HEADER-END")
        if end != -1:
            code = code[end + len("% TEMPLATE-HEADER-END") :].lstrip("\r\n")
    return f"```matlab\n{code}\n```"


def inject_in_text(text: str, repo_root: Path) -> str:
    def repl(match: re.Match) -> str:
        rel_path = match.group("path").strip()
        return load_template(repo_root, rel_path)

    return TAG_RE.sub(repl, text)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Replace <template/...> tags in wiki pages with the content of the referenced scripts."
    )
    parser.add_argument(
        "--wiki",
        type=str,
        default=None,
        help="Path to the wiki folder (default: <repo>.wiki sibling).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print files that would change without writing.",
    )
    args = parser.parse_args()

    repo_root = resolve_repo_root()
    wiki_dir = Path(args.wiki) if args.wiki else resolve_wiki_dir(repo_root)

    if not wiki_dir.exists():
        raise SystemExit(f"Wiki folder not found: {wiki_dir}")

    changed = []
    for md_path in wiki_dir.glob("*.md"):
        original = md_path.read_text(encoding="utf-8")
        updated = inject_in_text(original, repo_root)
        if updated != original:
            changed.append(md_path)
            if not args.dry_run:
                md_path.write_text(updated, encoding="utf-8")

    if args.dry_run:
        for p in changed:
            print(f"would update: {p}")
    else:
        for p in changed:
            print(f"updated: {p}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
