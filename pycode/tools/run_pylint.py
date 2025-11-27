#!/usr/bin/env python3
"""Run pylint for a given memilio package and write the jsonextended report."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args(pycode_root: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run pylint for a specific memilio pycode package."
    )
    parser.add_argument(
        "--package-dir",
        help="Directory name under pycode/ (e.g. memilio-epidata).",
    )
    args = parser.parse_args()

    # if package_dir is given, use it directly
    if args.package_dir:
        return args

    # If no package dir is provided, try to infer it from the current working directory
    # Used, when the script is called from inside pycode/<package>/...
    try:
        relative = Path.cwd().resolve().relative_to(pycode_root)
    except ValueError:
        parser.error(
            "--package-dir is required when not running from inside pycode/<package>"
        )

    if not relative.parts:
        parser.error(
            "--package-dir is required when not running from inside pycode/<package>"
        )

    args.package_dir = relative.parts[0]
    return args


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    pycode_root = repo_root / "pycode"
    args = parse_args(pycode_root)
    package_dir = pycode_root / args.package_dir

    if not package_dir.is_dir():
        sys.stderr.write(f"Package directory not found: {package_dir}\n")
        return 1

    rcfile = pycode_root / "pylintrc"
    if not rcfile.is_file():
        sys.stderr.write(f"Cannot find pylintrc at {rcfile}\n")
        return 1

    build_dir = package_dir / "build_pylint"
    build_dir.mkdir(exist_ok=True)
    output_file = build_dir / "pylint_extended.json"

    cmd = [
        sys.executable,
        "-m",
        "pylint",
        "--rcfile",
        str(rcfile),
        "--load-plugins",
        "pylint_json2html",
        "--output-format=jsonextended",
        "memilio/",
    ]

    with output_file.open("w", encoding="utf-8") as stream:
        result = subprocess.run(
            cmd,
            cwd=package_dir,
            stdout=stream,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

    if result.stderr:
        sys.stderr.write(result.stderr)

    if result.returncode:
        sys.stderr.write(
            f"pylint exited with {result.returncode}; "
            f"see {output_file} for details.\n"
        )

    # Always exit 0 so CI can upload and publish the report even when pylint fails.
    return 0


if __name__ == "__main__":
    sys.exit(main())
