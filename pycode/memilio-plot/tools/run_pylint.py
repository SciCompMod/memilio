#!/usr/bin/env python3
"""Generate Pylint reports in jsonextended format for the plot package."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    project_dir = Path(__file__).resolve().parent.parent
    repo_root = project_dir.parent
    build_dir = project_dir / "build_pylint"
    build_dir.mkdir(exist_ok=True)
    output_file = build_dir / "pylint_extended.json"

    cmd = [
        sys.executable,
        "-m",
        "pylint",
        "--rcfile",
        str(repo_root / "pylintrc"),
        "--load-plugins",
        "pylint_json2html",
        "--output-format=jsonextended",
        "memilio/",
    ]

    with output_file.open("w", encoding="utf-8") as stream:
        result = subprocess.run(
            cmd,
            cwd=project_dir,
            stdout=stream,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

    if result.stderr:
        sys.stderr.write(result.stderr)

    if result.returncode:
        sys.stderr.write(
            "pylint reported issues; see build_pylint/pylint_extended.json for details.\n"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
