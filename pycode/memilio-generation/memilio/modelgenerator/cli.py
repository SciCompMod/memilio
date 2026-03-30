#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Henrik Zunker
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################

"""
Command-line interface for the model generator.

Usage
-----
::

    memilio-modelgenerator path/to/model.yaml [--output-dir DIR] [--preview]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .generator import Generator
from .validator import ValidationError


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="memilio-modelgenerator",
        description="Generate MEmilio C++ model files and pybind11 bindings from a YAML config.",
    )
    parser.add_argument(
        "config",
        metavar="CONFIG",
        help="Path to the YAML model configuration file.",
    )
    parser.add_argument(
        "--output-dir",
        metavar="DIR",
        default=None,
        help=(
            "Root directory of the MEmilio repository where files are written. "
            "Defaults to the directory two levels above this package "
            "(i.e. the repository root when installed in editable mode)."
        ),
    )
    parser.add_argument(
        "--preview",
        action="store_true",
        help="Print all generated file contents instead of writing them to disk.",
    )

    args = parser.parse_args(argv)

    try:
        if args.config.endswith(".toml"):
            gen = Generator.from_toml(args.config)
        else:
            gen = Generator.from_yaml(args.config)
    except FileNotFoundError:
        print(f"ERROR: config file not found: {args.config}", file=sys.stderr)
        return 1
    except ValidationError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    if args.preview:
        output_dir = Path(args.output_dir) if args.output_dir else Path(
            __file__).resolve().parents[4]
        separator = "=" * 72
        for rel_path, content in gen.render().items():
            print(f"\n{separator}")
            print(f"  NEW FILE: {rel_path}")
            print(separator)
            print(content)
        for rel_path, content in gen.render_patches(output_dir).items():
            if content is None:
                print(f"\n{separator}")
                print(f"  PATCH (already present – no change): {rel_path}")
                print(separator)
            else:
                print(f"\n{separator}")
                print(f"  PATCH: {rel_path}")
                print(separator)
                print(content)
        return 0

    # Output directory
    if args.output_dir is not None:
        output_dir = Path(args.output_dir)
    else:
        # Go up from this file to the repo root
        output_dir = Path(__file__).resolve().parents[4]

    print(f"Writing model files to: {output_dir}")
    gen.write(output_dir)
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
