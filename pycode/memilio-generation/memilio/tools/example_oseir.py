#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
#
# Authors: Maximilian Betz
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
Example for the ode seir model.
"""
import argparse
import sys
import os

if sys.version_info >= (3, 9):
    # For python 3.9 and newer
    import importlib.resources as importlib_resources
else:
    # For older python versions
    import importlib_resources

from memilio.generation import Generator, Scanner, ScannerConfig


def run_memilio_generation(print_ast=False):
    # Define ScannerConfig and initialize Scanner
    pkg = importlib_resources.files("memilio.generation")
    with importlib_resources.as_file(pkg.joinpath('../tools/config.json')) as path:
        with open(path) as file:
            conf = ScannerConfig.schema().loads(file.read(), many=True)[0]
    scanner = Scanner(conf)

    # Extract results of Scanner into a intermed_repr
    intermed_repr = scanner.extract_results()

    # Generate code
    generator = Generator()
    generator.create_substitutions(intermed_repr)
    generator.generate_files(intermed_repr)

    # Print the abstract syntax tree to a file
    if (print_ast):
        scanner.output_ast_file()


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        'memilio_generation',
        description='Simple example demonstrating the memilio-generation package.')
    arg_parser.add_argument('-p', '--print_ast',
                            action='store_const', const=True, default=False)
    args = arg_parser.parse_args()
    run_memilio_generation(**args.__dict__)
