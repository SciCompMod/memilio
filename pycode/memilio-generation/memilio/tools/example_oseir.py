#############################################################################
# Copyright (C) 2020-2025 MEmilio
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
import time
import matplotlib.pyplot as plt
import statistics
import json
import os
from clang import cindex
from clang.cindex import CursorKind

if sys.version_info >= (3, 9):
    # For python 3.9 and newer
    import importlib.resources as importlib_resources
else:
    # For older python versions
    import importlib_resources

from memilio.generation import Generator, Scanner, ScannerConfig, AST, ASTHandler
from memilio.generation.graph_visualization import Visualization


def run_memilio_generation(print_ast=False):
    """

    :param print_ast:  (Default value = False)

    """
    # Define ScannerConfig and initialize Scanner
    pkg = importlib_resources.files("memilio.generation")
    with importlib_resources.as_file(pkg.joinpath('../tools/config.json')) as path:
        with open(path) as file:
            conf = ScannerConfig.schema().loads(file.read(), many=True)[0]

    file_path = os.path.dirname(os.path.abspath(__file__))

    conf.source_file = os.path.abspath(os.path.join(
        file_path, "..", "..", "..", "..", "cpp", "models", "ode_secirvvs", "model.cpp"))
    conf.target_folder = file_path

    scanner = Scanner(conf)
    asts = ASTHandler(conf)
    asts.add_source_file("model.cpp")

    asts.handle_ast_creation()

    aviz = Visualization()

    ast = asts.get_ast_by_id(0)

    intermed_repr = scanner.extract_results(ast.root_cursor)

    generator = Generator()
    generator.create_substitutions(intermed_repr)
    generator.generate_files(intermed_repr)

    # Print the abstract syntax tree to a file
    if (print_ast):
        aviz.output_ast_formatted(
            asts.get_ast_by_id(0), ast.get_node_by_index(0))
        # aviz.output_ast_terminal(asts.get_ast_by_id(0), ast.get_node_by_index(2))
        # aviz.output_ast_png(ast.get_node_by_index(1), 1)


if __name__ == "__main__":
    start_front = time.time()
    arg_parser = argparse.ArgumentParser(
        'memilio_generation',
        description='Simple example demonstrating the memilio-generation package.')
    arg_parser.add_argument('-p', '--print_ast',
                            action='store_const', const=True, default=False)
    args = arg_parser.parse_args()
    run_memilio_generation(**args.__dict__)
