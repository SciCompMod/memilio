#############################################################################
# Copyright (C) 2020-2026 MEmilio
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
import time
import os

from memilio.generation import Generator, Scanner, get_config, ASTHandler, Visualization

def run_memilio_generation(print_ast=False):
    """

    :param print_ast:  (Default value = False)
    """
    # Get ScannerConfig and initialize Scanner
    conf = get_config()

    file_path = os.path.dirname(os.path.abspath(__file__))

    asts = ASTHandler(conf)
    asts.add_source_file(os.path.abspath(os.path.join(
        file_path, "..", "..", "..", "..", "cpp", "models", "ode_seir", "model.cpp")))
    asts.handle_ast_creation()
    ast = asts.get_ast(0)
    scanner = Scanner(conf)

    intermed_repr = scanner.extract_results(ast.root_cursor)

    generator = Generator()
    generator.create_substitutions(intermed_repr)
    generator.generate_files(intermed_repr)

    if (print_ast):
        Visualization.output_ast_formatted(
            ast, ast.get_node_by_index(0))
        # Visualization.output_ast_terminal(ast, ast.get_node_by_index(3))
        # Visualization.output_ast_png(ast.get_node_by_index(0), 3)


if __name__ == "__main__":
    start_front = time.time()
    arg_parser = argparse.ArgumentParser(
        'memilio_generation',
        description='Simple example demonstrating the memilio-generation package.')
    arg_parser.add_argument('-p', '--print_ast',
                            action='store_const', const=True, default=False)
    args = arg_parser.parse_args()
    run_memilio_generation(**args.__dict__)
