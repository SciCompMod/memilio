#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Maximilian Betz, Daniel Richter
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

import os
import logging
from typing import Callable
from graphviz import Digraph
from clang.cindex import Cursor
from memilio.generation.ast import AST


class Visualization:
    """! Class for plotting the abstract syntax tree in different formats.
    """
    @staticmethod
    def output_ast_terminal(ast: AST, cursor: Cursor) -> None:
        """! Output the abstract syntax tree to terminal.

        @param ast: ast object from AST class.
        @param cursor: The current node of the AST as a cursor object from libclang.
        """

        def terminal_writer(level: int, cursor_label: str) -> None:
            print(indent(level) + cursor_label)

        _output_cursor_and_children(cursor, ast, terminal_writer)

        logging.info("AST-Terminal written.")

    @staticmethod
    def output_ast_png(cursor: Cursor, max_depth: int, output_file_name: str = 'ast_graph') -> None:
        """! Output the abstract syntax tree to a .png. Set the starting node and the max depth.

        To save the abstract syntax tree as an png with a starting node and a depth u cann use the following command

        Example command: aviz.output_ast_png(ast.get_node_by_index(1), 2)

        aviz -> instance of the Visualization class.

        ast -> instance of the AST class.

        .get_node_by_index -> get a specific node by id (use .output_ast_formatted to see node ids)

        The number 2 is a example for the depth the graph will show

        @param cursor: The current node of the AST as a cursor object from libclang.
        @param max_depth: Maximal depth the graph displays.
        """

        graph = Digraph(format='png')

        _output_cursor_and_children_graphviz_digraph(
            cursor, graph, max_depth, 0)

        graph.render(filename=output_file_name, view=False)

        output_path = os.path.abspath(f"{output_file_name}.png")
        logging.info(f"AST-png written to {output_path}")

    @staticmethod
    def output_ast_formatted(ast: AST, cursor: Cursor, output_file_name: str = 'ast_formatted.txt') -> None:
        """! Output the abstract syntax tree to a file.

        @param ast: ast object from AST class.
        @param cursor: The current node of the AST as a cursor object from libclang.
        """

        with open(output_file_name, 'w') as f:
            def file_writer(level: int, cursor_label: str) -> None:
                f.write(indent(level) + cursor_label + newline())
            _output_cursor_and_children(cursor, ast,  file_writer)

        output_path = os.path.abspath(f"{output_file_name}")
        logging.info(f"AST-formatted written to {output_path}")


def indent(level: int) -> str:
    """! Create an indentation based on the level.
    """
    return '│   ' * level + '├── '


def newline() -> str:
    """! Create a new line.
    """
    return '\n'


def _output_cursor_and_children(cursor: Cursor, ast: AST, writer: Callable[[int, str], None], level: int = 0) -> None:
    """!Generic function to output the cursor and its children with a specified writer.

    @param cursor: The current node of the AST as a libclang cursor object.
    @param ast: AST object from the AST class.
    @param writer: Function that takes `level` and `cursor_label` and handles output.
    @param level: The current depth in the AST for indentation purposes.
    """

    cursor_id = ast.get_node_id(cursor)

    cursor_kind = f"<CursorKind.{cursor.kind.name}>"
    file_path = cursor.location.file.name if cursor.location.file else ""
    line_number = cursor.location.line if cursor.location.file else ""

    if cursor.spelling:
        cursor_label = (f'ID:{cursor_id} {cursor.spelling} '
                        f'{cursor_kind}   '
                        f'{file_path}:{line_number}')
    else:
        cursor_label = f'ID:{cursor_id} {cursor_kind} {file_path}:{line_number}'

    writer(level, cursor_label)

    for child in cursor.get_children():
        _output_cursor_and_children(
            child, ast, writer, level + 1)


def _output_cursor_and_children_graphviz_digraph(cursor: Cursor, graph: Digraph, max_d: int, current_d: int, parent_node: str = None) -> None:
    """! Output the cursor and its children as a graph using Graphviz.

    @param cursor: The current node of the AST as a Cursor object from libclang.
    @param graph: Graphviz Digraph object where the nodes and edges will be added.
    @param max_d: Maximal depth.
    @param current_d: Current depth.
    @param parent_node: Name of the parent node in the graph (None for the root node).
    """

    if current_d > max_d:
        return

    node_label = f"{cursor.kind.name}{newline()}({cursor.spelling})" if cursor.spelling else cursor.kind.name

    current_node = f"{cursor.kind.name}_{cursor.hash}"

    graph.node(current_node, label=node_label)

    if parent_node:
        graph.edge(parent_node, current_node)

    if cursor.kind.is_reference():
        referenced_label = f"ref_to_{cursor.referenced.kind.name}{newline()}({cursor.referenced.spelling})"
        referenced_node = f"ref_{cursor.referenced.hash}"
        graph.node(referenced_node, label=referenced_label)
        graph.edge(current_node, referenced_node)

    for child in cursor.get_children():
        _output_cursor_and_children_graphviz_digraph(
            child, graph, max_d, current_d + 1, current_node)
