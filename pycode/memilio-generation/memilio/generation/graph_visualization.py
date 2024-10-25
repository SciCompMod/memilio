#############################################################################
# Copyright (C) 2020-2024 MEmilio
#
# Authors: Daniel Richter
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
from typing import TextIO
from graphviz import Digraph
from clang.cindex import Cursor
from memilio.generation import Scanner


class ASTViz:
    """
    Class for plotting the abstract-syntax-tree in different formats
    """

    @staticmethod
    def output_ast_terminal(ast_cursor: Cursor) -> None:
        """
        Output the abstract syntax tree to terminal. Not formatted.
        """

        _output_cursor_and_children_print(ast_cursor)

        print(f'AST-terminal written.')

    @staticmethod
    def output_ast_png(ast_cursor: Cursor, max_depth: int, ) -> None:
        """
        Output the abstract syntax tree to a .png. Set the starting node and the max depth.
        """

        graph = Digraph(format='png')

        _output_cursor_and_children_graphviz_digraph(
            ast_cursor, graph, max_depth, 0)

        output_file = 'ast_graph_png_limit'

        graph.render(output_file, view=False)

        output_path = os.path.abspath(f"{output_file}.png")
        print(f'AST-png written to {output_path}')

    @staticmethod
    def output_ast_formatted(ast_cursor: Cursor) -> None:
        """
        Output the abstract syntax tree to a file. Formatted.
        """

        with open('output_ast_format.txt', 'w') as f:
            _output_cursor_and_children_text(ast_cursor, f)

        print("AST-format written to " + str(os.path.abspath(f.name)))


def indent2(level: int) -> str:
    """
    Create an indentation based on the level.
    """
    return '│   ' * level + '├── '


def _output_cursor_and_children_text(cursor: Cursor, f: TextIO, level: int = 0, cursor_id: int = 0) -> int:
    """
    Output of the cursor and its children in text format, with highlighting for folder, spelling and child type.

    @param cursor: The current node of the AST as a cursor object from libclang.
    @param f: Open file object for output.
    @param level: The current depth in the AST for indentation purposes.
    @param cursor_id: A unique ID to identify each cursor.
    """

    cursor_kind = f"<CursorKind.{cursor.kind.name}>"
    file_path = cursor.location.file.name if cursor.location.file else ""

    if cursor.spelling:
        cursor_label = (f'ID={cursor_id} {cursor.spelling} '
                        f'{cursor_kind}   '
                        f'{file_path}')
    else:
        cursor_label = f'ID={cursor_id} {cursor_kind} [{file_path}]'

    f.write(indent2(level) + cursor_label + '\n')

    for child in cursor.get_children():
        cursor_id = _output_cursor_and_children_text(
            child, f, level + 1, cursor_id + 1)
    return cursor_id


def _output_cursor_and_children_print(cursor: Cursor, level: int = 0, cursor_id: int = 0) -> int:
    """
    Prints the current cursor and its child elements in text format,
    with highlighting for folder, case and node type.

    @param cursor: The current node of the AST as a libclang cursor object.
    @param f: Open file object for output.
    @param level: The current depth in the AST for indentation purposes.
    @param cursor_id: A unique ID to identify each cursor.
    @return: The current cursor ID for subsequent nodes.
    """

    cursor_id = Scanner.get_node_id(cursor) - 1

    cursor_kind = f"<CursorKind.{cursor.kind.name}>"
    file_path = cursor.location.file.name if cursor.location.file else ""
    cursor_label = (
        f'ID={cursor_id} {cursor.spelling} {cursor_kind} {file_path}'
        if cursor.spelling else
        f'ID={cursor_id} {cursor_kind} {file_path}'
    )

    print(indent2(level) + cursor_label)

    for child in cursor.get_children():
        cursor_id = _output_cursor_and_children_print(
            child, level + 1, cursor_id + 1)
    return cursor_id


def _output_cursor_and_children_graphviz_digraph(cursor: Cursor, graph: Digraph, max_d: int, current_d: int, parent_node: str = None) -> None:
    """
    Output the cursor and its children as a graph using Graphviz.

    @param cursor: The current node of the AST as a Cursor object from libclang.
    @param graph: Graphviz Digraph object where the nodes and edges will be added.
    @param parent_node: Name of the parent node in the graph (None for the root node).
    """
    if current_d > max_d:
        return

    node_label = f"{cursor.kind.name}\n({cursor.spelling})" if cursor.spelling else cursor.kind.name

    current_node = f"{cursor.kind.name}_{cursor.hash}"

    graph.node(current_node, label=node_label)

    if parent_node:
        graph.edge(parent_node, current_node)

    if cursor.kind.is_reference():
        referenced_label = f"ref_to_{cursor.referenced.kind.name}\n({cursor.referenced.spelling})"
        referenced_node = f"ref_{cursor.referenced.hash}"
        graph.node(referenced_node, label=referenced_label)
        graph.edge(current_node, referenced_node)

    for child in cursor.get_children():
        _output_cursor_and_children_graphviz_digraph(
            child, graph, max_d, current_d + 1, current_node)
