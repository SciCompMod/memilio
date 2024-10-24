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
import subprocess
from typing import Any, List, TextIO
from graphviz import Digraph
from clang.cindex import Cursor
import tempfile
import logging

from memilio.generation import scanner, utility, intermediate_representation, scanner_config


class ASTViz:

    def __init__(self):
        ast_cursor = Cursor()
        self.node_counter = 0

    def output_ast_terminal(self, ast_cursor: Cursor) -> None:
        """
        Output the abstract syntax tree to terminal. Not formatted.
        """
        _output_cursor_and_children_print(ast_cursor)

    def output_ast_file(self, ast_cursor: Cursor) -> None:
        """
        Output the abstract syntax tree to file. Not formatted.
        """
        with open('output_ast.txt', 'a') as f:
            _output_cursor_and_children_file(ast_cursor, f)
            print('AST-file written to ' + str(os.path.abspath(f.name)))

    # def output_ast_png_all(self, ast_cursor: Cursor) -> None:
    #     """
    #     Output the abstract syntax tree as a graph using Graphviz and save it to a file.
    #     """

    #     graph = Digraph(format='png')

    #     _output_cursor_and_children_graphviz_digraph(ast_cursor, graph)

    #     output_file = 'ast_graph'

    #     graph.render(output_file, view=False)

    #     output_path = os.path.abspath(f"{output_file}.png")
    #     print(f'AST png written to {output_path}')

    def output_ast_png_limit(self, ast_cursor: Cursor, max_depth: int, ) -> None:
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

    # def output_ast_digraph(self, ast_cursor: Cursor) -> None:
    #     """
    #     Output the abstract syntax tree as a graph using Graphviz and save it to a file.
    #     """

    #     graph = Digraph(format='dot')

    #     _output_cursor_and_children_graphviz_digraph(ast_cursor, graph)

    #     output_file = 'ast_graph'

    #     graph.save(output_file + '.dot')

    #     # Ausgabe der Pfad-Informationen
    #     output_path = os.path.abspath(f"{output_file}.dot")
    #     print(f'AST digraph written to {output_path}')

    def output_ast_formatted(self, ast_cursor: Cursor) -> None:
        """
        Output the abstract syntax tree to a file. Formatted.
        """

        with open('output_ast_format.txt', 'w') as f:
            _output_cursor_and_children_text(ast_cursor, f)

        print("AST format written to " + str(os.path.abspath(f.name)))


def indent2(level: int) -> str:
    """
    Create an indentation based on the level.
    """
    return '│   ' * level + '├── '


def _output_cursor_and_children_text(cursor: Cursor, f: TextIO, level: int = 0, node_counter: int = 0) -> None:
    """
    Output of the cursor and its children in text format, with highlighting for folder, spelling and child type.

    @param cursor: Der aktuelle Knoten des AST als Cursor-Objekt von libclang.
    @param f: Offenes Dateiobjekt für die Ausgabe.
    @param level: Die aktuelle Tiefe im AST für Einrückungszwecke.
    """

    node_counter += 1
    cursor_id = node_counter

    cursor_kind = f"<CursorKind.{cursor.kind.name}>"
    file_path = cursor.location.file.name if cursor.location.file else ""
    cursor_label = f'ID={cursor_id} {cursor.spelling} {cursor_kind} {file_path}' if cursor.spelling else f'{cursor_kind} {file_path}'

    if cursor.spelling:
        cursor_label = (f'ID={cursor_id}{cursor.spelling} '
                        f'{cursor_kind}   '
                        f'{file_path}')
    else:
        cursor_label = f'{cursor_kind} [{file_path}]'

    f.write(indent2(level) + cursor_label + '\n')

    for child in cursor.get_children():
        _output_cursor_and_children_text(child, f, level + 1)


def _output_cursor_and_children_graphviz_digraph(cursor: Cursor, graph: Digraph, max_d: int, current_d: int, parent_node: str = None) -> None:
    """
    Output the cursor and its children as a graph using Graphviz.

    @param cursor: The current node of the AST as a Cursor object from libclang.
    @param graph: Graphviz Digraph object where the nodes and edges will be added.
    @param parent_node: Name of the parent node in the graph (None for the root node).
    """
    if current_d > max_d:
        return
    # Define a label for the current node in the graph
    node_label = f"{cursor.kind.name}\n({cursor.spelling})" if cursor.spelling else cursor.kind.name
    # Unique node ID using kind and hash
    current_node = f"{cursor.kind.name}_{cursor.hash}"

    # Add the current node to the graph
    graph.node(current_node, label=node_label)

    # If there is a parent node, create an edge from the parent to the current node
    if parent_node:
        graph.edge(parent_node, current_node)

    # Check if the cursor is a reference, and add a reference node if so
    if cursor.kind.is_reference():
        referenced_label = f"ref_to_{cursor.referenced.kind.name}\n({cursor.referenced.spelling})"
        referenced_node = f"ref_{cursor.referenced.hash}"
        graph.node(referenced_node, label=referenced_label)
        graph.edge(current_node, referenced_node)

    # Recurse for children of this cursor
    for child in cursor.get_children():
        _output_cursor_and_children_graphviz_digraph(
            child, graph, max_d, current_d + 1, current_node)


def _output_cursor_and_children_file(
        cursor: Cursor, f: TextIO, spaces: int = 0) -> None:
    """ 
    Output this cursor and its children with minimal formatting to a file.

    @param cursor Represents the current node of the AST as an Cursor object from libclang.
    @param f Open file object for output.
    @param spaces Number of spaces.
    """
    output_cursor_file(cursor, f, spaces)
    if cursor.kind.is_reference():
        f.write(indent(spaces) + 'reference to:\n')
        output_cursor_file(cursor.referenced, f, spaces+1)

    # Recurse for children of this cursor
    has_children = False
    for c in cursor.get_children():
        if not has_children:
            f.write(indent(spaces) + '{\n')
            has_children = True
        _output_cursor_and_children_file(c, f, spaces+1)

    if has_children:
        f.write(indent(spaces) + '}\n')


def output_cursor_file(cursor: Cursor, f: TextIO, spaces: int) -> None:
    """ 
    Low level cursor output to a file.

    @param cursor Represents the current node of the AST as an Cursor object from libclang.
    @param f Open file object for output.
    @param spaces Number of spaces.
    """
    spelling = ''
    displayname = ''

    if cursor.spelling:
        spelling = cursor.spelling
    if cursor.displayname:
        displayname = cursor.displayname
    kind = cursor.kind

    f.write(indent(spaces) + spelling + ' <' + str(kind) + '> ')
    if cursor.location.file:
        f.write(cursor.location.file.name + '\n')
    f.write(indent(spaces+1) + '"' + displayname + '"\n')


def indent(spaces: int) -> str:
    """ 
    Indentation string for pretty-printing.

    @param spaces Number of spaces.
    """
    return '  '*spaces


def _output_cursor_and_children_print(cursor: Cursor, spaces: int = 0) -> None:
    """ 
    Output this cursor and its children with minimal formatting to the terminal.

    @param cursor Represents the current node of the AST as an Cursor object from libclang.
    @param spaces [Default = 0] Number of spaces.
    """

    _output_cursor_print(cursor, spaces)
    if cursor.kind.is_reference():
        print(indent2(spaces) + 'reference to:')
        _output_cursor_print(cursor.referenced, spaces+1)

    # Recurse for children of this cursor
    has_children = False
    for c in cursor.get_children():
        if not has_children:
            print(indent2(spaces) + '{')
            has_children = True
        _output_cursor_and_children_print(c, spaces+1)

    if has_children:
        print(indent2(spaces) + '}')


def _output_cursor_print(cursor: Cursor, spaces: int) -> None:
    """ 
    Low level cursor output to the terminal.

    @param cursor Represents the current node of the AST as an Cursor object from libclang.
    @param spaces Number of spaces.
    """
    spelling = ''
    displayname = ''

    if cursor.spelling:
        spelling = cursor.spelling
    if cursor.displayname:
        displayname = cursor.displayname
    kind = cursor.kind

    print(indent(spaces) + spelling, '<' + str(kind) + '>')
    print(indent(spaces+1) + '"' + displayname + '"')
