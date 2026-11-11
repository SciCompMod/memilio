#############################################################################
# Copyright (C) 2020-2026 MEmilio
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

import subprocess
import tempfile
import logging
import os
from typing import TYPE_CHECKING, Optional
from typing_extensions import Self

from clang.cindex import Cursor, TranslationUnit, Index, CompilationDatabase


class AST:
    """ Handles AST creation and ID assignment.
    Provides functionality to generate the AST from a source file using Clang,
    assign unique IDs to AST nodes, and retrieve nodes or their IDs.
    """

    def __init__(self: Self, source_file: str) -> None:
        """ Basic constructor of the AST class.
        """

        self.source_file = source_file
        self.cursor_id = -1
        self.id_to_val = dict()
        self.val_to_id = dict()
        self.tempfile_path: Optional[str] = None
        self.translation_unit: Optional[TranslationUnit] = None

    def process_ast_file(self):
        """ Process AST file and create the translation_unit.
        """
        # ToDo: error handling for tempfile_path
        idx = Index.create()
        self.translation_unit = idx.read(self.tempfile_path)
        os.remove(self.tempfile_path)
        self.tempfile_path = None
        self._assing_ast_with_ids(self.translation_unit.cursor)

    def create_ast_for_pickle(self, path_to_database: str):
        """ Create the abstract syntax tree and run external processes outside the pickling context.
        This method sets up the necessary arguments, invokes Clang to generate the AST,
        and stores the resulting AST in a temporary file.
        """

        self.cursor_id = -1
        self.id_to_val.clear()
        self.val_to_id.clear()

        file_args = self.get_file_args(path_to_database)

        clang_cmd = self.prepare_clang_command(file_args)

        clang_cmd_result = self.run_clang_command(clang_cmd)

        ast_file_path = self.create_temp_ast_file(clang_cmd_result)

        logging.info(f"AST file created at: {ast_file_path}")

    def get_file_args(self, path_to_database: str):
        """ Extract arguments from the compilation database, excluding unwanted ones.

        :returns: List of relevant file arguments.
        """
        unwanted_arguments = [
            '-Wno-unknown-warning', "--driver-mode=g++", "-O3", "-Werror", "-Wshadow"
        ]
        file_args = []

        # ToDo: Add error handling if compile comands dont exist or the path_to_database is wrong
        dirname = path_to_database
        compdb = CompilationDatabase.fromDirectory(dirname)
        commands = compdb.getCompileCommands(self.source_file)
        for command in commands:
            for argument in command.arguments:
                if argument not in unwanted_arguments:
                    file_args.append(argument)

        return file_args[1:-4]

    def prepare_clang_command(self, file_args):
        """ Prepare the Clang command without executing it.

        :param file_args: List of file arguments for the Clang command.
        :returns: List of command-line arguments for Clang.
        """
        clang_cmd = [
            "clang", self.source_file,
            "-std=c++20", '-emit-ast', '-o', '-']
        clang_cmd.extend(file_args)
        return clang_cmd

    def run_clang_command(self, clang_cmd):
        """ Execute the Clang command.

        :param clang_cmd: List of command-line arguments for Clang.
        :returns: Result of the Clang command execution.
        """
        try:
            clang_cmd_result = subprocess.run(
                clang_cmd, stdout=subprocess.PIPE)
            clang_cmd_result.check_returncode()
            return clang_cmd_result
        except subprocess.CalledProcessError as e:
            logging.error(f"Clang failed: {e}")
            raise RuntimeError("Clang AST generation failed.")

    def create_temp_ast_file(self, clang_cmd_result):
        """ Create a temporary file and write the Clang AST output to it.

        :param clang_cmd_result: Result of the Clang command execution.
        :returns: Path to the temporary AST file.
        """

        with tempfile.NamedTemporaryFile(delete=False) as ast_file:
            ast_file.write(clang_cmd_result.stdout)
            ast_file.flush()
            self.tempfile_path = ast_file.name

    def _assing_ast_with_ids(self, cursor) -> None:
        """ Traverse the AST and assign a unique ID to each node.
        """
        self.cursor_id += 1
        id = self.cursor_id
        self.id_to_val[id] = cursor

        if cursor.hash in self.val_to_id.keys():
            self.val_to_id[cursor.hash].append(id)
        else:
            self.val_to_id[cursor.hash] = [id]

        logging.info(
            f"Node {cursor.spelling or cursor.kind} assigned ID {id}")

        for child in cursor.get_children():
            self._assing_ast_with_ids(child)

    @property
    def root_cursor(self):
        """ Returns the root cursor of the AST."""
        return self.translation_unit.cursor

    def get_node_id(self, cursor: Cursor) -> int:
        """ Returns the id of the current node.
        Extracts the key from the current cursor from the dictonary id_to_val

        :param cursor: The current node of the AST as a cursor object from libclang.
        :returns: The ID of the current node.

        """
        for cursor_id in self.val_to_id[cursor.hash]:

            if self.id_to_val[cursor_id] == cursor:

                return cursor_id
        raise IndexError(f"Cursor {cursor} is out of bounds.")

    def get_node_by_index(self, index: int) -> Cursor:
        """ Returns the node at the specified index position.

        :param index: Node_id from the ast.
        :returns: The node at the specified index.
        """

        if index < 0 or index >= len(self.id_to_val):
            raise IndexError(f"Index {index} is out of bounds.")
        return self.id_to_val[index]
