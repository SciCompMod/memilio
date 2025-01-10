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
"""
@file ast.py
@brief Create the ast and assign ids. Get ids and nodes. 
"""
import subprocess
import tempfile
import logging
from clang.cindex import Cursor, TranslationUnit, Index, CompilationDatabase
from typing import TYPE_CHECKING
from memilio.generation import utility


if TYPE_CHECKING:
    from memilio.generation import ScannerConfig

from typing_extensions import Self


class AST:
    """! Create the ast and assign ids.
    Functions for getting nodes and node ids.
    """

    def __init__(self: Self, conf: "ScannerConfig") -> None:
        self.config = conf
        self.cursor_id = -1
        self.id_to_val = dict()
        self.val_to_id = dict()
        self.cursor = None
        self.translation_unit = self.create_ast()

    def create_ast(self: Self) -> TranslationUnit:
        """! Create an abstract syntax tree for the main model.cpp file with a corresponding CompilationDatabase.
        A compile_commands.json is required (automatically generated in the build process).
        """
        self.cursor_id = -1
        self.id_to_val.clear()
        self.val_to_id.clear()

        idx = Index.create()

        file_args = []

        unwanted_arguments = [
            '-Wno-unknown-warning', "--driver-mode=g++", "-O3", "-Werror", "-Wshadow"
        ]

        dirname = utility.try_get_compilation_database_path(
            self.config.skbuild_path_to_database)
        compdb = CompilationDatabase.fromDirectory(dirname)
        commands = compdb.getCompileCommands(self.config.source_file)
        for command in commands:
            for argument in command.arguments:
                if argument not in unwanted_arguments:
                    file_args.append(argument)
        file_args = file_args[1:-4]

        clang_cmd = [
            "clang-14", self.config.source_file,
            "-std=c++17", '-emit-ast', '-o', '-']
        clang_cmd.extend(file_args)

        try:
            clang_cmd_result = subprocess.run(
                clang_cmd, stdout=subprocess.PIPE)
            clang_cmd_result.check_returncode()
        except subprocess.CalledProcessError as e:
            # Capture standard error and output
            logging.error(
                f"Clang failed with return code {e.returncode}. Error: {clang_cmd_result.stderr.decode()}")
            raise RuntimeError(
                f"Clang AST generation failed. See error log for details.")

        # Since `clang.Index.read` expects a file path, write generated abstract syntax tree to a
        # temporary named file. This file will be automatically deleted when closed.
        with tempfile.NamedTemporaryFile() as ast_file:
            ast_file.write(clang_cmd_result.stdout)
            translation_unit = idx.read(ast_file.name)

        self._assing_ast_with_ids(translation_unit.cursor)

        logging.info("AST generation completed successfully.")

        return translation_unit

    def _assing_ast_with_ids(self, cursor: Cursor) -> None:
        """! Traverse the AST and assign a unique ID to each node during traversal.

        @param cursor: The current node (Cursor) in the AST to traverse.
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
        return self.translation_unit.cursor

    def get_node_id(self, cursor: Cursor) -> int:
        """! Returns the id of the current node.

        Extracts the key from the current cursor from the dictonary id_to_val

        @param cursor: The current node of the AST as a cursor object from libclang.
        """
        for cursor_id in self.val_to_id[cursor.hash]:

            if self.id_to_val[cursor_id] == cursor:

                return cursor_id
        raise IndexError(f"Cursor {cursor} is out of bounds.")

    def get_node_by_index(self, index: int) -> Cursor:
        """! Returns the node at the specified index position.

        @param index: Node_id from the AST.
        """

        if index < 0 or index >= len(self.id_to_val):
            raise IndexError(f"Index {index} is out of bounds.")
        return self.id_to_val[index]
