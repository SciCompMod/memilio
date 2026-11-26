#############################################################################
# Copyright (C) 2020-2025 MEmilio
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
##############################################################################

import os
import logging
import multiprocessing
import subprocess
from typing_extensions import Self

from memilio.generation.scanner_config import get_config
from memilio.generation.ast import AST

from clang.cindex import Config
    

class ASTHandler:
    """ ASTHandler class for creating ASTs.
    """

    def __init__(self: Self) -> None:
        """ Basic constructor of the ASTHandler class
        """
        try_set_libclang_path(
            get_config().libclang_library_path)
        self.ast_list = []
        self.source_files = []

    def handle_ast_creation(self: Self) -> None:
        """ Handle the creation of the ASTs. Single or Parallel.
        """
        if len(self.source_files) == 1:
            self.single_creation()
        else:
            self.parallel_creation()

    def add_source_file(self: Self, path_to_source_file: str, path_to_database: str = "") -> None:
        """ Add source file to list.

        :param path_to_source_file: Absolute path to source file.
        :param path_to_database: Package resource for the database needed to build the source file.
                                 Default uses the one shipped with building MEmilio.
        """
        if (os.path.isfile(path_to_source_file)):
            self.source_files.append((path_to_source_file, path_to_database))
            logging.info(f"Adding source file: {path_to_source_file}")

    def parallel_creation(self: Self) -> None:
        """ Create ASTs for all source files using multiprocessing.

        This method distributes the AST creation process across multiple CPU cores. 
        Each source file is processed independently using the `create_ast` function, 
        and the resulting ASTs are stored in `self.ast_list`.
        """
        with multiprocessing.Pool() as pool:
            try:

                results = pool.starmap(
                    create_ast, [(database, source_file) for source_file, database in self.source_files])

                for ast in results:
                    try:
                        ast.process_ast_file()
                        self.ast_list.append(ast)

                    except Exception as e:
                        logging.info(
                            f"Error processing AST file {ast.source_file}: {e}")
            except Exception as e:
                logging.info(f"Error in parallel AST creation: {e}")

    def single_creation(self: Self) -> None:
        """ Create AST for a single source file.
        """
        ast = create_ast(self.source_files[0][1], self.source_files[0][0])
        ast.process_ast_file()
        self.ast_list.append(ast)

    def get_ast(self: Self, idx: int) -> AST:
        """ Get AST out of ast_list.

        :param id: idx of the AST.
        :returns: AST object.
        """
        return self.ast_list[idx]


def create_ast(path_to_database: str, source_file: str) -> AST:
    """ Create AST for a file.

    :param path_to_database: Path to compile commands
    :param source_file: Sourcefile path
    :returns: The AST object.
    """
    ast = AST(source_file)
    ast.create_ast_for_pickle(path_to_database)

    return ast

def try_set_libclang_path(path: str) -> None:
    """Try to set the file path for the libclang library.
    If its already set, the returned Exception gets caught and discarded.
    If the given path string is empty or None, the path is determined with a call on the terminal.

    :param path: Path to the library files of libclang. Can be an empty string.
    """
    # Check if path was set in config. If not, try to get it with cmd.
    if (not path):
        clang_cmd = [get_config().clang_executable, '-print-file-name=']
        clang_cmd_result = subprocess.check_output(clang_cmd)
        path, dirname = os.path.split(clang_cmd_result)
        while ("llvm" not in str(dirname)):
            path, dirname = os.path.split(path)
            if str(path) == '':
                print("Error: define libclang path in config.json")
                raise (Exception)
        path = os.path.join(path.decode(), dirname.decode(),
                            "lib", "libclang.so")
    try:
        if not Config.loaded:
            Config.set_library_file(path)

    except Exception as e:

        if str(e) != "library file must be set before using any other functionalities in libclang.":
            raise (e)