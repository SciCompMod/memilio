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


from typing import TypeVar
from memilio.generation import ScannerConfig, AST
from clang.cindex import TranslationUnit, Index
import os
import logging
import copy
import multiprocessing

T = TypeVar("T", bound="ASTHandler")


class ASTHandler:
    """ ASTHandler class for creating ASTs.
    """

    def __init__(self: T, conf: ScannerConfig) -> None:
        """ Basic constructor of the ASTHandler class

        :param conf: ScannerConfig dataclass with the configurations.
        """
        self.conf = conf
        self.ast_list = []
        self.source_files = []

    def handle_ast_creation(self: T) -> None:
        """ Handle the creation of the ASTs. Single or Parallel.
        """
        if len(self.source_files) == 1:
            self.single_creation()
        else:
            self.parallel_creation()

    def add_source_file(self: T, source_file_str: str) -> None:
        """ Add source file to list.

        :param source_file_str: String that represents another file.
        """
        if (os.path.isfile(source_file_str)):
            self.source_files.append(source_file_str)
            logging.info(f"Adding source file: {source_file_str}")

    def parallel_creation(self: T) -> None:
        """ Create ASTs for all source files using multiprocessing.

        This method distributes the AST creation process across multiple CPU cores. 
        Each source file is processed independently using the `create_ast` function, 
        and the resulting ASTs are stored in `self.ast_list`.
        """
        with multiprocessing.Pool() as pool:
            try:

                results = pool.starmap(
                    create_ast, [(self.conf, source_file) for source_file in self.source_files])

                for ast in results:
                    try:
                        ast.process_ast_file()
                        self.ast_list.append(ast)

                    except Exception as e:
                        logging.info(
                            f"Error processing AST file {ast.source_file}: {e}")
            except Exception as e:
                logging.info(f"Error in parallel AST creation: {e}")

    def single_creation(self: T) -> None:
        """ Create AST for a single source file.

        :param conf: ScannerConfig dataclass with the configurations.
        :returns: TranslationUnit object.
        """
        ast = create_ast(self.conf, self.source_files[0])
        ast.process_ast_file()
        self.ast_list.append(ast)

    def get_ast_by_id(self: T, id: int) -> AST:
        """ Get AST by id.

        :param id: ID of the AST.
        :returns: AST object.
        """
        return self.ast_list[id]


def create_ast(conf_template: ScannerConfig, src: str) -> AST:
    """ Create AST for a file.

    :param conf_template: Base ScannerConfig object to be copied and modified.
    :param src: Sourcefile path
    :returns: The AST object.
    """
    ast = AST(conf_template, src)
    ast.create_ast_for_pickle()

    return ast
