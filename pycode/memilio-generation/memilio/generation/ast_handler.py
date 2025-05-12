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

    # Todo: clean __init__()
    def __init__(self: T, conf: ScannerConfig) -> None:
        """ Basic constructor of the ASTHandler class

        :param conf: ScannerConfig dataclass with the configurations.
        """
        self.conf = conf
        self.ast_list = []
        self.source_files = []
        self.source_files.append(self.conf.source_file)
        self.handle_ast_creation()

    def handle_ast_creation(self: T) -> None:
        """ Handle the creation of the ASTs. Single or Parallel.
        """
        if len(self.source_files) == 1:
            self.single_creation(self.conf)
        else:
            self.parallel_creation(self.conf)

    def add_source_file(self: T, source_file_str: str) -> None:
        """ Add source file to list. No need for the full filepath. 
        Just filename with ending.

        :param source_file_str: String that represents another file.
        """
        source_file = self.conf.source_file
        model_folder = os.path.dirname(source_file)
        added_file = os.path.join(
            model_folder, source_file_str)

        if (os.path.isfile(added_file)):
            self.source_files.append(added_file)
            logging.info(f"Adding source file: {added_file}")

    @staticmethod
    def process_ast_file(ast_path: str, ast: AST) -> TranslationUnit:
        """ Process AST file and create the translation_unit.

        :param ast_path: Represents the path to the generated AST file created by the clang process.
        :param ast: AST instance.
        :returns: TranslationUnit object.
        """
        idx = Index.create()
        tu = idx.read(ast_path)
        ast._assing_ast_with_ids(tu.cursor)
        ast.set_translation_unit(tu)
        os.remove(ast_path)
        return tu

    def parallel_creation(self: T, conf: ScannerConfig) -> None:
        """ Create ASTs for all source files using multiprocessing.

        This method distributes the AST creation process across multiple CPU cores. 
        Each source file is processed independently using the `create_ast` function, 
        and the resulting ASTs are stored in `self.ast_list`.

        :param conf: ScannerConfig dataclass with the configurations.
        """
        with multiprocessing.Pool() as pool:
            try:

                results = pool.starmap(
                    create_ast, [(conf, source_file) for source_file in self.source_files])

                for ast_path, ast in results:
                    try:
                        self.process_ast_file(
                            ast_path, ast)

                        self.add_ast_to_list(ast)

                    except Exception as e:
                        logging.info(
                            f"Error processing AST file {ast_path}: {e}")
            except Exception as e:
                logging.info(f"Error in parallel AST creation: {e}")

    def single_creation(self: T, conf: ScannerConfig) -> TranslationUnit:
        """ Create AST for a single source file.

        :param conf: ScannerConfig dataclass with the configurations.
        :returns: TranslationUnit object.
        """
        ast_path, ast = create_ast(conf, self.source_files[0])
        translation_unit = self.process_ast_file(ast_path, ast)
        self.add_ast_to_list(ast)
        return translation_unit

    def add_ast_to_list(self: T, ast: AST) -> None:
        """ Add AST to list.

        :param ast: AST instance.
        """
        self.ast_list.append(ast)

    def remove_ast_from_list(self: T, ast: AST) -> None:
        """ Remove AST from list.

        :param ast: AST instance.
        """
        self.ast_list.remove(ast)

    def get_ast_list(self: T) -> list:
        """ Get list of ASTs.

        :returns: List of AST objects.
        """
        return self.ast_list

    def get_ast_by_id(self: T, id: int) -> AST:
        """ Get AST by id.

        :param id: ID of the AST.
        :returns: AST object.
        """
        return self.ast_list[id]


def create_ast(conf_template: ScannerConfig, src: str) -> tuple[str, AST]:
    """ Create AST for a file.

    :param conf_template: Base ScannerConfig object to be copied and modified.
    :param src: Sourcefile path
    :returns: A tuple containing the path to the temporary AST file and the AST object.
    """
    conf = copy.deepcopy(conf_template)
    conf.source_file = src
    ast = AST(conf)
    ast_path = ast.create_ast_for_pickle()
    return ast_path, ast
