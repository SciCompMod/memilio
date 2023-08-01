#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
@file scanner.py
@brief Analyze the model and extract the needed information. Information get passed to the IntermediateRepresenation.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from typing import TYPE_CHECKING, Any, Callable
from warnings import catch_warnings

from clang.cindex import *
from typing_extensions import Self

from memilio.generation import IntermediateRepresentation, utility

if TYPE_CHECKING:
    from memilio.generation import ScannerConfig


class Scanner:
    """
    Analyze the model and extract the needed information.
    """

    def __init__(self: Self, conf: ScannerConfig) -> None:
        """
        Basic Constructor of Scanner class.

        @param conf ScannerConfig dataclass with the configurations.
        """
        self.config = conf
        utility.try_set_libclang_path(
            self.config.optional.get("libclang_library_path"))
        self.ast = None
        self.create_ast()

    def create_ast(self: Self) -> None:
        """
        Create an abstract syntax tree for the main model.cpp file with a corresponding CompilationDatabase. 
        A compile_commands.json is required (automatically generated in the build process).
        """
        idx = Index.create()

        # Create the cmd arguments
        file_args = []

        dirname = utility.try_get_compilation_database_path(
            self.config.skbuild_path_to_database)
        compdb = CompilationDatabase.fromDirectory(dirname)
        commands = compdb.getCompileCommands(self.config.source_file)
        for command in commands:
            for argument in command.arguments:
                if (argument != '-Wno-unknown-warning' and
                        argument != "--driver-mode=g++" and argument != "-O3"):
                    file_args.append(argument)
        file_args = file_args[1:-4]
        clang_cmd = [
            "clang-14", self.config.source_file,
            "-std=c++17", '-emit-ast', '-o', '-']
        clang_cmd.extend(file_args)

        clang_cmd_result = subprocess.run(clang_cmd, stdout=subprocess.PIPE)
        clang_cmd_result.check_returncode()

        # Since `clang.Index.read` expects a file path, write generated abstract syntax tree to a
        # temporary named file. This file will be automatically deleted when closed.
        with tempfile.NamedTemporaryFile() as ast_file:
            ast_file.write(clang_cmd_result.stdout)
            self.ast = idx.read(ast_file.name)

    def extract_results(self: Self) -> IntermediateRepresentation:
        """
        Extract the information of the abstract syntax tree and save them in the dataclass intermed_repr.
        Call find_node to visit all nodes of abstract syntax tree and finalize to finish the extraction.

        @return Information extracted from the model saved as an IntermediateRepresentation. 
        """
        intermed_repr = IntermediateRepresentation()
        utility.output_cursor_print(self.ast.cursor, 1)
        self.find_node(self.ast.cursor, intermed_repr)
        self.finalize(intermed_repr)
        return intermed_repr

    def find_node(self: Self, node: Cursor,
                  intermed_repr: IntermediateRepresentation, namespace: str = "") -> None:
        """
        Recursively walk over every node of an abstract syntax tree. Save the namespace the node is in.
        Call check_node_kind for extracting information from the nodes.

        @param node Represents the current node of the abstract syntax tree as a Cursor object from libclang.
        @param intermed_repr Dataclass used for saving the extracted model features.
        @param namespace [Default = ""] Namespace of the current node.
        """
        if node.kind == CursorKind.NAMESPACE:
            namespace = (namespace + node.spelling + "::")
        elif namespace == self.config.namespace:
            self.switch_node_kind(node.kind)(node, intermed_repr)

        for n in node.get_children():
            self.find_node(n, intermed_repr, namespace)

    def switch_node_kind(self: Self, kind: CursorKind) -> Callable[[Any,
                                                                   IntermediateRepresentation],
                                                                   None]:
        """
        Dictionary to map CursorKind to methods. Works like a switch.

        @param Underlying kind of the current node.
        @return Appropriate method for the given kind.
        """
        switch = {
            CursorKind.ENUM_DECL: self.check_enum,
            CursorKind.ENUM_CONSTANT_DECL: self.check_enum_const,
            CursorKind.CLASS_DECL: self.check_class,
            CursorKind.CLASS_TEMPLATE: self.check_class,
            CursorKind.CXX_BASE_SPECIFIER: self.check_base_specifier,
            CursorKind.CONSTRUCTOR: self.check_constructor,
            CursorKind.STRUCT_DECL: self.check_struct,
            CursorKind.TYPE_ALIAS_DECL: self.check_type_alias
        }
        return switch.get(kind, lambda *args: None)

    def check_enum(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """
        Inspect the nodes of kind ENUM_DECL and write needed information into intermed_repr.
        Information: Name of Enum

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        if node.spelling.strip() != "":  # alternative self.folder in node.location.file.name:
            intermed_repr.enum_populations[node.spelling] = []

    def check_enum_const(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """
        Inspect the nodes of kind ENUM_CONSTANT_DECL and write needed information into intermed_repr.
        Information: Keys of an Enum

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        if node.semantic_parent.spelling in intermed_repr.enum_populations.keys():
            key = node.semantic_parent.spelling
            intermed_repr.enum_populations[key].append(node.spelling)

    def check_class(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """
        Inspect the nodes of kind CLASS_DECL and write information 
        (model_class, model_base, simulation_class, parameterset_wrapper) into intermed_repr.

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        if node.spelling == self.config.model_class:
            intermed_repr.model_class = node.spelling
            self.check_model_base(node, intermed_repr)
            self.check_model_includes(node, intermed_repr)
            if self.config.optional.get("age_group"):
                self.check_age_group(node, intermed_repr)
        elif (self.config.optional.get("simulation_class")
              and node.spelling == self.config.optional.get("simulation_class")):
            intermed_repr.simulation_class = node.spelling
        elif (self.config.optional.get("parameterset_wrapper") and self.config.namespace
                + self.config.parameterset in [base.spelling for base in node.get_children()]):
            intermed_repr.parameterset_wrapper = node.spelling

    def check_model_base(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """
        Helper function to retreive the model base.

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        for base in node.get_children():
            if base.kind != CursorKind.CXX_BASE_SPECIFIER:
                continue
            base_type = base.get_definition().type
            intermed_repr.model_base = utility.get_base_class_string(base_type)

    def check_base_specifier(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """ 
        Not used yet.
        Inspect nodes which represent base specifier.
        For now this is handled by the parent node, which represents the class.
        """
        pass

    def check_model_includes(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """
        Helper function to retrieve the model specific includes needed for pybind.

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        filepath = node.location.file.name
        filepaths = filepath.split("../")
        model_has_analyze_results = False

        intermed_repr.include_list.append(filepaths[1])
        intermed_repr.include_list.append(
            filepaths[1].replace("model.h", "") + "infection_state.h")

        # Iterate throught files in the directory of the model and check for files
        for file in os.listdir(filepaths[0]):
            if file == "parameter_space.h":
                intermed_repr.include_list.append(filepaths[1].replace(
                    "model.h", "") +
                    "parameter_space.h")
            elif file == "analyze_result.h":
                model_has_analyze_results = True

        if model_has_analyze_results:
            intermed_repr.include_list.append(
                filepaths[1].replace("model.h", "") + "analyze_result.h")
        else:
            intermed_repr.include_list.append("memilio/data/analyze_result.h")

    def check_age_group(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """
        Inspect the nodes of kind CLASS_DECL with the name defined in 
        config.age_group and write needed information into intermed_repr.
        Information: age_group

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        for base in node.get_children():
            if base.kind != CursorKind.CXX_BASE_SPECIFIER:
                continue
            for base_template_arg in base.get_children():
                if (base_template_arg.kind == CursorKind.TYPE_REF
                        and "AgeGroup" in base_template_arg.spelling):
                    for child in base_template_arg.get_definition().get_children():
                        if child.kind == CursorKind.CXX_BASE_SPECIFIER:
                            intermed_repr.age_group["base"] = child.get_definition(
                            ).type.spelling
                        elif child.kind == CursorKind.CONSTRUCTOR:
                            intermed_repr.age_group["init"] = [
                                arg.spelling
                                for arg in child.type.argument_types()]

    def check_constructor(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """
        Inspect the nodes of kind CONSTRUCTOR and write needed information into intermed_repr.
        Information: intermed_repr.init

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        if node.spelling == intermed_repr.model_class:
            init = {"type": [], "name": []}
            for arg in node.get_arguments():
                tokens = []
                for token in arg.get_tokens():
                    tokens.append(token.spelling)
                init["type"].append(" ".join(tokens[:-1]))
                init["name"].append(tokens[-1])
            intermed_repr.model_init.append(init)

    def check_type_alias(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """
        Inspect the nodes of kind TYPE_ALIAS_DECL and write needed information into intermed_repr.
        Information: intermed_repr.parameterset

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        if node.spelling == self.config.parameterset:
            intermed_repr.parameterset = node.spelling

    def check_struct(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """ Not used yet."""
        pass

    def finalize(self: Self, intermed_repr: IntermediateRepresentation) -> None:
        """
        Finalize the IntermediateRepresenation as last step of the Scanner.
        Write needed information from config into intermed_repr,
        delet unnecesary enums and check for missing model features.

        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        # remove unnecesary enum
        population_groups = []
        for value in intermed_repr.model_base[1:]:
            if "Population" in value[0]:
                population_groups = [pop[0].split(
                    "::")[-1] for pop in value[1:]]
        intermed_repr.population_groups = population_groups
        new_enum = {}
        for key in intermed_repr.enum_populations:
            if key in population_groups:
                new_enum[key] = intermed_repr.enum_populations[key]
        intermed_repr.enum_populations = new_enum

        # pass information from config
        intermed_repr.set_attribute("namespace", self.config.namespace)
        intermed_repr.set_attribute(
            "python_module_name", self.config.python_module_name)
        intermed_repr.set_attribute("target_folder", self.config.target_folder)
        intermed_repr.set_attribute(
            "python_generation_module_path", self.config.python_generation_module_path)

        # check for missing data
        intermed_repr.check_complete_data(self.config.optional)

    def output_ast(self: Self) -> None:
        """
        Output the abstract syntax tree to terminal.
        """
        utility.output_cursor_and_children(self.ast.cursor)

    def output_ast_file(self: Self) -> None:
        """
        Output the abstract syntax tree to file.
        """
        with open('output_ast.txt', 'a') as f:
            utility.output_cursor_and_children_file(self.ast.cursor, f)
            print('AST written to ' + str(os.path.abspath(f.name)))
