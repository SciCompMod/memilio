#############################################################################
# Copyright (C) 2020-2024 MEmilio
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
from typing import TYPE_CHECKING, Any, Callable

from clang.cindex import *
from typing_extensions import Self

from memilio.generation import IntermediateRepresentation, utility


if TYPE_CHECKING:
    from memilio.generation import ScannerConfig


class Scanner:
    """! Analyze the model and extract the needed information.
    """

    def __init__(self: Self, conf: ScannerConfig) -> None:
        """
        Basic Constructor of Scanner class.

        @param conf ScannerConfig dataclass with the configurations.
        """
        self.config = conf
        utility.try_set_libclang_path(
            self.config.optional.get("libclang_library_path"))

    def extract_results(self: Self, root_cursor: Cursor) -> IntermediateRepresentation:
        """! Extract the information of the abstract syntax tree and save them in the dataclass intermed_repr.
        Call find_node to visit all nodes of abstract syntax tree and finalize to finish the extraction.

        @param root_cursor Represents the root node of the abstract syntax tree as a Cursor object from libclang.
        @return Information extracted from the model saved as an IntermediateRepresentation.
        """
        intermed_repr = IntermediateRepresentation()
        self.find_node(root_cursor, intermed_repr)
        self.finalize(intermed_repr)
        return intermed_repr

    def find_node(self: Self, node: Cursor,
                  intermed_repr: IntermediateRepresentation, namespace: str = "") -> None:
        """! Recursively walk over every node of an abstract syntax tree. Save the namespace the node is in.
        Call check_node_kind for extracting information from the nodes.

        @param node Represents the current node of the abstract syntax tree as a Cursor object from libclang.
        @param intermed_repr Dataclass used for saving the extracted model features.
        @param namespace [Default = ""] Namespace of the current node.
        """
        if node.spelling == "draw_sample":
            intermed_repr.has_draw_sample = True
        if node.kind == CursorKind.NAMESPACE:
            # Alternative
            # namespace += node.spelling + "::"
            namespace = (namespace + node.spelling + "::")
        elif namespace == self.config.namespace:
            self.switch_node_kind(node.kind)(node, intermed_repr)

        for n in node.get_children():
            self.find_node(n, intermed_repr, namespace)

    def switch_node_kind(self: Self, kind: CursorKind) -> Callable[[Any,
                                                                   IntermediateRepresentation],
                                                                   None]:
        """! Dictionary to map CursorKind to methods. Works like a switch.

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
            CursorKind.TYPE_ALIAS_DECL: self.check_type_alias,
            CursorKind.TYPE_ALIAS_TEMPLATE_DECL: self.check_type_alias


        }
        return switch.get(kind, lambda *args: None)

    def check_enum(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """! Inspect the nodes of kind ENUM_DECL and write needed information into intermed_repr.
        Information: Name of Enum

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        if node.spelling.strip() != "":  # alternative self.folder in node.location.file.name:
            intermed_repr.enum_populations[node.spelling] = []
            # print(intermed_repr.enum_populations)

    def check_enum_const(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """! Inspect the nodes of kind ENUM_CONSTANT_DECL and write needed information into intermed_repr.
        Information: Keys of an Enum

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        if node.semantic_parent.spelling in intermed_repr.enum_populations.keys():
            key = node.semantic_parent.spelling
            intermed_repr.enum_populations[key].append(node.spelling)
            # print(intermed_repr.enum_populations)
        # Alternative
        # enum_name = node.semantic_parent.spelling
        # if enum_name in intermed_repr.enum_populations:
        #     intermed_repr.enum_populations[enum_name].append(node.spelling)

    def check_class(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """! Inspect the nodes of kind CLASS_DECL and write information
        (model_class, model_base, simulation_class, parameterset_wrapper) into intermed_repr.

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        # print(node.spelling)
        if node.spelling == self.config.model_class:
            intermed_repr.model_class = node.spelling
            self.check_model_base(node, intermed_repr)
            self.check_model_includes(node, intermed_repr)
            self.check_age_group(node, intermed_repr)
        # elif (self.config.optional.get("simulation_class")
        #       and node.spelling == self.config.optional.get("simulation_class")):
        #     intermed_repr.simulation_class = node.spelling
        elif (node.spelling == "Simulation"):
            intermed_repr.simulation = True

        # Alternative
        # elif self.config.optional.get("simulation_class") == node.spelling:
        #     intermed_repr.simulation_class = node.spelling
        elif (intermed_repr.has_age_group and self.config.parameterset + "<FP>" in [base.spelling for base in node.get_children()]):
            intermed_repr.parameterset_wrapper = node.spelling

    def check_model_base(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """! Helper function to retreive the model base.

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """

        for base in node.get_children():

            if base.kind != CursorKind.CXX_BASE_SPECIFIER:
                continue

            base_name = base.spelling
            base_type = base.type

            if "FlowModel" in base_name:
                intermed_repr.is_flowmodel = True

            if "CompartmentalModel" in base_name and "mio" in node.semantic_parent.spelling:
                intermed_repr.is_compartmentalmodel = True

            intermed_repr.model_base.append(
                utility.get_base_class_string(base_type))

            self.check_model_base(base.referenced, intermed_repr)
            # print(f"Base type: {intermed_repr.model_base}")

    def check_base_specifier(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """! Not used yet.
        Inspect nodes which represent base specifier.
        For now this is handled by the parent node, which represents the class.
        """
        pass

    def check_model_includes(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """! Helper function to retrieve the model specific includes needed for pybind.

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
        """! Inspect the nodes of kind CLASS_DECL with the name defined in 
        config.age_group and write needed information into intermed_repr.
        Information: age_group

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        for base in node.get_children():
            if base.kind != CursorKind.CXX_BASE_SPECIFIER:
                continue
            for base_template_arg in base.get_children():
                if (base_template_arg.kind == CursorKind.TYPE_REF and "AgeGroup" in base_template_arg.spelling):
                    intermed_repr.has_age_group = True
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
        """! Inspect the nodes of kind CONSTRUCTOR and write needed information into intermed_repr.
        Information: intermed_repr.init

        @param node Current node represented as a Cursor object.
        @param intermed_repr Dataclass used for saving the extracted model features.
        """

        # wieder alter version, template klammern rauswerfen und dann vergleichen
        if intermed_repr.model_class == "Model":
            if node.spelling.startswith("Model") and ('<' in node.spelling and '>' in node.spelling):
                init = {"type": [], "name": []}
                for arg in node.get_arguments():
                    tokens = []
                    for token in arg.get_tokens():
                        tokens.append(token.spelling)
                    # init["type"].append(arg.type.spelling)
                    init["type"].append(" ".join(tokens[:-1]))
                    init["name"].append(tokens[-1])
                intermed_repr.model_init.append(init)

    def check_type_alias(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """! Inspect the nodes of kind TYPE_ALIAS_DECL and write needed information into intermed_repr.
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
        """! Finalize the IntermediateRepresenation as last step of the Scanner.
        Write needed information from config into intermed_repr,
        delet unnecesary enums and check for missing model features.

        @param intermed_repr Dataclass used for saving the extracted model features.
        """
        # remove unnecesary enum
        population_groups = []
        for value in intermed_repr.model_base[0:]:
            if "Populations" in value[0]:
                population_groups = [pop[0].split(
                    "::")[-1] for pop in value[1:]]
        intermed_repr.population_groups = population_groups
        new_enum = {}
        for key in intermed_repr.enum_populations:

            if key in population_groups:
                new_enum[key] = intermed_repr.enum_populations[key]
        intermed_repr.enum_populations = new_enum  # ist leer
        # pass information from config
        intermed_repr.set_attribute("namespace", self.config.namespace)
        intermed_repr.set_attribute(
            "python_module_name", self.config.python_module_name)
        intermed_repr.set_attribute("target_folder", self.config.target_folder)
        intermed_repr.set_attribute(
            "python_generation_module_path", self.config.python_generation_module_path)

        # check for missing data
        intermed_repr.check_complete_data(self.config.optional)
