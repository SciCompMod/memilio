#############################################################################
# Copyright (C) 2020-2025 MEmilio
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
:strong:`scanner.py`
Analyze the model and extract the needed information. Information get passed to the IntermediateRepresenation.
"""
from __future__ import annotations

import os
from typing import TYPE_CHECKING, Any, Callable

from clang.cindex import *
from typing_extensions import Self

from memilio.generation import IntermediateRepresentation, utility
from memilio.generation.default_generation_dict import default_dict, general_bindings_dict


if TYPE_CHECKING:
    from memilio.generation import ScannerConfig


class Scanner:
    """ Analyze the model and extract the needed information."""

    def __init__(self: Self, conf: ScannerConfig) -> None:
        """
        Basic Constructor of Scanner class.

        :param conf: ScannerConfig dataclass with the configurations.
        """
        self.config = conf
        utility.try_set_libclang_path(
            self.config.libclang_library_path)
        source_file = self.config.source_file
        from_folder = os.path.basename(os.path.dirname(source_file))
        self.python_module_name = "o" + from_folder.split("_")[1]
        self.namespace = default_dict["mio"] + \
            "::" + self.python_module_name + "::"
        self.handled_bindings = set()

    def extract_results(self: Self, root_cursor: Cursor) -> IntermediateRepresentation:
        """ Extract the information of the abstract syntax tree and save them in the dataclass intermed_repr.
        Call find_node to visit all nodes of abstract syntax tree and finalize to finish the extraction.

        :param root_cursor: Represents the root node of the abstract syntax tree as a Cursor object from libclang.
        :param self: Self:
        :returns: Information extracted from the model saved as an IntermediateRepresentation.

        """
        if self.config.model_class != default_dict["model"]:
            raise AssertionError("set a model name")

        intermed_repr = IntermediateRepresentation()

        self.find_node(root_cursor, intermed_repr, general_bindings_dict)
        self.finalize(intermed_repr)
        self.check_parameter_space(intermed_repr)
        return intermed_repr

    def search_binding_target(self: Self, binding_list: list[dict], node: Cursor,
                              intermed_repr: IntermediateRepresentation, namespace: str) -> list:
        """ Search for a binding target in the current node.
        If the node matches the binding, set the information in intermed_repr.
        :param binding_list: List of dictionaries containing the bindings.
        :param node: Represents the current node of the abstract syntax tree as a Cursor object from libclang.
        :param intermed_repr: Dataclass used for saving the extracted model features.
        :param namespace: Namespace of the current node."""

        for binding in binding_list:
            kind_name = binding.get("cursorkind")  # name of the CursorKind
            # name of the function or class that is going to be bound
            binding_name = binding.get("name")
            # type of the binding, e.g. "class", "function"
            binding_type = binding.get("type")

            if not kind_name or not binding_name:
                continue

            try:
                expected_kind = getattr(CursorKind, kind_name)
            except AttributeError:
                continue

            function_arguments = self.get_function_arguments(
                node)

            key = (binding_name, kind_name, tuple(function_arguments.get(
                "arg_types", [])), tuple(function_arguments.get("arg_names", [])))

            if key in self.handled_bindings:
                continue

            if node.kind == expected_kind and node.spelling == binding_name:
                self.set_info(intermed_repr, node, binding_type,
                              namespace, function_arguments)
                self.handled_bindings.add(key)

        return intermed_repr.found_bindings

    def set_info(self: Self, intermed_repr: IntermediateRepresentation, node: Cursor, binding_type, namespace: str, function_arguments: dict) -> None:
        """ Set the information of the current node into the intermed_repr.

        :param intermed_repr: Dataclass used for saving the extracted model features.
        :param node: Represents the current node of the abstract syntax tree as a Cursor object from libclang.
        :param arg_types: List of argument types of the current node.
        :param arg_names: List of argument names of the current node.
        :param binding_type: Type of the binding, e.g. "class", "function", "extern_function".
        :param namespace: Namespace of the current node.
        """

        found_info = {
            "type": binding_type,
            "name": node.spelling,
            "kind": node.kind.name,
            "namespace": namespace,
            "return_type": node.result_type.spelling if node.kind.is_declaration() else "",
            "arg_types": function_arguments.get("arg_types", []),
            "arg_names": function_arguments.get("arg_names", []),
            "parent_name": function_arguments.get("parent_name", ""),
            "is_const": function_arguments.get("is_const", False),
            "is_member": function_arguments.get("is_member", False)
        }

        intermed_repr.found_bindings.append(found_info)

    def get_function_arguments(self, node: Cursor) -> dict:
        """ Get the argument types and names of a function node.
        :param node: Represents the current node of the abstract syntax tree as a Cursor object from libclang.
        :returns: A tuple containing a list of argument types and a list of argument names."""
        arg_types = []
        arg_names = []

        for child in node.get_children():
            if child.kind == CursorKind.PARM_DECL:

                arg_types.append(child.type.spelling)

                arg_names.append(child.spelling)

        is_const = node.is_const_method() if hasattr(node, 'is_const_method') else False
        is_member, parent_name = self.is_member(node)

        return {
            "arg_types": arg_types,
            "arg_names": arg_names,
            "parent_name": parent_name,
            "is_const": is_const,
            "is_member": is_member
        }

    def is_member(self, node: Cursor) -> tuple[bool, str]:
        """Check if the node is a member function (also for function templates).
        :param node: Represents the current node of the abstract syntax tree as a Cursor object from libclang.
        :returns: True if the node is a member function, False otherwise.
        """

        if node.kind not in (
            CursorKind.CXX_METHOD,
            CursorKind.FUNCTION_TEMPLATE,
            CursorKind.CONSTRUCTOR,
            CursorKind.DESTRUCTOR,
        ):
            return False, ""

        parent = node.semantic_parent

        if parent is None:
            return False, ""

        if parent.kind in (
            CursorKind.CLASS_DECL,
            CursorKind.STRUCT_DECL,
            CursorKind.CLASS_TEMPLATE,
        ):
            parent_name = parent.spelling
            return True, parent_name

        return False, ""

    def check_parameter_space(self: Self, intermed_repr: IntermediateRepresentation) -> None:
        """! Checks for parameter_space.cpp in the model folder and set has_draw_sample

        @param intermed_repr: Dataclass used for saving the extracted model features.
        """
        source_file = self.config.source_file
        model_folder = os.path.dirname(source_file)
        parameter_space_file = os.path.join(
            model_folder, default_dict["parameterspacefile"])

        if (os.path.isfile(parameter_space_file)):
            intermed_repr.has_draw_sample = True

    def find_node(self: Self, node: Cursor,
                  intermed_repr: IntermediateRepresentation, binding_list: list[dict], namespace: str = "") -> None:
        """ Recursively walk over every node of an abstract syntax tree. Save the namespace the node is in.
        Call check_node_kind for extracting information from the nodes.

        :param node: Represents the current node of the abstract syntax tree as a Cursor object from libclang.
        :param intermed_repr: Dataclass used for saving the extracted model features.
        :param namespace: Default = ""] Namespace of the current node.
        :param self: Self: 

        """
        if node.kind == CursorKind.NAMESPACE:
            namespace = (namespace + node.spelling + "::")

        elif namespace == self.namespace:
            self.search_binding_target(
                binding_list, node, intermed_repr, namespace)
            self.switch_node_kind(node.kind)(node, intermed_repr)

        for n in node.get_children():
            self.find_node(n, intermed_repr, binding_list, namespace)

    def switch_node_kind(self: Self, kind: CursorKind) -> Callable[[Any,
                                                                   IntermediateRepresentation],
                                                                   None]:
        """ Dictionary to map CursorKind to methods. Works like a switch.

        :param Underlying: kind of the current node.
        :param self: Self: 
        :param kind: CursorKind: 
        :returns: Appropriate method for the given kind.

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
        """ Inspect the nodes of kind ENUM_DECL and write needed information into intermed_repr.
        Information: Name of Enum

        :param node: Current node represented as a Cursor object.
        :param intermed_repr: Dataclass used for saving the extracted model features.
        :param self: Self: 

        """
        if node.spelling.strip() != default_dict["emptystring"]:
            intermed_repr.enum_populations[node.spelling] = []

    def check_enum_const(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """ Inspect the nodes of kind ENUM_CONSTANT_DECL and write needed information into intermed_repr.
        Information: Keys of an Enum

        :param node: Current node represented as a Cursor object.
        :param intermed_repr: Dataclass used for saving the extracted model features.

        """
        if node.semantic_parent.spelling in intermed_repr.enum_populations.keys():
            key = node.semantic_parent.spelling
            intermed_repr.enum_populations[key].append(node.spelling)

    def check_class(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """! Inspect the nodes of kind CLASS_DECL and write information
        (model_class, model_base, simulation_class, parameterset_wrapper) into intermed_repr.

        :param node: Current node represented as a Cursor object.
        :param intermed_repr: Dataclass used for saving the extracted model features.

        """
        if node.spelling == self.config.model_class:
            intermed_repr.model_class = node.spelling
            self.check_model_base(node, intermed_repr)
            self.check_model_includes(node, intermed_repr)
            self.check_age_group(node, intermed_repr)

        elif (node.spelling == default_dict["simulation"]):
            intermed_repr.simulation = True

        elif (intermed_repr.has_age_group and self.config.parameterset + "<FP>" in [base.spelling for base in node.get_children()]):
            intermed_repr.parameterset_wrapper = node.spelling

    def check_model_base(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """ Helper function to retreive the model base.

        :param node: Current node represented as a Cursor object.
        :param intermed_repr: Dataclass used for saving the extracted model features.

        """

        for base in node.get_children():

            if base.kind != CursorKind.CXX_BASE_SPECIFIER:
                continue

            base_name = base.spelling
            base_type = base.type

            if default_dict["flowmodel"] in base_name:
                intermed_repr.is_flowmodel = True

            if default_dict["compartmentalmodel"] in base_name and default_dict["mio"] in node.semantic_parent.spelling:
                intermed_repr.is_compartmentalmodel = True

            intermed_repr.model_base.append(
                utility.get_base_class_string(base_type))

            self.check_model_base(base.referenced, intermed_repr)

    def check_base_specifier(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """ Not used yet.
        Inspect nodes which represent base specifier.
        For now this is handled by the parent node, which represents the class.

        :param self: Self: 
        :param node: Cursor: 
        :param intermed_repr: IntermediateRepresentation: 

        """
        pass

    def check_model_includes(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """ Helper function to retrieve the model specific includes needed for pybind.

        :param node: Current node represented as a Cursor object.
        :param intermed_repr: Dataclass used for saving the extracted model features.
        :param self: Self: 

        """
        filepath = node.location.file.name
        filepaths = filepath.split("../")

        model_has_analyze_results = False

        intermed_repr.include_list.append(filepaths[1])

        intermed_repr.include_list.append(
            filepaths[1].replace(default_dict["modelfile"], default_dict["emptystring"]) + default_dict["infectionstatefile"])

        for file in os.listdir(filepaths[0]):
            if file == default_dict["parameterspacefile"]:
                intermed_repr.include_list.append(filepaths[1].replace(
                    default_dict["modelfile"], default_dict["emptystring"]) +
                    default_dict["parameterspacefile"])
            elif file == default_dict["analyzeresultfile"]:
                model_has_analyze_results = True

        if model_has_analyze_results:
            intermed_repr.include_list.append(
                filepaths[1].replace(default_dict["modelfile"], default_dict["emptystring"]) + default_dict["analyzeresultfile"])
        else:
            intermed_repr.include_list.append("memilio/data/analyze_result.h")

    def check_age_group(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """! Inspect the nodes of kind CLASS_DECL with the name defined in
        config.age_group and write needed information into intermed_repr.
        Information: age_group

        :param node: Current node represented as a Cursor object.
        :param intermed_repr: Dataclass used for saving the extracted model features.
        :param self: Self: 
        """
        for base in node.get_children():
            if base.kind != CursorKind.CXX_BASE_SPECIFIER:
                continue
            for base_template_arg in base.get_children():
                if (base_template_arg.kind == CursorKind.TYPE_REF and default_dict["agegroup"] in base_template_arg.spelling):
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
        """ Inspect the nodes of kind CONSTRUCTOR and write needed information into intermed_repr.
        Information: intermed_repr.init

        :param node: Current node represented as a Cursor object.
        :param intermed_repr: Dataclass used for saving the extracted model features.
        :param self: Self: 

        """
        if intermed_repr.model_class == default_dict["model"]:
            if node.spelling.startswith(default_dict["model"]) and ('<' in node.spelling and '>' in node.spelling):
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
        """ Inspect the nodes of kind TYPE_ALIAS_DECL and write needed information into intermed_repr.
        Information: intermed_repr.parameterset

        :param node: Current node represented as a Cursor object.
        :param intermed_repr: Dataclass used for saving the extracted model features.
        :param self: Self: 

        """
        if node.spelling == self.config.parameterset:
            intermed_repr.parameterset = node.spelling

    def check_struct(
        self: Self, node: Cursor,
            intermed_repr: IntermediateRepresentation) -> None:
        """Not used yet.

        :param self: Self: 
        :param node: Cursor: 
        :param intermed_repr: IntermediateRepresentation: 

        """
        pass

    def finalize(self: Self, intermed_repr: IntermediateRepresentation) -> None:
        """ Finalize the IntermediateRepresenation as last step of the Scanner.
        Write needed information from config into intermed_repr,
        delet unnecesary enums and check for missing model features.

        :param intermed_repr: Dataclass used for saving the extracted model features.
        :param self: Self: 
        """

        population_groups = []
        for value in intermed_repr.model_base[0:]:
            if default_dict["flowmodel"] in value[0].strip():
                start = value[0].find("Populations<")
                end = value[0].find(">", start)

                if start != -1 and end != -1:
                    populations_part = value[0][start + 11:end]
                    population_groups = [
                        part.strip(" <>").split("::")[-1]
                        for part in populations_part.split(",")
                    ]

        intermed_repr.population_groups = population_groups

        new_enum = {}
        for key in intermed_repr.enum_populations:

            if key in population_groups:

                new_enum[key] = intermed_repr.enum_populations[key]

                intermed_repr.enum_populations = new_enum

        intermed_repr.set_attribute(
            "namespace", self.namespace)
        intermed_repr.set_attribute(
            "python_module_name", self.python_module_name)
        intermed_repr.set_attribute("target_folder", self.config.target_folder)
        intermed_repr.set_attribute(
            "python_generation_module_path", self.config.python_generation_module_path)
