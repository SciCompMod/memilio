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
@file generator.py
@brief Generate the model specific python bindings code with the information given by the IntermediateRepresantation.
"""
from __future__ import annotations

import os
import string
from typing import TYPE_CHECKING

from typing_extensions import Self

from memilio.generation.template import template_string as StringTemplates

if TYPE_CHECKING:
    from memilio.generation import IntermediateRepresentation


class Generator:
    """! Generates the model specific python bindings code with the information given by the IntermediateRepresentation.
    """

    def __init__(self: Self) -> None:
        self.substitutions_py = {}
        self.substitutions_cpp = {}

    def create_substitutions(
            self: Self, intermed_repr: IntermediateRepresentation) -> None:
        """! Create the substitutions needed to generate the bindings.
        Divided into substitutions for the Python- and C++-file. 
        Uses the string template methods from the template folder.

        @param intermed_repr: Dataclass holding the model features.
        """

        self.substitutions_py = {
            "python_module_name": intermed_repr.python_module_name
        }

        if len(intermed_repr.model_base) > 0:
            model_base_templates = ", ".join(
                entry[0] for entry in intermed_repr.model_base if len(entry) > 0
            )
        else:
            raise IndexError("model_base is empty. No base classes found.")

        self.substitutions_cpp = {
            "namespace": intermed_repr.namespace,
            "model_class_name": intermed_repr.model_class,
            "model_base": intermed_repr.model_base[0],
            "model_base_templates": model_base_templates,
            "python_module_name": intermed_repr.python_module_name,
            "parameterset": intermed_repr.parameterset,
            "includes": StringTemplates.includes(intermed_repr),
            "pretty_name_function": StringTemplates.pretty_name_function(intermed_repr),
            "population_enums": StringTemplates.population_enums(intermed_repr),
            "model_init": StringTemplates.model_init(intermed_repr),
            "population": StringTemplates.population(intermed_repr),
            "parameterset_indexing": StringTemplates.parameterset_indexing(intermed_repr),
            "parameterset_wrapper": StringTemplates.parameterset_wrapper(intermed_repr),
            "simulation": StringTemplates.simulation(intermed_repr),
            "ScalarType": StringTemplates.ScalarType(intermed_repr),
            "draw_sample": StringTemplates.draw_sample(intermed_repr),
            "simulation_vector_definition": StringTemplates.simulation_vector_definition(intermed_repr)
        }

    def generate_files(
            self: Self, intermed_repr: IntermediateRepresentation) -> None:
        """! Generate the python bindings to the C++ code.
        Template files for python and cpp from the template folder are used 
        and the identifiers substituted with the corresponding substitutions.

        @param intermed_repr Dataclass holding the model features.
        """
        with open(os.path.join(intermed_repr.python_generation_module_path,
                               "memilio/generation/template/template_py.txt")) as t:
            template_py = string.Template(t.read())
        with open(os.path.join(intermed_repr.python_generation_module_path,
                               "memilio/generation/template/template_cpp.txt")) as t:
            template_cpp = string.Template(t.read())

        output_py = template_py.safe_substitute(**self.substitutions_py)
        output_cpp = template_cpp.safe_substitute(**self.substitutions_cpp)

        py_filename = intermed_repr.python_module_name + ".py"
        cpp_filename = intermed_repr.python_module_name + ".cpp"
        with open(os.path.join(intermed_repr.target_folder, py_filename), "w") as output:
            output.write(output_py)
        with open(os.path.join(intermed_repr.target_folder, cpp_filename), "w") as output:
            output.write(output_cpp)
