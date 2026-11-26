#############################################################################
# Copyright (C) 2020-2026 MEmilio
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
:strong:`generator.py`
Generate the model specific python bindings code with the information given by the IntermediateRepresantation.
"""
from __future__ import annotations

import os
import string
from typing import TYPE_CHECKING
from importlib.resources import files

from typing_extensions import Self

from memilio.generation.template import template_string as StringTemplates

if TYPE_CHECKING:
    from memilio.generation.intermediate_representation import IntermediateRepresentation


class Generator:
    """ Generates the model specific python bindings code with the information given by the IntermediateRepresentation.
    """

    def __init__(self: Self) -> None:
        self.substitutions_py = {}
        self.substitutions_cpp = {}

    def create_substitutions(
            self: Self, intermed_representation: IntermediateRepresentation) -> None:
        """ Create the substitutions needed to generate the bindings.
        Divided into substitutions for the Python- and C++-file. 
        Uses the string template methods from the template folder.

        :param intermed_repr: Dataclass holding the model features.

        """
        intermed_repr = intermed_representation.found_model_specific_bindings

        if len(intermed_repr.model_base) > 0:
            model_base_templates = ", ".join(
                entry[0] for entry in intermed_repr.model_base if len(entry) > 0
            )
        else:
            raise IndexError("model_base is empty. No base classes found.")

        self.substitutions_cpp = {
            **StringTemplates.handle_all_model_specific_bindings(intermed_representation),
            "handle_all_bindings": StringTemplates.handle_all_bindings(intermed_representation)
        }
        self.substitutions_py = {
            "python_module_name": intermed_repr.python_module_name
        }

    def generate_files(
            self: Self, intermed_repr: IntermediateRepresentation, target_folder: str = "") -> None:
        """ Generate the python bindings to the C++ code.
        Template files for python and cpp from the template folder are used 
        and the identifiers substituted with the corresponding substitutions.

        :param intermed_repr: Dataclass holding the model features.
        :param target_folder: Output folder for generated files.
        """
        text = (files("memilio.generation") / "template" / "template_py.txt").read_text(encoding="utf-8")
        template_py = string.Template(text)

        text = (files("memilio.generation") / "template" / "template_cpp.txt").read_text(encoding="utf-8")
        template_cpp = string.Template(text)

        output_py = template_py.safe_substitute(**self.substitutions_py)
        output_cpp = template_cpp.safe_substitute(**self.substitutions_cpp)

        py_filename = intermed_repr.python_module_name + ".py"
        cpp_filename = intermed_repr.python_module_name + ".cpp"
        with open(os.path.join(target_folder, py_filename), "w") as output:
            output.write(output_py)
        with open(os.path.join(target_folder, cpp_filename), "w") as output:
            output.write(output_cpp)
