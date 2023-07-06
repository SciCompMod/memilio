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
@file template_string.py
@brief Generate small pieces of the target code as strings.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from memilio.generation import IntermediateRepresentation


def includes(intermed_repr: IntermediateRepresentation) -> str:
    """
    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    substitution_string = (
        "//Includes from pymio\n"
        "#include \"pybind_util.h\"\n"
        "#include \"utils/custom_index_array.h\"\n"
        "#include \"utils/parameter_set.h\"\n"
        "#include \"utils/index.h\"\n"
    )

    if "CompartmentalModel" in intermed_repr.model_base[0]:
        substitution_string += (
            "#include \"compartments/simulation.h\"\n"
            "#include \"compartments/compartmentalmodel.h\"\n"
            "#include \"epidemiology/populations.h\"\n"
        )

    if intermed_repr.simulation_class is not None:
        substitution_string += (
            "#include \"mobility/graph_simulation.h\"\n"
            "#include \"mobility/meta_mobility_instant.h\"\n"
        )

    substitution_string += "\n//Includes from Memilio\n"
    for inlcude in intermed_repr.include_list:
        substitution_string += "#include \"" + inlcude + "\"\n"

    substitution_string += "\n#include \"pybind11/pybind11.h\"\n"
    if intermed_repr.simulation_class is not None:
        substitution_string += (
            "#include \"pybind11/stl_bind.h\"\n"
            "#include \"Eigen/Core\"\n"
            "#include <vector>\n"
        )
    return substitution_string


def pretty_name_function(intermed_repr: IntermediateRepresentation) -> str:
    """
    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    substitution_string = (
        "namespace pymio"
        "{"
        "\n"
        "//specialization of pretty_name\n"
    )

    for key in intermed_repr.population_groups:
        if key == "AgeGroup":
            substitution_string += (
                "template <>\n"
                "std::string pretty_name<mio::{enum_class}>()\n"
                "{{\n"
                "   return \"{enum_class}\";\n"
                "}}\n"
                "\n"
            ).format(
                enum_class=key
            )
        else:
            substitution_string += (
                "template <>\n"
                "std::string pretty_name<{namespace}{enum_class}>()\n"
                "{{\n"
                "   return \"{enum_class}\";\n"
                "}}\n"
                "\n"
            ).format(
                namespace=intermed_repr.namespace,
                enum_class=key
            )

    return substitution_string + "} // namespace pymio\n"


def population_enums(intermed_repr: IntermediateRepresentation) -> str:
    """
    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    substitution_string = ""
    for key, values in intermed_repr.enum_populations.items():
        substitution_string += (
            "pymio::iterable_enum<{namespace}{enum_class}>(m, \"{enum_class}\")\n\t"
        ).format(
            namespace=intermed_repr.namespace,
            enum_class=key
        )

        for value in values[:-1]:
            substitution_string += (
                "    .value(\"{comp_class}\", {namespace}{enum_class}::{comp_class})\n\t"
            ).format(
                namespace=intermed_repr.namespace,
                comp_class=value,
                enum_class=key
            )
        substitution_string = substitution_string.rstrip() + ";\n\n"
    return substitution_string


def population(intermed_repr: IntermediateRepresentation) -> str:
    """
    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    for value in intermed_repr.model_base[1:]:
        if "Population" in value[0]:
            return value[0]


def model_init(intermed_repr: IntermediateRepresentation) -> str:
    """
    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    substitution_string = ""
    for init in intermed_repr.model_init:

        if len(init["type"]) > 1:
            continue
        elif len(init["type"]) == 0:
            substitution_string += "    .def(py::init<>())\n\t"
        else:
            substitution_string += (
                "    .def(py::init<{type}>(), py::arg(\"{name}\"))\n\t"
            ).format(
                type=init["type"][0],
                name=init["name"][0]
            )
    return substitution_string.rstrip() + ";\n"


def parameterset_indexing(intermed_repr: IntermediateRepresentation) -> str:
    """
    Generate the code for the AgeGroup class.
    Not used by every model.

    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    if not intermed_repr.parameterset_wrapper:
        return ""

    return (
        "pymio::bind_CustomIndexArray<mio::UncertainValue, mio::AgeGroup>(m, \"AgeGroupArray\");\n"
    )


def parameterset_wrapper(intermed_repr: IntermediateRepresentation) -> str:
    """
    Generate the code for the parameterset_wrapper needed when using age groups.
    Not used by every model.

    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    if not intermed_repr.parameterset_wrapper:
        return ""

    return (
        "py::class_<{namespace}{parameterset_wrapper}, {namespace}{parameterset}>(m, \"{parameterset_wrapper}\")\n"
        "\t.def(py::init<mio::AgeGroup>())\n"
        "\t.def(\"check_constraints\", &{namespace}{parameterset_wrapper}::check_constraints)\n"
        "\t.def(\"apply_constraints\", &{namespace}{parameterset_wrapper}::apply_constraints);\n"
    ).format(
        namespace=intermed_repr.namespace,
        parameterset=intermed_repr.parameterset,
        parameterset_wrapper=intermed_repr.parameterset_wrapper
    )


def age_group(intermed_repr: IntermediateRepresentation) -> str:
    """
    Generate the code for the AgeGroup class.
    Not used by every model.

    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    if not intermed_repr.age_group:
        return ""

    return (
        "py::class_<mio::AgeGroup, {base}>(m, \"AgeGroup\").def(py::init<{init}>());"
    ).format(
        namespace=intermed_repr.namespace,
        base=intermed_repr.age_group["base"],
        init=intermed_repr.age_group["init"][0]
    )


def simulation(intermed_repr: IntermediateRepresentation) -> str:
    """
    Generate the code for the Simulation class.
    Not used by every model.

    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    if intermed_repr.simulation_class is None or (
            not intermed_repr.simulation_class.strip()):
        return ""

    return (
        "pymio::bind_Simulation<{namespace}{simulation_class}<>>(m, \"{simulation_class}\");\n"
    ).format(
        namespace=intermed_repr.namespace,
        simulation_class=intermed_repr.simulation_class
    )


def simulation_graph(intermed_repr: IntermediateRepresentation) -> str:
    """
    Generate the code of the classes for graph simulations.
    Not used by every model.

    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    if intermed_repr.simulation_class is None or (
            not intermed_repr.simulation_class.strip()):
        return ""

    return (
        "pymio::bind_ModelNode<{namespace}{model_class}>(m, \"ModelNode\");\n\t"
        "pymio::bind_SimulationNode<{namespace}{simulation_class}<>>(m, \"SimulationNode\");\n\t"
        "pymio::bind_ModelGraph<{namespace}{model_class}>(m, \"ModelGraph\");\n\t"
        "pymio::bind_MigrationGraph<{namespace}{simulation_class}<>>(m, \"MigrationGraph\");\n\t"
        "pymio::bind_GraphSimulation<mio::Graph<mio::SimulationNode<{namespace}{simulation_class}<>>, mio::MigrationEdge>>(m, \"MigrationSimulation\");\n\t"
    ).format(
        namespace=intermed_repr.namespace,
        model_class=intermed_repr.model_class,
        simulation_class=intermed_repr.simulation_class
    )


def simulation_vector_definition(
        intermed_repr: IntermediateRepresentation) -> str:
    """
    Generate the code for vector definition.
    Not used by every model.

    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    if intermed_repr.simulation_class is None or (
            not intermed_repr.simulation_class.strip()):
        return ""

    return (
        "PYBIND11_MAKE_OPAQUE(std::vector<mio::Graph<mio::SimulationNode<{namespace}{simulation_class}<>>, mio::MigrationEdge>>);\n"
    ).format(
        namespace=intermed_repr.namespace,
        simulation_class=intermed_repr.simulation_class
    )
