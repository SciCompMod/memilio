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
@file template_string.py
@brief Generate small pieces of the target code as strings.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from memilio.generation import IntermediateRepresentation


def ScalarType(intermed_repr: IntermediateRepresentation) -> str:
    """! Set the datatype for the bindings via. intermediate_representation.
    @return string from intermediate_representation
    """
    scalartypestr = intermed_repr.scalartype
    return scalartypestr


def includes(intermed_repr: IntermediateRepresentation) -> str:
    """! Fills in the Includes for the binding
    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    substitution_string = (
        "//Includes from pymio\n"
        "#include \"pybind_util.h\"\n"
    )

    if intermed_repr.is_flowmodel:
        substitution_string += (
            "#include \"compartments/flow_simulation.h\"\n"
        )

    if intermed_repr.is_compartmentalmodel:
        substitution_string += (
            "#include \"compartments/simulation.h\"\n"
            "#include \"compartments/compartmentalmodel.h\"\n"
            "#include \"epidemiology/populations.h\"\n"
        )

    if intermed_repr.has_age_group:
        substitution_string += (
            "#include \"epidemiology/age_group.h\"\n"
            "#include \"utils/parameter_set.h\"\n"
        )

    substitution_string += (
        "#include \"utils/custom_index_array.h\"\n"
        "#include \"utils/index.h\"\n"
        "#include \"mobility/graph_simulation.h\"\n"
        "#include \"mobility/metapopulation_mobility_instant.h\"\n"
        "#include \"io/mobility_io.h\"\n"
        "#include \"io/result_io.h\"\n"
    )

    substitution_string += "\n//Includes from MEmilio\n"
    for inlcude in intermed_repr.include_list:
        substitution_string += "#include \"" + inlcude + "\"\n"

    substitution_string += (
        "#include \"memilio/compartments/parameter_studies.h\"\n"
        "#include \"memilio/data/analyze_result.h\"\n"
    )

    substitution_string += "\n#include \"pybind11/pybind11.h\"\n"
    if intermed_repr.simulation is True:
        substitution_string += (
            "#include \"pybind11/stl_bind.h\"\n"
            "#include \"Eigen/Core\"\n"
            "#include <vector>\n"
        )
    return substitution_string


def pretty_name_function(intermed_repr: IntermediateRepresentation) -> str:
    """ ! pretty_name_function 
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
        if key == "FP":
            continue

        if key == "AgeGroup":
            if intermed_repr.has_age_group is False:
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


def age_group(intermed_repr: IntermediateRepresentation) -> str:
    """! Generate the code for the AgeGroup class.
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


def population_enums(intermed_repr: IntermediateRepresentation) -> str:
    """! Sets the population enums as part of the bindings.

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
    for value in intermed_repr.model_base[0:]:
        if "Population" in value[0]:
            return "Populations"


def draw_sample(intermed_repr: IntermediateRepresentation) -> str:
    """! Sets the draw_sample function as part of the bindings.

    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.  
    """
    if not intermed_repr.has_draw_sample:
        return ""
    return (

        'm.def(\n\t\t"draw_sample",\n\t\t[]({namespace}Model<' +
        ScalarType(intermed_repr) +
        '>& model) {{\n\t\t\treturn {namespace}draw_sample(model);\n\t}},\n\tpy::arg("model"));\n'
    ).format(
        namespace=intermed_repr.namespace
    )


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
    """! Generate the code for the AgeGroup class.
    Not used by every model.

    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    if not intermed_repr.has_age_group:
        return ""

    return (
        "pymio::bind_CustomIndexArray<mio::UncertainValue, mio::AgeGroup>(m, \"AgeGroupArray\");\n"
    )


def parameterset_wrapper(intermed_repr: IntermediateRepresentation) -> str:
    """! Generate the code for the parameterset_wrapper needed when using age groups.
    Not used by every model.

    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    if intermed_repr.has_age_group is False:
        return ""

    return (
        "py::class_<{namespace}Parameters<"+ScalarType(
            intermed_repr)+">, pymio::EnablePickling::Required, {namespace}{parameterset}<"+ScalarType(intermed_repr)+">>(m, \"Parameters\");\n"
        "\t\t.def(py::init<mio::AgeGroup>());\n"
        "\t\t.def(\"check_constraints\", &{namespace}Parameters<" +
        ScalarType(intermed_repr)+">::check_constraints);\n"
        "\t\t.def(\"apply_constraints\", &{namespace}Parameters<" +
        ScalarType(intermed_repr)+">::apply_constraints);\n"
    ).format(
        namespace=intermed_repr.namespace,
        parameterset=intermed_repr.parameterset,
        parameterset_wrapper=intermed_repr.parameterset_wrapper
    )


def simulation(intermed_repr: IntermediateRepresentation) -> str:
    """! Generate the code for the Simulation class.
    Not used by every model.

    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    namespace = intermed_repr.namespace

    if intermed_repr.simulation is True:

        simulation = "&" + namespace + "simulate<" + \
            ScalarType(intermed_repr)+">"
        flow_simulation = "&" + namespace + \
            "flow_simulation<"+ScalarType(intermed_repr)+">"
        bind_simulation_class_string = "" + namespace + "Simulation"
        bind_flowsimulation_class_string = "" + namespace + "FlowSimulation"
        bind_simulation_string = ""
    else:

        simulation = "&mio::simulate<"+ScalarType(intermed_repr) + \
            "," + namespace + "Model<"+ScalarType(intermed_repr)+">>"
        flow_simulation = "&mio::flow_simulation<" + \
            namespace + "Model<"+ScalarType(intermed_repr)+">>"
        bind_simulation_class_string = "mio::Simulation"
        bind_flowsimulation_class_string = "mio::FlowSimulation"
        bind_simulation_string = "" + \
            ScalarType(intermed_repr)+"," + namespace + \
            "Model<"+ScalarType(intermed_repr)+">"

    sub_string = ""

    sub_string += (
        "pymio::bind_Simulation<{b_sim}<{sub}>>(m, \"Simulation\");\n\n\t"

        'm.def(\n\t\t"simulate", {sim},\n\t\t'
        'py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("model"), py::arg("integrator") = py::none(),\n\t\t'
        '"Simulate a {namespace} from t0 to tmax."\n\t\t'
        ');\n\n\t'

        "pymio::bind_ModelNode<{namespace}Model<" +
        ScalarType(intermed_repr)+">>(m, \"ModelNode\");\n\t"
        "pymio::bind_SimulationNode<{namespace}Simulation<>>(m, \"SimulationNode\");\n\t"
        "pymio::bind_ModelGraph<{namespace}Model<" +
        ScalarType(intermed_repr)+">>(m, \"ModelGraph\");\n\t"
        "pymio::bind_MobilityGraph<{namespace}Simulation<>>(m, \"MobilityGraph\");\n\t"
        "pymio::bind_GraphSimulation<mio::Graph<mio::SimulationNode<{b_sim}<>>, mio::MobilityEdge<" +
        ScalarType(intermed_repr)+">>>(m, \"MobilitySimulation\");\n\n\t"

    ).format(
        namespace=intermed_repr.namespace,
        sim=simulation,
        b_sim=bind_simulation_class_string,
        sub=bind_simulation_string
    )

    if intermed_repr.is_flowmodel is True:

        sub_string += (
            "pymio::bind_Flow_Simulation<{b_sim}<"+ScalarType(intermed_repr) +
            ", mio::FlowSimulation<{sub}>>>(m, \"FlowSimulation\");\n\n\t"

            'm.def(\n\t\t"simulate_flows", {flow_sim},\n\t\t'
            'py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("model"), py::arg("integrator") = py::none(),\n\t\t'
            '"Simulate a {namespace} with flows from t0 to tmax."\n\t\t'
            ');\n\n\t'

            "pymio::bind_SimulationNode<{namespace}FlowSimulation<>>(m, \"SimulationNode\");\n\t"
            "pymio::bind_MobilityGraph<{namespace}FlowSimulation<>>(m, \"MobilityGraph\");\n\t"
            "pymio::bind_GraphSimulation<mio::Graph<mio::SimulationNode<{b_flowsim}<>>, mio::MobilityEdge<" +
            ScalarType(intermed_repr)+">>>(m, \"MobilitySimulation\");\n\t"

        ).format(
            namespace=intermed_repr.namespace,
            flow_sim=flow_simulation,
            b_sim=bind_simulation_class_string,
            b_flowsim=bind_flowsimulation_class_string,
            sub=bind_simulation_string
        )
    return sub_string


def simulation_vector_definition(
        intermed_repr: IntermediateRepresentation) -> str:
    """! Generate the code for vector definition.
    Not used by every model.

    @param intermed_repr Dataclass holding the model features.
    @return Formatted string representing a part of the bindings.
    """
    namespace = intermed_repr.namespace
    sub_string = ""
    if intermed_repr.simulation is True:
        sim_class = "" + namespace + "Simulation"
        flowsim_class = "" + namespace + "FlowSimulation"
    else:
        sim_class = "mio::Simulation"
        flowsim_class = "mio::FlowSimualtion"

    sub_string += (
        "PYBIND11_MAKE_OPAQUE(std::vector<mio::Graph<mio::SimulationNode<{simulation_class}<>>, mio::MobilityEdge<"+ScalarType(
            intermed_repr)+">>);\n"
    ).format(
        simulation_class=sim_class
    )

    if intermed_repr.is_flowmodel:
        sub_string += (
            "PYBIND11_MAKE_OPAQUE(std::vector<mio::Graph<mio::SimulationNode<{flowsimulation_class}<>>, mio::MobilityEdge<"+ScalarType(
                intermed_repr)+">>);\n"
        ).format(
            flowsimulation_class=flowsim_class
        )
    return sub_string
