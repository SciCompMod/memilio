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


def handle_all_bindings(intermed_repr: IntermediateRepresentation) -> str:
    result = ""
    print(intermed_repr.found_bindings)
    for bindings in intermed_repr.found_bindings:
        if bindings.get("type") == "class":
            result += bind_classes(intermed_repr, bindings)
        elif bindings.get("type") == "function":
            result += bind_functions(intermed_repr, bindings)
    return result.strip()


def bind_classes(intermed_repr: IntermediateRepresentation, bindings: dict) -> str:

    class_block = ""
    scalarType = ScalarType(intermed_repr)
    class_namespace = bindings.get("namespace", "")
    class_name = bindings.get("name", "")
    base_classes = bindings.get("base_classes", "")

    template_params = f"{class_namespace}{class_name}<{scalarType}>"
    if base_classes:

        base_template = ", ".join(
            f"{class_namespace}{base_class.replace('FP', scalarType)}" for base_class in ensure_list(base_classes))
        class_block += f'\tpy::class_<{template_params}, {base_template}>(m, "{class_name}")\n'
    else:
        class_block += f'\tpy::class_<{template_params}>(m, "{class_name}")\n'

    init_params = bindings.get("init", [])
    class_block += generate_constructor(init_params)

    methods = bindings.get("methods", [])
    for method in methods:
        class_block += (f'\t' + bind_functions(intermed_repr, method))
    return class_block + '\t;\n'


def bind_functions(intermed_repr: IntermediateRepresentation, bindings: dict) -> str:
    """ Bind the functions found in the intermediate representation.
    :param intermed_repr: Dataclass holding the model features.
    :return: Formatted string representing a part of the bindings."""

    result = ""
    name = bindings.get("name", "")
    namespace = bindings.get("namespace", "")
    arg_types = ensure_list(bindings.get("arg_types", []))
    arg_names = ensure_list(bindings.get("arg_names", []))
    call_args = ", ".join(arg_names)
    py_args = ", ".join(f'py::arg("{n}")' for n in arg_names)
    is_member = bindings.get("is_member", False)
    parent_name = bindings.get("parent_name", "")
    scalartype = ScalarType(intermed_repr)
    prefix = set_prefix(bindings)

    if is_overloaded_function(bindings, intermed_repr.found_bindings):

        result += overloaded_function_bindings(
            name, namespace, scalartype, arg_types, py_args, bindings, prefix)

    elif needs_lambda_binding(bindings):

        lambda_args = generate_lambda_args(
            namespace, scalartype, arg_types, arg_names)

        result += lambda_binding(name, namespace, lambda_args,
                                 call_args, py_args, prefix)
    else:
        result += direct_binding(
            name, namespace, scalartype, py_args, prefix, bool(is_member), parent_name)

    return result


def generate_constructor(init_params: list) -> str:
    """ Generate the constructor for a class.
    :param init_params: List of initialization parameters for the constructor.
    :return: Formatted string representing the constructor.
    """
    if not init_params:
        return "\t\t.def(py::init<>())\n"

    types_list = ", ".join(param["type"] for param in init_params)
    args_list = ", ".join(
        f'py::arg("{param["name"]}")' for param in init_params)

    return f'\t\t.def(py::init<{types_list}>(), {args_list})\n'


def generate_lambda_args(namespace: str, scalartype: str, arg_types: list, arg_names: list) -> str:
    """ Generate the arguments for a lambda function.
    :param namespace: Namespace of the function.
    :param scalartype: Scalar type used in the function.
    :param arg_types: List of argument types for the function.
    :param arg_names: List of argument names for the function.
    :return: Formatted string representing the lambda arguments.
    """

    return ", ".join(
        f"{namespace}{t.replace('FP', scalartype)} {n}"
        for t, n in zip(arg_types, arg_names)
    )


def set_prefix(bindings: dict) -> str:
    """ Set the prefix for the bindings.
    :param bindings: A dictionary representing a function or method node from the AST.
    :return: Formatted string representing a part of the bindings.
    """
    is_member = bindings.get("is_member", False)
    if not is_member:
        prefix = "m"
    else:
        prefix = ""

    return prefix


def needs_lambda_binding(bindings: dict) -> bool:
    """
    Determines whether a lambda binding is necessary for the given function or method.
    :param bindings: A dictionary representing a function or method node from the AST.
    :return: True if a lambda binding should be generated, False otherwise.
    """

    arg_types = ensure_list(bindings.get("arg_types", []))
    if arg_types and "Model" in arg_types[0]:
        return True

    return False


def lambda_binding(name: str, namespace: str, lambda_args: str, call_args: str, py_args: str, prefix: str) -> str:
    """ Generate a lambda binding for a function.
    :param name: Name of the function.
    :param namespace: Namespace of the function.
    :param lambda_args: Arguments for the lambda function.
    :param call_args: Arguments to call the function.
    :param py_args: Python arguments for the function.
    :return: Formatted string representing a part of the bindings."""

    return (
        f'\t{prefix}.def(\n'
        f'\t\t"{name}",\n'
        f'\t\t[]({lambda_args}) {{\n'
        f'\t\treturn {namespace}{name}({call_args});\n'
        f'\t\t}},\n'
        f'\t\t{py_args}\n'
        f'\t);\n\n'
    )


def direct_binding(name: str, namespace: str, scalartype: str, py_args: str, prefix: str, is_member: bool, parent_name: str) -> str:
    """ Generate a direct binding for a function.
    :param name: Name of the function.
    :param namespace: Namespace of the function.
    :param scalartype: Scalar type used in the function.
    :param py_args: Python arguments for the function.
    :return: Formatted string representing a part of the bindings.
    """
    if not is_member:
        return (
            f'\t{prefix}.def('
            f'\n\t\t"{name}",'
            f'\n\t\t&{namespace}{name}<{scalartype}>,'
            f'\n\t\t{py_args}\n\t);\n\n'
        )
    else:
        return (
            f'\t{prefix}.def("{name}", &{namespace}{parent_name}<{scalartype}>::{name})\n'

        )


def is_overloaded_function(binding: dict, found_bindings: list) -> bool:
    """ Check if a function binding is overloaded.
    :param binding: The binding dictionary containing function details.
    :param found_bindings: List of all found bindings.
    :return: True if the function is overloaded, False otherwise."""

    name = binding.get("name", "")
    namespace = binding.get("namespace", "")
    count = sum(
        1 for b in found_bindings
        if b.get("name", "") == name and b.get("namespace", "") == namespace
    )
    return count > 1


def overloaded_function_bindings(name: str, namespace: str, scalartype: str, arg_types: list, py_args: str, bindings: dict, prefix: str) -> str:
    """ Generate bindings for overloaded functions.
    :param name: Name of the function.
    :param namespace: Namespace of the function.
    :param scalartype: Scalar type used in the function.
    :param arg_types: List of argument types for the function.
    :param py_args: Python arguments for the function.
    :param bindings: The binding dictionary containing function details.
    :return: Formatted string representing a part of the bindings."""

    formatted_arg_types = [
        t.replace("FP", scalartype) for t in arg_types]
    return_type = bindings.get(
        "return_type", "void").replace("FP", scalartype)

    is_const = bindings.get("is_const", False)

    func_ptr = f"{return_type} ({namespace}*)({', '.join(formatted_arg_types)})"

    if is_const:
        func_ptr += " const"

    return (
        f'\t{prefix}.def(\n'
        f'\t\t"{name}",\n'
        # overload_cast helps generate the correct function pointer type from the template arguments
        f'\t\tpy::overload_cast<{", ".join(formatted_arg_types)}>(static_cast<{func_ptr}>(&{namespace}{name})),'
        f'\t\t{"py::const_," if is_const else ""}\n'
        f'\t\t{py_args}\n'
        f'\t);\n\n'
    )


def ensure_list(value) -> list:
    """ Ensure the value is a list. If not, wrap it in a list.
    :param value: The input value, possibly a list.
    :return: A list containing the value(s).
    """
    if isinstance(value, list):
        return value
    elif value is None:
        return []
    else:
        return [value]


def ScalarType(intermed_repr: IntermediateRepresentation) -> str:
    """ Set the datatype for the bindings via. intermediate_representation.
    @return string from intermediate_representation
    """
    scalartypestr = intermed_repr.scalartype
    return scalartypestr


def includes(intermed_repr: IntermediateRepresentation) -> str:
    """ Fills in the Includes for the binding
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
    """ Sets the draw_sample function as part of the bindings.

    :param intermed_repr: Dataclass holding the model features.
    :return: Formatted string representing a part of the bindings.
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

    :param intermed_repr: Dataclass holding the model features.
    :param intermed_repr: IntermediateRepresentation:
    :returns: Formatted string representing a part of the bindings.

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

    :param intermed_repr: Dataclass holding the model features.
    :param intermed_repr: IntermediateRepresentation:
    :returns: Formatted string representing a part of the bindings.

    """
    if not intermed_repr.has_age_group:
        return ""

    return (
        "pymio::bind_CustomIndexArray<mio::UncertainValue, mio::AgeGroup>(m, \"AgeGroupArray\");\n"
    )


def parameterset_wrapper(intermed_repr: IntermediateRepresentation) -> str:
    """! Generate the code for the parameterset_wrapper needed when using age groups.
    Not used by every model.

    :param intermed_repr: Dataclass holding the model features.
    :param intermed_repr: IntermediateRepresentation:
    :returns: Formatted string representing a part of the bindings.

    """
    if intermed_repr.has_age_group is False:
        return ""

    return (
        "py::class_<{namespace}Parameters<"+ScalarType(
            intermed_repr)+">, pymio::EnablePickling::Required, {namespace}{parameterset}<"+ScalarType(intermed_repr)+">>(m, \"Parameters\")\n"
        "\t\t.def(py::init<mio::AgeGroup>())\n"
        "\t\t.def(\"check_constraints\", &{namespace}Parameters<" +
        ScalarType(intermed_repr)+">::check_constraints)\n"
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

    :param intermed_repr: Dataclass holding the model features.
    :param intermed_repr: IntermediateRepresentation:
    :returns: Formatted string representing a part of the bindings.

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

    :param intermed_repr: Dataclass holding the model features.
    :param intermed_repr: IntermediateRepresentation:
    :returns: Formatted string representing a part of the bindings.

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
