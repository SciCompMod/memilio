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
@file template_string.py
@brief Generate small pieces of the target code as strings.
"""
from __future__ import annotations

from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from memilio.generation import IntermediateRepresentation
    from memilio.generation.info_types import binding_type_info, method_type_info


def handle_all_bindings(intermed_repr: IntermediateRepresentation) -> str:
    """Generate the bindings for all found classes and functions in the intermediate representation.

    :param intermed_repr: Dataclass holding the model features.
    """
    result = ""
    for bindings in intermed_repr.found_bindings:
        if bindings.type == "class":
            if bindings.cursorkind == "CLASS_TEMPLATE" and needs_template(bindings, intermed_repr):
                result += template_class_wrapper(intermed_repr, bindings)
            else:
                result += bind_classes(intermed_repr, bindings)

        elif bindings.type == "function":
            result += bind_functions(intermed_repr, bindings)
    return result.strip()


def needs_template(bindings: binding_type_info, intermed_repr: IntermediateRepresentation) -> bool:
    """ Check if a template wrapper is needed for a class.

    :param bindings: A dataclass entry representing a class node from the AST.
    :param intermed_repr: Dataclass holding the model features.
    :return: True if a template wrapper should be generated, False otherwise.
    """
    # Extract the globally defined scalar type
    scalartype = ScalarType(intermed_repr)

    template_args = set(bindings.template_args)

    # If there are no template parameters, it is not a template class, so no template wrapper is needed
    if not bindings.template_args:
        return False

    template_used = False

    # If a base class depends on a template parameter → template is required
    for base in ensure_list(bindings.base_classes):
        if any(arg in base for arg in template_args):
            template_used = True
            break

    # If any method uses template parameters → template is required
    if not template_used:
        for method in bindings.methods:
            arg_types = ensure_list(method.arg_types)

            for types in arg_types:
                if any(arg in str(types) for arg in template_args):
                    template_used = True
                    break

    if not template_used:
        return False

    # If a concrete scalar type exists, templates are already resolved, so no template wrapper is needed
    if scalartype is not None:
        return False

    return False


def template_class_wrapper(intermed_repr: IntermediateRepresentation, bindings: binding_type_info) -> str:
    """ Generate the template class wrapper for the bindings.

    :param intermed_repr: Dataclass holding the model features.
    :param bindings: A dataclass entry representing a class node from the AST.
    :return: Formatted string representing a part of the bindings.
    """

    template_class_block = ""
    template_args = bindings.template_args
    template_string = ", ".join(
        f"typename {arg}" for arg in template_args)
    class_name = bindings.name
    if template_string:
        full_template = f"class {class_name}, {template_string}"
    else:
        full_template = f"class {class_name}"

    template_class_block += (

        f'template <{full_template}>\n\t'
        f'void bind_{class_name}(py::module& m, const std::string& name)\n\t'
        f'{{\n\t')
    template_class_block += f'\t' + bind_classes(intermed_repr, bindings)
    template_class_block += (
        '\t}\n')

    return template_class_block


def bind_classes(intermed_repr: IntermediateRepresentation, bindings: binding_type_info) -> str:
    """ Bind the classes found in the intermediate representation.

    :param intermed_repr: Dataclass holding the model features.
    :param bindings: A dataclass entry representing a class node from the AST.
    :return: Formatted string representing a part of the bindings.
    """

    class_block = ""
    scalarType = ScalarType(intermed_repr)
    class_namespace = bindings.namespace
    class_name = bindings.name
    base_classes = bindings.base_classes

    template_params = f"{class_namespace}{class_name}<{scalarType}>"

    if bindings.cursorkind == "CLASS_TEMPLATE":
        if base_classes:
            base_template = ", ".join(
                f"{class_namespace}{base_class.replace('FP', scalarType)}" for base_class in ensure_list(base_classes))
            class_block += (
                f'py::class_<{template_params}, {base_template}>(m, name.c_str())\n')
        else:
            class_block += (
                f'py::class_<{template_params}>(m, name.c_str())\n')
    else:
        if base_classes:
            base_template = ", ".join(
                f"{class_namespace}{base_class.replace('FP', scalarType)}" for base_class in ensure_list(base_classes))
            class_block += f'\tpy::class_<{template_params}, {base_template}>(m, "{class_name}")\n'
        else:
            class_block += f'\tpy::class_<{template_params}>(m, "{class_name}")\n'

    init_params = bindings.init
    class_block += generate_constructor(init_params)

    methods = bindings.methods
    for method in methods:
        class_block += (f'\t' + bind_functions(intermed_repr,
                        method, class_name))
    return class_block + '\t;\n\n'


def bind_functions(intermed_repr: IntermediateRepresentation, bindings: binding_type_info, class_name: str = "") -> str:
    """ Bind the functions found in the intermediate representation.

    :param intermed_repr: Dataclass holding the model features.
    :param bindings: A dataclass entry representing a function or method node from the AST
    :param class_name: Name of the class if the function is a member function.
    :return: Formatted string representing a part of the bindings."""

    result = ""
    name = bindings.name
    namespace = bindings.namespace
    arg_types = ensure_list(bindings.arg_types)
    arg_names = ensure_list(bindings.arg_names)
    py_args = ", ".join(f'py::arg("{n}")' for n in arg_names)
    is_member = bindings.is_member
    parent_name = bindings.parent_name
    scalartype = ScalarType(intermed_repr)
    prefix = set_prefix(bindings)

    if is_overloaded_function(bindings, intermed_repr.found_bindings):

        result += overloaded_function_bindings(
            name, namespace, scalartype, arg_types, py_args, bindings, prefix)

    elif needs_lambda_binding(bindings):

        result += lambda_binding(name, namespace, scalartype, class_name)
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


def set_prefix(bindings: binding_type_info) -> str:
    """ Set the prefix for the bindings.

    :param bindings: A dictionary representing a function or method node from the AST.
    :return: Formatted string representing a part of the bindings.
    """
    return "m" if not bindings.is_member else ""


def needs_lambda_binding(bindings: binding_type_info) -> bool:
    """ Determines whether a lambda binding is necessary for the given function or method.

    :param bindings: A dataclass entry representing a function or method node from the AST.
    :return: True if a lambda binding should be generated, False otherwise.
    """
    name = bindings.name
    return bool(bindings.is_member and (name.startswith("get") or name.startswith("set")))


def lambda_binding(name: str, namespace: str, scalartype: str, class_name: str = "") -> str:
    """ Generate a lambda binding for a function.

    :param name: Name of the function.
    :param namespace: Namespace of the function.
    :param scalartype: Scalar type used in the function.
    :param class_name: Name of the class if the function is a member function.
    :return: Formatted string representing a part of the bindings."""

    qualified_class = f"{namespace}{class_name}<{scalartype}>"
    property_name = name[4:]
    getter = (
        f'\t\t[]({qualified_class} const& self) -> {scalartype} {{\n'
        f'\t\t\treturn self.{name}();\n'
        f'\t\t\t}}'
    )

    if name.startswith("set"):
        setter = (
            f',\n'
            f'\t\t[]({qualified_class}& self, {scalartype} value) {{\n'
            f'\t\t\tself.{name}(value);\n'
            f'\t\t\t}}'
        )
    else:
        setter = ''

    return (
        f'\t.def_property(\n'
        f'\t\t\t"{property_name}",\n'
        f'\t{getter}{setter})\n'
    )


def direct_binding(name: str, namespace: str, scalartype: str, py_args: str, prefix: str, is_member: bool, parent_name: str) -> str:
    """ Generate a direct binding for a function.

    :param name: Name of the function.
    :param namespace: Namespace of the function.
    :param scalartype: Scalar type used in the function.
    :param py_args: Python arguments for the function.
    :param prefix: Prefix for the binding (e.g., "m" for module-level functions).
    :param is_member: Boolean indicating whether the function is a member function.
    :param parent_name: Name of the parent class if the function is a member function.
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


def is_overloaded_function(binding: binding_type_info, found_bindings: list[binding_type_info]) -> bool:
    """ Check if a function binding is overloaded.

    :param binding: A dataclass entry representing a function or method node from the AST.
    :param found_bindings: List of all found bindings.
    :return: True if the function is overloaded, False otherwise."""

    count = sum(
        1 for b in found_bindings
        if b.name == binding.name and b.namespace == binding.namespace
    )
    return count > 1


def overloaded_function_bindings(name: str, namespace: str, scalartype: str, arg_types: list, py_args: str, bindings: binding_type_info, prefix: str) -> str:
    """ Generate bindings for overloaded functions.

    :param name: Name of the function.
    :param namespace: Namespace of the function.
    :param scalartype: Scalar type used in the function.
    :param arg_types: List of argument types for the function.
    :param py_args: Python arguments for the function.
    :param bindings: A dataclass entry representing a overloaded function or method node from the AST.
    :param prefix: Prefix for the binding (e.g., "m" for module-level functions).
    :return: Formatted string representing a part of the bindings."""

    formatted_arg_types = [
        t.replace("FP", scalartype) for t in arg_types]
    return_type = bindings.get(
        "return_type", "void").replace("FP", scalartype)

    is_const = bindings.is_const

    func_ptr = f"{return_type} ({namespace}*)({', '.join(formatted_arg_types)})"

    if is_const:
        func_ptr += " const"

    return (
        f'\t{prefix}.def(\n'
        f'\t\t"{name}",\n'
        f'\t\tpy::overload_cast<{", ".join(formatted_arg_types)}>(&{namespace}{name})),'
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

    :param intermed_repr: Dataclass holding the model features.
    :return: String from intermediate_representation
    """
    scalartypestr = intermed_repr.scalartype
    return scalartypestr


def includes(intermed_repr: IntermediateRepresentation) -> str:
    """ Fills in the Includes for the binding
    :param intermed_repr: Dataclass holding the model features.
    :return: Formatted string representing a part of the bindings.
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
            "#include \"compartments/compartmental_model.h\"\n"
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
    """  pretty_name_function
    :param intermed_repr: Dataclass holding the model features.
    :return: Formatted string representing a part of the bindings.
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
    """ Generate the code for the AgeGroup class.
    Not used by every model.

    :param intermed_repr: Dataclass holding the model features.
    :return: Formatted string representing a part of the bindings.
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
    """ Sets the population enums as part of the bindings.

    :param intermed_repr: Dataclass holding the model features.
    :return: Formatted string representing a part of the bindings.
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
    """ Sets the population as part of the bindings.

    :param intermed_repr: Dataclass holding the model features.
    :return: Formatted string representing a part of the bindings.

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
    """  Generate the code for the model initialization in the bindings.
    Not used by every model.

    :param intermed_repr: Dataclass holding the model features.
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
    """ Generate the code for the AgeGroup class.
    Not used by every model.

    :param intermed_repr: Dataclass holding the model features.
    :returns: Formatted string representing a part of the bindings.

    """
    if not intermed_repr.has_age_group:
        return ""

    return (
        "pymio::bind_CustomIndexArray<mio::UncertainValue, mio::AgeGroup>(m, \"AgeGroupArray\");\n"
    )


def parameterset_wrapper(intermed_repr: IntermediateRepresentation) -> str:
    """ Generate the code for the parameterset_wrapper needed when using age groups.
    Not used by every model.

    :param intermed_repr: Dataclass holding the model features.
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
    """ Generate the code for the Simulation class.
    Not used by every model.

    :param intermed_repr: Dataclass holding the model features.
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
    """ Generate the code for vector definition.
    Not used by every model.

    :param intermed_repr: Dataclass holding the model features.
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
