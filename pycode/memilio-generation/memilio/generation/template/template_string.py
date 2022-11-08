import string
from memilio.generation import IntermediateRepresentation

def includes(intermed_repr: IntermediateRepresentation) -> str:
    substition_string = (
        "#include \"pybind_util.h\"\n"
        "//Includes from pymio\n"
        "//Includes for the model\n"
        "//Includes from Memilio\n"
        "//Optional Includes: \"pybind11/stl_bind.h\", \"Eigen/Core\"\n"
    )
    if intermed_repr.simulation_class is not None:
        substition_string += "#include <vector>\n"
    return substition_string

def pretty_name_function(intermed_repr: IntermediateRepresentation) -> str:
    substition_string = (
        "namespace pymio"
        "{"
        "\n"
        "//specialization of pretty_name\n"
    )

    for key in intermed_repr.population_groups:
        if key == "AgeGroup":
            substition_string += (
                "template <>\n"
                "std::string pretty_name<mio::{enum_class}>()\n"
                "{{\n"
                "   return \"{enum_class}\";\n"
                "}}\n"
                "\n"
            ).format(
                enum_class = key
            )
        else:
            substition_string += (
                "template <>\n"
                "std::string pretty_name<{namespace}{enum_class}>()\n"
                "{{\n"
                "   return \"{enum_class}\";\n"
                "}}\n"
                "\n"
            ).format(
                namespace = intermed_repr.namespace,
                enum_class = key
            )

    return substition_string + "} // namespace pymio\n"

def population_enums(intermed_repr: IntermediateRepresentation) -> str:
    substition_string = ""
    for key, values in intermed_repr.enum_populations.items():
        substition_string += (
            "pymio::iterable_enum<{namespace}{enum_class}>(m, \"{enum_class}\")\n\t"
            ).format(
                namespace = intermed_repr.namespace,
                enum_class=key
            )

        for value in values[:-1]:
            substition_string += (
                "    .value(\"{comp_class}\", {namespace}{enum_class}::{comp_class})\n\t"
            ).format(
                namespace = intermed_repr.namespace,
                comp_class = value,
                enum_class = key
            )
        substition_string = substition_string.rstrip() + ";\n\n"
    return substition_string

def population(intermed_repr: IntermediateRepresentation) -> str:
    for value in intermed_repr.model_base[1:]:
        if "Population" in value[0]:
            return value[0]

def model_init(intermed_repr: IntermediateRepresentation) -> str:
    substition_string = ""
    for init in intermed_repr.model_init:

        if len(init["type"]) > 1:
            continue
        elif len(init["type"]) == 0:
            substition_string += "    .def(py::init<>())\n\t"
        else:
            substition_string += (
                "    .def(py::init<{type}>(), py::arg(\"{name}\"))\n\t"
            ).format(
                type = init["type"][0],
                name = init["name"][0]
            )
    return substition_string.rstrip() + ";\n"

def parameterset_indexing(intermed_repr: IntermediateRepresentation) -> str:
    if not intermed_repr.parameterset_wrapper:
        return ""

    return (
        "pymio::bind_CustomIndexArray<mio::UncertainValue, mio::AgeGroup>(m, \"AgeGroupArray\");\n"
    )

def parameterset_wrapper(intermed_repr: IntermediateRepresentation) -> str:
    if not intermed_repr.parameterset_wrapper:
        return ""

    return (
        "py::class_<{namespace}{parameterset_wrapper}, {namespace}{parameterset}>(m, \"{parameterset_wrapper}\")\n"
        "\t.def(py::init<mio::AgeGroup>())\n"
        "\t.def(\"check_constraints\", &{namespace}{parameterset_wrapper}::check_constraints)\n"
        "\t.def(\"apply_constraints\", &{namespace}{parameterset_wrapper}::apply_constraints);\n"
    ).format(
        namespace = intermed_repr.namespace,
        parameterset = intermed_repr.parameterset,
        parameterset_wrapper = intermed_repr.parameterset_wrapper
    )

def age_group(intermed_repr: IntermediateRepresentation) -> str:
    if not intermed_repr.age_group:
        return ""

    return (
        "py::class_<mio::AgeGroup, {base}>(m, \"AgeGroup\").def(py::init<{init}>());"
    ).format(
        namespace = intermed_repr.namespace,
        base = intermed_repr.age_group["base"],
        init = intermed_repr.age_group["init"][0]
    )

def simulation(intermed_repr: IntermediateRepresentation) -> str:
    if intermed_repr.simulation_class is None or (not intermed_repr.simulation_class.strip()):
        return ""

    return (
        "pymio::bind_Simulation<{namespace}{simulation_class}<>>(m, \"{simulation_class}\");\n"
    ).format(
        namespace = intermed_repr.namespace,
        simulation_class = intermed_repr.simulation_class
    )

def simulation_graph(intermed_repr: IntermediateRepresentation) -> str:
    if intermed_repr.simulation_class is None or (not intermed_repr.simulation_class.strip()):
        return ""

    return (
        "pymio::bind_ModelNode<{namespace}{model_class}>(m, \"ModelNode\");\n\t"
        "pymio::bind_SimulationNode<{namespace}{simulation_class}<>>(m, \"SimulationNode\");\n\t"
        "pymio::bind_ModelGraph<{namespace}{model_class}>(m, \"ModelGraph\");\n\t"
        "pymio::bind_MigrationGraph<{namespace}{simulation_class}<>>(m, \"MigrationGraph\");\n\t"
        "pymio::bind_GraphSimulation<mio::Graph<mio::SimulationNode<{namespace}{simulation_class}<>>, mio::MigrationEdge>>(m, \"MigrationSimulation\");\n\t"
    ).format(
        namespace = intermed_repr.namespace,
        model_class = intermed_repr.model_class,
        simulation_class = intermed_repr.simulation_class
    )

def simulation_vector_definition(intermed_repr: IntermediateRepresentation) -> str:
    if intermed_repr.simulation_class is None or (not intermed_repr.simulation_class.strip()):
        return ""

    return (
        "PYBIND11_MAKE_OPAQUE(std::vector<mio::Graph<mio::SimulationNode<{namespace}{simulation_class}>, mio::MigrationEdge>>);\n"
    ).format(
        namespace = intermed_repr.namespace,
        simulation_class = intermed_repr.simulation_class
    )