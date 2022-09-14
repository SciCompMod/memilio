import string

def include_string(intermed_repr):
    str = (
        "#include \"pybind_util.h\"\n"
        "//Includes from pymio\n"
        "//Includes for the model\n"
        "//Includes from Memilio\n"
        "//Optional Includes: \"pybind11/stl_bind.h\", \"Eigen/Core\"\n"
    )
    if intermed_repr.simulation_class is not None:
        str += "#include <vector>\n"
    return str

def pretty_name_string(intermed_repr):
    str = (
        "namespace pymio"
        "{"
        "\n"
        "//specialization of pretty_name\n"
    )

    for key in intermed_repr.population_groups:
        if key == "AgeGroup":
            str += (
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
            str += (
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

    return str + "} // namespace pymio\n"

def enum_populations(intermed_repr):
    str = ""
    for key, values in intermed_repr.enum_populations.items():
        str += (
            "pymio::iterable_enum<{namespace}{enum_class}>(m, \"{enum_class}\")\n\t"
            ).format(
                namespace = intermed_repr.namespace,
                enum_class=key
            )

        for value in values[:-1]:
            str += (
                "    .value(\"{comp_class}\", {namespace}{enum_class}::{comp_class})\n\t"
            ).format(
                namespace = intermed_repr.namespace,
                comp_class = value,
                enum_class = key
            )
        str = str.rstrip() + ";\n\n"
    return str

def population_string(intermed_repr):
    for value in intermed_repr.model_base[1:]:
        if "Population" in value[0]:
            return value[0]

def model_init(intermed_repr):
    str = ""
    for init in intermed_repr.model_init:

        if len(init["type"]) > 1:
            continue
        elif len(init["type"]) == 0:
            str += "    .def(py::init<>())\n\t"
        else:
            str += (
                "    .def(py::init<{type}>(), py::arg(\"{name}\"))\n\t"
            ).format(
                type = init["type"][0],
                name = init["name"][0]
            )
    return str.rstrip() + ";\n"

def parameterset_indexing(intermed_repr):
    if not intermed_repr.parameterset_wrapper:
        return ""
    return (
        "pymio::bind_CustomIndexArray<mio::UncertainValue, mio::AgeGroup>(m, \"AgeGroupArray\");\n"
    )

def parameterset_wrapper(intermed_repr):
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

def age_group(intermed_repr):
    if not intermed_repr.age_group:
        return ""
    return (
        "py::class_<mio::AgeGroup, {base}>(m, \"AgeGroup\").def(py::init<{init}>());"
    ).format(
        namespace = intermed_repr.namespace,
        base = intermed_repr.age_group["base"],
        init = intermed_repr.age_group["init"][0]
    )

def simulation(intermed_repr):
    if intermed_repr.simulation_class is None or (not intermed_repr.simulation_class.strip()):
        return ""
    return (
        "pymio::bind_Simulation<{namespace}{simulation_class}<>>(m, \"{simulation_class}\");\n"
    ).format(
        namespace = intermed_repr.namespace,
        simulation_class = intermed_repr.simulation_class
    )

def simulation_graph(intermed_repr):
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

def simulation_vector_definition(intermed_repr):
    if intermed_repr.simulation_class is None or (not intermed_repr.simulation_class.strip()):
        return ""
    return (
        "PYBIND11_MAKE_OPAQUE(std::vector<mio::Graph<mio::SimulationNode<{namespace}{simulation_class}>, mio::MigrationEdge>>);\n"
    ).format(
        namespace = intermed_repr.namespace,
        simulation_class = intermed_repr.simulation_class
    )