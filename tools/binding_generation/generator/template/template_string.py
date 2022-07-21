import string

def include_string(model):
   return (
        "#include \"templates.h\"\n"
        "#include \"ode_seir/model.h\"\n"
        "#include \"Eigen/Core\"\n"
        "#include \"pybind11/stl_bind.h\"\n"
        "#include <vector>\n"
    ).format(
    )

def pretty_name_string(model):
    str = (
        "namespace pymio"
        "{"
        "\n"
        "//specialization of pretty_name\n"
    )

    for key in model.population_groups:

        str += (
            "template <>\n"
            "std::string pretty_name<{namespace}{enum_class}>()\n"
            "{{\n"
            "   return \"{enum_class}\";\n"
            "}}\n"
            "\n"
        ).format(
            namespace = model.namespace,
            enum_class = key
        )

    return str + "} // namespace pymio\n"

def enum_populations(model):
    str = ""
    for key, values in model.enum_populations.items():
        str += (
            "pymio::iterable_enum<{namespace}{enum_class}>(m, \"{enum_class}\")\n\t"
            ).format(
                namespace = model.namespace,
                enum_class=key
            )

        for value in values:
            str += (
                "    .value(\"{comp_class}\", {namespace}{comp_class})\n\t"
            ).format(
                namespace = model.namespace,
                comp_class = value
            )
        str = str.rstrip() + ";\n\n"
    return str

def population_string(model):
    for value in model.model_base[1:]:
        if "Population" in value[0]:
            return value[0]

def model_init(model):
    str = ""
    for init in model.model_init:

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

def agegroup(model):
    if not model.age_group:
        return ""
    return (
        "py::class_<mio::AgeGroup, {base}>(m, \"AgeGroup\").def(py::init<{init}>());"
    ).format(
        base = model.age_group["base"],
        init = model.age_group["init"][0]
    )

def simulation(model):
    if model.simulation_class is None or (not model.simulation_class.strip()):
        return ""
    return (
        "pymio::bind_Simulation<{namespace}{simulation_class}<>>(m, \"{simulation_class}\");\n"
    ).format(
        namespace = model.namespace,
        simulation_class = model.simulation_class
    )

def simulation_graph(model):
    if model.simulation_class is None or (not model.simulation_class.strip()):
        return ""
    return (
        "pymio::bind_ModelNode<{namespace}{model_class}>(m, \"ModelNode\");\n\t"
        "pymio::bind_SimulationNode<{namespace}{simulation_class}<>>(m, \"SimulationNode\");\n\t"
        "pymio::bind_ModelGraph<{namespace}{model_class}>(m, \"ModelGraph\");\n\t"
        "pymio::bind_MigrationGraph<{namespace}{simulation_class}<>>(m, \"MigrationGraph\");\n\t"
        "pymio::bind_GraphSimulation<mio::Graph<mio::SimulationNode<{namespace}{simulation_class}<>>, mio::MigrationEdge>>(m, \"MigrationSimulation\");\n\t"
    ).format(
        namespace = model.namespace,
        model_class = model.model_class,
        simulation_class = model.simulation_class
    )