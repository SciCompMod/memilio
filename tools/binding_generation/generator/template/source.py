from generator.utility import join

def print_source(model, file=None):
    s = source_string(model)
    if file is None:
        print(s)
    else:
        with open(file, 'w') as file:
            file.write(s)
    return

def source_string(model):
    return (
        "/*\n"
        "* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)\n"
        "*\n"
        "* Authors: Martin Siggel, Daniel Abele, Martin J. Kuehn, Jan Kleinert\n"
        "*\n"
        "* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>\n"
        "*\n"
        "* Licensed under the Apache License, Version 2.0 (the \"License\");\n"
        "* you may not use this file except in compliance with the License.\n"
        "* You may obtain a copy of the License at\n"
        "*\n"
        "*     http://www.apache.org/licenses/LICENSE-2.0\n"
        "*\n"
        "* Unless required by applicable law or agreed to in writing, software\n"
        "* distributed under the License is distributed on an \"AS IS\" BASIS,\n"
        "* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.\n"
        "* See the License for the specific language governing permissions and\n"
        "* limitations under the License.\n"
        "*/\n"
        "\n"
        "{includes}"
        "\n"
        "{pretty_name}"
        "\n"
        "PYBIND11_MODULE(_simulation_{pymio_name}, m)\n"
        "{{\n"
        "{enum_infection_state}"
        "\n"
        "{indexing}"
        "\n"
        "pymio::bind_ParameterSet<{namespace}ParametersBase>(m, \"ParametersBase\");\n"
        "\n"
        "pymio::bind_Population{population_template_arguments}(m, \"Population\");\n"
        "\n"
        "using Populations = mio::Populations{population_template_arguments};\n"
        "pymio::bind_CompartmentalModel<Populations, {namespace}ParametersBase>(m, \"ModelBase\");\n"
        "py::class_<{namespace}Model, mio::CompartmentalModel<Populations, {namespace}ParametersBase>>(m, \"Model\")\n"
        "   .def({model_init});\n"
        "\n"
        "m.def(\n"
        "\"simulate\",\n"
        "    [](double t0, double tmax, double dt, const {namespace}Model& model) {{\n"
        "        return mio::simulate(t0, tmax, dt, model);\n"
        "    }},\n"
        "    \"Simulates a {name} from t0 to tmax.\", py::arg(\"t0\"), py::arg(\"tmax\"), py::arg(\"dt\"), py::arg(\"model\"));\n"
        "\n"
        "m.attr(\"__version__\") = \"dev\";\n"
        "}}\n"
    ).format(
        includes = include_string(model),
        pretty_name = pretty_name_string(model),
        name = model.name,
        pymio_name = model.pymio_name,
        namespace = model.namespace,
        indexing = indexing_binds_string(model),
        population_template_arguments = population_string(model),
        enum_infection_state = enum_infection(model),
        model_init = model_init(model)
    )

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
    for group in model.population_groups:
        str += (
            "template <>\n"
            "std::string pretty_name<{namespace}{group}>()\n"
            "{{\n"
            "    return \"{group}\";\n"
            "}}\n"
            "\n"
        ).format(
            namespace=model.namespace,
            group = group
        )

    return str + "} // namespace pymio\n"

def population_string(model):
    str = "<"

    for group in model.population_groups:
        str += (
            "{group}, "
        ).format(
            group = group
        )
    return str[:-2] + ">"

def indexing_binds_string(model):
    return (
        "pymio::bind_Index<{namespace}InfectionState>(m, \"Index_InfectionState\");\n"
        "pymio::bind_CustomIndexArray<mio::UncertainValue, {namespace}InfectionState>(m, \"PopulationArray\");\n"
    ).format(
        namespace = model.namespace
    )

def model_init(model):

    if model.init[0].strip():
        return (
            "py::init<{type}>(), py::arg(\"{name}\")"
        ).format(
            type = model.init[0].split(" ")[0],
            name = model.init[0].split(" ")[1]
        )

    return "py::init<>()"

def enum_infection(model):
    str = (
        "pymio::iterable_enum<{namespace}InfectionState>(m, \"InfectionState\")\n"
        ).format(
            namespace=model.namespace
        )

    for comp in model.compartments:
        str += (
            "    .value(\"{comp}\", {namespace}InfectionState::{comp})\n"
        ).format(
            comp = comp,
            namespace=model.namespace
        )
    return str.rstrip() + ";\n"
