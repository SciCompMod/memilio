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
        "{enum_population_states}"
        "\n"
        "{indexing}"
        "\n"
        "pymio::bind_ParameterSet<{namespace}ParametersBase>(m, \"ParametersBase\");\n"
        "\n"
        "pymio::bind_Population{population_template_arguments}(m, \"Population\");\n"
        "\n"
        "using Populations = mio::Populations{population_template_arguments};\n"
        "pymio::bind_CompartmentalModel<Populations, {namespace}ParametersBase>(m, \"ModelBase\");\n"
        "py::class_<{model_class}, mio::CompartmentalModel<Populations, {namespace}ParametersBase>>(m, \"{model_cpp_name}\")\n"
        "{model_init}"
        "\n"
        "m.def(\n"
        "\"simulate\",\n"
        "    [](double t0, double tmax, double dt, const {model_class}& model) {{\n"
        "        return mio::simulate(t0, tmax, dt, model);\n"
        "    }},\n"
        "    \"Simulates a {model_cpp_name} from t0 to tmax.\", py::arg(\"t0\"), py::arg(\"tmax\"), py::arg(\"dt\"), py::arg(\"model\"));\n"
        "\n"
        "m.attr(\"__version__\") = \"dev\";\n"
        "}}\n"
    ).format(
        includes = include_string(model),
        pretty_name = pretty_name_string(model),
        model_class = model.name,
        model_cpp_name = model.name.split("::")[-1],
        pymio_name = model.pymio_name,
        namespace = ("::".join(model.name.split("::")[:-1]) + "::"),
        indexing = indexing_binds_string(model),
        population_template_arguments = population_string(model),
        enum_population_states = enum(model),
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
    for key in model.enum_dict:

        if key.split("::")[-1] not in model.population_groups:
            continue

        str += (
            "template <>\n"
            "std::string pretty_name<{enum_class}>()\n"
            "{{\n"
            "    return \"{enum_name}\";\n"
            "}}\n"
            "\n"
        ).format(
            enum_class=key,
            enum_name=key.split("::")[-1]
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
    for key in model.enum_dict:

        if key.split("::")[-1] not in model.population_groups:
            continue

        return (
            "pymio::bind_Index<{enum_class}>(m, \"Index_{enum_name}\");\n"
            "pymio::bind_CustomIndexArray<mio::UncertainValue, {enum_class}>(m, \"PopulationArray\");\n"
        ).format(
            enum_class=key,
            enum_name=key.split("::")[-1]
        )

def model_init(model):
    str = ""
    for init in model.init:

        if len(init["type"]) > 1:
            continue
        elif len(init["type"]) == 0:
            str += "   .def(py::init<>())\n"
        else:
            str += (
                "   .def(py::init<{type}>(), py::arg(\"{name}\"))\n"
            ).format(
                type = init["type"][0],
                name = init["name"][0]
            )
    return str.rstrip() + ";\n"

def enum(model):
    str = ""
    for key, values in model.enum_dict.items():

        if key.split("::")[-1] not in model.population_groups:
            continue

        str += (
            "pymio::iterable_enum<{enum_class}>(m, \"{enum_name}\")\n"
            ).format(
                enum_class=key,
                enum_name=key.split("::")[-1]
            )

        for value in values:
            str += (
                "    .value(\"{comp_name}\", {comp_class})\n"
            ).format(
                comp_class = value,
                comp_name = value.split("::")[-1]
            )
        str = str.rstrip() + ";\n\n"
    return str
