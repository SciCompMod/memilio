import string

class Generator:
    def __init__(self):
        self.substitutions_py = {}
        self.substitutions_cpp = {}

    def create_substitutions(self, model):
        # create the substitutions with a given model

        self.substitutions_py = {}
        self.substitutions_py["model_name"] = model.pymio_name

        self.substitutions_cpp = {}
        self.substitutions_cpp["model_class"] = model.name
        self.substitutions_cpp["model_cpp_name"] = model.name.split("::")[-1]
        self.substitutions_cpp["pymio_name"] = model.pymio_name
        self.substitutions_cpp["simulation_class"] = model.simulation_name
        self.substitutions_cpp["namespace"] = ("::".join(model.name.split("::")[:-1]) + "::")
        self.substitutions_cpp["includes"] = self.include_string(model)
        self.substitutions_cpp["pretty_name"] = self.pretty_name_string(model)
        self.substitutions_cpp["population_template_arguments"] = self.population_string(model)
        self.substitutions_cpp["indexing"] = self.indexing_binds_string(model)
        self.substitutions_cpp["model_init"] = self.model_init(model)
        self.substitutions_cpp["enum_population_states"] = self.enum(model)

    def generate_files(self, model):
        with open("./generator/template/template_py.txt") as t:
            template_py = string.Template(t.read())
        with open("./generator/template/template_cpp.txt") as t:
            template_cpp = string.Template(t.read())
        
        output_py = template_py.safe_substitute(**self.substitutions_py)
        output_cpp = template_cpp.safe_substitute(**self.substitutions_cpp)

        py_filename = model.pymio_name + ".py"
        cpp_filename = model.pymio_name + ".cpp"
        with open(py_filename, "w") as output:
            output.write(output_py)
        with open(cpp_filename, "w") as output:
            output.write(output_cpp)

    # helperfunctions to create strings for substitution
    @staticmethod
    def include_string(model):
        return (
            "#include \"templates.h\"\n"
            "#include \"ode_seir/model.h\"\n"
            "#include \"Eigen/Core\"\n"
            "#include \"pybind11/stl_bind.h\"\n"
            "#include <vector>\n"
        ).format(
        )

    @staticmethod
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

    @staticmethod
    def population_string(model):
        str = "<"

        for group in model.population_groups:
            str += (
                "{group}, "
            ).format(
                group = group
            )
        return str[:-2] + ">"

    @staticmethod
    def indexing_binds_string(model):
        for key in model.enum_dict:

            if key.split("::")[-1] not in model.population_groups:
                continue

            return (
                "pymio::bind_Index<{enum_class}>(m, \"Index_{enum_name}\");\n\t"
                "pymio::bind_CustomIndexArray<mio::UncertainValue, {enum_class}>(m, \"PopulationArray\");\n"
            ).format(
                enum_class=key,
                enum_name=key.split("::")[-1]
            )

    @staticmethod
    def model_init(model):
        str = ""
        for init in model.init:

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

    @staticmethod
    def enum(model):
        str = ""
        for key, values in model.enum_dict.items():

            if key.split("::")[-1] not in model.population_groups:
                continue

            str += (
                "pymio::iterable_enum<{enum_class}>(m, \"{enum_name}\")\n\t"
                ).format(
                    enum_class=key,
                    enum_name=key.split("::")[-1]
                )

            for value in values:
                str += (
                    "    .value(\"{comp_name}\", {comp_class})\n\t"
                ).format(
                    comp_class = value,
                    comp_name = value.split("::")[-1]
                )
            str = str.rstrip() + ";\n\n"
        return str