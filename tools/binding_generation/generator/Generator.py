from statistics import mode
import string
import generator.template.template_string as StringTemplates

class Generator:
    def __init__(self):
        self.substitutions_py = {}
        self.substitutions_cpp = {}

    def create_substitutions(self, model):
        # create the substitutions with a given model

        self.substitutions_py = {}
        self.substitutions_py["python_module_name"] = model.python_module_name

        self.substitutions_cpp = {}
    
        self.substitutions_cpp["namespace"] = model.namespace
        self.substitutions_cpp["model_class"] = model.model_class
        self.substitutions_cpp["model_base"] = model.model_base[0]
        self.substitutions_cpp["model_base_templates"] = model.model_base[1][0] + ", " + model.model_base[2][0]
        self.substitutions_cpp["python_module_name"] = model.python_module_name
        
        self.substitutions_cpp["includes"] = StringTemplates.include_string(model)
        self.substitutions_cpp["pretty_name"] = StringTemplates.pretty_name_string(model)
        self.substitutions_cpp["enum_populations"] = StringTemplates.enum_populations(model)
        self.substitutions_cpp["model_init"] = StringTemplates.model_init(model)
        self.substitutions_cpp["population"] = StringTemplates.population_string(model)

        # optionals
        self.substitutions_cpp["parameterset_wrapper"] = None
        self.substitutions_cpp["parameterset_indexing"] = None
        self.substitutions_cpp["age_group"] = StringTemplates.agegroup(model)

        self.substitutions_cpp["simulation"] = StringTemplates.simulation_graph(model)
        self.substitutions_cpp["simulation_graph"] = StringTemplates.simulation_graph(model)


    def generate_files(self, model):
        with open("./generator/template/template_py.txt") as t:
            template_py = string.Template(t.read())
        with open("./generator/template/template_cpp.txt") as t:
            template_cpp = string.Template(t.read())
        
        output_py = template_py.safe_substitute(**self.substitutions_py)
        output_cpp = template_cpp.safe_substitute(**self.substitutions_cpp)

        py_filename = model.python_module_name + ".py"
        cpp_filename = model.python_module_name + ".cpp"
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
    
    @staticmethod
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

    @staticmethod
    def population_string(model):
        for value in model.model_base[1:]:
            if "Population" in value[0]:
                return value[0]

    @staticmethod
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

    @staticmethod
    def agegroup(model):
        if not model.age_group:
            return ""
        return (
            "py::class_<mio::AgeGroup, {base}>(m, \"AgeGroup\").def(py::init<{init}>());"
        ).format(
            base = model.age_group["base"],
            init = model.age_group["init"][0]
        )
    
    @staticmethod
    def simulation(model):
        if model.simulation_class is None or (not model.simulation_class.strip()):
            return ""
        return (
            "pymio::bind_Simulation<{namespace}{simulation_class}<>>(m, \"{simulation_class}\");\n"
        ).format(
            namespace = model.namespace,
            simulation_class = model.simulation_class
        )

    @staticmethod
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