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