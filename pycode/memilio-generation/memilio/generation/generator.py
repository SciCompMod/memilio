import string
import os
import memilio.generation.template.template_string as StringTemplates

class Generator:
    def __init__(self):
        self.substitutions_py = {}
        self.substitutions_cpp = {}

    def create_substitutions(self, intermed_repr):
        # create the substitutions with a given intermed_repr

        self.substitutions_py = {}
        self.substitutions_py["python_module_name"] = intermed_repr.python_module_name

        self.substitutions_cpp = {}
    
        self.substitutions_cpp["namespace"] = intermed_repr.namespace
        self.substitutions_cpp["model_class"] = intermed_repr.model_class
        self.substitutions_cpp["model_base"] = intermed_repr.model_base[0]
        self.substitutions_cpp["model_base_templates"] = intermed_repr.model_base[1][0] + ", " + intermed_repr.model_base[2][0] + ", " + intermed_repr.model_base[3][0]
        self.substitutions_cpp["python_module_name"] = intermed_repr.python_module_name
        self.substitutions_cpp["parameterset"] = intermed_repr.parameterset
        
        self.substitutions_cpp["includes"] = StringTemplates.include_string(intermed_repr)
        self.substitutions_cpp["pretty_name"] = StringTemplates.pretty_name_string(intermed_repr)
        self.substitutions_cpp["enum_populations"] = StringTemplates.enum_populations(intermed_repr)
        self.substitutions_cpp["model_init"] = StringTemplates.model_init(intermed_repr)
        self.substitutions_cpp["population"] = StringTemplates.population_string(intermed_repr)

        # optionals
        self.substitutions_cpp["parameterset_indexing"] = StringTemplates.parameterset_indexing(intermed_repr)
        self.substitutions_cpp["parameterset_wrapper"] = StringTemplates.parameterset_wrapper(intermed_repr)
        self.substitutions_cpp["age_group"] = StringTemplates.age_group(intermed_repr)

        self.substitutions_cpp["simulation"] = StringTemplates.simulation(intermed_repr)
        self.substitutions_cpp["simulation_graph"] = StringTemplates.simulation_graph(intermed_repr)
        self.substitutions_cpp["simulation_vector_definition"] = StringTemplates.simulation_vector_definition(intermed_repr)


    def generate_files(self, intermed_repr):
        with open(os.path.join(intermed_repr.project_path + "/pycode/memilio-generation/memilio/generation/template/template_py.txt")) as t:
            template_py = string.Template(t.read())
        with open(os.path.join(intermed_repr.project_path + "/pycode/memilio-generation/memilio/generation/template/template_cpp.txt")) as t:
            template_cpp = string.Template(t.read())
        
        output_py = template_py.safe_substitute(**self.substitutions_py)
        output_cpp = template_cpp.safe_substitute(**self.substitutions_cpp)

        py_filename = "/" + intermed_repr.python_module_name + ".py"
        cpp_filename = "/" + intermed_repr.python_module_name + ".cpp"
        with open(os.path.join(intermed_repr.target_folder + py_filename), "w") as output:
            output.write(output_py)
        with open(os.path.join(intermed_repr.target_folder + cpp_filename), "w") as output:
            output.write(output_cpp)
