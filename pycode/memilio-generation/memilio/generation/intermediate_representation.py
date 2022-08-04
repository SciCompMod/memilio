from dataclasses import dataclass, field

@dataclass
class IntermediateRepresentation:
    namespace           : str   = None
    model_class         : str   = None
    python_module_name  : str   = None
    parameterset        : str   = None
    parameterset_wrapper: str   = None
    simulation_class    : str   = None
    project_path        : str   = None
    target_folder       : str   = None
    enum_populations    : dict  = field(default_factory=dict)
    model_init          : list  = field(default_factory=list)
    model_base          : list  = field(default_factory=list)
    population_groups   : list  = field(default_factory=list)
    age_group           : dict  = field(default_factory=dict)

    def set_attribute(self, attribute_name, value):
        self.__setattr__(attribute_name, value)
    
    def check_complete_data(self, optional):
        assert(self.model_class != None), "set a model name"
        assert(self.namespace != None), "set a model name_space"
        assert(self.parameterset != None), "set a parameterset"
        if optional.get("parameterset_wrapper"):
            assert(self.parameterset_wrapper != None), "No Parameterset_Wrapper found. If None is used in this model set parameterset_wrapper to false in config.json"
        if optional.get("age_group"):
            assert(self.age_group != {}), "No AgeGroup found. If None is used in this model set age_group to false in config.json"
