from dataclasses import dataclass, field

@dataclass
class Model:
    namespace           : str   = None
    model_class         : str   = None
    python_module_name  : str   = None
    parameterset        : str   = None
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
