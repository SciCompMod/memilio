from dataclasses import dataclass, field

@dataclass
class Model:
    name                : str   = None
    python_module_name  : str   = None
    init                : list  = field(default_factory=list)
    namespace           : str   = None
    enum_dict           : dict  = field(default_factory=dict)
    population_groups   : list  = field(default_factory=list)
    parameterset        : str   = None
    compartments        : list  = field(default_factory=list)
    simulation_name     : str   = None

    def set_attribute(self, attribute_name, value):
        self.__setattr__(attribute_name, value)
