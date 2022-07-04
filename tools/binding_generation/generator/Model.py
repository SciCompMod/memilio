import generator.utility
from sympy import symbols as Symbols

class Model:
    def __init__(self):
        self.name               = None
        self.pymio_name         = None
        self.init               = []
        self.namespace          = None
        self.enum_dict          = {}
        self.population_groups  = []
        self.compartments       = []

    def set_name(self, name):
        err = "model name must be a string"
        if type(name) is not str: raise TypeError(err)
        self.name = name
        return
        
    def set_namespace(self, namespace):
        err = "namespace must be a string"
        if type(namespace) is not str: raise TypeError(err)
        self.namespace = namespace
        return

    def __str__(self):
        return(
            "class Model()\n"
            "name: {name}, enum_dict: {enum}\n"
            "init: {init}, population_groups: {pop_group}\n"
            "compartments: {comp}"
            ).format(
                name = self.name,
                init = self.init,
                enum = self.enum_dict,
                pop_group = self.population_groups,
                comp = self.compartments
            )

    def finalize(self, conf):
        if conf.model_name is None:
            self.pymio_name = self.name
        self.pymio_name = conf.model_name
        #assert(self.name != None), "set a model name using generator.name()"
        #assert(self.namespace!= None), "set a model name using generator.namespace()"
        #assert(len(self.compartments) != 0), "add compartments using generator.compartments"