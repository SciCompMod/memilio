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
        """
        Print all attributes of object
        """
        out = ""
        for key, value in self.__dict__.items():
            out += (key + ": " + str(value) + "\n")
        return out

    def finalize(self, conf):
        """
        Finalize the input of the model. Call before using model to generate code.
        """
        self.pymio_name = conf.python_module_name
        if self.pymio_name is None:
            self.pymio_name = self.name.split("::")[-1]
        
        self.simulation_name = conf.simulation_name
        if self.simulation_name is None:
            self.simulation_name = "mio::Simulation<" + self.name + ">"

        #assert(self.name != None), "set a model name using generator.name()"
        #assert(self.namespace!= None), "set a model name using generator.namespace()"
        #assert(len(self.compartments) != 0), "add compartments using generator.compartments"