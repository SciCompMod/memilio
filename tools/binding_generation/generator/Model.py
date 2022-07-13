
class Model:
    def __init__(self):
        self.name               = None
        self.python_module_name = None
        self.init               = []
        self.namespace          = None
        self.enum_dict          = {}
        self.population_groups  = []
        self.compartments       = []
        self.simulation_name    = None

    def set_attribute(self, attribute_name, value):
        self.__setattr__(attribute_name, value)

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