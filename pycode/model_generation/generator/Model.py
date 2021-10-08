import generator.utility
import generator.ModelManager
from sympy import symbols as Symbols

class Model:
    def __init__(self, num):
        self.num                = num
        self.name               = None
        self.namespace          = None
        self.compartment_name   = "Compartments"
        self.population_name    = "Populations"
        self.parameterset_name  = "Parameters"
        self.parameters         = []
        self.compartments       = []
        self.equations          = []
        self.variables          = []

    def add_parameter(self, name, c_type, default_value):
        generator.ModelManager.set_active(self)
        symbol = Symbols(name)
        self.parameters += [generator.utility.Parameter(name, c_type, default_value, symbol)]
        return symbol

    def set_compartments(self, *args):
        assert(len(args) > 0)
        err = "compartments for this model are already set and cannot be overwritten"
        if len(self.compartments) > 0: raise RuntimeError(err)
        generator.ModelManager.set_active(self)
        names = []
        for arg in args:
            names += arg.split()
        symbols = Symbols(names)
        if len(names) == 1:
            self.compartments = [generator.utility.Variable(names[0], symbols)]
        else:
            self.compartments = [generator.utility.Variable(n, s) for n, s in zip(names, symbols)]
        return symbols

    def set_name(self, name):
        generator.ModelManager.set_active(self)
        err = "model name must be a string"
        if type(name) is not str: raise TypeError(err)
        self.name = name
        return
        
    def set_namespace(self, namespace):
        generator.ModelManager.set_active(self)
        err = "namespace must be a string"
        if type(namespace) is not str: raise TypeError(err)
        self.namespace = namespace
        return

    def add_equation(self, lhs, op, rhs):
        generator.ModelManager.set_active(self)
        self.equations += [generator.utility.assign(lhs, rhs, op)]
        return

    def add_variable(self, *args):
        assert(len(args) > 0)
        generator.ModelManager.set_active(self)
        names = []
        for arg in args:
            names += arg.split()
        symbols = Symbols(names)
        self.variables += [generator.utility.Variable(n, s) for n, s in zip(names, symbols)]
        return symbols if len(names) > 1 else symbols[0]

    def finalize(self):
        assert(self.name != None), "set a model name using generator.name()"
        assert(self.namespace!= None), "set a model name using generator.namespace()"
        assert(len(self.compartments) != 0), "add compartments using generator.compartments"
        if len(self.parameters) == 0:
            self.parameters = [generator.utility.Parameter("dummy", "char", 0, None)]