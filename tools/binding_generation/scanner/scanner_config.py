class ScannerConfig:
    def __init__(self):
        self.python_module_name         = None
        self.folder             = None
        self.simulation_name    = None

    def set_python_module_name(self, python_module_name):
        """
        Set the name of the python module.
        """
        err = "module name must be a string"
        if type(python_module_name) is not str: raise TypeError(err)
        self.python_module_name = python_module_name
        return
        
    def set_folder(self, folder):
        """
        Set the folder with the source code for the cpp model.
        """
        err = "namespace must be a string"
        if type(folder) is not str: raise TypeError(err)
        self.folder = folder
        return
    
    def set_simulation_name(self, simulation_name):
        """
        Simulation name is the name of the class used for simulating a model with its namespace.
        """
        err = "simulation name must be a string"
        if type(simulation_name) is not str: raise TypeError(err)
        self.simulation_name = simulation_name
        return
