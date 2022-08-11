import os
from dataclasses import dataclass, field
from dataclasses_json import dataclass_json

@dataclass_json
@dataclass
class ScannerConfig:
    source_file             : str
    path_database           : str
    namespace               : str
    optional                : dict = field(default_factory=dict)
    model_class             : str  = field(init = False)
    parameterset            : str  = field(init = False)
    project_path            : str  = field(init = False)
    target_folder           : str  = field(init = False)

    def __post_init__(self):
        # Predefined Variables
        #self.model_class            = "SecirModel"
        #self.parameterset           = "SecirParamsBase"
        self.model_class            = "Model"
        self.parameterset           = "Parameters"
        # Get the project path. If this is used outside of Memilio it needs to be changed.
        project_directory_path      = os.path.dirname(os.path.abspath(__file__))
        self.project_path           = project_directory_path.split('pycode')[0].rstrip("/")
        self.target_folder          = self.project_path + self.optional.get("target_folder", "")
