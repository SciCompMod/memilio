from dataclasses import dataclass, field
from dataclasses_json import dataclass_json
from subprocess import check_output

from pyparsing import Optional

@dataclass_json
@dataclass
class ScannerConfig:
    libclang_library_path   : str
    source_file             : str
    path_database           : str
    namespace               : str
    optional                : dict = field(default_factory=dict)
    model_name              : str  = field(init = False)
    main_file               : str  = field(init = False)
    parameterset_name       : str  = field(init = False)
    project_path            : str  = field(init = False)

    def __post_init__(self):
        # Predefined Variables
        self.model_name         = "Model"
        self.parameterset_name  = "Parameters"
        self.project_path       = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1]



