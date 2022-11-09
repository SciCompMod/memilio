#############################################################################
# Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
#
# Authors: Maximilian Betz
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
"""
@file scanner_config.py
@brief Dataclass to import the configurations from the config.json.
"""
import os
from dataclasses import dataclass, field
from typing_extensions import Self
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

    def __post_init__(self: Self) -> None:
        # Predefined Variables
        self.model_class            = "Model"
        self.parameterset           = "Parameters"
        # Get the project path. If this is used outside of Memilio it needs to be changed.
        project_directory_path      = os.path.dirname(os.path.abspath(__file__))
        self.project_path           = project_directory_path.split('/pycode')[0]
        self.target_folder          = self.project_path + self.optional.get("target_folder", "")
