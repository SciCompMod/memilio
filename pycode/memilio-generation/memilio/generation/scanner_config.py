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
from dataclasses_json import dataclass_json
from typing_extensions import Self


@dataclass_json
@dataclass
class ScannerConfig:
    """
    Provide configurations from JSON-file in Python as dataclass.

    Attributes (and config.json parameters):
        source_file: Path to the main file of the mode e.g., model.cpp
        path_database: Path to the folder of the compile_commands.json
        namespace: C++ namespace of the model class
        python_module_name: Individual name for binded python module
        optional: List with optional arguments
            libclang_library_path: Path to the local libclang library, if empty string path gets obtained with terminal command
            simulation_class: Name of simulation class, if not used set as empty string
            age_group: Boolean defining if model uses age groups
            parameterset_wrapper": Boelean defining if model uses wrapper for parameterset
            target_folder: Target folder for generated files
    """
    source_file: str
    path_database: str
    namespace: str
    python_module_name: str
    optional: dict = field(default_factory=dict)
    model_class: str = field(init=False)
    parameterset: str = field(init=False)
    project_path: str = field(init=False)
    target_folder: str = field(init=False)

    def __post_init__(self: Self) -> None:
        """
        Additional data for the generation, that is not set in the config.json. 
        Typically these should not be changed by the user. __post_init__() is automatically called after the Constructor.
        """
        # Predefined Variables
        self.model_class = "Model"
        self.parameterset = "Parameters"

        # Get the project path. If this is used outside of Memilio it needs to be changed.
        project_directory_path = os.path.dirname(os.path.abspath(__file__))
        self.project_path = project_directory_path.split('/pycode')[0]
        self.target_folder = self.project_path + \
            self.optional.get("target_folder", "")
