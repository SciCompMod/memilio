#############################################################################
# Copyright (C) 2020-2026 MEmilio
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
:strong:`scanner_config.py`
Dataclass to import the configurations from the config.json.
"""
from dataclasses import dataclass
from dataclasses_json import dataclass_json
import importlib.resources as importlib_resources

@dataclass_json
@dataclass
class ScannerConfig:
    """Provide configurations from JSON-file in Python as dataclass.

    Attributes (and config.json parameters):
        python_generation_module_path: Path to the ``pyproject.toml`` of the generation module
        path_to_database: Path to compile_commands.json
        libclang_library_path: Path to the local libclang library. If the string is empty, the program tries to retrieve the path from a terminal command.
    """

    python_generation_module_path: str
    path_to_database: str
    libclang_library_path: str

def get_config():
    pkg = importlib_resources.files("memilio.generation")
    with importlib_resources.as_file(
            pkg.joinpath("generation/config.json")) as path:
        with open(path) as file:
            conf = ScannerConfig.schema().loads(file.read(), many=True)[0]
    return conf

config = get_config()
