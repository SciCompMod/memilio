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
import shutil
from dataclasses import dataclass
from dataclasses_json import dataclass_json
from importlib.resources import files

@dataclass_json
@dataclass
class ScannerConfig:
    """Provide configurations from JSON-file in Python as dataclass.

    Attributes (and config.json parameters):
        libclang_library_path: Path to the local libclang library. If the string is empty, the program tries to retrieve the path from a terminal command.
        clang_executable: Command for clang executable.
    """
    
    libclang_library_path: str
    clang_executable: str

    def __post_init__(self) -> None:
        # Check libclang_library_path correctly set
        if self.libclang_library_path.endswith("-NOTFOUND"):
            # if not found derive the library from clang_executable
            self.libclang_library_path = ""

        if not self.clang_executable or self.clang_executable.endswith("-NOTFOUND"):
            raise ValueError(
                "No clang executable configured. Set MEMILIO_CLANG_VERSION in CMake "
                "or clang_executable in config.json.")
        if shutil.which(self.clang_executable) is None:
            raise FileNotFoundError(
                f"clang executable not found: {self.clang_executable!r}")
        
def get_config():
    text = files("memilio.generation").joinpath("config.json").read_text(encoding="utf-8")
    conf = ScannerConfig.from_json(text)
    return conf

config = get_config()
