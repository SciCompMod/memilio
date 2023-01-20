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

import json
import os
import subprocess
import tempfile
import unittest

import importlib_resources
from memilio.generation import Generator, Scanner, ScannerConfig


class TestOseirGeneration(unittest.TestCase):

    # Get a file object with write permission.
    here = os.path.dirname(os.path.abspath(__file__))
    project_path = here.split('/pycode')[0]

    # Load expected results for oseir generation.
    with open(os.path.join(here, "test_data/test_oseir.py.txt")) as expected:
        expected_test_oseir_py = expected.read()

    with open(os.path.join(here, "test_data/test_oseir.cpp.txt")) as expected:
        expected_test_oseir_cpp = expected.read()

    # Create a temporary directory
    test_dir = tempfile.TemporaryDirectory(dir=project_path)

    # load config.json
    pkg = importlib_resources.files("memilio")
    with importlib_resources.as_file(pkg.joinpath("tools/config.json")) as path:
        with open(path) as file:
            loaded_config_json = json.load(file)

    def setUp(self):
        config_json = {
            "source_file": self.loaded_config_json[0]['source_file'],
            "namespace": "mio::oseir::",
            "python_module_name": "test_oseir",
            "skbuild_path_to_database": self.loaded_config_json[0]['skbuild_path_to_database'],
            "python_generation_module_path": self.loaded_config_json[0]['python_generation_module_path'],
            "target_folder": self.test_dir.name,
            "optional": {
                "libclang_library_path": self.loaded_config_json[0]['optional']['libclang_library_path'],
                "simulation_name": "",
                "age_group": False,
                "parameterset_wrapper": True
            }
        }

        conf = ScannerConfig.from_dict(config_json)
        self.scanner = Scanner(conf)

    def test_clean_oseir(self):
        irdata = self.scanner.extract_results()

        generator = Generator()
        generator.create_substitutions(irdata)
        generator.generate_files(irdata)

        with open(os.path.join(irdata.target_folder, "test_oseir.py")) as result:
            self.assertEqual(result.read(), self.expected_test_oseir_py)
        with open(os.path.join(irdata.target_folder, "test_oseir.cpp")) as result:
            self.assertEqual(result.read(), self.expected_test_oseir_cpp)

    def test_wrong_model_name(self):
        self.scanner.config.model_class = "wrong_name"
        with self.assertRaises(AssertionError) as error:
            self.scanner.extract_results()

        error_message = "set a model name"
        self.assertEqual(str(error.exception), error_message)


if __name__ == '__main__':
    unittest.main()
