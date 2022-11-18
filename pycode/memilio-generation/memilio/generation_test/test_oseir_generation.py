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

from memilio.generation import Scanner, ScannerConfig, Generator
import unittest
import subprocess
import tempfile

import os


class TestOseirGeneration(unittest.TestCase):

    # Get a file object with write permission.
    here = os.path.dirname(os.path.abspath(__file__))
    project_path = here.split('/pycode')[0]

    # Load expected results for oseir generation.
    with open(os.path.join(here + "/test_data/test_oseir.py.txt"), 'r') as expected:
        expected_test_oseir_py = expected.read()

    with open(os.path.join(here + "/test_data/test_oseir.cpp.txt"), 'r') as expected:
        expected_test_oseir_cpp = expected.read()

    # Create a temporary directory
    test_dir = tempfile.TemporaryDirectory(dir=project_path)

    # Create compile_commands
    build_path = os.path.join(test_dir.name + '/build')
    source_path = os.path.join(project_path + '/cpp')
    clang_cmd = ["cmake", '-S', source_path, '-B',
                 build_path, '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON']
    subprocess.run(clang_cmd)

    def setUp(self):
        path_database = self.build_path.split('memilio')[-1]
        config_json = {
            "source_file": "/cpp/models/ode_seir/model.cpp",
            "path_database": path_database,
            "namespace": "mio::oseir::",
            "python_module_name": "test_oseir",
            "optional": {
                "libclang_library_path": "",
                "simulation_name": "",
                "age_group": False
            }
        }

        conf = ScannerConfig.from_dict(config_json)
        conf.project_path = self.project_path
        conf.target_folder = os.path.join(self.project_path + path_database)
        self.scanner = Scanner(conf)

    def test_clean_oseir(self):
        irdata = self.scanner.extract_results()

        generator = Generator()
        generator.create_substitutions(irdata)
        generator.generate_files(irdata)

        with open(os.path.join(irdata.target_folder + "/test_oseir.py")) as result:
            self.assertEqual(result.read(), self.expected_test_oseir_py)
        with open(os.path.join(irdata.target_folder + "/test_oseir.cpp")) as result:
            self.assertEqual(result.read(), self.expected_test_oseir_cpp)

    def test_wrong_model_name(self):
        self.scanner.config.model_class = "wrong_name"
        with self.assertRaises(AssertionError) as error:
            self.scanner.extract_results()

        error_message = "set a model name"
        self.assertEqual(str(error.exception), error_message)


if __name__ == '__main__':
    unittest.main()
