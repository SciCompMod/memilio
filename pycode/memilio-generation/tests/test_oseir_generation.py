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

import os
import sys
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from memilio.generation import Generator, Scanner, ScannerConfig, AST

if sys.version_info >= (3, 9):
    # For python 3.9 and newer
    import importlib.resources as importlib_resources
else:
    # For older python versions
    import importlib_resources


class TestOseirGeneration(unittest.TestCase):
    """ """
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
    build_dir = tempfile.TemporaryDirectory(dir=test_dir.name)

    # run cmake
    cmake_cmd = ["cmake", os.path.join(
        project_path, "pycode/memilio-generation"), "-Wno-dev"]

    cmake_cmd_result = subprocess.run(
        cmake_cmd, stdout=subprocess.PIPE, cwd=build_dir.name)
    cmake_cmd_result.check_returncode()

    @patch('memilio.generation.scanner.utility.try_get_compilation_database_path')
    def setUp(self, try_get_compilation_database_path_mock):
        """

        :param try_get_compilation_database_path_mock: 

        """
        try_get_compilation_database_path_mock.return_value = self.build_dir.name
        config_json = {

            "python_generation_module_path": self.project_path + "/pycode/memilio-generation",
            "skbuild_path_to_database": "",
            "libclang_library_path": ""

        }

        conf = ScannerConfig.from_dict(config_json)

        conf.source_file = self.project_path + \
            "/cpp/models/ode_seir/model.cpp"

        # Could be any target folder
        conf.target_folder = self.test_dir.name
        self.scanner = Scanner(conf)
        self.scanner.python_module_name = "test_oseir"
        self.ast = AST(conf)

    @unittest.skip("Skip test until resolved")
    def test_clean_oseir(self):
        """ """
        irdata = self.scanner.extract_results(self.ast.root_cursor)

        generator = Generator()
        generator.create_substitutions(irdata)
        generator.generate_files(irdata)

        self.maxDiff = None

        with open(os.path.join(irdata.target_folder, "test_oseir.py")) as result:
            self.assertEqual(result.read(), self.expected_test_oseir_py)
        with open(os.path.join(irdata.target_folder, "test_oseir.cpp")) as result:
            self.assertEqual(result.read(), self.expected_test_oseir_cpp)

    def test_wrong_model_name(self):
        """ """
        self.scanner.config.model_class = "wrong_name"
        with self.assertRaises(AssertionError) as error:
            self.scanner.extract_results(self.ast.root_cursor)

        error_message = "set a model name"
        self.assertEqual(str(error.exception), error_message)


if __name__ == '__main__':
    unittest.main()
