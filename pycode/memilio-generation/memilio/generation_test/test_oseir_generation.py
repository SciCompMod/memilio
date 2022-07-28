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

class TestOseirGeneration(unittest.TestCase):

    def setUp(self):
        config_json =   {
                        "libclang_library_path": "/tools/modulesystem/spack-22.1/opt/spack/linux-ubuntu20.04-haswell/gcc-9.4.0/llvm-13.0.0-e4gcooqwzwwiiprr7u7nnorvgrqqhl63/lib/libclang.so",
                        "source_file": "../../cpp/models/ode_seir/model.cpp",
                        "path_database": "/build",
                        "namespace": "mio::oseir::",
                        "optional": {
                            "python_module_name": "generated_oseir",
                            "simulation_name": "",
                            "age_group": False,
                            "target_folder": "/pycode/memilio-generation/memilio/generation_test"
                            }
                        }
        conf = ScannerConfig.from_dict(config_json)
        self.scanner = Scanner(conf)
        
    def test_clean_oseir(self):
        model = self.scanner.extract_results()

        generator = Generator()
        generator.create_substitutions(model)
        generator.generate_files(model)

        with open(model.target_folder + "/generated_oseir.py") as result:
            with open(model.target_folder + "/test_oseir.py.txt") as expected:
                self.assertMultiLineEqual(result.read(), expected.read())
        with open(model.target_folder + "/generated_oseir.cpp") as result:
            with open(model.target_folder + "/test_oseir.cpp.txt") as expected:
                self.assertMultiLineEqual(result.read(), expected.read())

    def test_wrong_model_name(self):
        self.scanner.config.model_class = "wrong_name"
        with self.assertRaises(AssertionError) as error:
            self.scanner.extract_results()
    
        error_message = "set a model name"
        self.assertEqual(str(error.exception), error_message)



if __name__ == '__main__':
    unittest.main()
