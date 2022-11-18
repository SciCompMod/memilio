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
Example for the ode seir model.
"""
from memilio.generation import ScannerConfig, Generator, Scanner
import os


here = os.path.dirname(os.path.abspath(__file__))

# Define ScannerConfig and initialize Scanner
with open(os.path.join(here + '/config.json'), 'r') as file:
    conf = ScannerConfig.schema().loads(file.read(), many=True)[0]
scanner = Scanner(conf)

# Extract results of Scanner into a intermed_repr
intermed_repr = scanner.extract_results()

# Generate code
generator = Generator()
generator.create_substitutions(intermed_repr)
generator.generate_files(intermed_repr)

# scanner.output_ast_file()
