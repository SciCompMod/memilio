
from sympy import symbols as Symbols
from sympy.printing import ccode
from generator.utility import join, assign

def print_python_file(model, file=None):
    s = python_file_string(model)
    if file is None:
        print(s)
    else:
        with open(file, 'w') as file:
            file.write(s)
    return

def python_file_string(model):
    return (
        "#############################################################################\n"
        "# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)\n"
        "#\n"
        "# Authors: Daniel Abele\n"
        "#\n"
        "# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>\n"
        "#\n"
        "# Licensed under the Apache License, Version 2.0 (the \"License\");\n"        
        "# you may not use this file except in compliance with the License.\n"
        "# You may obtain a copy of the License at\n"
        "#\n"
        "#     http://www.apache.org/licenses/LICENSE-2.0\n"
        "#\n"
        "# Unless required by applicable law or agreed to in writing, software\n"
        "# distributed under the License is distributed on an \"AS IS\" BASIS,\n"
        "# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.\n"
        "# See the License for the specific language governing permissions and\n"
        "# limitations under the License.\n"
        "#############################################################################\n"
        "\n"
        "\"\"\"\n"
        "Python bindings for MEmilio {model_name} model.\n"
        "\"\"\"\n"
        "\n"
        "from memilio._simulation_{model_name} import *\n"
    ).format(
        model_name      = model.pymio_name
    )