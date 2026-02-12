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
import importlib.util
import os
import shutil
import subprocess
import sys

pyproject_content_protected_module = """[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "memilio-stubs"
version = "0.1"
description = "Stubs for the memilio.simulation package."
requires-python = ">=3.8"
dependencies = []

[tool.setuptools.packages.find]
where = ["."]
include = ["memilio-stubs"]

[tool.setuptools.package-data]
"memilio-stubs" = ["simulation/*.pyi"]
"""

# Define if the generated stubs of mypy should be configured
# For testing generate stubs with flag_configure_generated_stubs=False and run mypy on codebase
flag_configure_generated_stubs = True

if __name__ == "__main__":

    python_interpreter = sys.executable

    # Check for needed packages. If it fails either pacakge is not installed or the wrong python interpreter is detected.
    # For the latter try setting python_interpreter with full path
    if importlib.util.find_spec('mypy') is None:
        print('pybind11_stubgen is not installed')
        exit()
    if importlib.util.find_spec('memilio.simulation') is None:
        print('memilio.simulation is not installed')
        exit()

    file_path = os.path.dirname(os.path.abspath(__file__))
    package_dir = os.path.abspath(os.path.join(
        file_path, "../../memilio-simulation-stubs"))

    # delete stubs if they already exist
    try:
        shutil.rmtree(os.path.join(package_dir, "memilio-stubs"))
    except:
        pass
    # create folders, if they do not exist
    try:
        os.makedirs(package_dir)
    except:
        pass

    # memilio-stubs/simulation module needs same structure as memilio/simulation
    # test --include-docstrings, --doc-dir PATH for better docs
    subprocess.check_call(
        ['stubgen', '--include-docstrings', '-o', package_dir, '-p',
         'memilio.simulation'])

    # TODO:
    #   - fix numpy.float64[m, 1] to
    #       1. numpy.float64[M, Literal(1)] with M = typing.TypeVar("M", bound=int) at top of file
    #       2. or numpy.ndarray[numpy.float64[m, 1]] to numpy.ndarray[tuple[M, typing.Literal[1]], numpy.dtype[numpy.float64]]
    #   - additionaly for numpy.ndarray[numpy.float64[m, n], flags.f_contiguous] n needs to be changed similar to m
    #   - change following error:
    #        def __eq__(self, other: object) -> bool:
    #            if not isinstance(other, Index_SimulationDay):
    #                return NotImplemented
    #            return <logic to compare two Index_SimulationDay instances>
    #   - Default values are not defined

    if flag_configure_generated_stubs:

        # get all model modules from memilio.simulation
        # if package structure changes this needs to be adjusted
        models = [f.rstrip(".py") for f in os.listdir(os.path.join(
            file_path, "../memilio/simulation")) if f.endswith(".py")]
        models.remove("__init__")
        stub_files_path = os.path.join(package_dir, "memilio/simulation")

        # read .pyi from protected simulation module
        protected_module_file_path = os.path.join(
            stub_files_path, "_simulation.pyi")
        with open(protected_module_file_path, encoding='utf-8') as file:
            content_protected_module = file.read()

        # Replace the protected namespaces with the desired ones
        content_protected_module = content_protected_module.replace(
            "simulation._simulation", "simulation")

        # read .pyi from pure python simulation file
        module_file_path = os.path.join(stub_files_path, "__init__.pyi")
        with open(module_file_path, encoding='utf-8') as file:
            content_module = file.read()

        content_module = content_module.replace(
            "from memilio.simulation._simulation import *" + os.linesep, "")

        import_submodules_string = ["from memilio.simulation import ("]
        import_submodules_string += [(f"    {model} as {model},")
                                     for model in models]
        import_submodules_string.append(")")
        import_submodules_string = (os.linesep).join(import_submodules_string)

        content_protected_module = import_submodules_string + os.linesep + \
            content_protected_module + os.linesep + content_module

        # Write the modified content_protected_module back to the file
        with open(module_file_path, 'w', encoding='utf-8') as file:
            file.write(content_protected_module)

        # remove the protected module file
        os.remove(protected_module_file_path)

        for model in models:

            # read .pyi from protected model module
            protected_module_file_path = os.path.join(
                stub_files_path, "_simulation_" + model + ".pyi")
            with open(protected_module_file_path, encoding='utf-8') as file:
                content_protected_module = file.read()

            # Replace the protected namespaces with the desired ones
            content_protected_module = content_protected_module.replace(
                "_simulation_" + model, model)
            content_protected_module = content_protected_module.replace(
                "simulation._simulation", "simulation")

            # read .pyi from pure python code
            module_file_path = os.path.join(stub_files_path, model + ".pyi")
            with open(module_file_path, encoding='utf-8') as file:
                content_module = file.read()

            content_module = content_module.replace(
                "from memilio.simulation._simulation_" + model + " import *" + os.linesep, "")

            content_protected_module = content_protected_module + os.linesep + content_module

            # Write the modified content_protected_module back to the file
            with open(module_file_path, 'w', encoding='utf-8') as file:
                file.write(content_protected_module)

            # remove the protected module file
            os.remove(protected_module_file_path)

    # rename directory memilio to memilio-stubs
    shutil.move(os.path.join(package_dir, "memilio"),
                os.path.join(package_dir, "memilio-stubs"))

    # create pyproject.toml and install package
    with open(os.path.join(package_dir, "pyproject.toml"), "w", encoding="utf-8") as pyproject_file:
        pyproject_file.write(pyproject_content_protected_module)
    subprocess.check_call(
        [python_interpreter, '-m', 'pip', 'install', package_dir])
