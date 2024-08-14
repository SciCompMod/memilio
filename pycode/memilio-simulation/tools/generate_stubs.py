import os
import subprocess
import sys
import importlib.util
import shutil

setup_content = f"""
from setuptools import setup, find_packages

setup(
    name='memilio-stubs',
    version='0.1',
    packages=['memilio-stubs'],
    package_data={{
        'memilio-stubs': ['simulation/*.pyi'],
    }},
)
"""

if __name__ == "__main__":

    python_interpreter = sys.executable

    # Check for needed packages. If it fails either pacakge is not installed or the wrong python interpreter is detected.
    # For the latter try setting python_interpreter with full path
    if importlib.util.find_spec('pybind11_stubgen') is None:
        print('pybind11_stubgen is not installed')
        exit()
    if importlib.util.find_spec('memilio.simulation') is None:
        print('memilio.simulation is not installed')
        exit()

    file_path = os.path.dirname(os.path.abspath(__file__))
    package_dir = os.path.abspath(os.path.join(
        file_path, "../../memilio-simulation-stubs"))

    # create folders, if they do not exist
    try:
        os.makedirs(package_dir)
    except:
        pass

    # delete stubs if they already exist
    try:
        shutil.rmtree(os.path.join(package_dir, "memilio-stubs"))
    except:
        pass

    # memilio-stubs/simulation module needs same structure as memilio/simulation
    # can change between [--numpy-array-wrap-with-annotated|--numpy-array-use-type-var|--numpy-array-remove-parameters|] for differend representations of numpy arrays
    subprocess.check_call(
        [python_interpreter, '-m', 'pybind11_stubgen', '--ignore-all-errors', '--root-suffix=-stubs', '--numpy-array-use-type-var', '-o', package_dir, 'memilio._simulation'])

    # rename all instances of memilio._simulation to memilio.simulation in .pyi files
    for root, _, files in os.walk(os.path.join(package_dir, "memilio-stubs/_simulation")):
        for file in files:
            if file.endswith(".pyi"):  # Process only Python files
                file_path = os.path.join(root, file)
                with open(file_path, encoding='utf-8') as file:
                    content = file.read()

                # Replace the old namespace with the new one
                new_content = content.replace(
                    "memilio._simulation", "memilio.simulation")

                # Write the modified content back to the file
                with open(file_path, 'w', encoding='utf-8') as file:
                    file.write(new_content)

    # rename directory _simulation to simulation
    shutil.move(os.path.join(package_dir, "memilio-stubs/_simulation"),
                os.path.join(package_dir, "memilio-stubs/simulation"))

    # create setup.py and install package
    with open(os.path.join(package_dir, "setup.py"), "w") as setup_file:
        setup_file.write(setup_content)
    subprocess.check_call(
        [python_interpreter, '-m', 'pip', 'install', package_dir])
