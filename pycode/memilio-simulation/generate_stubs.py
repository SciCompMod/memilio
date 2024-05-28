import os
import subprocess
import sys
import importlib.util

# all current models
models = ["osir", "oseir", "secir", "osecirvvs", "abm"]
setup_content = f"""
from setuptools import setup, find_packages

setup(
    name='memilio-stubs',
    version='0.1',
    packages=['memilio-stubs'],
    package_data={{
        'memilio-stubs/simulation': ['*.pyi'],
    }},
)
"""

if __name__ == "__main__":

    python_interpreter = sys.executable

    # Check for needed packages. If it fails either pacakge is not installed or the wron python interpreter is detected.
    # For later try setting python_interpreter with full path
    if importlib.util.find_spec('pybind11_stubgen') is None:
        print('pybind11_stubgen is not installed')
        exit()
    if importlib.util.find_spec('memilio.simulation') is None:
        print('memilio.simulation is not installed')
        exit()

    file_path = os.path.dirname(os.path.abspath(__file__))
    package_dir = os.path.abspath(os.path.join(
        file_path, "../memilio-simulation-stubs"))
    output_dir = os.path.join(package_dir, "memilio-stubs/simulation")
    output_module_dir = os.path.join(output_dir, 'memilio')

    # create folders, if they do not exist
    try:
        os.makedirs(output_dir)
    except:
        pass

    # generate stubs and moce them into correct folder with right name
    # memilio-stubs/simulation module needs same structure as memilio/simulation
    subprocess.check_call(
        [python_interpreter, '-m', 'pybind11_stubgen', '-o', output_dir, 'memilio._simulation'])
    os.rename(os.path.join(output_module_dir, '_simulation.pyi'),
              os.path.join(output_dir, '__init__.pyi'))

    for model in models:
        module_name = "memilio._simulation_" + model
        subprocess.check_call(
            [python_interpreter, '-m', 'pybind11_stubgen', '-o', output_dir, module_name])
        os.rename(os.path.join(output_module_dir, '_simulation_' + model + '.pyi'),
                  os.path.join(output_dir, model + '.pyi'))

    os.rmdir(output_module_dir)

    # create setup.py and install package
    with open(os.path.join(package_dir, "setup.py"), "w") as setup_file:
        setup_file.write(setup_content)
    subprocess.check_call(
        [python_interpreter, '-m', 'pip', 'install', package_dir])
