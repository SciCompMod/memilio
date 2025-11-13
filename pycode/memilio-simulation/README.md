# MEmilio Python Bindings

This package contains Python bindings for the MEmilio C++ library. It enables setting up and running simulations from Python code.

## Installation

This project is configured via ``pyproject.toml`` and is built with [scikit-build-core](https://scikit-build-core.readthedocs.io). CMake and Ninja must be available on the system. The package uses the [Pybind11 C++ library](https://pybind11.readthedocs.io) to create the bindings.

To install the package, use the command (from the directory containing ``pyproject.toml``)

```bash
pip install .
```

This builds the C++ library and C++ Python extension module and copies everything required to your site-packages.

All the requirements of the [C++ library](../../cpp/README.md) must be met in order to build and use the python bindings. A virtual environment is recommended. 

CMake is executed internally by scikit-build-core. All the options provided by the CMake configuration of the C++ library are available when building the Python extension as well. Additionally, the CMake configuration for the bindings provide the following CMake options:

- MEMILIO_USE_BUNDLED_PYBIND11: ON or OFF, default ON. If ON, downloads Pybind11 automatically from a repository during CMake configuration. If OFF, Pybind11 needs to be installed on the system.

When building the bindings, CMake options can be forwarded with configuration settings, e.g.

```bash
pip install . --config-settings=cmake.args="-DCMAKE_BUILD_TYPE=Debug" --config-settings=cmake.args="-DMEMILIO_USE_BUNDLED_PYBIND11=OFF"
```

Alternatively, edit the `CMakeCache.txt` in the directory created by scikit-build-core.

## Development

For developement of the cpp bindings use

```bash
pip install -e .[dev]
```

This command allows you to work on the code without having to reinstall the package after a change. Note that this only works for changes to Python code. If C++ code is modified, the install command has to be repeated every time. The command also installs all additional dependencies required for development and maintenance. 

The [bindings](memilio/simulation/bindings/) folder contains all the C++ code to expose the MEmilio library to Python and follows its structure, expect for the models that are bundled in a subfolder.

In order to add a new model, the following steps need to be taken:
1. Add new bindings including a model file in the [models](memilio/simulation/bindings/models/) folder that defines the new module
2. Add the new module to the building process by modifying [CMakeLists.txt](CMakeLists.txt)
3. Add a python module file similar to the other models, e.g. [osir.py](memilio/simulation/osir.py) (needed for the structure of the python package) and modify getter function in [__init__.py](memilio/simulation/__init__.py)
4. Write new tests and examples

The bindings can also be expanded with pure Python code by expanding the existing python files or adding new ones underneath the [simulation](memilio/simulation/) folder. 

## Stubs

A stub file is a file containing a skeleton of the public interface of that Python module including classes, variables, functions and their types. They help by enabling autocompletes and type annotations. `mypy stubgen` is used to generate the stub files for the MEmilio Python Bindings and provide them as a separate stubs-only package.

For installing stubs you first need to install our package `memilio.simulation` (not in editable mode -e) and the external dependency `mypy` for your python interpreter. Then run [generate_stubs.py](tools/generate_stubs.py) to generate the stubs-only package and install it as `memilio-stubs`, e.g. from the [current folder](.)

```bash
python ./tools/generate_stubs.py
```

## Usage

The package provides the following modules:

- `memilio.simulation`: core simulation framework and utilities, corresponds to the framework in `cpp/memilio`.
- `memilio.simulation.osir`: ODE SIR model and simulation with demographic resolution, corresponds to the model in `cpp/models/ode_sir`.
- `memilio.simulation.oseir`: ODE SEIR model and simulation with demographic resolution, corresponds to the model in `cpp/models/ode_seir`.
- `memilio.simulation.osecir`: ODE SECIHURD model and simulation with demographic and geographic resolution, corresponds to the model in `cpp/models/ode_secir`.
- `memilio.simulation.osecirvvs`: Extended ODE SECIHURD model to include vaccinations and multi-layered immunity, among other enhancements, with the ability to integrate new disease variants. Simulation includes demographic and geographic resolution, corresponds to the model in `cpp/models/ode_secirvvs`.
`memilio.simulation.abm`: Agent based model and simulation, corresponds to the model in `cpp/models/abm` (python model is incomplete and work in progress).

Detailed documentation under construction. See the scripts in the [examples](../examples/simulation/) directory for more information.

## Testing

The package provides a test suite in `memilio/simulation_test`. To run the tests, simply run the following command

```bash
python -m unittest
```

Note that these tests do not cover every case of the C++ library, they are only intended to test the binding code. To verify correctness of the C++ library itself, build and run the [C++ unit tests](../../cpp/README.md).
