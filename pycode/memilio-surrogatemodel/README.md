MEmilio Surrogate Model Package
=======================
This package contains machine learning models that make predictions based on the memilio ODE models. ODE simulations are used for data generation. The goal is to create a powerful tool that predicts the spread of Covid19 faster than an ODE simulation while aiming for high precision. 
 
## Installation

Use the provided `setup.py` script to build the bindings and install the package. The script requires CMake and the Scikit-Build packages. Both are installed by the script if not available on the system. The package uses the [Pybind11 C++ library](https://pybind11.readthedocs.io) to create the bindings.

To install the package, use the command (from the directory containing `setup.py`)

```bash
pip install .
```

This builds the C++ library and C++ Python extension module and copies everything required to your site-packages. 

For developement of code use

```bash
pip install -e .[dev]
``` 

## Usage

The package provides the following modules:

- `models`: models for different specific tasks
   Currently we have the following models: 
   - `ode_secir_simple`: a very simplified model with only one population and one age group.
   - `ode_secir_age_groups`: extension of the simple model by considering multiple age groups. In this model one contact pattern reduction (damping) is considered. 
   - `ode_multiple_damping`: this model consideres multiple dampings and makes predictions for a longer time period. 

     each model folder contains the following files: 
     - `data_generation`: data is generated from memilio secir simulation.
     - `model`: training and evaluation of the model. 
     - `different_networks`: multiple network architectures are saved in this file.


- `tests`: this file contains all tests 

## Testing 
The package provides a test suite in `memilio/tests`. To run the tests, simply run the following command.

```bash
python -m unittest
```
