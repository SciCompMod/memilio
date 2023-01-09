MEmilio Surrogate Model Package
=======================
This package contains machine learning based surrogate models that make predictions based on the MEmilio simulation models. Currently there are only surrogate models for ODE models. These simulations of these equation-based models are used for data generation. The goal is to create a powerful tool that predicts the dynamics faster than a simulation of an expert model, e.g., a metapopulation or agent-based model while still having acceptable errors with respect to the original simulations.
 
## Installation

Use the provided `setup.py` script install the package. 
To install the package, use the command (from the directory containing `setup.py`)

```bash
pip install .
```

For developement of code use

```bash
pip install -e .[dev]
``` 

Since we are running simulations to generate the data, the MEmilio `memilio-simulation` package (https://github.com/DLR-SC/memilio/tree/main/pycode/memilio-simulation) also needs to be installed.
## Usage
The package currently provides the following modules:

- `models`: models for different specific tasks
   Currently we have the following models: 
   - `ode_secir_simple`: A simple model allowing for asymptomatic as well as symptomatic infection not stratified by age groups.

     Each model folder contains the following files: 
     - `data_generation`: data generated from expert model simulation.
     - `model`: training and evaluation of the model. 
     - `network_architectures`: multiple network architectures are saved in this file.


- `tests`: this file contains all tests 

## Testing 
The package provides a test suite in `memilio/surrogatemodel_test`. To run the tests, simply run the following command.

```bash
python -m unittest
```
