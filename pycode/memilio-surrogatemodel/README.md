MEmilio Surrogate Model Package
=======================
This package contains machine learning based surrogate models that make predictions based on the MEmilio simulation models. Currently there are only surrogate models for equation-basesd models. These ODE simulations are used for data generation. The goal is to create a powerful tool that predicts the spread of Covid19 faster than an ODE simulation while aiming for high precision. 
 
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

Since we are running simulations to generate the data, the MEmilio `simulation <https://github.com/DLR-SC/memilio/blob/main/pycode/memilio-simulation/setup.py>`_ package also needs to be installed.
## Usage
The package currently provides the following modules:

- `models`: models for different specific tasks
   Currently we have the following models: 
   - `ode_secir_simple`: a very simplified model with only one population and one age group.

     Each model folder contains the following files: 
     - `data_generation`: data is generated from memilio secir simulation.
     - `model`: training and evaluation of the model. 
     - `different_networks`: multiple network architectures are saved in this file.


- `tests`: this file contains all tests 

## Testing 
The package provides a test suite in `memilio/surrogatemodel_test`. To run the tests, simply run the following command.

```bash
python -m unittest
```
