MEmilio Inference Package
=======================
This package contains...
 
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

Since we are running simulations to generate the data, the MEmilio `memilio-simulation` package (https://github.com/SciCompMod/memilio/tree/main/pycode/memilio-simulation) also needs to be installed.

## Usage
...

## Testing 
The package provides a test suite in `memilio/inference_test`. To run the tests, simply run the following command.

```bash
python -m unittest
```
