# MEmilio Automatic Code Generation of Python Bindings

This package contains Python bindings for the MEmilio C++ library. 
It enables the automatic generation of a part of the [Python Bindings](../memilio-simulation/README.md). It focuses on the model specific binding files e.g., for the seir model the files `oseir.cpp` and `oseir.py`. 

## Prerequisites

The package uses the [Clang C++ library](https://clang.llvm.org/) and the [LibClang Python library](https://libclang.readthedocs.io/en/latest/index.html) to analyze the C++ code of a model. Both need to be installed and share the same version.

To use the package you need a usable model including a compilation database (compile_commands.json). For models of the C++ library the compilation database is generated when building the [C++ Library](../../cpp/README.md).

## Installation

Use the provided `setup.py` script to build and install the package. To install the package, use the command (from the directory containing `setup.py`)

```bash
pip install .
```

## Usage

The package provides an example script on how to use it in `memilio/tools`. The example uses the ode_seir model from the [C++ Library](../../cpp/models/ode_seir/README.md).

Before running the example you have to do multiple steps of setup:
- Change `memilio/tools/config.json`
- Change `memilio/generation/scanner_config.py`

Run the example with the command (path according to the current folder):

```bash
python memilio/tools/seir.py 
```

When working on a new model you can copy the example script and add an additional segment to the config.json. The setup works similar to the example. Additionaly you can print the AST of your model into a file (Usefull for development/debugging).

## Testing

The package provides a test suite in `memilio/generation_test`. To run the tests, simply run the following command

```bash
python -m unittest
```

## Development

When implementing new model features you can follow these steps:
- Add necessary configurations to [config file](./memilio/tools/config.json) and add corresponding attributes to the [ScannerConfig class](./memilio/generation/scanner_config.py).
- Find the nodes in the AST with the features you want to implement (use method Scanner.output_ast_file()).
- Add the extraction of those features. Therefore you need to change the "check_*"-methods corresponding to the CursorKind of your nodes in the [Scanner class](./memilio/generation/scanner.py). If there is no corresponding "check_*"-method you need to write a new one and add it to the switch-method (scanner.switch_node_kind()).
- Extend the [IntermediateRepresentation](./memilio/generation/intermediate_representation.py) for the new model features.
- Adjust the [cpp-template](./memilio//generation/template/template_cpp.txt) and the [string-template-methods](./memilio/generation/template/template_string.py). If needed use new identifiers and write new string-template-methods for them.
- Adjust the substitution dictionaries in the [Generator class](./memilio/generation/generator.py).
- Write new/Adjust script in the [tool folder](./memilio/tools/) for the model an try to run.
- Update [tests](./memilio/generation_test/)