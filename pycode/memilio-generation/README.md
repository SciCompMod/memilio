# MEmilio Automatic Code Generation of Python Bindings

This package contains Python bindings generating code for the MEmilio C++ library. 
It enables the automatic generation of a part of the [Python Bindings](../memilio-simulation/README.md) that is common across multiple models. For a particular example, see the SEIR model with its files `oseir.cpp` and `oseir.py`.

This generating software was developed as a part of the Bachelor thesis [Automatische Codegenerierung für nutzerfreundliche mathematisch-epidemiologische Modelle](https://elib.dlr.de/190367/). The following figure from Chapter 5 outlines the workflow of the generator. Blue boxes represent parts of the code generator and orange ones the input and output. Rectangular boxes contain classes with logic, the rest represent data.

[<img src="generator_workflow.png" width="50%"/>](generator_workflow.png "Workflow of the code generator")

## Prerequisites

The package uses the [Clang C++ library](https://clang.llvm.org/) and the [LibClang Python library](https://libclang.readthedocs.io/en/latest/index.html) to analyze the C++ code of a model. Both need to be installed and share the same version.

## Installation

Use the provided `setup.py` script to build and install the package. To install the package, use the command (from the directory containing `setup.py`)

```bash
pip install -e .[dev]
```

During this step the package creates a compilation database (compile_commands.json) for the models of the C++ library[C++ Library](../../cpp/README.md).

Afterwards you can use the following command to make a full installation of the package (not necessary for usage)

```bash
pip install .
```

## Usage

The package provides an example script on how to use it in `memilio/tools`. The example uses the ode_seir model from the [C++ Library](../../cpp/models/ode_seir/README.md).

Before running the example you have to do these steps of setup:
- Change [config.json.txt](./memilio/tools/config.json.txt). You can find a documentation in the ScannerConfig class.
- Check if the parameters set in __post_init__() of the [ScannerConfig class](./memilio/generation/scanner_config.py) match with the cpp-class names.

Example:
After processing as described in the previous paragraph, run the example with the command (path according to the current folder):

```bash
python memilio/tools/example_oseir.py 
```

When working on a new model you can copy the example script and add an additional segment to the config.json.txt. The setup works similar to the example. Additionaly you can print the AST of your model into a file (Usefull for development/debugging).

## Testing

The package provides a test suite in `memilio/generation_test`. To run the tests, simply run the following command:

```bash
python -m unittest
```

## Development

When implementing new model features you can follow these steps:
- Add necessary configurations to [config txt-file](./memilio/tools/config.json.txt) and add corresponding attributes to the [ScannerConfig class](./memilio/generation/scanner_config.py).
- For the features you want to implement, find the nodes in the abstract syntax tree (AST) (use method Scanner.output_ast_file(); see the example in tools/).
- Add the extraction of those features. Therefore you need to change the "check_..."-methods corresponding to the CursorKind of your nodes in the [Scanner class](./memilio/generation/scanner.py). If there is no corresponding "check_..."-method you need to write a new one and add it to the switch-method (scanner.switch_node_kind()).
- Extend the [IntermediateRepresentation](./memilio/generation/intermediate_representation.py) for the new model features.
- Adjust the [cpp-template](./memilio//generation/template/template_ode_cpp.txt) and the [string-template-methods](./memilio/generation/template/template_ode_string.py). If needed, use new identifiers and write new string-template-methods for them.
- Adjust the substitution dictionaries in the [Generator class](./memilio/generation/generator.py).
- Write new/Adjust script in the [tool folder](./memilio/tools/) for the model and try to run.
- Update [tests](./memilio/generation_test/).