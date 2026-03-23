# PyGen usage of the generation and visualization

A guide for the generation package, including installation instructions and simple examples, can be found in [MEmilio Automatic Code Generation of Python Bindings](/pycode/memilio-generation/README.md).

This README provides guidance on how to use the generation package and the visualization features. The explanation follows the example script in [Example](/pycode/memilio-generation/memilio/tools/example_oseir.py). This script allows you to select an ODE model you want to bind, based on the model.cpp from the [Models folder](/cpp/models/..) provided through the corresponding file path

## Requirements

The package uses the [Clang C++ library](https://clang.llvm.org/) and the [LibClang Python library](https://libclang.readthedocs.io/en/latest/index.html) to analyze the C++ code of a model. Both need to be installed and share the same version.


## Usage

Before running the example script, configure the following settings:

`conf.source_file`: Specifies the source file to be bound.

`conf.target_folder`: Defines the folder where the generated Python file will be stored.


Example Configuration:
If you want to use the ode_secirvvs model the config should be like this -> `conf.source_file = os.path.abspath(os.path.join(file_path, "..", "..", "..", "..", "cpp", "models", "ode_secirvvs", "model.cpp"))` 
This descibes the filepath to your model.cpp outgoing from the tools folder.

To use a different model, replace ode_secirvvs with the desired model's name in the file path.

Adding Additional Source Files:
In addition to the main source file, it is possible to include further source files in the binding process.
For this purpose, the `asts´ instance of the [ASTHandler class](/pycode/memilio-generation/memilio/generation/ast_handler.py) provides the method `add_source_file´:
```python
asts.add_source_file("name_of_the_file.cpp")
```

To set up a target folder, specify the desired output directory for the generated bindings here:
`conf.target_folder = file_path`
If no custom path is specified and `file_path`is used as the target folder, the generated files will be saved in
the tools directory by default, as this is the working directory of the skript.

After setting up the source file and target folder, set the path to the example script in the terminal:
If the terminal shows the memilio package, just give the path with: 
```bash
cd pycode/memilio-generation/memilio/tools
```

After the path is set, run the example script with the following command:

```bash
python example_oseir.py
```

## Visualization

The package contains a [Visualization class](/pycode/memilio-generation/memilio/generation/graph_visualization.py) to display the generated AST.
This allows for:

- Printing the AST in the terminal

- Saving it as a graph (PNG format)

- Formatting it in a text file

You can print the AST in the example script with the aviz instance of the Visualisation() class.


Example:
`aviz.output_ast_formatted(asts.get_ast_by_id(0), ast.get_node_by_index(1))` displays the second node of the first AST and its children in a file called ast_formatted.txt. 
The root node `.get_node_by_index(0)` displays the whole AST.

`aviz.output_ast_terminal(asts.get_ast_by_id(0), ast.get_node_by_index(1))` displays the second node of the AST and its children in terminal.

The first argument of the statements specify the given AST. The second argument specifies the node and its children that you want to display.

`aviz.output_ast_png(ast.get_node_by_index(2), 2)` displays the third node of the AST and its children with a depth of 2 as a png file. The second argument of the statement specifies the depth up to which the function displays child nodes. This means that any nodes beyond the specified depth (e.g., all nodes at level 3 and beyond if the depth is set to 2) will not be shown.
The output of the AST as a PNG should be restricted, since not many nodes can be displayed at once in a PNG. So only subtrees with a restricted depth should be printed.

The visualization is specified under the `if (print_ast):` statement in the example script.
It is enabled when running the example with the following command:


```bash
python example_oseir.py -p
```





