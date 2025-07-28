MEmilio Generation Package
==========================

This package provides an automatic code generator for python bindings of the MEmilio C++ library. 
It enables the automatic generation of a part of the :doc:`Python Bindings <memilio_simulation>` that is common across multiple models. 
For a particular example, see the SEIR model with its files `oseir.cpp` and `oseir.py`.

This generating software was developed as a part of the Bachelor thesis `Automatische Codegenerierung für nutzerfreundliche mathematisch-epidemiologische Modelle <https://elib.dlr.de/190367/>`_. 
The following figure from Chapter 5 outlines the workflow of the generator. Blue boxes represent parts of the code generator and orange ones the input and output. Rectangular boxes contain classes with logic, the rest represent data.

.. image:: ../../../pycode/memilio-generation/generator_workflow.png
   :alt: tikzGeneratorWorkflow

Dependencies
------------

The package uses the `Clang C++ library <https://clang.llvm.org/>`_ and the `LibClang Python library <https://libclang.readthedocs.io/en/latest/index.html>`_ to analyze the C++ code of a model. Both need to be installed and share the same version.

Required python packages:

* scikit-build
* dataclasses
* dataclasses_json
* graphviz
* importlib-resources

Usage
-----

During the installation the package creates a compilation database (compile_commands.json) for the models of the `C++ Library <https://github.com/SciCompMod/memilio/blob/main/cpp/>`_.

The package provides an example script on how to use it in `memilio/tools`. The example uses the ode_seir model.

Before running the example you have to do these steps of setup:

* Set the source file path in `example_oseir.py <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/tools/example_oseir.py>`_ under ``conf.source_file`` to the path of the C++ model you want to generate bindings for.
* Set the target folder path in `example_oseir.py <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/tools/example_oseir.py>`_ under ``conf.target_folder`` to the path where you want the generated bindings to be. 


Example:
After processing as described in the previous paragraph, run the example with the command (adjust the path according to your current folder):

.. code-block:: console 

    python memilio/tools/example_oseir.py 



Visualization
-------------

The package contains a `Visualization class <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/generation/graph_visualization.py>`_  to display the generated AST.
This allows you to visualize the abstract syntax tree (AST) of the C++ model.

This allows for:

* Printing the AST in the terminal.
* Saving the AST as a PDF file.
* Formatting the AST in a text file.

Example:
``aviz.output_ast_formatted(ast, ast.get_node_by_index(1))`` displays the second node of the AST and its children in a file called ast_formatted.txt. 
With the root node ``.get_node_by_index(0)`` you can display the whole AST.

``aviz.output_ast_terminal(ast, ast.get_node_by_index(1))`` displays the second node of the AST and its children in terminal.

The first argument of the statements specifies the given AST. The second argument specifies the node and its children that you want to display.

``aviz.output_ast_png(ast.get_node_by_index(2), 2)`` displays the third node of the AST and its children with a depth of 2 as a png file. 
The second argument of the statement specifies the depth up to which the function displays child nodes. 
This means that any nodes beyond the specified depth (e.g., all nodes at level 3 and beyond if the depth is set to 2) will not be shown.

Notice that the visualization as a png-file should not print the whole AST, as it is not possible to display the whole AST in a single image.

Testing
-------

The package provides a test suite in `memilio/generation_test <https://github.com/SciCompMod/memilio/tree/main/pycode/memilio-generation/memilio/generation_test>`_. 
To run the tests, simply use the following command:

.. code-block:: console 

    python -m unittest


Development
-----------

When implementing new model features you can follow these steps:

* Add necessary configurations to `config.json.txt <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/tools/config.json.txt/>`_ and add corresponding attributes to the ``ScannerConfig``.
* For the features you want to implement, find the nodes in the abstract syntax tree (AST) (use method Scanner.output_ast_file(); see the example in tools/).
* Add the extraction of those features. Therefore you need to change the "check_..."-methods corresponding to the ``CursorKind`` of your nodes in the ``Scanner``. If there is no corresponding "check\_..."-method you need to write a new one and add it to the switch-method (``scanner.switch_node_kind()``).
* Extend the ``IntermediateRepresentation`` for the new model features.
* Adjust the `cpp-template <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/generation/template/template_cpp.txt>`_ and the `string-template-methods <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/generation/template/template_string.py>`_. If needed, use new identifiers and write new string-template-methods for them.
* Adjust the substitution dictionaries in the ``Generator``.
* Write new/Adjust scripts in the `tool folder <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/tools/>`_ for the model and try to run.
* Update tests.