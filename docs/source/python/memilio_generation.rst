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
* Change `config.json.txt <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/tools/config.json.txt/>`_.
* Check if the parameters set in __post_init__() of the [ScannerConfig class](./memilio/generation/scanner_config.py) match with the cpp-class names.

Example:
After processing as described in the previous paragraph, run the example with the command (path according to the current folder):

.. code-block:: console 

    python memilio/tools/example_oseir.py 


When working on a new model you can copy the example script and add an additional segment to the config.json.txt. The setup works similar to the example. Additionaly you can print the AST of your model into a file (Usefull for development/debugging).

Testing
-------

The package provides a test suite in `memilio/generation_test <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-surrogatemodel/memilio/generation_test>`_. 
To run the tests, simply run the following command:

.. code-block:: console 

    python -m unittest


Development
-----------

When implementing new model features you can follow these steps:
* Add necessary configurations to `config.json.txt <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/tools/config.json.txt/>`_ and add corresponding attributes to the ``ScannerConfig``.
* For the features you want to implement, find the nodes in the abstract syntax tree (AST) (use method Scanner.output_ast_file(); see the example in tools/).
* Add the extraction of those features. Therefore you need to change the "check_..."-methods corresponding to the ``CursorKind`` of your nodes in the ``Scanner``. If there is no corresponding "check_..."-method you need to write a new one and add it to the switch-method (scanner.switch_node_kind()).
* Extend the ``IntermediateRepresentation`` for the new model features.
* Adjust the `cpp-template <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/generation/template/template_ode_cpp.txt>`_ and the `string-template-methods <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/generation/template/template_ode_string.py>`_. If needed, use new identifiers and write new string-template-methods for them.
* Adjust the substitution dictionaries in the ``Generator``.
* Write new/Adjust script in the `tool folder <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/tools/>`_ for the model and try to run.
* Update tests.