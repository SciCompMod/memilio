How to: Common patterns and Guidlines
==========================================

We want to look at some common patterns developers will encounter when extending the ``memilio.simulation`` package and explain how to handle them. 
This will include patterns specific to the MEmilio C++ library and development guidelines on connecting general C++ functionalities.

MEmilio binding namespace pymio
--------------------------------

MEmilio simulation provides different functionalities for extending the basic bindings in the namespace ``pymio``. 


Custom class binding functor
-----------------------------

As part of ``pymio`` a custom binding functor ``pymio::bind_class`` extends the pybind11 ``pybind11::class_`` to create bindings for C++ classes. 
It should always be used instead of the pybind11 functor. Through this interface additional functionalities can be added to the bindings of C++ classes specific for MEmilio.

Currently, the custom binding functor has a flag as template parameter defining the strategy with which classes are exposed for serialization. The three strategies defined in the enum ``EnablePickling`` are:

* ``Never``: Pickling support is disabled
* ``Required``: Pickling support is required, and an error is raised if not available
* ``IfAvailable``: Pickling support is enabled if the class matches SFINEA-clause of ``serialize_internal``

Therefore, for the following example class

.. code-block:: c++

    class MyClass {
        ...
    };

the bindings with enabled serialization based on SFINEA-clauses could look like

.. code-block:: c++

    PYBIND11_MODULE(module_name, m)
    {
        pymio::bind_class<MyClass, pymio::EnablePickling::IfAvailable>(m, "MyClass")
    }

Preprocessor macros
--------------------

##TODO

Template classes/functions
--------------------------

C++ introduces template programming. When exposing templated C++ code to Python, Python cannot work with templates, because templates are a compile-time C++ feature, 
but Python is dynamically typed and runtime-based.

- introduce templates with example

- explain the problem for bindings, ie. python can only connect to 

So, in your documentation, you're referring to the process of choosing and naming specific instantiations of the C++ templates to expose them to Python.

.. Yes, the process you're referring to in C++ is generally known as template instantiation, and more specifically explicit template instantiation or template specialization when you are selecting or defining specific versions of a templated class or function.

IOResult as return
------------------

- mio usage of IOResult for handling return values of io functions
- short example
- IOResult is not directly exposed, instead generally using lambda function calls to unpack it before returning the result to python
- example code of binding