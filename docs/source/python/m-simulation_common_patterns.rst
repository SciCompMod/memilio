How to: Common patterns and guidlines
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
        pymio::bind_class<MyClass, pymio::EnablePickling::IfAvailable>(m, "MyClass");
    }

Configuration Macros
--------------------

Similar to the MEmilio C++ library, MEmilio Simulation uses compile-time configuration options through CMake to enable or disable specific features. 
The correct preprocessor conditionals need to be added around the bindings, if those features are needed. For more information read :doc:`Configuration Macros <../cpp/cpp_makros_preprocessing>`.

Template classes/functions
--------------------------

C++ introduces template programming. Python cannot work with templates, because templates are a compile-time C++ feature, 
but Python is dynamically typed and runtime-based. Therefore, when exposing templated C++ code to Python only specific template instantiations
can be connected, i.e. each template specification that should be exposed to Python needs to be binded explitily with its own name.

To reduce the redundancy while programming MEmilio often uses templated binding methods handling a single template object. Those can be defined in the ``pymio`` namespace and
used in the module definition code. We will follow with a short example by defining a simple templated class.

.. code-block:: c++

    template <class Template>
    class MyTemplatedClass {
        ...
    };

To implement the interface with pybind11, we can define a function taking the same template arguments as the template class.

.. code-block:: c++

    namespace pymio {

    template <class Template>
    void bind_MyTemplatedClass(pybind11::module_& m, std::string const& name)
    {
    bind_class<MyTemplatedClass<Template>, pymio::EnablePickling::IfAvailable>(m, name.c_str())
        .def(...)
        ...;
    }

    } // namespace pymio

Now, with a single call to ``bind_MyTemplatedClass`` a specific instantiations is exposed. Let´s consider we have three different instantiations we know of, then the module definition could look like:

.. code-block:: c++

    PYBIND11_MODULE(module_name, m)
    {
        pymio::bind_MyTemplatedClass<TemplateType1>(m, "MyTemplatedClass1");
        pymio::bind_MyTemplatedClass<TemplateType2>(m, "MyTemplatedClass2");
        pymio::bind_MyTemplatedClass<TemplateType3>(m, "MyTemplatedClass3");

    }

This can be extended to include more definitions for the templated class or adding the pickling flag as a template aswell, so they can differ for each instantiation.
Additionally, the binding function could return the ``pybind11::class_`` object, if some instantiations need custom interfaces which are not generalizable.
For an example of a MEmilio class using this approach look at the ``CompartmentalModel`` class bindings defined in `compartmental_model.h <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-simulation/memilio/simulation/bindings/compartments/compartmental_model.h>`_.

IOResult as return
------------------

MEmilio uses ``mio::IOResult`` for handling return values and errors of io functions. An example is the ``save_result`` function in `result_io.cpp <https://github.com/SciCompMod/memilio/blob/main//cpp/memilio/io/result_io.cpp>`_ with the following function interface.

.. code-block:: c++

    mio::IOResult<void> save_result(const std::vector<TimeSeries<double>>& results, const std::vector<int>& ids, int num_groups,
                               const std::string& filename) 
    {
        ...;
        return mio::success();
    }


``mio::IOResult`` is not exposed to Python, instead the needed functions get exposed by wrapping them inside a lambda function.
The return object is handled through ``pymio::check_and_throw``, which either throws a Python error or returns the internal object. 

.. code-block:: c++

    m.def(
        "save_result",
        [](...) {
            auto result = mio::save_result(...);
            return pymio::check_and_throw(result);
        },
        py::return_value_policy::move
    );
