Input/Output functionalities
============================

MEmilio provides a set of input/output (IO) functionalities to read and write data from and to files in different formats. IO functionality has been harmonized as much as possible across different models but due to the inherent differences between models, IO functionality is to some extent model-specific; for details see the descriptions below and the respective model (group) documentations.

Aggregated output
-----------------

All aggregated models save their results in a :code:`mio::TimeSeries` that contains all subpopulations for every time 
point, see :doc:`data_types` for a description of the :code:`mio::TimeSeries` data type. The time series can be printed 
to the terminal with the :code:`TimeSeries::print_table()` function and exported to a CSV file using :code:`TimeSeries::export_csv`. 
As e.g. for the ODE-based models adaptive step size methods are used, the time series will not be necessarily available 
on equidistant time points or days. To obtain a  :code:`mio::TimeSeries` interpolated on days or user-defined time points, 
the :code:`mio::interpolate_simulation_result()` can be used.

This document describes utilities for reading and writing data from and to files in different formats, in cases where
``TimeSeries::print_table()`` is not enough. The main sections are:

- The :ref:`serialization framework<serialization>`, that can be used to define the structure of data without using a specific file format.
  There are implementations of the framework for different formats. The framework is described in detail below, also
  see the `serialization example <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/serialize.cpp>`__.
  
- The :ref:`command line interface<command line>`, that can be used to set (and get) values of a ``ParameterSet``.

- The :ref:`History class<history>` for logging almost arbitrary information. It can be thought of as a generalization of a results
  ``TimeSeries``, and is mainly used for the ABM (and currently only implemented for it). 

There also exist other IO modules that are implemented for specific models.

- We have HDF5 support classes for C++, e.g. for reading of mobility matrix files in the graph-based metapopulation model, see :doc:`graph_metapop`.

.. _command line:

The command line interface
--------------------------

We provide a function ``mio::command_line_interface`` in the header ``memilio/io/cli.h``, that can be used to write to
or read from a parameter set. It can take parameters from command line arguments (i.e. the content of ``argv`` in the
main function), and assign them to or get them from a ``mio::ParameterSet``. A small example can be seen in
`cpp/examples/cli.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/cli.cpp>`_.

The command line interface (CLI) provides some non-parameter options listed below.

====================== =====================================
Name  (Alias)          Description
====================== =====================================
``--help`` (``-h``)    Shows the basic usage of the CLI, and lists each parameter by name, as well as any alias and
                       description. Takes priority before all other options and exits the program.
``--print_option``     Can be used with a (space-separated) list of parameter names or aliases (both without dashes) to
                       print the current values of each parameter to the terminal. This shows the correct JSON format
                       used by the parameters. Exits after use.
``--read_from_json``   Allows reading parameters from a file instead of the command line. Both parameter names and
                       aliases can be used, for example:

                       .. code-block::

                          {"<ParameterName>" : <value>, "<ParameterAlias>" : <value> }

``--write_to_json``    Writes *all* parameters with their current values to a specified file.
====================== =====================================

In general, an option is defined as a string, which consists either of two dashes followed by a name (e.g. ``--help``),
or a single dash followed by an alias (e.g. ``-h``). Apart from the built-in options, the names each refer to a
parameter that can be set.

To set the value of a parameter from the command line, first type the corresponding parameter option (see ``--help``),
followed by the value that should be assigned (reference ``--print_option``). Values are given as a JSON value
corresponding to the Type of the parameter. Note that some characters may need to be escaped or quoted. For example, the
JSON string ``"some string"`` must be entered as ``\\"some string\\"`` or ``'"some string"'``.

.. _history:

Working with the History object
-------------------------------

The History object provides a way to save data throughout the simulation process. It offers an interface where users can
define the data to be saved from a given object using Loggers and the method of saving it using ``Writer``\s. Afterward, the
user can access this data from the History object and manipulate it. For a basic ``Logger`` use case, refer to
`this example <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/history.cpp>`__. For an example demonstrating using a ``Logger`` in the ABM
`this example <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/abm_history_object.cpp>`_.

Loggers
~~~~~~~

The ``Logger`` struct is a tool for logging data from a given object. Each user-implemented ``Logger`` must have a ``Type``
and implement two functions: ``Type log(const T&)`` and ``bool should_log(const T&)``. The input ``T`` for these
functions is the same as the one given to the ``History`` member-function ``History::log``, e.g. ``Model&`` in 

- ``Type``: Return Type of ``log``.

- ``log``: This function determines which data from the input ``T`` is saved. It must have the same return Type ``Type``
  as the Loggers Type ``Type``.

- ``should_log``: This function must return a boolean to determine if data should be logged and can use the input ``T``
  for this, e.g. if ``T`` fulfills some criteria.

Users can derive their Loggers from ``LogOnce`` or ``LogAlways`` to use a predefined ``should_log`` function.
``LogOnce`` logs only at the first call of ``Logger::log()``, while ``LogAlways`` logs every time ``log`` is called.
All implemented Loggers must be default constructible/destructible. For user-defined examples in the ABM see
`this file <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/common_abm_loggers.h>`_.

.. code-block:: cpp

    struct LoggerExample { /* : public LogOnce/LogAlways if one wants to derive the should_log from these. */
        using Type = /* type of the record */;
        /* Below, T must be replaced by the type T from History::log(t). */
        Type log(const T& t) 
        {
            return /* something of type Type */;
        }
        bool should_log(const T& t) 
        {
              /* Determine whether log and add_record should be called by History::log(t). */
              return /* true or false */;
        }
    };

Writers
~~~~~~~

The ``Writer`` struct defines how to store the logged data from one or more implemented ``Loggers``. Each
user-implemented ``Writer`` must have a ``Data`` Type and implement the
``template <class Logger> static void add_record(const typename Logger::Type& t, Data& data)`` function.

- ``Data``: This is some kind of container that stores the data returned by the Loggers. For example, this can be a
  ``TimeSeries`` or depend on the Loggers (like ``std::tuple<std::vector<Logger::Type>...>``).

- ``add_record``: This manipulates the passed Data member of the ``History`` class to store the value ``t`` returned by
  the Loggers. It is used whenever ``History::log`` is called and ``Logger::should_log`` is true.

A predefined universal ``Writer`` called ``DataWriterToMemory`` is already implemented in `history.h <https://github.com/SciCompMod/memilio/blob/main/cpp/memilio/io/history.h>`__.
This stores the data from the loggers in a tuple of vectors every time the ``Logger`` is called. Another ``Writer`` named
``TimeSeriesWriter`` can be found in `this file <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/common_abm_loggers.h>`_, which saves 
Timeseries. The according ``Logger`` has to have a suitable return type.

.. code-block:: cpp

    template <class... Loggers>
    struct DataWriterExample {
        using Data = /* Container for the stored data of the Loggers */;
        template <class Logger>
        static void add_record(const typename Logger::Type& t, Data& data)
        {
              /* Manipulation of data to store the value t returned by the Loggers */;
        }
    };

History
~~~~~~~

The ``History`` class manages the ``Writer``\s and Loggers and provides an interface to log data. Currently it is only available in the ABM. It is templated on one
``Writer`` and several suitable and unique ``Logger``\s. To use the Writer to log something, the ``History`` provides the
function ``void log(const T& t)`` to call the ``add_record`` function of the ``Writer`` if the ``Logger`` function
``should_log`` returns true.

To access the data from the ``History`` class after logging, we provide the function ``get_log`` to access all records.
For this, the lifetime of the ``History`` has to be as long as one wants to have access to the data, e.g., a history
should not be constructed in the function it is called in when data is needed later.

To access data from a specific ``Logger``, one can use ``std::get<x>`` where x is the position of the ``Logger`` in the template
argument list of the ``History`` object. Refer to `this example <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/history.cpp>`__ for a simple
implementation of a history object and `this full ABM example <https://github.com/SciCompMod/memilio/blob/main/cpp/simulations/abm.cpp>`__ for a more 
of the History object with several History objects in use.

As mentioned, if multiple ``Writer``\s have to be used simultaneously, a separate History object is needed for each Writer.
For a use case of this, refer to `the ABM Simulation advance function <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/simulation.h>`__

.. _serialization:

The serialization framework
---------------------------

Serialization is the process of converting a data structure or object into a different format that can be stored or
transmitted. In this section we will show you how to make use of and implement MEmilio's serialization feature, as
well as explaining concepts, error handling, and extension of the feature to new types and formats.
Our guiding example will be a humble struct ``Foo``:

.. code-block:: cpp

   struct Foo {
     int i;
   };


Using serialization
~~~~~~~~~~~~~~~~~~~

In the next sections we will explain how to implement serialization (both for types and formats), here we quickly show
how to use it once it already is implemented for a type. In the following examples, we serialize (write) ``Foo`` to a
file in Json format, then deserialize (read) the Json again.

.. code-block:: cpp

   Foo foo{5};
   mio::IOResult<void> io_result = mio::write_json("path/to/foo.json", foo);

.. code-block:: cpp

   mio::IOResult<Foo> io_result = mio::read_json("path/to/foo.json", mio::Tag<Foo>{});
   if (io_result) {
     Foo foo = io_result.value();
   }

There is also support for a binary format. If you want to use a format directly instead of writing to a file, use the
``serialize_json``/``deserialize_json`` and ``serialize_binary``/``deserialize_binary`` functions.

Main functions and types
~~~~~~~~~~~~~~~~~~~~~~~~

- **Functions serialize and deserialize**:
  Main entry points to the framework to write and read values, respectively. The functions expect an `IOContext`
  (see :ref:`Concepts<concepts>` below) that stores the serialized data. (De-)serialization can be customized by providing a
  (de-)serialize_internal overload or a (de-)serialize member function for the type. See the section 
  :ref:`Adding a new data type to be serialized<adding new serialization>` or the documentation for ``serialize`` and ``deserialize``.
- **IOStatus and IOResult**:
  Used for error handling, see section :ref:`Error Handling<error handling>` below.

Default serialization
~~~~~~~~~~~~~~~~~~~~~

Before we get into the details of the framework, this feature provides an easy and convenient alternative to
implementing the serialize and deserialize functions. To give an example:

.. code-block:: cpp

   struct Foo {
     int i;
     auto default_serialize() {
       return Members("Foo").add("i", i);
     }
   };
   
Additional class members are added by repeated ``add`` calls, e.g. ``return Members("Foo").add("i", i).add("j", j)``,
where the first argument is a (descriptive) name and the second is a class member.

The default serialization is intentionally less flexible than the serialize and deserialize functions
(which will be explained later) and has additional requirements:

- The class must be default constructible.

  - If there is a default constructor that is *private*, it can still be used by marking the struct ``DefaultFactory``
    as a friend. For the example above, the line ``friend DefaultFactory<Foo>;`` would be added to the struct
    definition.
    
  - Alternatively, you may provide a specialization of the struct ``DefaultFactory``. For more details, view the
    struct's documentation.

- Every class member must be added to ``Members`` exactly once, and the provided names must be unique.

  - The members must be passed directly, like in the example. No copies, accessors, dereferencing, etc.

  - It is recommended, but not required, to add member variables to ``Members`` in the same order they are declared in
    the class, using the variables' names or something very similar. 

- Every class member itself must be serializable, deserializable and assignable.

This feature is primarily meant to make data classes easy to (de)serialize, avoiding some repetition that is necessary
when writing both a serialize and deserialize function. It can, however, be used for any class that should be
serialized in its entirety, and that does not need to make any decisions or computations while doing so. For example,
default serialization cannot be used if your class has optional members or values, or if one of its members is stored
as a pointer.

As to the feature set, default-serialization only supports the ``add_element`` and ``expect_element`` operations defined
in the :ref:`Concepts<concepts>` section, where each operation's arguments are provided through the ``add`` function. Note that the
value provided to ``add`` is also used to assign a value during deserialization, hence the class members must be used
directly in the function (i.e. as a non-const lvalue reference).

.. _concepts:

Concepts
~~~~~~~~

1. **IOContext**

   Stores data that describes serialized objects of any type in some unspecified format and provides structured
   access to the data for deserialization. Implementations of this concept may store the data in any format
   they want including binary. The data may also be written directly to disk. The context also keeps track
   of errors. An IOContext object ``io`` allows the following operations:

   - ``io.create_object("Type")``:
       Returns an IOObject for the type called ``"Type"``. The IOObject (see below) allows adding data that describes
       the object to be serialized. The function must return something that can be assigned to a local
       variable, e.g., a temporary or copyable function. IOObject may store references to the context internally,
       so the lifetime of the local IOObject may not exceed the lifetime of the IOContext that created it.
   - ``io.expect_object("Type")``:
       Returns an IOObject for the type called ``"Type"``. The IOObject (see below) provides access to the data needed
       for deserialization.
   - ``io.flags()``:
       Returns the flags that determine the behavior of serialization; see IOFlags.
   - ``io.error()``:
       Returns an ``IOStatus`` object to check if there were any errors during serialization. Usually it is not necessary to
       check this manually but can be used to report the error faster and avoid expensive operations that would be
       wasted anyway.
   - ``io.set_error(s)`` with some ``IOStatus`` object:
       Stores an error that was generated outside of the IOContext, e.g., if a value that was deserialized is outside an
       allowed range.

2. **IOObject**

   Gives structured access to serialized data. During serialization, data can be added with ``add_...`` operations.
   During deserialization, data can be retrieved with ``expect_...`` operations. Data must be retrieved in the same
   order as it was added since, e.g., binary format does not allow lookup by key. The following operations are supported
   for an IOObject ``obj``:

   - ``obj.add_element("Name", t)``:
     Stores an object ``t`` in the IOObject under the key "Name". If ``t`` is of basic type (i.e., int, string),
     IOObject is expected to handle it directly. Otherwise, the object uses ``mio::serialize`` to get the data for ``t``.
   - ``obj.add_list("Name", b, e)``:
     Stores the elements in the range represented by iterators ``b`` and ``e`` under the key "Name". The individual
     elements are not named. The elements are either handled directly by the IOObject or using ``mio::serialize`` just
     like ``add_element``.
   - ``obj.add_optional("Name", p)``:
     Stores the element pointed to by pointer ``p`` under the key "Name". The pointer may be null. Otherwise identical
     to add_element.
   - ``obj.expect_element("Name", Tag<T>{})``:
     If an object of type T can be found under the key "Name" and can be deserialized, returns the object. Otherwise
     returns an error. Analogously to serialization, the IOObject is expected to handle basic types directly and use
     ``mio::deserialize`` otherwise.
   - ``obj.expect_list("Name", Tag<T>{})``:
     If a list of objects of type T can be found under the key "Name" and can be deserialized, returns a range that can
     be iterated over. Otherwise returns an error.
   - ``obj.expect_optional("Name", Tag<T>{})``:
     Returns ``boost::optional<T>`` if an optional value of type T can be found under the key "Name". The optional may
     contain a value or it may be empty. Otherwise returns an error. Note that for some formats a wrong key is
     indistinguishable from an empty optional, so make sure to provide the correct key.

.. _error handling:

Error handling
~~~~~~~~~~~~~~

Errors are handled by returning error codes. The type ``IOStatus`` contains an error code and an optional string with
additional information. The type ``IOResult`` contains either a value or an ``IOStatus`` that describes an error. Operations
that can fail return an ``IOResult<T>`` where T is the type of the value that is produced by the operation if it is
successful. Except where necessary because of dependencies, the MEmilio framework does neither throw nor catch any exceptions.
IOContext and IOObject implementations are expected to store errors. During serialization, ``add_...`` operations fail
without returning errors, but the error is stored in the IOObject and subsequent calls are usually no-ops. During
deserialization, the values produced must usually be used or inspected, so ``expect_...`` operations return an IOResult.
The ``apply`` utility function provides a simple way to inspect the result of multiple ``expect_...`` operations and use
the values if all are successful. See the documentation of ``IOStatus``, ``IOResult`` and ``apply`` below for more
details.

.. _adding new serialization:
Adding a new data type to be serialized
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Serialization of a new type T can be customized by providing *either* member functions ``serialize`` and ``deserialize``
*or* free functions ``serialize_internal`` and ``deserialize_internal``.

The ``void serialize(IOContext& io)`` member function takes an IOContext and uses ``create_object`` and ``add_...``
operations to add data. The static ``IOResult<T> deserialize(IOContext& io)`` member function takes an IOContext and
uses ``expect_...`` operations to retrieve the data. The ``apply`` utility function can be used to inspect the result of
the ``expect_...`` operations and construct the object of type T.
E.g.:

.. code-block:: cpp

    struct Foo {
      int i;
      template<class IOContext>
      void serialize(IOContext& io) {
        auto obj = io.create_object("Foo");
        obj.add_element("i", i);
      }
      template<class IOContext>
      static IOResult<Foo> deserialize(IOContext& io) {
        auto obj = io.expect_object("Foo");
        auto i_result = obj.expect_element("i", mio::Tag<int>{});
        return mio::apply(io, [](auto&& i) { return Foo{i}; }, i_result);
      }
    };

The free functions ``serialize_internal`` and ``deserialize_internal`` must be found with argument-dependent lookup
(ADL). They can be used if no member function should or can be added to the type. See the code in `memilio/io/io.h <https://memilio.readthedocs.io/en/latest/api/program_listing_file__home_docs_checkouts_readthedocs.org_user_builds_memilio_checkouts_latest_cpp_memilio_io_io.h.html>`_
for examples where this was done for, e.g., Eigen3 matrices and STL containers.

Adding a new format
~~~~~~~~~~~~~~~~~~~

Implement concepts IOContext and IOObject that provide the operations listed above. Your implementation should handle
all built-in types as well as ``std::string``. It may handle other types (e.g., STL containers) as well if it can do so
more efficiently than the provided general free functions.

.. _epidemiological_data:

Epidemiological Data Integration
--------------------------------

For equation-based models, MEmilio provides direct integration with real-world epidemiological data through the :doc:`memilio-epidata <../python/memilio_epidata>` package. This allows models to be initialized and calibrated with actual 
population, case, hospitalization, and ICU data from various sources.

MEmilio provides specialized functions in the parameter I/O modules to read and process this epidemiological data 
for model initialization, such as population data, confirmed cases, vaccination data, and ICU data.

Additionally, the integration supports different administrative levels such as districts, counties, and states. 
Therefore, MEmilio provides high-level functions that combine multiple data sources:

.. code-block:: cpp

    // Read all input data for Germany
    template <class Model>
    IOResult<void> read_input_data_germany(std::vector<Model>& model, 
                                           Date date,
                                           const std::vector<double>& scaling_factor_inf,
                                           double scaling_factor_icu,
                                           const std::string& pydata_dir);
    
    // Read input data for specific counties
    template <class Model>
    IOResult<void> read_input_data_county(std::vector<Model>& model, 
                                          Date date,
                                          const std::vector<int>& county,
                                          const std::vector<double>& scaling_factor_inf,
                                          double scaling_factor_icu,
                                          const std::string& pydata_dir);

These functions automatically handle:

- Reading population data from census files
- Loading case data from RKI sources
- Integrating ICU data from DIVI register
- Applying scaling factors for data adjustment
- Setting initial conditions in epidemiological models
