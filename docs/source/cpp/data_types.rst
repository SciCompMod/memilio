Common data types
-----------------

Here we want to explain the intention behind the non-standard data types that are used throughout MEmilio.

FP
~~

The template :code:`FP` and the derived types like :code:`UncertainValue<FP>` in these examples are commonly used throughout MEmilio.
:code:`FP` is a floating point type, usually :code:`double`. An :code:`UncertainValue<FP>` holds a value of type
:code:`FP` as well as (optionally) a distribution to sample new values from, e.g. for a parameter study.

.. dropdown:: :fa:`gears` Expert's settings

    The type ``mio::AgeGroup`` is a typesafe ``size_t``, meaning an integer that cannot be confused with other integer
    types. So assignment, addition, etc. only works with another ``mio::AgeGroup``, not ``size_t`` or another integer
    type. This is useful for function interfaces or indexing, as it makes it (nearly) impossible to mix up, e.g., age
    groups with infection states. Check out ``mio::Index`` if you want to learn more.

    The type ``mio::Populations`` is an extension of a ``mio::CustomIndexArray``, which is a template type that manages
    a flat array. Its main purpose is to allow multidimensional indexing into this array, using typesafe indices like
    a ``mio::Index`` or a ``enum class``.