Common data types
-----------------

The follwing list expalins the non-standard data types that are used throughout MEmilio.

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Data type name
     - Description
   * - :code:`FP`
     - A floating point type. Usually :code:`double` is used, but for instane in the optimization using optimal control :code:`FP` is equal to :code:`Ipopt::Number`, see models/oseair and `examples/ode_seair_optimization.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_seair_optimization.cpp>`_.
   * - :code:`UncertainValue`
     - This data type describes a value sampled from a given distribution. The value is intialized with a given :code:`FP` and can be (re)sampled with the :code:`draw_sample()` function.
   * - :code:`AgeGroup`
     - A typesafe ``size_t``, i.e. an integer that cannot be confused with other integer types so operations like assignment, addition etc. only work with other :code:`AgeGroup`s.
   * - :code:`Region`
     - A typesafe ``size_t``.
   * - :code:`CustomIndexArray`
     - An array whose values can be accesses by a multi-index. This datatype is for example used in the parameter :code:`mio::abm::TimeExposedToNoSymptoms` making it dependent on :code:`mio::abm::VirusVariant` and :code:`mio::AgeGroup`. Its values can then be set for a specific :code:`virus_variant` and :code:`age_group` using :code`model.parameters.template get<mio::abm::TimeInfectedSevereToCritical>()[{virus_variant, age_group}]`.
   * - :code:`Populations`
     - Is a :code:`mio::CustomIndexArray` with :code:`mio::UncertainValue<FP>` as values.
