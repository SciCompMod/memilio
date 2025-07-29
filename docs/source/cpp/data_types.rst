Common data types
-----------------

The following list explains the nonstandard data types that are used throughout MEmilio.

.. list-table::
   :header-rows: 1
   :widths: 20 60

   * - Data type name
     - Description
   * - :code:`FP`
     - A floating point type. Usually :code:`double` is used, but for instance in the optimization using optimal control :code:`FP` is equal to :code:`Ipopt::Number`, see :doc:`models/oseair` and `examples/ode_seair_optimization.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_seair_optimization.cpp>`_.
   * - :code:`UncertainValue`
     - This data type describes a value sampled from a given distribution. The value is initialized with a given :code:`FP` and can be (re)sampled with the :code:`draw_sample()` function.
   * - :code:`AgeGroup`
     - A typesafe ``size_t``, i.e. an integer that cannot be confused with other integer types, so operations like assignment, addition etc. only work with other :code:`AgeGroup`\s. Derived from :code:`mio::Index`.
   * - :code:`Region`
     - A typesafe ``size_t``, derived from :code:`mio::Index`.
   * - :code:`CustomIndexArray`
     - A contiguous array whose values can be accessed by a multi-index, e.g. a :code:`mio::Index` with one or more categories. This datatype is, for example, used in the parameter :code:`mio::abm::TimeExposedToNoSymptoms` making it dependent on :code:`mio::abm::VirusVariant` and :code:`mio::AgeGroup`. Its values can then be set for a specific :code:`virus_variant` and :code:`age_group` using :code:`model.parameters.template get<mio::abm::TimeInfectedSevereToCritical>()[{virus_variant, age_group}]`.
   * - :code:`Populations`
     - Is a :code:`mio::CustomIndexArray` with :code:`mio::UncertainValue<FP>` as values. Adds some convenient functions like :code:`get_group_total`.
   * - :code:`TimeSeries`
     - Stores vectors of values at time points. Each time point has a vector of values of the same size with operations like adding time points, retrieving values, exporting to CSV, etc. It's also used for storing and analyzing simulation results over time.
   * - :code:`Graph`
     - A generic graph structure that represents a network of nodes connected by edges. Each node and edge can have associated properties. The graph is used to model geographical regions connected by mobility patterns (e.g., commuting), where each node is represented by its own epidemiological model.
   * - :code:`Node`
     - Represents a node in a graph with a unique ID and associated properties. 
   * - :code:`Edge`
     - Represents a directed connection between two nodes in a graph with associated properties.
   * - :code:`EdgeBase`, :code:`InEdgeBase`, :code:`OutEdgeBase`
     - Base classes for Edge that define start and end node indices for connections in the graph.
   * - :code:`SimulationNode`
     - Represents a simulation in one node of a mobility graph. Contains a simulation model of any type and keeps track of the last state and time point.
   * - :code:`MobilityCoefficients`
     - Time-dependent mobility coefficients used to model how populations move between nodes in a graph.
   * - :code:`MobilityCoefficientGroup`
     - A collection of time-dependent mobility coefficients that differentiate between various sources of mobility.
   * - :code:`MobilityParameters`
     - Parameters that influence mobility between nodes, including coefficients and dynamic nonpharmaceutical interventions (NPIs).
   * - :code:`MobilityEdge`
     - Represents mobility between two nodes in a graph. Handles the movement of populations between nodes, tracks mobile populations, and applies mobility returns according to epidemiological models.
   * - :code:`ContactMatrix`
     - Time-dependent contact frequencies between groups, derived from ``DampingMatrixExpression``. Models how the contact rates between different age groups change over time due to interventions.
   * - :code:`ContactMatrixGroup`
     - A collection of contact matrices that represent different contexts (e.g., home, school, work) whose sum is the total number of contacts, derived from ``DampingMatrixExpressionGroup``.
   * - :code:`DampingMatrixExpression`
     - Represents a coefficient-wise matrix expression :math:`B - D \odot (B - M)`, where :math:`B` is a baseline matrix, :math:`M` is a minimum matrix, :math:`D` is a time-dependent complex damping factor, and :math:`\odot` is element wise multiplication. Used as the base for time-dependent contact matrices.
   * - :code:`DampingMatrixExpressionGroup`
     - Represents a collection of ``DampingMatrixExpression``\s that are summed up. Used for representing multiple sources of contacts or mobility.
