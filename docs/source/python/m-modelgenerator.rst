Model Generator
===============

The model generator is part of the ``memilio-generation`` package and provides a high-level way to create new
compartmental ODE models for MEmilio from a simple configuration file. Instead of writing C++ code by hand, you
describe your model in a YAML or TOML file and the generator produces all required source files automatically.
C++ knowledge is not required to use the generator, but you can of course edit the generated C++ code afterwards
if you want to add custom features. Additionally, a Python example application is generated that you can run
immediately after generation is done.

.. note::

    Here, you start with a model specification and get C++ source files and Python bindings as output.
    If the already have a C++ model and want to generate Python bindings for it, you can use the 
    **bindings generator**, see :doc:`m-generation` for details.

Overview
--------

Given a configuration file, the generator produces the following files:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Output file
     - Description
   * - ``cpp/models/<prefix>/infection_state.h``
     - C++ enum ``InfectionState`` with all compartments
   * - ``cpp/models/<prefix>/parameters.h``
     - Parameter structs, ``ParametersBase``, and ``Parameters`` class with constraint checks
   * - ``cpp/models/<prefix>/model.h``
     - ``Model`` class with the ``get_flows()`` implementation
   * - ``cpp/models/<prefix>/model.cpp``
     - Minimal translation unit (includes ``model.h``)
   * - ``cpp/models/<prefix>/CMakeLists.txt``
     - CMake library target for the new model
   * - ``pycode/memilio-simulation/memilio/simulation/bindings/models/<prefix>.cpp``
     - pybind11 module transferring the model to Python
   * - ``pycode/examples/simulation/<prefix>_simple.py``
     - Ready-to-run Python simulation example

In addition, the two existing CMakeLists files
``cpp/CMakeLists.txt`` and ``pycode/memilio-simulation/CMakeLists.txt``
are **patched in-place** to register the new model.

Configuration file format
--------------------------

Both YAML and TOML are supported. YAML is recommended because it does not require quotes around
string values, making the files easier to read and write.

.. note::

    In TOML, all string values must be enclosed in quotes as this is a
    hard requirement of the TOML specification and cannot be avoided.

The configuration file has four sections:

model
~~~~~

Metadata about the model.

.. code-block:: yaml

    model:
      name: SEIR          # Human-readable name used in comments and doc-strings
      namespace: oseir    # Inner C++ namespace  ->  mio::oseir
      prefix: ode_seir    # Folder name and CMake target prefix

infection_states
~~~~~~~~~~~~~~~~

A list of compartment names. At least two are required and all names must be unique.
A ``Count`` auxiliary compartment is added automatically at the end of the C++ enum.

.. code-block:: yaml

    infection_states:
      - Susceptible
      - Exposed
      - Infected
      - Recovered

parameters
~~~~~~~~~~

A list of model parameters. Each parameter entry has the following fields:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Field
     - Required
     - Description
   * - ``name``
     - yes
     - C++ struct name, e.g. ``TransmissionProbabilityOnContact``
   * - ``description``
     - yes
     - Short description used in the Doxygen comment
   * - ``type``
     - yes
     - ``probability`` (scalar in [0,1]), ``time`` (positive duration in days), or ``custom``
   * - ``default``
     - yes
     - Default value passed to ``get_default()``
   * - ``per_age_group``
     - no
     - ``true`` (default) -> stored as ``CustomIndexArray<UncertainValue, AgeGroup>``; ``false`` -> plain ``UncertainValue``
   * - ``bounds``
     - no
     - ``[lower, upper]`` - use ``null`` for an open bound. Inferred from ``type`` if omitted.

**Built-in types and their bounds:**

- ``probability``: constraint check enforces ``[0.0, 1.0]``
- ``time``: constraint check enforces ``[0.1, ∞)`` with a tolerance warning to prevent near zero values.
- ``custom``: no automatic constraint check is generated

.. note::

    When at least one ``infection`` transition is present, a ``ContactPatterns`` parameter is
    added to the model **automatically**,  you do not need to declare it in the ``parameters``
    list. It stores the age-group contact matrix (``UncertainContactMatrix``) and is used by
    every infection transition to compute the force of infection.
    In the generated Python example and in your own scripts, set it up like this:

    .. code-block:: python

        model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones((num_groups, num_groups))
        model.parameters.ContactPatterns.cont_freq_mat[0].minimum  = np.zeros((num_groups, num_groups))

.. code-block:: yaml

    parameters:
      - name: TransmissionProbabilityOnContact
        description: probability of getting infected from a contact
        type: probability
        default: 1.0
        per_age_group: true
        bounds: [0.0, 1.0]

      - name: TimeExposed
        description: the latent time in day unit
        type: time
        default: 5.2
        per_age_group: true
        bounds: [0.1, null]

transitions
~~~~~~~~~~~

A list of compartment flows. Each transition has the following fields:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Field
     - Required
     - Description
   * - ``from``
     - yes
     - Source compartment (must be in ``infection_states``)
   * - ``to``
     - yes
     - Target compartment (must be in ``infection_states``, must differ from ``from``)
   * - ``type``
     - yes
     - ``infection``, ``linear``, or ``custom``
   * - ``parameter``
     - for ``infection`` and ``linear``
     - Name of the driving parameter (must be in ``parameters``)
   * - ``infectious_state``
     - for ``infection``
     - Compartment whose population drives the force of infection (e.g. ``Infected``)
   * - ``custom_formula``
     - no
     - Optional hint placed in a ``TODO`` comment in the generated code

**Transition types:**

``infection``
    Force-of-infection flow. Generates a double loop over contact age groups using the
    ``ContactPatterns`` contact matrix. The ``ContactPatterns`` parameter is added to the
    model automatically when at least one infection transition is present.

    .. math::

       \dot{S}_i \leftarrow -\sum_j c_{ij} \cdot \phi \cdot \frac{I_j}{N_j} \cdot S_i

    where :math:`c_{ij}` is the contact rate between age groups *i* and *j*,
    :math:`\phi` is the transmission probability, and :math:`N_j` is the total
    population of age group *j*.

``linear``
    Simple outflow proportional to the compartment size:

    .. math::

       \dot{X}_i \leftarrow -\frac{1}{\tau_i} \cdot X_i

    where :math:`\tau_i` is the time parameter for age group *i*.

``custom``
    A placeholder is inserted into ``get_flows()`` with a ``TODO`` comment.
    If ``custom_formula`` is provided, it is shown as a hint next to the placeholder.
    **The generated code will not compile until you fill in the expression.**

.. code-block:: yaml

    transitions:
      - from: Susceptible
        to: Exposed
        type: infection
        parameter: TransmissionProbabilityOnContact
        infectious_state: Infected

      - from: Exposed
        to: Infected
        type: linear
        parameter: TimeExposed

      - from: Infected
        to: Recovered
        type: linear
        parameter: TimeInfected

Complete example: SEIR model
-----------------------------

The following YAML file fully specifies an SEIR model:

.. code-block:: yaml

    model:
      name: SEIR
      namespace: oseir
      prefix: ode_seir

    infection_states:
      - Susceptible
      - Exposed
      - Infected
      - Recovered

    parameters:
      - name: TransmissionProbabilityOnContact
        description: probability of getting infected from a contact
        type: probability
        default: 1.0
        per_age_group: true
        bounds: [0.0, 1.0]

      - name: TimeExposed
        description: the latent time in day unit
        type: time
        default: 5.2
        per_age_group: true
        bounds: [0.1, null]

      - name: TimeInfected
        description: the infectious time in day unit
        type: time
        default: 6.0
        per_age_group: true
        bounds: [0.1, null]

    transitions:
      - from: Susceptible
        to: Exposed
        type: infection
        parameter: TransmissionProbabilityOnContact
        infectious_state: Infected

      - from: Exposed
        to: Infected
        type: linear
        parameter: TimeExposed

      - from: Infected
        to: Recovered
        type: linear
        parameter: TimeInfected

More example configurations (including an SEIRD model with a ``custom`` transition and a TOML
version of the SEIR model) can be found in
`pycode/examples/modelgenerator/ <https://github.com/SciCompMod/memilio/blob/main/pycode/examples/modelgenerator/>`_.

Usage
-----

Command-line interface
~~~~~~~~~~~~~~~~~~~~~~

The generator is installed as the command ``memilio-modelgenerator``:

.. code-block:: console

    # Write all files into the MEmilio repository root
    memilio-modelgenerator path/to/seir.yaml --output-dir /path/to/memilio

    # Preview the generated files without writing them
    memilio-modelgenerator path/to/seir.yaml --preview

    # TOML input works the same way
    memilio-modelgenerator path/to/seir.toml --output-dir /path/to/memilio

    # Overwrite an existing model directory (see warning below)
    memilio-modelgenerator path/to/seir.yaml --output-dir /path/to/memilio --force

Python API
~~~~~~~~~~

.. code-block:: python

    from memilio.modelgenerator import Generator

    # Load from YAML
    gen = Generator.from_yaml("seir.yaml")

    # Load from TOML
    gen = Generator.from_toml("seir.toml")

    # Load from a dict (useful in scripts or tests)
    gen = Generator.from_dict(raw_dict)

    # Render all files to a dict {relative_path: content}
    files = gen.render()

    # Write all files and patch existing CMakeLists
    gen.write("/path/to/memilio")

    # Overwrite an existing model directory (see warning below)
    gen.write("/path/to/memilio", overwrite=True)

.. warning::

    The generator refuses to write into a model directory that already exists
    (``cpp/models/<prefix>/``) unless ``overwrite=True`` (Python API) or ``--force``
    (CLI) is passed explicitly.
    This guard is intentional: ``prefix`` and ``namespace`` must be unique across the
    whole MEmilio repository. Using the same values as an existing model (e.g.
    ``prefix: ode_seir``) would replace a handwritten C++ source file of
    that model with generated ones.

After generation
~~~~~~~~~~~~~~~~

1. **Fill in custom transitions** (if any): open the generated ``model.h`` and replace the
   ``/* YOUR EXPRESSION HERE */`` placeholder with the actual expression before compiling.

2. **Compile the model** by building the MEmilio C++ library as usual (CMake).
   The patched ``cpp/CMakeLists.txt`` picks up the new model directory automatically.

3. **Install the Python bindings** by reinstalling ``memilio-simulation``:

   .. code-block:: console

       pip install -e pycode/memilio-simulation

4. **Run the generated example**:

   .. code-block:: console

       python pycode/examples/simulation/<prefix>_simple.py

Validation
----------

The generator validates the configuration before any code is produced.
All errors are collected and reported together.

Common validation errors:

* Missing or empty ``model``, ``infection_states``, ``parameters``, or ``transitions`` section
* Fewer than 2 infection states, or duplicate state names
* Parameter ``type`` is not one of ``probability``, ``time``, ``custom``
* ``parameter`` or ``infectious_state`` in a transition references an unknown name
* A transition has the same ``from`` and ``to`` state (self-loop)

Development
-----------

Adding a new transition type or template feature:

1. Add the new type constant to ``TransitionType`` in
   `schema.py <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/modelgenerator/schema.py>`_.
2. Add validation logic to
   `validator.py <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/modelgenerator/validator.py>`_.
3. Update the relevant Jinja2 templates under
   `templates/ <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/memilio/modelgenerator/templates/>`_.
4. Update the tests in
   `tests/test_modelgenerator.py <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-generation/tests/test_modelgenerator.py>`_.
