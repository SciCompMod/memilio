MEmilio Generation
===================

The ``memilio-generation`` package contains two independent tools, each with its own dependencies and installation requirements:

* **Model Generator:** generates a C++ compartmental model with Python bindings from a YAML/TOML configuration file. See :doc:`m-modelgenerator`.
* **Bindings Generator:** automatically generates Python bindings from existing C++ model source files using libclang. See :doc:`m-bindingsgenerator`.

.. toctree::
   :maxdepth: 1
   :hidden:

   m-modelgenerator
   m-bindingsgenerator

Installation
------------

.. code-block:: console

    # Model Generator only (no libclang required)
    pip install -e pycode/memilio-generation

    # Model Generator + Bindings Generator (installs libclang)
    pip install -e "pycode/memilio-generation[bindings]"

The Model Generator only requires ``jinja2``, ``pyyaml``, and (on Python < 3.11) ``tomli``.
The Bindings Generator additionally needs ``libclang==18.1.1`` and the
``python3.x-dev`` system libraries (see :doc:`m-bindingsgenerator` for details).
