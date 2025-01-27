Usage
=====

.. _installation:

Installation
------------

To use MEmilio, first install it using ...

.. code-block:: console

   c++ init

C++ Tutorial
----------------

TBD

.. autofunction:: lumache.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import memilio.epidata import progress_indicator


Other examples can be found in the :doc:`models` page.


