FAQ
====

.. _model-faq:

What do you mean when talking about a *model* (implementation)?
--------------------------------------------------------------------

That's a very good question as *model* is used as short hand for different model abstraction or realization levels. The reason for the complexity in MEmilio is the modular structure and templated design which allows to transfer abstraction layers or abstract models to many different applications. With this, novel applications can be built with small overhead, benefitting from prior set-ups and scenarios, eventually realizing swift reaction to novel threats by emerging pathogens.
 
The question can hopefully be explained best from the end-user perspective. A particular example as implemented or used by an end-user and , e.g., returning concrete numbers of new infections on different simulation days is a *particular* model *realization* of a *templated* model such as ODE-SEIR in the C++ model folder. This templated model does not yet define the stratification into a precise number of age groups. Furthermore, while it provides default values for parameters, the model realization should be parametrized for the particular use case as the templated model is considered to be generic for many different use cases or pathogens with the same transmission pathway. In addition, a model *realization* can be a combination of two or more model realizations. For instance, we can use a graph structure to model the mobility and use different model realizations of ODE-SEIR models in different nodes of the graph. We then have one realization of the graph for modeling the mobility and two model realizations of ODE-SEIR type for two different regions. For an example, see here: `https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secir_graph.cpp`_

If it will not be clear from the context, feel free to reach out directly.

I found a bug. What should I do?
--------------------------------

Please open an issue on our `GitHub <https://github.com/SciCompMod/memilio/issues>`_ page. Please try to describe the bug as detailed as possible, including the steps to reproduce it. 

I solved a bug. What should I do?
---------------------------------

Please open a pull request on our `GitHub <https://github.com/SciCompMod/memilio/pulls>`_ page. If there is not yet a bug report, please open an issue first and reference it in your pull request. For more information, please refer to :doc:`development`.