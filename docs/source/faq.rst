FAQ
====

What do you mean when talking about a model implementation?
--------------------------------------------------------------------

That's a very good question as a model for us can mean both a modelling framework as well as an actual (example) implementation. 
We seperate between the abstract model, which e.g. could be a SEIR-ODE-Model and would be put into a model
folder (in the C++ code), and the implemented model which fills the abstract model with parameter values, time spans etc. However, the model
framework is also implemented, so it's also an implemented model... Hopefully it will be clear from the context what exactly
we are talking about. 

I found a bug. What should I do?
--------------------------------

Please open an issue on our `GitHub <https://github.com/SciCompMod/memilio/issues>`_ page. Please try to describe the bug as detailed as possible, including the steps to reproduce it. 

I solved a bug. What should I do?
---------------------------------

Please open a pull request on our `GitHub <https://github.com/SciCompMod/memilio/pulls>`_ page. If there is not yet a bug report, please open an issue first and reference it in your pull request. For more information, please refer to :doc:`development`.