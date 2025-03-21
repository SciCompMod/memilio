ODE model creation
==================

Model definition
----------------

Before implementing a model in MEmilio, we need to do a some math, in particular, define an initial value problem
given by a system of ordinary differential equations. For example we consider a SIRD model given by

.. math::  
    \begin{aligned}
        S'(t) & = -\rho\phi\ \frac{S(t)*I(t)}{N_{\perp D}} \\
        I'(t) & = \rho\phi\ \frac{S(t)*I(t)}{N_{\perp D}} - \frac{\mu_R + \mu_D}{T_I}I(t) \\
        R'(t) & = \frac{\mu_R}{T_I}I(t) \\
        D'(t) & = \frac{\mu_D}{T_I}I(t) \\
    \end{aligned}

and some initial values for :math:`t=0`. Here :math:`N_{\perp D} := S(t) + I(t) + R(t)`.

This type of model is called compartmental model, because the model population is represented by discrete infection
states **S** usceptible, **I** nfected, **R** ecovered, **D** ead, also called compartments.

How to define an ODE model
--------------------------
