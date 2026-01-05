# Optimization of dynamic NPIs for the SECIR Graph Model

This directory contains a C++ implementation to optimize dynamic non-pharmaceutical interventions (NPIs) for a graph SECIR epidemiological model.
The focus is on optimizing damping-based interventions while respecting path and terminal constraints over a simulation horizon.

The simulation is based on the results of `optimization_model/results_run0.h`.
This file supplies the initial population data at time t=60 and we optimize further dynamic NPIs up to t=120.

The strength and costs of NPIs count towards the objective function, when the incidence is above the dynamic NPI threshold.
It is possible to define multiple dynamic NPIs for a single treshold, though in this example we just use a single damping for treshold_A = 50, threshold_B = 250 based on physical distancing.

It is important to choose the integrator resolution = 20 since otherwise the gradient information is not accurate enough.
