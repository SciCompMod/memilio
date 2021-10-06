# MEmilio C++ modelling framework #

This directory contains the framework used to build all our epidemiological models. 

Subdirectories:
- utils: Generic utility functions, i.e. not specific to epidemiological programms.
- io: Utility for writing and reading data to and from files in different formats.
- math: Math utility.
- epidemiology: Classes and functions that are useful in any type of epidemiological model.
- compartments: Classes and functions that can be used to create compartment (i.e. SIR-type) models.
- mobility: Model-agnostic loose coupling to model mobility between different model instances that represent e.g. geographic regions.