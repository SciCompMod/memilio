# Graph Neural Network (GNN) Surrogate Models

This module implements Graph Neural Network-based surrogate models for epidemiological simulations, specifically designed to accelerate and enhance pandemic response modeling.

## Overview

The GNN surrogate models are based on the research presented in:

**Schmidt A, Zunker H, Heinlein A, Kühn MJ. (2025).** *Graph Neural Network Surrogates to leverage Mechanistic Expert Knowledge towards Reliable and Immediate Pandemic Response*. Submitted for publication.  
https://doi.org/10.48550/arXiv.2411.06500

This implementation leverages the underlying [ODE SECIR model](https://memilio.readthedocs.io/en/latest/cpp/models/osecir.html)  and applies Graph Neural Networks to create fast, reliable surrogate models that can be used for immediate pandemic response scenarios. The models are stratified by age groups and incorporate spatial connectivity through mobility data.

## Module Structure

The GNN module consists of the following components:

- **`data_generation.py`**: Generates training and evaluation data for GNN surrogate models by simulating epidemiological scenarios using the mechanistic SECIR model. Handles parameter sampling, compartment initialization, damping factors, and mobility connections between regions.

- **`network_architectures.py`**: Defines various GNN architectures using different layer types (e.g., Graph Convolutional Networks, Graph Attention Networks). Provides functionality to configure network depth, width, activation functions, and preprocessing transformations.

- **`evaluate_and_train.py`**: Implements training and evaluation pipelines for GNN surrogate models. Loads generated data, trains models with specified hyperparameters, evaluates performance metrics, and saves trained model weights.

- **`grid_search.py`**: Provides hyperparameter optimization through systematic grid search over network architectures, training configurations, and model parameters to identify optimal GNN configurations for epidemiological forecasting.

- **`GNN_utils.py`**: Contains utility functions for data preprocessing (e.g., removing confirmed compartments, scaling data), and building graphs or baseline contact matrices used throughout the GNN workflow.

## Documentation

Comprehensive documentation for the GNN surrogate models, including tutorials and usage examples, is available in our [documentation](https://memilio.readthedocs.io/en/latest/python/m-surrogate.html).
