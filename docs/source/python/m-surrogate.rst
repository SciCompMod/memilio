MEmilio Surrogate Model
========================

MEmilio Surrogate Model contains machine learning based surrogate models that make predictions based on the MEmilio simulation models. 
Currently there are only surrogate models for ODE-type models. The simulations of these models are used for data generation. 
The goal is to create a powerful tool that predicts the infection dynamics faster than a simulation of an expert model, 
e.g., a metapopulation or agent-based model while still having acceptable errors with respect to the original simulations.
 
The package can be found in `pycode/memilio-surrogatemodel <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-surrogatemodel>`_.

For more details, we refer to: Schmidt A, Zunker H, Heinlein A, Kühn MJ. (2025). *Graph Neural Network Surrogates to leverage Mechanistic Expert Knowledge towards Reliable and Immediate Pandemic Response*. Submitted for publication. `arXiv:2411.06500 <https://arxiv.org/abs/2411.06500>`_

Dependencies
------------

Required python packages:

- pandas >= 1.2.2
- numpy >= 1.22, !=1.25.*
- tensorflow
- matplotlib
- scikit-learn
- progress

Since we are running simulations to generate the data, the MEmilio :doc:`memilio-simulation <m-simulation>` package also needs to be installed.

Usage
-----

The package currently provides the following modules:

- `models`: models for different tasks
   Currently we have the following models: 
   - `ode_secir_simple`: A simple model allowing for asymptomatic as well as symptomatic infection states not stratified by age groups.
   - `ode_secir_groups`: A model allowing for asymptomatic as well as symptomatic infection states stratified by age groups and including one damping.

     Each model folder contains the following files: 
     - `data_generation`: Data generation from expert model simulation outputs.
     - `model`: Training and evaluation of the model. 
     - `network_architectures`: Contains multiple network architectures.
     - `grid_search`: Utilities for hyperparameter optimization.
     - `hyperparameter_tuning`: Scripts for tuning model hyperparameters.


- `tests`: This file contains all tests.

ODE-SECIR Simple Model
----------------------

The `ode_secir_simple` module provides surrogate models for the basic ODE-SECIR epidemiological model. This model is not stratified by age groups and simulates disease progression through the following compartments:

- **S**: Susceptible
- **E**: Exposed
- **C**: Infected (asymptomatic/pre-symptomatic)
- **I**: Infected (symptomatic)
- **R**: Recovered
- **H**: Hospitalized (severe cases)
- **U**: ICU (critical cases)
- **D**: Dead

For more details on the model structure and parameters, we refer to the ODE-SECIR model documentation.

Data Generation
~~~~~~~~~~~~~~~

The `data_generation.py` module provides functionality to generate training data for the surrogate models by running multiple simulations of the basic ODE-SECIR model with randomized initial conditions. The data generation process involves:

.. code-block:: python

    # Generate dataset with 10,000 simulation runs
    # Each with 5 days of input data and 30 days of prediction horizon
    data = generate_data(
        num_runs=10000,
        path=path_data,
        input_width=5,
        label_width=30,
        normalize=True,
        save_data=True
    )

The data generation process can be summarized as follows:

1. Randomly initializes the model parameters and initial compartment populations.
2. Runs the ODE-SECIR simulation using the C++ backend via Python bindings.
3. Applies logarithmic normalization to improve training stability.
4. Splits each time series into input and label segments.
5. Saves the dataset as a pickle file for later use.

Network Architectures
~~~~~~~~~~~~~~~~~~~~~

The `network_architectures.py` module provides different neural network architectures for time series prediction:

1. **MLP (Multi-Layer Perceptron)**:

   - Simple feedforward networks that take flattened time series as input
   - Available in both single-output and multi-output variants
   
2. **LSTM (Long Short-Term Memory)**:

   - Recurrent neural networks specialized for sequence modeling
   - Can process variable-length time series while maintaining temporal information
   
3. **CNN (Convolutional Neural Network)**:

   - Uses 1D convolutions to detect patterns in time series data
   - Particularly efficient for capturing local temporal patterns

Model Training and Evaluation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `model.py` module provides functionality for:

1. **Preparing data**:

   - Splitting data into training, validation, and test sets
   - Processing data for different model architectures (classic vs. time series)
   
2. **Model training**:

   - Initializing models with customizable hyperparameters
   - Training with early stopping and customizable loss functions
   
3. **Evaluation**:

   - Computing error metrics (MAE, MAPE) across compartments
   - Visualizing predictions versus ground truth

Example usage:

.. code-block:: python

    # Define model and training parameters
    model_parameters = (label_width, num_outputs, hidden_layers, 
                       neurons_per_layer, activation, modelname)
    training_parameters = (early_stop, max_epochs, loss, optimizer, metrics)
    
    # Initialize and train model
    model = initialize_model(model_parameters)
    history = network_fit(model, modeltype, training_parameters, path_data)
    
    # Plot results
    plot_compartment_prediction_model(test_inputs, test_labels, 
                                     modeltype, model, 'InfectedSymptoms')

Hyperparameter Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `grid_search.py` and `hyperparameter_tuning.py` modules provide tools for systematic hyperparameter optimization:

1. **Cross-validation**:

   - K-fold cross-validation to prevent overfitting
   - Evaluation of multiple model architectures and training configurations

2. **Grid search**:

   - Systematic exploration of hyperparameter space
   - Tracking and storage of performance metrics
   
3. **Result analysis**:

   - Visualization of hyperparameter importance
   - Selection of optimal model configurations

SECIR Groups Model
------------------

To be added...
