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
- numpy >= 1.22, < 1.25
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



Graph Neural Network (GNN) Surrogate Models
--------------------------------------------

The Graph Neural Network (GNN) module provides advanced surrogate models that leverage spatial connectivity and age-stratified epidemiological dynamics. These models are designed for immediate and reliable pandemic response by combining mechanistic expert knowledge with machine learning efficiency.

Overview and Scientific Foundation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GNN surrogate models are based on the research presented in:

|Graph_Neural_Network_Surrogates|

The implementation leverages the mechanistic ODE-SECIR model (see :doc:`ODE-SECIR documentation <../models/ode_secir>`) as the underlying expert model, using Python bindings to the C++ backend for efficient simulation during data generation.

Module Structure
~~~~~~~~~~~~~~~~

The GNN module is located in `pycode/memilio-surrogatemodel/memilio/surrogatemodel/GNN <https://github.com/SciCompMod/memilio/tree/main/pycode/memilio-surrogatemodel/memilio/surrogatemodel/GNN>`_ and consists of:

- **data_generation.py**: Generates training and evaluation data by simulating epidemiological scenarios with the mechanistic SECIR model
- **network_architectures.py**: Defines various GNN architectures (GCN, GAT, GIN) with configurable layers and preprocessing
- **evaluate_and_train.py**: Implements training and evaluation pipelines for GNN models
- **grid_search.py**: Provides hyperparameter optimization through systematic grid search
- **GNN_utils.py**: Contains utility functions for data preprocessing, graph construction, and population data handling

Data Generation
~~~~~~~~~~~~~~~

The data generation process in ``data_generation.py`` creates graph-structured training data through mechanistic simulations:

.. code-block:: python

    from memilio.surrogatemodel.GNN import data_generation
    
    # Generate training dataset
    dataset = data_generation.generate_dataset(
        num_runs=1000,                    # Number of simulation scenarios
        num_days=30,                      # Simulation horizon
        num_age_groups=6,                 # Age stratification
        data_dir='path/to/contact_data',  # Contact matrices location
        mobility_dir='path/to/mobility',  # Mobility data location
        save_path='gnn_training_data.pickle'
    )

**Data Generation Workflow:**

1. **Parameter Sampling**: Randomly sample epidemiological parameters (transmission rates, incubation periods, recovery rates) from predefined distributions to create diverse scenarios.

2. **Compartment Initialization**: Initialize epidemic compartments for each age group in each region based on realistic demographic data. Compartments are initialized using shared base factors.

3. **Mobility Graph Construction**: Build a spatial graph where:
   
   - Nodes represent geographic regions (e.g., German counties)
   - Edges represent mobility connections with weights from commuting data
   - Node features include age-stratified population sizes

4. **Contact Matrix Configuration**: Load and configure baseline contact matrices for different location types (home, school, work, other) stratified by age groups.

5. **Damping Application**: Apply time-varying dampings to contact matrices to simulate NPIs:
   
   - Multiple damping periods with random start days
   - Location-specific damping factors (e.g., stronger school closures, moderate workplace restrictions)
   - Realistic parameter ranges based on observed intervention strengths

6. **Simulation Execution**: Run the mechanistic ODE-SECIR model using MEmilio's C++ backend through Python bindings to generate the dataset.

7. **Data Processing**: Transform simulation results into graph-structured format:
   
   - Extract compartment time series for each node (region) and age group
   - Apply logarithmic transformation for numerical stability
   - Store graph topology, node features, and temporal sequences

Network Architectures
~~~~~~~~~~~~~~~~~~~~~

The ``network_architectures.py`` module provides flexible GNN model construction for different layer types.

.. code-block:: python

    from memilio.surrogatemodel.GNN import network_architectures
    
    # Define GNN architecture
    model_config = {
        'layer_type': 'GCN',           # GNN layer type
        'num_layers': 3,                # Network depth
        'hidden_dim': 64,               # Hidden layer dimensions
        'activation': 'relu',           # Activation function
        'dropout_rate': 0.2,            # Dropout for regularization
        'use_batch_norm': True,         # Batch normalization
        'aggregation': 'mean',          # Neighborhood aggregation method
    }
    
    # Build model
    model = network_architectures.build_gnn_model(
        config=model_config,
        input_shape=(num_timesteps, num_features),
        output_dim=num_compartments * num_age_groups
    )


Training and Evaluation
~~~~~~~~~~~~~~~~~~~~~~~

The ``evaluate_and_train.py`` module provides the training functionality:

.. code-block:: python

    from memilio.surrogatemodel.GNN import evaluate_and_train
    
    # Load training data
    with open('gnn_training_data.pickle', 'rb') as f:
        dataset = pickle.load(f)
    
    # Define training configuration
    training_config = {
        'epochs': 100,
        'batch_size': 32,
        'learning_rate': 0.001,
        'optimizer': 'adam',
        'loss_function': 'mse',
        'early_stopping_patience': 10,
        'validation_split': 0.2
    }
    
    # Train model
    history = evaluate_and_train.train_gnn_model(
        model=model,
        dataset=dataset,
        config=training_config,
        save_weights='best_gnn_model.h5'
    )
    
    # Evaluate on test set
    metrics = evaluate_and_train.evaluate_model(
        model=model,
        test_data=test_dataset,
        metrics=['mae', 'mape', 'r2']
    )

**Training Features:**

1. **Mini-batch Training**: Graph batching for efficient training on large datasets
2. **Custom Loss Functions**: MSE, MAE, MAPE, or custom compartment-weighted losses
3. **Early Stopping**: Monitors validation loss to prevent overfitting
4. **Learning Rate Scheduling**: Adaptive learning rate reduction on plateaus
5. **Save Best Weights**: Saves best model weights based on validation performance

**Evaluation Metrics:**

- **Mean Absolute Error (MAE)**: Average absolute prediction error per compartment
- **Mean Absolute Percentage Error (MAPE)**: Mean absolute error as percentage
- **R² Score**: Coefficient of determination for prediction quality

**Data Splitting:**

- **Training Set (70%)**: For model parameter optimization
- **Validation Set (15%)**: For hyperparameter tuning and early stopping
- **Test Set (15%)**: For final performance evaluation

Hyperparameter Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``grid_search.py`` module enables systematic exploration of hyperparameter space:

.. code-block:: python

    from memilio.surrogatemodel.GNN import grid_search
    
    # Define search space
    param_grid = {
        'layer_type': ['GCN', 'GAT', 'GIN'],
        'num_layers': [2, 3, 4, 5],
        'hidden_dim': [32, 64, 128, 256],
        'learning_rate': [0.001, 0.0005, 0.0001],
        'dropout_rate': [0.0, 0.1, 0.2, 0.3],
        'batch_size': [16, 32, 64],
        'activation': ['relu', 'elu', 'tanh']
    }
    
    # Run grid search with cross-validation
    results = grid_search.run_hyperparameter_search(
        param_grid=param_grid,
        data_path='gnn_training_data.pickle',
        cv_folds=5,
        metric='mae',
        save_results='grid_search_results.csv'
    )
    
    # Analyze best configuration
    best_config = grid_search.get_best_configuration(results)
    print(f"Best configuration: {best_config}")

Utility Functions
~~~~~~~~~~~~~~~~~

The ``GNN_utils.py`` module provides essential helper functions used throughout the GNN workflow:

**Data Preprocessing:**

.. code-block:: python

    from memilio.surrogatemodel.GNN import GNN_utils
    
    # Remove confirmed compartments (simplify model)
    simplified_data = GNN_utils.remove_confirmed_compartments(
        dataset_entries=dataset,
        num_groups=6
    )
    
    # Apply logarithmic scaling
    scaled_inputs, scaled_labels = GNN_utils.scale_data(
        data=dataset,
        transform=True
    )

**Graph Construction:**

.. code-block:: python

    # Create mobility graph from commuting data
    graph = GNN_utils.create_mobility_graph(
        mobility_dir='path/to/mobility',
        num_regions=401,            # German counties
        county_ids=county_list,
        models=models_per_region    # SECIR models for each region
    )
    
    # Get baseline contact matrix
    contact_matrix = GNN_utils.get_baseline_contact_matrix(
        data_dir='path/to/contact_matrices'
    )

Practical Usage Example
~~~~~~~~~~~~~~~~~~~~~~~

Here is a complete example workflow from data generation to model evaluation:

.. code-block:: python

    import pickle
    from pathlib import Path
    from memilio.surrogatemodel.GNN import (
        data_generation, 
        network_architectures, 
        evaluate_and_train
    )
    
    # Step 1: Generate training data
    print("Generating training data...")
    dataset = data_generation.generate_dataset(
        num_runs=5000,
        num_days=30,
        num_age_groups=6,
        data_dir='/path/to/memilio/data/Germany',
        mobility_dir='/path/to/mobility_data',
        save_path='gnn_dataset_5000.pickle'
    )
    
    # Step 2: Define and build GNN model
    print("Building GNN model...")
    model_config = {
        'layer_type': 'GCN',
        'num_layers': 4,
        'hidden_dim': 128,
        'activation': 'relu',
        'dropout_rate': 0.2,
        'use_batch_norm': True
    }
    
    model = network_architectures.build_gnn_model(
        config=model_config,
        input_shape=(1, 48),  # 6 age groups × 8 compartments
        output_dim=48         # Predict all compartments
    )
    
    # Step 3: Train the model
    print("Training model...")
    training_config = {
        'epochs': 200,
        'batch_size': 32,
        'learning_rate': 0.001,
        'optimizer': 'adam',
        'loss_function': 'mae',
        'early_stopping_patience': 20,
        'validation_split': 0.2
    }
    
    history = evaluate_and_train.train_gnn_model(
        model=model,
        dataset=dataset,
        config=training_config,
        save_weights='gnn_weights_best.h5'
    )
    
    # Step 4: Evaluate on test data
    print("Evaluating model...")
    test_metrics = evaluate_and_train.evaluate_model(
        model=model,
        test_data='gnn_test_data.pickle',
        metrics=['mae', 'mape', 'r2']
    )
    
    # Print results
    print(f"Test MAE: {test_metrics['mae']:.4f}")
    print(f"Test MAPE: {test_metrics['mape']:.2f}%")
    print(f"Test R²: {test_metrics['r2']:.4f}")
    
    # Step 5: Make predictions on new scenarios
    with open('new_scenario.pickle', 'rb') as f:
        new_data = pickle.load(f)
    
    predictions = model.predict(new_data)
    print(f"Predictions shape: {predictions.shape}")

**GPU Acceleration:**

- TensorFlow automatically uses GPU when available
- Spektral layers are optimized for GPU execution
- Training time can be heavily reduced with appropriate GPU hardware

Additional Resources
~~~~~~~~~~~~~~~~~~~~

**Code and Examples:**

- `GNN Module <https://github.com/SciCompMod/memilio/tree/main/pycode/memilio-surrogatemodel/memilio/surrogatemodel/GNN>`_
- `GNN README <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-surrogatemodel/memilio/surrogatemodel/GNN/README.md>`_
- `Test Scripts <https://github.com/SciCompMod/memilio/tree/main/pycode/memilio-surrogatemodel/memilio/surrogatemodel_test>`_

**Related Documentation:**

- :doc:`ODE-SECIR Model <../models/ode_secir>`
- :doc:`MEmilio Simulation Package <m-simulation>`
- :doc:`Python Bindings <python_bindings>`
