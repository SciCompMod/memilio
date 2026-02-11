How To: Coupling MEmilio with BayesFlow
========================================

Simulation-based inference explanation

Through pymio the data simulator can be incorporated into the learning process enabling shared interfaces 
and even further enable online training, i.e. generating training data on the spot during the learning phase.

Most of the important includes for BayesFlow:

.. code-block:: python

    import numpy as np

    # ensure the backend is set
    import os
    if "KERAS_BACKEND" not in os.environ:
        # set this to "torch", "tensorflow", or "jax"
        os.environ["KERAS_BACKEND"] = "tensorflow"

    import keras
    import bayesflow as bf
    import pandas as pd

.. code-block:: python

    def prior():
        """Generates a random draw from the joint prior."""

        lambd = RNG.lognormal(mean=np.log(0.4), sigma=0.5)
        mu = RNG.lognormal(mean=np.log(1 / 8), sigma=0.2)
        D = RNG.lognormal(mean=np.log(8), sigma=0.2)
        I0 = RNG.gamma(shape=2, scale=20)
        psi = RNG.exponential(5)
        return {"lambd": lambd, "mu": mu, "D": D, "I0": I0, "psi": psi}


Define the simulator function with the MEmilio python model. We will use a simple ODE SIR model for this example.

.. code-block:: python

    import memilio.simulation as mio
    import memilio.simulation.osir as osir

    def simulate_sir_model(lambd: float, mu: float, D: float, I0: float, scale_I: float, phi_i: float, N: int = 83e6, T: int = 91, sim_lag: int = 15, eps: float = 1e-5) -> npt.NDArray[np.float64]:
        """Performs a forward simulation from the stationary SIR model given a random draw from the prior."""

        # Define boundaries of parameters
        I0 = max(1, np.round(I0))
        D = int(round(D))

        # Configure model
        num_groups = 1
        sir_model = osir.Model(num_groups)
        A0 = mio.AgeGroup(0)

        # sim_lag is the maximum number of days for the delay
        # we simulate sim_lag days before t_0 to have values for t_0 - D
        t_max = T + sim_lag

        # Initial conditions
        sir_model.populations[A0, osir.InfectionState.Infected] = I0
        sir_model.populations[A0, osir.InfectionState.Recovered] = 0
        sir_model.populations.set_difference_from_total(
            (A0, osir.InfectionState.Susceptible), N)

        # Initialize Parameters
        sir_model.parameters.TimeInfected[A0] = mu
        sir_model.parameters.TransmissionProbabilityOnContact[A0] = lambd

        # Check logical constraints to parameters, should we really use apply_constraints?!
        sir_model.apply_constraints()

        # Reported new cases
        I_data = np.zeros(T)
        fs_i = np.zeros(T)

        # Run Simulation
        integrator = mio.RKIntegratorCore(dt_max=1)
        (result, flows) = simulate_flows(
            0, t_max, 1, sir_model, integrator)

        # interpolate results
        flows = osir.interpolate_simulation_result(flows)

        # consistency check, None values get deleted in configure input
        try:
            assert flows.get_last_time() == t_max
        except AssertionError as e:
            print('Invalid value simulated...return nan')
            return np.stack(([np.nan] * T, )).T

        # Adding new cases with delay D
        # Note, we assume the same delay
        shifted_t0 = sim_lag-D_i+1
        I_data = np.diff(flows.as_ndarray()[
                         1, shifted_t0:shifted_t0+T+1])
        I_data = np.clip(I_data, 10 ** -14, N)

        # Compute lags
        fs_i = (1-f_i)*(1 -
                        np.abs(np.sin((np.pi/7) * np.arange(0, T, 1) - 0.5*phi_i)))

        # Compute weekly modulation
        I_data = (1-fs_i) * I_data

        # check for negative values
        try:
            scale = np.sqrt(I_data)*scale_I
            assert np.all(scale >= 0)
        except AssertionError as e:
            print('Invalid value simulated...return nan')
            return np.stack(([np.nan] * T, )).T

        # Add noise
        I_data = stats.t(df=4, loc=I_data, scale=np.sqrt(I_data)*scale_I).rvs()

        # bound all negative values to 0
        I_data = np.clip(I_data, 10 ** -14, N)
        return dict(cases=np.stack((I_data, )).T)


.. code-block:: python

    simulator = bf.make_simulator([prior, stationary_SIR])

    adapter = (
        bf.adapters.Adapter()
        .convert_dtype("float64", "float32")
        .as_time_series("cases")
        .concatenate(["lambd", "mu", "D", "I0", "scale_I", "phi_I"], into="inference_variables")
        .rename("cases", "summary_variables")
        # since all our variables are non-negative (zero or larger), the next call transforms them
        # to the unconstrained real space and can be back-transformed under the hood
        .log(["inference_variables", "summary_variables"], p1=True)
    )

.. code-block:: python

    summary_network = bf.networks.TimeSeriesNetwork(summary_dim=4)
    inference_network = bf.networks.CouplingFlow()

.. code-block:: python

    workflow = bf.BasicWorkflow(
        simulator=simulator,
        adapter=adapter,
        inference_network=inference_network,
        summary_network=summary_network,
    )

.. code-block:: python

    history = workflow.fit_online(epochs=100, batch_size=64)


Load data, first need to download them using epidata

.. code-block:: python
    
    def load_observation_data(date_data_begin: datetime.date, T: int, data_path: str) -> np.ndarray:
        """Helper function to load cumulative cases and transform them to new cases."""

        # Use correct corona data based on the model (either reporting or reference date)
        confirmed_cases_json = data_path
        confirmed_cases = pd.read_json(confirmed_cases_json)
        confirmed_cases = confirmed_cases.set_index('Date')

        date_data_end = date_data_begin + datetime.timedelta(T)
        cases_obs = np.array(
            confirmed_cases.loc[date_data_begin:date_data_end]
        ).flatten()
        new_cases_obs = np.diff(cases_obs)
        return new_cases_obs