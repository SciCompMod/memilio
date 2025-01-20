import argparse
import os

import numpy as np
from memilio.simulation import AgeGroup, SimulationDay, Date
from memilio.simulation.osecirvvs import InfectionState, Model, simulate
from memilio.simulation.osecirvvs import set_vaccination_data, read_vaccination_data, set_divi_data, set_confirmed_cases_data, set_population_data, set_vaccination_data_from_entries


def run_ode_secirvvs_data_example(data_dir, show_plot=True):

    t0 = 0
    tmax = 30
    dt = 0.1
    num_groups = 6
    regions = [1001]
    scaling_factor_inf = [1.0]
    scaling_factor_icu = 1.0
    set_deaths = False

    # Initialize the model
    model = [Model(num_groups)]  # A list of models for each region

    # Set up real-world data for the simulation
    vaccination_data_path = os.path.join(
        data_dir, "vacc_county_ageinf_ma7.json")
    divi_data_path = os.path.join(data_dir, "county_divi_ma7.json")
    confirmed_cases_path = os.path.join(
        data_dir, "cases_all_county_age_ma7.json")
    population_data_path = os.path.join(
        data_dir, "county_current_population.json")

    vacc_data = read_vaccination_data(vaccination_data_path)

    date = Date(2022, 1, 1)
    set_vaccination_data_from_entries(
        model, vacc_data, date, regions, 10)
    # set_divi_data(model, divi_data_path, regions, date, scaling_factor_icu)
    # set_confirmed_cases_data(model, confirmed_cases_path,
    #                          regions, date, scaling_factor_inf, set_deaths)
    # set_population_data(model, population_data_path,
    #                     confirmed_cases_path, regions, date)

    for m in model:
        m.apply_constraints()

    # Run the simulation
    result = simulate(t0, tmax, dt, model[0])


if __name__ == "__main__":
    data_dir = "/localdata1/code_2024/memilio/data/pydata/Germany/"
    run_ode_secirvvs_data_example(data_dir)
