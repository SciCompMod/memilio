#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Henrik Zunker
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
"""
Example: For a given country, this example loads the age-structured contact
matrix, resolves the total population, builds a simple age-resolved ODE SEIR
model, and runs a simulation.
"""

import io

import numpy as np
import pandas as pd
import requests

from memilio.epidata.getContactData import (load_contact_matrix)
from memilio.simulation import AgeGroup, ContactMatrix
from memilio.simulation.oseir import InfectionState as State
from memilio.simulation.oseir import (Model, interpolate_simulation_result,
                                      simulate)

POPULATION_URL = (
    "https://raw.githubusercontent.com/kieshaprem/synthetic-contact-matrices/"
    "master/generate_synthetic_matrices/input/pop/popage_total2020.csv"
)


def get_population_by_age(country_name: str):
    """
    Loads population data for a specific country from a GitHub source (based on UN data).
    Returns the population in 5-year steps (Unit: number of people).

    Source: https://github.com/kieshaprem/synthetic-contact-matrices
    Note: Data is originally in thousands, converted here to absolute numbers.
    """

    try:
        # download the file and save only in variable
        response = requests.get(POPULATION_URL, timeout=10)
        response.raise_for_status()
        df = pd.read_csv(io.StringIO(response.text))

        # clean column names
        df.columns = df.columns.str.strip()

        # Filter by country (case-insensitive)
        country_col = "Region, subregion, country or area *"
        row = df[df[country_col].str.lower() == country_name.lower()]

        if row.empty:
            raise ValueError(f"Country '{country_name}' not found in data.")

        # Extract relevant columns
        age_cols = [c for c in df.columns if c.startswith('age')]
        pop_data = row.iloc[0][age_cols].astype(float)

        # Convert from thousands to absolute numbers
        pop_data = pop_data * 1000

        # --- Aggregation for Memilio / Prem Matrices (usually up to 75+) ---
        # Prem matrices in memilio often end at 75+.
        # The CSV goes up to 100. We sum everything from 75 upwards.

        # Columns up to 70 (0-4, ..., 70-74)
        cols_up_to_70 = [f'age{i}' for i in range(0, 75, 5)]

        # Columns from 75 (75-79, ..., 100+)
        cols_75_plus = [
            f'age{i}' for i in range(75, 105, 5)
            if f'age{i}' in pop_data.index]

        final_pop = []

        # 1. Keep groups up to 70-74
        for col in cols_up_to_70:
            final_pop.append(pop_data[col])

        # 2. Sum groups from 75+
        sum_75_plus = pop_data[cols_75_plus].sum()
        final_pop.append(sum_75_plus)

        return final_pop

    except Exception as e:
        raise RuntimeError(f"Error loading population data: {e}")


def build_country_seir_model(
        contact_matrix: np.ndarray,
        population_by_age: list,
        transmission_probability: float = 0.06,
        exposed_share: float = 1e-5,
        infected_share: float = 5e-6):
    """
    Build a simple age-resolved ODE SEIR model.
    """
    contact_matrix = np.asarray(contact_matrix, dtype=float)
    num_groups = contact_matrix.shape[0]

    model = Model(num_groups)
    contacts = ContactMatrix(contact_matrix)
    contacts.minimum = np.zeros_like(contact_matrix)
    model.parameters.ContactPatterns.cont_freq_mat[0] = contacts

    # Use actual population distribution
    group_pop = np.array(population_by_age)

    if len(group_pop) != num_groups:
        # If dimensions mismatch (e.g. contact matrix has different age groups),
        # we might need to adjust. For now, we assume they match (16 groups for Prem 75+).
        # If not, we fallback to uniform distribution of total sum (not ideal but safe).
        print(
            f"Warning: Population groups ({len(group_pop)}) do not match contact matrix groups ({num_groups}). Using uniform distribution.")
        total_pop = np.sum(group_pop)
        group_weights = np.full(num_groups, 1.0 / num_groups)
        group_pop = group_weights * total_pop

    exposed_init = np.maximum(1.0, group_pop * exposed_share)
    infected_init = np.maximum(1.0, group_pop * infected_share)

    # Set parameters and initial conditions equal for all age groups.
    for idx in range(num_groups):
        age_group = AgeGroup(idx)
        model.parameters.TimeExposed[age_group] = 5.2
        model.parameters.TimeInfected[age_group] = 6.0
        model.parameters.TransmissionProbabilityOnContact[age_group] = (
            transmission_probability)

        model.populations[age_group, State.Exposed] = exposed_init[idx]
        model.populations[age_group, State.Infected] = infected_init[idx]
        model.populations[age_group, State.Recovered] = 0.0
        model.populations.set_difference_from_group_total_AgeGroup(
            (age_group, State.Susceptible), group_pop[idx])

    model.check_constraints()
    return model


def simulate_country_seir(
        country: str,
        days: float = 120.0,
        dt: float = 0.25,
        transmission_probability: float = 0.06,
        exposed_share: float = 1e-5,
        infected_share: float = 5e-6,
        interpolate: bool = True):
    """
    Load contact matrix, fetch population, build the ODE SEIR model, and run
    a simulation.
    """
    contacts = load_contact_matrix(country, reduce_to_rki_groups=False)
    population = get_population_by_age(country)
    model = build_country_seir_model(
        contacts.values,
        population_by_age=population,
        transmission_probability=transmission_probability,
        exposed_share=exposed_share,
        infected_share=infected_share)

    result = simulate(t0=0.0, tmax=days, dt=dt, model=model)
    if interpolate:
        result = interpolate_simulation_result(result)

    return result


def run_demo(country: str,
             days: float = 120.0,
             dt: float = 0.25,
             transmission_probability: float = 0.06,
             exposed_share: float = 1e-5,
             infected_share: float = 5e-6):
    """
    Run the SEIR simulation demo for a user defined country and parameters.
    """
    result = simulate_country_seir(
        country,
        days=days,
        dt=dt,
        transmission_probability=transmission_probability,
        exposed_share=exposed_share,
        infected_share=infected_share,
    )
    print(result.get_last_value())


if __name__ == "__main__":
    country = "Germany"
    run_demo(country=country)
