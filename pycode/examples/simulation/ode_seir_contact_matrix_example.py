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

OWID_POPULATION_URL = (
    "https://raw.githubusercontent.com/owid/covid-19-data/master/"
    "scripts/input/un/population_latest.csv"
)


def _normalize_country_name(country: str):
    """Return a case-insensitive key without whitespace or punctuation."""
    return "".join(ch for ch in country.casefold() if ch.isalnum())


def _download_population_table(url: str = OWID_POPULATION_URL):
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    return pd.read_csv(io.StringIO(response.text))


def _pick_population_from_table(df: pd.DataFrame,
                                country: str):
    # Identify candidate population columns provided by the OWID/UN table.
    pop_cols = [c
                for c
                in ["population", "Population", "pop", "PopTotal", "Value"]
                if c in df.columns]
    if not pop_cols:
        return None

    # Try multiple possible country-name columns to find the first match.
    country_cols = [
        c
        for c
        in
        ["entity", "Entity", "location", "Location", "country", "Country",
         "country_name", "Country Name"] if c in df.columns]
    key = _normalize_country_name(country)
    for c_col in country_cols:
        normalized = df[c_col].astype(str).map(_normalize_country_name)
        matches = df.loc[normalized == key]
        if matches.empty:
            continue
        # If multiple years are present, prefer the latest year.
        matches = matches.sort_values("year", ascending=False)
        pop = matches.iloc[0][pop_cols[0]]
        if pd.notna(pop):
            return int(pop)
    return None


def get_population_total(country: str):
    """
    get a countrys total population using online data.
    """
    try:
        df_pop = _download_population_table()
        pop_online = _pick_population_from_table(
            df_pop, country)
        if pop_online:
            return pop_online
    except Exception:
        pass
    raise RuntimeError(
        "No population found. Provide population_override or extend population_fallback.")


def build_country_seir_model(
        contact_matrix: np.ndarray,
        population_total: int,
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

    # Distribute population and initial exposed/infected uniformly across age groups.
    group_weights = np.full(num_groups, 1.0 / num_groups)
    group_pop = group_weights * population_total
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
    a simulation. Returns (result, contacts_df, population_int).
    """
    contacts = load_contact_matrix(country)
    population = get_population_total(country)
    model = build_country_seir_model(
        contacts.values,
        population_total=population,
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
    country = "China"
    run_demo(country=country)
