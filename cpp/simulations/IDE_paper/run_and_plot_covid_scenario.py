#############################################################################
# Copyright (C) 2020-2024 German Aerospace Center (DLR-SC)
#
# Authors: Anna Wendler
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
import subprocess
import os
import pandas as pd

from get_lognormal_parameters import get_lognormal_parameters
from plot_real_scenario import plot_daily_new_transmissions, plot_infectedsymptoms_deaths, plot_icu, get_scale_contacts, get_scale_confirmed_cases


def run_real_scenario(result_dir, data_dir, start_date, simulation_time, timestep, scale_contacts):

    year = start_date.split("-")[0]
    month = start_date.split("-")[1]
    day = start_date.split("-")[2]

    subprocess.call([f"./../../../build/bin/ide_covid_scenario", data_dir, result_dir,
                     f"{year}", f"{month}", f"{day}", f"{simulation_time}", f"{timestep}", f"{scale_contacts}"])


def contact_scaling(result_dir, reported_data_dir, start_date, simulation_time, timestep):
    scale_contacts = get_scale_contacts([os.path.join(
        result_dir, f"ide_{start_date}_{simulation_time}_{timestep}_flows")], reported_data_dir, pd.Timestamp(start_date), simulation_time)

    return scale_contacts


def confirmed_cases_scaling(result_dir, reported_data_dir, start_date, simulation_time, timestep):
    scale_confirmed_cases = get_scale_confirmed_cases([os.path.join(
        result_dir, f"ide_{start_date}_{simulation_time}_{timestep}_compartments")], reported_data_dir,  pd.Timestamp(start_date))

    return scale_confirmed_cases


def plot_covid_inspired_scenario(result_dir, data_dir, plot_dir, start_date, simulation_time, timestep):
    plot_daily_new_transmissions([os.path.join(result_dir, f"ode_{start_date}_{simulation_time}_{timestep}_flows"),
                                  os.path.join(result_dir, f"ide_{start_date}_{simulation_time}_{timestep}_flows")], data_dir,
                                 pd.Timestamp(start_date), simulation_time,
                                 fileending=f"{start_date}_{simulation_time}_{timestep}", save_dir=plot_dir)

    plot_infectedsymptoms_deaths([os.path.join(result_dir, f"ode_{start_date}_{simulation_time}_{timestep}_compartments"),
                                  os.path.join(result_dir, f"ide_{start_date}_{simulation_time}_{timestep}_compartments")],
                                 data_dir,
                                 pd.Timestamp(
                                     start_date), simulation_time,
                                 fileending=f"{start_date}_{simulation_time}_{timestep}",  save_dir=plot_dir)

    plot_icu([os.path.join(result_dir, f"ode_{start_date}_{simulation_time}_{timestep}_compartments"),
              os.path.join(result_dir, f"ide_{start_date}_{simulation_time}_{timestep}_compartments")],
             data_dir,
             pd.Timestamp(start_date), simulation_time,  fileending=f"{start_date}_{simulation_time}_{timestep}", save_dir=plot_dir)


def run_scenario(result_dir, data_dir, plot_dir, start_date, simulation_time, timestep):
    reported_data_dir = data_dir + "/pydata/Germany/"
    # First run the simulation with a contact scaling of 1.
    scale_contacts = 1.
    scale_confirmed_cases = 1.
    run_real_scenario(result_dir, data_dir, start_date, simulation_time,
                      timestep, scale_contacts)
    # Then determine contact scaling such that IDE results and RKI new infections match at first timestep.
    scale_contacts = contact_scaling(
        result_dir, reported_data_dir, start_date, simulation_time, timestep)
    scale_confirmed_cases = confirmed_cases_scaling(
        result_dir, reported_data_dir, start_date, simulation_time, timestep)
    print(scale_confirmed_cases)
    # Run simulation with new contact scaling.
    run_real_scenario(result_dir, data_dir, start_date, simulation_time,
                      timestep, scale_contacts)

    plot_covid_inspired_scenario(result_dir, reported_data_dir, plot_dir, start_date,
                                 simulation_time, timestep)


def october_scenario(timestep):
    start_date = '2020-10-1'
    simulation_time = 45

    # data_dir = "./data"
    # result_dir = f"./results/real/"
    # plot_dir = f"./plots/covid_scenario/{start_date}/"

    # Paths are valid if file is executed e.g. in memilio/cpp/simulations/IDE_paper.
    # Path where simulation results (generated with ide_real_scenario.cpp) are stored.
    result_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/simulation_results/covid_inspired_scenario/")

    # Path where data for contacts and reported data is stored.
    data_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/")

    # Path where plots will be stored.
    plot_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/plots/covid_inspired_scenario/")

    run_scenario(result_dir, data_dir, plot_dir, start_date, simulation_time,
                 timestep)


def main():
    # Folder in plots/real/save_folder where plots will be stored.

    timestep = "0.0100"
    october_scenario(timestep)


if __name__ == "__main__":
    main()
