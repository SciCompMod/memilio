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
from plot_real_scenario import plot_new_infections, plot_infectedsymptoms_deaths, plot_icu, get_scale_contacts, get_scale_confirmed_cases


def run_real_scenario(data_dir, save_dir, start_date, simulation_time, timestep,  scale_contacts):

    year = start_date.split("-")[0]
    month = start_date.split("-")[1]
    day = start_date.split("-")[2]

    subprocess.call([f"./build/bin/ide_covid_scenario", data_dir, save_dir,
                     f"{year}", f"{month}", f"{day}", f"{simulation_time}", f"{timestep}", f"{scale_contacts}"])


def contact_scaling(save_dir, start_date, simulation_time, timestep):
    scale_contacts = get_scale_contacts([os.path.join(
        save_dir, f"ide_{start_date}_{simulation_time}_{timestep}_flows")], pd.Timestamp(start_date), simulation_time)

    return scale_contacts


def confirmed_cases_scaling(save_dir, start_date, simulation_time, timestep):
    scale_confirmed_cases = get_scale_confirmed_cases([os.path.join(
        save_dir, f"ide_{start_date}_{simulation_time}_{timestep}_compartments")], pd.Timestamp(start_date))

    return scale_confirmed_cases


def plot_real_scenario(save_dir, plot_dir, start_date, simulation_time, timestep):
    plot_new_infections([os.path.join(save_dir, f"ode_{start_date}_{simulation_time}_{timestep}_flows"),
                        os.path.join(save_dir, f"ide_{start_date}_{simulation_time}_{timestep}_flows")],
                        pd.Timestamp(start_date), simulation_time,
                        fileending=f"{start_date}_{simulation_time}_{timestep}", save=True, save_dir=plot_dir)

    plot_infectedsymptoms_deaths([os.path.join(save_dir, f"ode_{start_date}_{simulation_time}_{timestep}_compartments"),
                                  os.path.join(save_dir, f"ide_{start_date}_{simulation_time}_{timestep}_compartments")],
                                 pd.Timestamp(
                                     start_date), simulation_time,
                                 fileending=f"{start_date}_{simulation_time}_{timestep}", save=True, save_dir=plot_dir)

    plot_icu([os.path.join(save_dir, f"ode_{start_date}_{simulation_time}_{timestep}_compartments"),
              os.path.join(save_dir, f"ide_{start_date}_{simulation_time}_{timestep}_compartments")],
             pd.Timestamp(start_date), simulation_time,  fileending=f"{start_date}_{simulation_time}_{timestep}", save=True, save_dir=plot_dir)


def run_scenario(data_dir, save_dir, plot_dir, start_date, simulation_time, timestep):
    # First run the simulation with a contact scaling of 1.
    scale_contacts = 1.
    scale_confirmed_cases = 1.
    run_real_scenario(data_dir, save_dir, start_date, simulation_time,
                      timestep, scale_contacts)
    # Then determine contact scaling such that IDE results and RKI new infections match at first timestep.
    scale_contacts = contact_scaling(
        save_dir, start_date, simulation_time, timestep)
    # scale_confirmed_cases=confirmed_cases_scaling(save_dir, start_date, simulation_time, timestep)
    print(scale_confirmed_cases)
    # Run simulation with new contact scaling.
    run_real_scenario(data_dir, save_dir, start_date, simulation_time,
                      timestep, scale_contacts)
    plot_real_scenario(save_dir, plot_dir, start_date,
                       simulation_time, timestep)


def october_scenario(timestep):
    start_date = '2020-10-1'
    simulation_time = 45

    data_dir = "./data"
    save_dir = f"./results/real/"
    plot_dir = f"./plots/covid_scenario/{start_date}/"

    run_scenario(data_dir, save_dir, plot_dir, start_date, simulation_time,
                 timestep)


def main():
    # Folder in plots/real/save_folder where plots will be stored.

    timestep = "0.0100"
    october_scenario(timestep)


if __name__ == "__main__":

    main()
