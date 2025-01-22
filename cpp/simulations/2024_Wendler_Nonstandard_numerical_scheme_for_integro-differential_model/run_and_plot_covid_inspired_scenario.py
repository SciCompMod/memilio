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
import math
import h5py
import pandas as pd
import numpy as np

from get_lognormal_parameters import get_lognormal_parameters
from plot_covid_inspired_scenario import load_data, plot_covid_inspired_scenario


def run_covid_inspired_scenario(result_dir, data_dir, start_date, simulation_time, timestep, scale_contacts):
    """ Run the Covid inspired scenarias defined in ide_covid_inspired_scenario.cpp with the here defined paths and
    parameters.

    @param[in] result_dir Directory where simiulation results are stored.
    @param[in] data_dir Directory where folders with data on contacts and reported data are located. 
    @param[in] start_date Start date of the simulations.
    @param[in] simulation_time Duration of the simulation.
    @param[in] timestep Time step used for the simulations. 
    @param[in] scale_contacts Scaling factor that is applied to the used contact matrix.
    """
    year = start_date.split("-")[0]
    month = start_date.split("-")[1]
    day = start_date.split("-")[2]

    subprocess.call([f"./../../../build/bin/ide_covid_inspired_scenario", data_dir, result_dir,
                     f"{year}", f"{month}", f"{day}", f"{simulation_time}", f"{timestep}", f"{scale_contacts}"])


def get_scale_contacts(files, reported_data_dir, start_date, simulation_time):
    """ Gets the scaling factor for the contacts so that the simulation results obtained with the IDE model match the
    reported number of daily new transmissions (as calculated in load_data()).

    @param[in] file Expects file with IDE simulation results for flows. 
    @param[in] data_dir Directory where RKI data is stored. 
    @param[in] start_date Start date of interest. 
    @param[in] simulation_time Number of days to be simulated.
    @returns scale_contacts Scaling factor for contacts to scale IDE results to RKI data.
    """

    datafile = os.path.join(reported_data_dir, "cases_all_germany.json")
    data_rki = load_data(datafile, start_date, simulation_time)

    # Load IDE data.
    for file in range(len(files)):
        h5file = h5py.File(str(files[file]) + '.h5', 'r')

        if (len(list(h5file.keys())) > 1):
            raise gd.DataError("File should contain one dataset.")
        if (len(list(h5file[list(h5file.keys())[0]].keys())) > 3):
            raise gd.DataError("Expected only one group.")

        data = h5file[list(h5file.keys())[0]]

        # As there should be only one Group, total is the simulation result
        total = data['Total'][:, :]

        dates = data['Time'][:]
        timestep = np.diff(dates)[0]

    # Date index referring to first time step of IDE simulation.
    date_idx = math.floor(-simulation_time / timestep)

    # Get daily new transmissions from IDE simulation results.
    new_transmissions_ide = total[date_idx, 0] / timestep
    print("")
    print(f"IDE new infections at time {dates[date_idx]}: ", new_transmissions_ide)

    # Get daily new transmissions from RKI data at first timestep.
    new_transmissions_rki = data_rki[data_rki["Date"] == start_date]["NewInfectionsDay"].values[0] + timestep * (data_rki[data_rki["Date"] == start_date + pd.DateOffset(
        days=1)]["NewInfectionsDay"].values[0] - data_rki[data_rki["Date"] == start_date]["NewInfectionsDay"].values[0])
    print(f"Expected new infections (based on RKI reports) at time {timestep}: {new_transmissions_rki}")
    print("")

    # Compute scaling factor for contacts.
    scale_contacts = new_transmissions_rki/new_transmissions_ide

    return scale_contacts

def run_scenario_with_adjusted_contact_scaling(result_dir, data_dir, plot_dir, start_date, simulation_time, timestep):
    """ Here, we run the scenario with a contact scaling that is adjusted to the data as reported by RKI.
    For this, we first run the scenario with a contact sclaing of 1. We then compare the IDE simulation results with
    the reported data at the first time step. Based on this, we compute the contact scaling such that the IDE 
    simulation results match the reported data. Finally, we run the simulation again with the adjusted contact scaling 
    and plot the results. 

    @param[in] result_dir Directory where simiulation results are stored.
    @param[in] data_dir Directory where folders with data on contacts and reported data are located. 
    @param[in] plot_dir Directory where plots will be stored.
    @param[in] start_date Start date of the simulations.
    @param[in] simulation_time Duration of the simulation.
    @param[in] timestep Time step used for the simulations. 
    """
    # First run the simulation with a contact scaling of 1.
    scale_contacts = 1.
    run_covid_inspired_scenario(result_dir, data_dir, start_date, simulation_time,
                      timestep, scale_contacts)
    
    # Then determine contact scaling such that IDE results and RKI new transmissions match at first timestep.
    reported_data_dir = data_dir + "/pydata/Germany/"

    scale_contacts = get_scale_contacts([os.path.join(
        result_dir, f"ide_{start_date}_{simulation_time}_{timestep}_flows")], reported_data_dir, pd.Timestamp(start_date), simulation_time)
    
    # Run simulation with new contact scaling.
    run_covid_inspired_scenario(result_dir, data_dir, start_date, simulation_time,
                      timestep, scale_contacts)

    plot_covid_inspired_scenario(result_dir, reported_data_dir, plot_dir, start_date,
                                 simulation_time, timestep)

def main():
    # Scenario starting on October 1st, 2020. 
    start_date = '2020-10-1'
    simulation_time = 45
    timestep = "0.0100"

    # Paths are valid if file is executed e.g. in memilio/cpp/simulations/IDE_paper.
    # Path where simulation results (generated with ide_covid_inspired_scenario.cpp) are stored.
    result_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/simulation_results/covid_inspired_scenario/")

    # Path where data for contacts and reported data is stored.
    data_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/")

    # Path where plots will be stored.
    plot_dir = os.path.join(os.path.dirname(
        __file__), "../../..", "data/plots/covid_inspired_scenario/")

    run_scenario_with_adjusted_contact_scaling(result_dir, data_dir, plot_dir, start_date, simulation_time,
                 timestep)


if __name__ == "__main__":
    main()
