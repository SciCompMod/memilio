import os

import numpy as np
import pandas as pd

from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd


def main():
    """! Main program entry."""

    arg_dict = gd.cli("commuter_official")

    directory = os.path.join(
        arg_dict['out_folder'].split('/pydata/')[0], 'Germany/')
    directory_mobility = os.path.join(directory, 'mobility/')
    directory_population = os.path.join(directory, 'pydata/')
    mobility_file = 'commuter_mobility_2022'
    population_file = 'county_current_population'

    mobility_matrix = pd.read_csv(
        os.path.join(directory_mobility + mobility_file + '.txt'),
        sep=' ', header=None)

    # get county and state IDs
    countyIDs = geoger.get_county_ids()
    stateIDs = geoger.get_state_ids()
    # get state ID to county ID map
    stateID_to_countyID = geoger.get_stateid_to_countyids_map()

    # iterate over state_to_county map and replace IDs by numbering 0, ..., n
    state_indices = []
    county_indices = []
    for state, counties in stateID_to_countyID.items():
        state_indices.append(stateIDs.index(state))
        county_indices.append(
            np.array([countyIDs.index(county) for county in counties]))

    mobility_matrix_nrw = mobility_matrix.loc[county_indices[4],
                                              county_indices[4]]

    gd.write_dataframe(
        mobility_matrix_nrw, directory_mobility, mobility_file + '_nrw', 'txt',
        param_dict={'sep': ' ', 'header': None, 'index': False})

    population = pd.read_json(os.path.join(
        directory_population + population_file + '.json'))
    population_nrw = population.loc[county_indices[4]]
    gd.write_dataframe(population_nrw, directory_population,
                       population_file + '_nrw', 'json')


if __name__ == "__main__":

    main()
