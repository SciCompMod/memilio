#############################################################################
# Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn
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
import os
import pandas as pd
import numpy as np

from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger
from memilio.epidata import getDataIntoPandasDataFrame as gd


def getMobilityFromFile(directory, mobility_file):
    """! Gets a mobility matrix that is written in a plain txt file
    under the given directory into a pandas data frame.

    @param directory Path to folder where data is read.
    @param mobility_file Mobility matrix file which has to be updated.
    @return Mobility matrix data frame.
    """
    mobility_matrix = gd.loadCsv(
        '', os.path.join(directory, mobility_file), extension='.txt',
        param_dict={'sep': ' ', 'header': None})

    return mobility_matrix


def createFederalStatesMobility(directory, mobility_file):
    """! Creates mobility matrices for German federal states based on
    county mobility. If mobility matrix dimensions are different from the
    number of German counties, nothing is done.

    @param directory Path to folder where data is read and written.
    @param mobility_file Mobility matrix file which has to be updated.
    @return State-aggregated mobility matrix if input matrix is sized correctly, zero otherwise.
    """
    mobility_matrix = getMobilityFromFile(directory, mobility_file)

    if (len(mobility_matrix.index) == len(geoger.get_county_ids())) and (len(mobility_matrix.columns) == len(geoger.get_county_ids())):
        # get county and state IDs
        countyIDs = geoger.get_county_ids()
        stateIDs = geoger.get_state_ids()
        # get state ID to county ID map
        stateID_to_countyID = geoger.get_stateid_to_countyids_map()
        # initialize federal states mobility matrix
        mobility_matrix_states = np.zeros(
            (len(stateID_to_countyID), len(stateID_to_countyID)))

        # iterate over state_to_county map and replace IDs by numbering 0, ..., n
        state_indices = []
        county_indices = []
        for state, counties in stateID_to_countyID.items():
            state_indices.append(stateIDs.index(state))
            county_indices.append(
                np.array([countyIDs.index(county) for county in counties]))

        state_to_county = dict(zip(state_indices, county_indices))
        # iterate over all states
        for state, counties in state_to_county.items():
            # iterate over all neighbors
            for state_neighb, counties_neighb in state_to_county.items():

                if state != state_neighb:

                    mobility_matrix_states[state, state_neighb] = mobility_matrix.iloc[counties, counties_neighb].sum(
                        axis=0).sum()
                    mobility_matrix_states[state_neighb, state] = mobility_matrix.iloc[counties_neighb, counties].sum(
                        axis=1).sum()

        mobility_matrix_states = pd.DataFrame(mobility_matrix_states)
        gd.write_dataframe(
            mobility_matrix_states, directory, mobility_file + '_states', 'txt',
            param_dict={'sep': ' ', 'header': None, 'index': False})
        return mobility_matrix_states

    else:
        return pd.DataFrame()


def updateMobility2022(directory, mobility_file):
    """! Merges rows and columns of Eisenach to Wartburgkreis which has
    become one single county by July 2021. If mobility matrix dimension is different
    from 401x401, nothing is done.

    @param directory Path to folder where data is read and written.
    @param mobility_file Mobility matrix file which has to be updated.
    @return Reduced mobility matrix or input mobility matrix if matrix was already reduced.
    """
    mobility_matrix = getMobilityFromFile(directory, mobility_file)

    if (len(mobility_matrix.index) == 401) and (len(mobility_matrix.columns) == 401):
        gd.write_dataframe(mobility_matrix, directory, mobility_file + '_dim401', 'txt',
            param_dict={'sep': ' ', 'header': None, 'index': False})
        # merge eisenach
        ids400 = geoger.get_county_ids()
        ids401 = geoger.get_county_ids(merge_eisenach=False)
        indices = [ids401.index(i) for i in ids400]
        mobility_matrix_new = mobility_matrix.iloc[indices, indices]
        idx_eisenach_old = ids401.index(geoger.CountyMerging[16063][1])
        idx_wartburg_new = ids400.index(geoger.CountyMerging[16063][0])
        mobility_matrix_new.iloc[idx_wartburg_new,
                                 :] += mobility_matrix.iloc[idx_eisenach_old, indices].values
        mobility_matrix_new.iloc[:, idx_wartburg_new] += mobility_matrix.iloc[indices,
                                                                              idx_eisenach_old].values

        gd.write_dataframe(mobility_matrix_new, directory, mobility_file, 'txt',
            param_dict={'sep': ' ', 'header': None, 'index': False})

        return mobility_matrix_new

    else:
        return mobility_matrix


def main():
    """! Main program entry."""

    path = os.path.join(os.getcwd(), 'data')
    directory = os.path.join(path, 'mobility/')

    # Merge Eisenach and Wartbugkreis in Input Data if need be
    updateMobility2022(directory, mobility_file='twitter_scaled_1252')
    updateMobility2022(directory, mobility_file='commuter_migration_scaled')
    # create federal states mobility matrix (not used in simulation for now)
    createFederalStatesMobility(directory, mobility_file='twitter_scaled_1252')
    createFederalStatesMobility(
        directory, mobility_file='commuter_migration_scaled')


if __name__ == "__main__":

    main()
