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

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd
from memilio.epidata import geoModificationGermany as geoger

def createFederalStatesMobility(mobility_file, file_format=dd.defaultDict['file_format'],
                         out_folder=dd.defaultDict['out_folder'].split('/pydata')[0]):
    """! Creates mobility matrices for German federal states based on
    county mobility.

    @param file_format File format which is used for writing the data. 
        Default defined in defaultDict.
    @param out_folder Path to folder where data is read and written.
        Default is data/mobility folder in root directory of MEmilio.
    """
    directory = out_folder
    directory = os.path.join(directory, 'mobility/')

    directory = directory.replace('memilio','memilio-fed-state')

    mobility_matrix = pd.read_csv(
        directory + mobility_file + '.txt', sep=' ', header=None)   

    if(len(mobility_matrix) == len(geoger.get_county_ids())):
        # get county and state IDs
        countyIDs = geoger.get_county_ids()
        stateIDs = geoger.get_state_ids()
        # get state ID to county ID map
        stateID_to_countyID = geoger.get_stateid_to_countyids_map()
        # initialize federal states mobility matrix
        mobility_matrix_states = np.zeros((len(stateID_to_countyID), len(stateID_to_countyID)))

        # iterate over state_to_county map and replace IDs by numbering 0, ..., n
        state_indices = []
        county_indices = []
        for state, counties in stateID_to_countyID.items():
            state_indices.append(stateIDs.index(state))
            county_indices.append(np.array([countyIDs.index(county) for county in counties]))

        state_to_county = dict(zip(state_indices, county_indices))
        # iterate over all states
        for state, counties in state_to_county.items():
            # iterate over all neighbors
            for state_neighb, counties_neighb in state_to_county.items():
                
                if state != state_neighb:

                    mobility_matrix_states[state, state_neighb] = mobility_matrix.iloc[counties, counties_neighb].sum(axis=0).sum()
                    mobility_matrix_states[state_neighb, state] = mobility_matrix.iloc[counties_neighb, counties].sum(axis=1).sum()



def updateMobility2022(mobility_file, file_format=dd.defaultDict['file_format'],
                         out_folder=dd.defaultDict['out_folder'].split(
    '/pydata')[0]
):
    """! Merges rows and columns of Eisenach to Wartburgkreis which has
    become one single county by July 2021. To be optimized for production code.

    @param mobility_file Mobility matrix file which has to be updated.
    @param file_format File format which is used for writing the data. 
        Default defined in defaultDict.
    @param out_folder Path to folder where data is read and written.
        Default is data/mobility folder in root directory of MEmilio.
    """

    directory = out_folder
    directory = os.path.join(directory, 'mobility/')

    directory = directory.replace('memilio','memilio-fed-state')

    mobility_matrix = pd.read_csv(
        directory + mobility_file + '.txt', sep=' ', header=None)

    if len(mobility_matrix) == 401:
        mobility_matrix.to_csv(
            directory + mobility_file + '_dim401.txt', sep=' ', header=None, index=False)
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
        # TODO: this file is generally to be optimized
        if abs(mobility_matrix_new.iloc[0:382, 0:382]-mobility_matrix.iloc[0:382, 0:382]).max().max() > 1e-10:
            print('Error...')
        mobility_matrix_new.to_csv(
            directory + mobility_file + '.txt', sep=' ', header=None, index=False)


def main():
    """! Main program entry."""

    # Merge Eisenach and Wartbugkreis in Input Data if need be
    updateMobility2022(mobility_file='twitter_scaled_1252')
    updateMobility2022(mobility_file='commuter_migration_scaled')
    # create federal states mobility matrix (not used in simulation for now)
    createFederalStatesMobility(mobility_file='twitter_scaled_1252')
    createFederalStatesMobility(mobility_file='commuter_migration_scaled')


if __name__ == "__main__":

    main()
