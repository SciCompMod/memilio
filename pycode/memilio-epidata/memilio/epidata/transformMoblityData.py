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

    if(len(mobility_matrix) == geoger.get_county_ids()):

        x=15

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
        twitter_new = mobility_matrix.iloc[indices, indices]
        idx_eisenach_old = ids401.index(geoger.CountyMerging[16063][1])
        idx_wartburg_new = ids400.index(geoger.CountyMerging[16063][0])
        twitter_new.iloc[idx_wartburg_new,
                         :] += mobility_matrix.iloc[idx_eisenach_old, indices].values
        twitter_new.iloc[:, idx_wartburg_new] += mobility_matrix.iloc[indices,
                                                              idx_eisenach_old].values
        # TODO: this file is generally to be optimized
        if abs(twitter_new.iloc[0:382, 0:382]-mobility_matrix.iloc[0:382, 0:382]).max().max() > 1e-10:
            print('Error...')
        twitter_new.to_csv(
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
