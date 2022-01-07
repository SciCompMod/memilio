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

import getDataIntoPandasDataFrame as gd
import defaultDict as dd
import geoModificationGermany as geoger


def transformTwitterData(file_format=dd.defaultDict['file_format'],
                         out_folder=dd.defaultDict['out_folder'].split(
    '/pydata')[0]
):
    """! Merges rows and columns of Eisenach to Wartburgkreis. To be 
    optimized for production code.

    @param file_format File format which is used for writing the data. 
        Default defined in defaultDict.
    @param out_folder Path to folder where data is written in folder 
        out_folder/Germany.
    """

    directory = out_folder
    directory = os.path.join(directory, 'mobility/')

    twitter = pd.read_csv(
        directory + 'twitter_scaled_1252.txt', sep=' ', header=None)

    if len(twitter) == 401:
        twitter.to_csv(
            directory + 'twitter_scaled_1252_dim401.txt', sep=' ', header=None, index=False)
        # merge eisenach
        ids400 = geoger.get_county_ids()
        ids401 = geoger.get_county_ids(merge_eisenach=False)
        indices = [ids401.index(i) for i in ids400]
        twitter_new = twitter.iloc[indices, indices]
        idx_eisenach_old = ids401.index(geoger.CountyMerging[16063][1])
        idx_wartburg_new = ids400.index(geoger.CountyMerging[16063][0])
        twitter_new.iloc[idx_wartburg_new,
                         :] += twitter.iloc[idx_eisenach_old, indices].values
        twitter_new.iloc[:, idx_wartburg_new] += twitter.iloc[indices,
                                                              idx_eisenach_old].values
        # TODO: this file is generally to be optimized
        if abs(twitter_new.iloc[0:382, 0:382]-twitter.iloc[0:382, 0:382]).max().max() > 1e-10:
            print('Error...')
        twitter_new.to_csv(
            directory + 'twitter_scaled_1252.txt', sep=' ', header=None, index=False)


def main():
    """! Main program entry."""

    transformTwitterData()


if __name__ == "__main__":

    main()
