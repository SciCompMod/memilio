# Script to extract one data row for every countx each to generate test
# data for test_epidata_get_Case_Data_read.py which is not to large

import pandas as pd
from memilio.epidata import defaultDict as dd
from memilio.epidata import getDataIntoPandasDataFrame as gd


# Put any CaseDataFull.json file or Update file has below here
file_name = "UpdateDataCases2021-07-08.json"

# read in
df = pd.read_json(file_name)

# empty dataframe to get every county jut once
df_neu = pd.DataFrame()

for k in dd.County.keys():

    df_k = df[df["IdLandkreis"] == k].reset_index()

    print(k, len(df_k))

    if len(df_k) > 0:

        df_append = df_k.iloc[[0]]
        df_neu = df_neu.append(df_append)

print("Number of counties in data:", len(df_neu))


gd.write_dataframe(df_neu, ".", "test_epidata_get_Case_Data_data_read", "json")
