import os
import csv
import pandas as pd
import numpy as np

from memilio.epidata import getDataIntoPandasDataFrame as gd
from memilio.epidata import defaultDict as dd

directory = '/home/wend_aa/memilio/data/pydata/Germany'

#############################################################################################################
# read old data for subcategories

df_npis_old_data = pd.read_csv(
    os.path.join(directory, 'kr_massnahmen_unterkategorien.csv'),
    sep=',')  # , nrows=numberofcities*1248

df_npis_old_data.rename(dd.GerEng, axis=1, inplace=True)

#############################################################################################################
# read new data for subcategories

codelist = ['m01a', 'm01b', 'm02a', 'm02b', 'm03', 'm04', 'm05', 'm06', 'm07', 'm08', 'm09',
            'm10', 'm11', 'm12', 'm13', 'm14', 'm15', 'm16', 'm17', 'm18', 'm19', 'm20', 'm21']
counter_codes = 0
for code in codelist:
    print(code)
    df_npis_per_code = pd.read_csv(
        os.path.join(directory,
                     f'kr_massn_unterkat_{code}.csv'),
        sep=',')

    # set some parameters for dataframe
    if counter_codes == 0:
        counties = np.sort(df_npis_per_code.ags5.unique())
        num_counties = len(df_npis_per_code.ags5.unique())

        # extract dates from data
        dates = df_npis_per_code.iloc[:int(
            df_npis_per_code.shape[0]/num_counties), 5]
        # rename dates so that they match dates from other npi dataframe
        dates_new = ['d' + date.replace('-', '') for date in dates]

        df_local = [pd.DataFrame() for i in range(num_counties)]

    # set df for all counties
    for i in range(0, num_counties):
        print(i)
        if counter_codes == 0:
            df_local[i] = pd.DataFrame(columns=list(
                df_npis_per_code.columns[0:5]) + ['code'] + dates_new)

        dummy_to_append = pd.DataFrame(columns=[
                                       'code'] + dates_new, data=df_npis_per_code[df_npis_per_code.ags5 == counties[i]].iloc[:, 6:].T.reset_index().values.copy())

        df_local[i] = pd.concat([df_local[i], dummy_to_append])

        if df_npis_per_code.iloc[i*len(dates):(i+1)*len(dates), 3].nunique() > 1:
            raise gd.DataError('Dates are not sorted as expected.')

        # Set first five columns so that they match old format of data frame (from kr_massnahmen_unterkategorien.csv)
        if counter_codes == len(codelist)-1:
            df_local[i].iloc[:, 0:5] = df_npis_per_code.iloc[i *
                                                             len(dates), 0:5].values

    counter_codes += 1

df_npis_new_data = pd.concat([df_local[i] for i in range(num_counties)])
df_npis_new_data.rename(dd.GerEng, axis=1, inplace=True)
df_npis_new_data['NPI_code'] = df_npis_new_data['NPI_code'].str.replace(
    'code_m', 'M')


#############################################################################################################
# compare dataframes

# check if all rows for code M22, M23 and M24 in df_npis_old_data are empty
codesnotused = ((df_npis_old_data[df_npis_old_data["NPI_code"].str.contains(
    "M22|M23|M24")].iloc[:, 6:] == -99).all() == True).all()
if codesnotused == True:
    print("Codes M22, M23 and M24 are not used in old data (as expected).")
else:
    print("Something wrong with data.")

# remove rows for codes M22, M23 and M24 from df_npis_old_data
df_npis_old_data = df_npis_old_data[~df_npis_old_data["NPI_code"].str.contains(
    "M22|M23|M24")].copy()

# check how many days are covered in each dataframe and adjust accordingly so that both dataframes have same size
# we already know that df_npis_new_data has more columns than df_npis_old_data
df_npis_new_data = df_npis_new_data.iloc[:, :len(df_npis_old_data.columns)]

# assert if frames are equal (except index and column '_id')

if (pd.testing.assert_frame_equal(df_npis_old_data.iloc[:, 1:].reset_index(drop=True), df_npis_new_data.iloc[:, 1:].reset_index(drop=True), check_dtype=False) == None):
    print('Data frames are equal.')
else:
    print('Data frames are not equal.')
