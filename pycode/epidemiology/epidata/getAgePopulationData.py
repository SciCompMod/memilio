# This script reads out the age specific population data from the 2011 census. In order to bring the data to a more
# recent level, the ratio of the total populations in 2019 and 2011 are computed and applied to all age groups in 2011.
# Both the original and the adjusted data are stored in json files.

import os
import pandas as pd
import numpy as np



path_counties = 'http://hpcagainstcorona.sc.bs.dlr.de/data/migration/'
path_reg_key = 'https://www.zensus2011.de/SharedDocs/Downloads/DE/Pressemitteilung/DemografischeGrunddaten/1A_EinwohnerzahlGeschlecht.xls?__blob=publicationFile&v=5'
path_zensus = 'https://opendata.arcgis.com/datasets/abad92e8eead46a4b0d252ee9438eb53_1.csv'

#read tables
counties = pd.read_excel(os.path.join(path_counties,'kreise_deu.xlsx'),sheet_name=1, header=3)
reg_key = pd.read_excel(path_reg_key, sheet_name='Tabelle_1A', header=12)
zensus = pd.read_csv(path_zensus)


#find region keys for census population data
key = np.zeros((len(zensus)))
for i in range(len(key)):
    for j in range(len(reg_key)):
        if zensus.Name.values[i] == reg_key['NAME'].values.astype(str)[j]:
            if zensus.EWZ.values[i] == round(reg_key['Zensus_EWZ'].values[j]*1000):
                key[i] = reg_key['AGS'].values[j]

unique, inds, count = np.unique(key, return_index=True, return_counts=True)

male = ['M_Unter_3', 'M_3_bis_5', 'M_6_bis_14', 'M_15_bis_17', 'M_18_bis_24',
        'M_25_bis_29', 'M_30_bis_39', 'M_40_bis_49', 'M_50_bis_64',
        'M_65_bis_74', 'M_75_und_aelter']
female = ['W_Unter_3', 'W_3_bis_5', 'W_6_bis_14', 'W_15_bis_17', 'W_18_bis_24',
        'W_25_bis_29', 'W_30_bis_39', 'W_40_bis_49', 'W_50_bis_64',
        'W_65_bis_74', 'W_75_und_aelter']
columns = ['Key', 'Total', '<3 years', '3-5 years', '6-14 years', '15-17 years', '18-24 years',
           '25-29 years', '30-39 years', '40-49 years', '50-64 years',
           '65-74 years', '>74 years']

# add male and female population data
data = np.zeros((len(inds), len(male)+2))
data[:,0] = key[inds].astype(int)
data[:,1] = zensus['EWZ'].values[inds].astype(int)
for i in range(len(male)):
    data[:, i+2] = zensus[male[i]].values[inds].astype(int) + zensus[female[i]].values[inds].astype(int)

# compute ratio of current and 2011 population data
ratio = np.ones(len(inds))
for i in range(len(inds)):
    for j in range(len(counties)-11):
        if not counties['Schlüssel-nummer'].isnull().values[j]:
            if data[i,0] == int(counties['Schlüssel-nummer'].values[j]):
                ratio[i] = counties['Bevölkerung2)'].values[j]/data[i, 1]

# adjust population data for all ages to current level
data_current = np.zeros(data.shape)
data_current[:, 0] = data[:, 0]
for i in range(len(data[0, :]) - 1):
    data_current[:, i + 1] = np.multiply(data[:, i + 1], ratio)


#create dataframe
df = pd.DataFrame(data.astype(int), columns=columns)
df_current = pd.DataFrame(np.round(data_current).astype(int), columns=columns)

#save dataframe
df.to_json('county_population.json')
df_current.to_json('county_current_population.json')