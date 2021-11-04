# TODO: to be either removed or substantially changed before merging...
import pandas as pd

mobi = 'mobility'

df1 = pd.read_csv('/home/kueh_mj/memilio/data/' + mobi + '/twitter_scaled_1252.txt', sep=' ', header=None)
print("Size twitter data " + str(len(df1)) + " x " + str(len(df1.columns)))
print(df1.iloc[0,1])
print('partially empty lines: ' + str(len(df1[df1.isnull().any(axis=1)])))

df2 = pd.read_csv('/home/kueh_mj/memilio/data/' + mobi + '/commuter_migration_scaled_2020.txt', sep=' ', header=None)
print("Size commuter data " + str(len(df2)) + " x " + str(len(df2.columns)))
print(df2.iloc[0,1])
print('partially empty lines: ' + str(len(df2[df2.isnull().any(axis=1)])))

df3 = pd.read_json('/home/kueh_mj/memilio/data/pydata/Germany/county_divi_ma7.json')
print("Size DIVI infection data " + str(len(df3)) + " x " + str(len(df3.columns)) + ", division of length by 400: " + str(len(df3)/400))
print(df3.columns)
print('partially empty lines: ' + str(len(df3[df3.isnull().any(axis=1)])))

df4 = pd.read_json('/home/kueh_mj/memilio/data/pydata/Germany/all_county_age_ma7_rki.json')
print("Size RKI infection data " + str(len(df4)) + " x " + str(len(df4.columns)) + ", division of length by 400: " + str(len(df4)/400))
print(df4.columns)
print('partially empty lines: ' + str(len(df4[df4.isnull().any(axis=1)])))

df5 = pd.read_json('/home/kueh_mj/memilio/data/pydata/Germany/county_current_population.json')
print("Size population data " + str(len(df5)) + " x " + str(len(df5.columns)) + ", division of length by 400: " + str(len(df5)/400))
print(df5.columns)
print('partially empty lines: ' + str(len(df5[df5.isnull().any(axis=1)])))

df6 = pd.read_json('/home/kueh_mj/memilio/data/pydata/Germany/all_county_ageinf_vacc_ma7.json')
print("Size RKI vaccination data " + str(len(df6)) + " x " + str(len(df6.columns)) + ", division of length by 400: " + str(len(df6)/400))
print(df6.columns)
print('partially empty lines: ' + str(len(df6[df6.isnull().any(axis=1)])))