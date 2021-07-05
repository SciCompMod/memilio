import pandas
import requests
import io
import numpy as np
import os

from epidemiology.epidata  import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from epidemiology.epidata import getPopulationData

# Downloads vaccine data from RKI
def download_vaccine_data():
    # get rki vaccine data
    url = 'https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Daten/Impfquotenmonitoring.xlsx?__blob=publicationFile'
    header = {'User-Agent': 'Mozilla/5.0'}
    r = requests.get(url, headers=header)
    with io.BytesIO(r.content) as fh:
        df = pandas.io.excel.ExcelFile(fh, engine='openpyxl')
        sheet_names = df.sheet_names
        vaccine = pandas.read_excel(df, sheet_name=sheet_names[1], header=[1])[2:18].drop(['Unnamed: 1', 'Unnamed: 10'], axis=1)
        name = sheet_names[2].replace('.', '_', 2).replace('.', '', 1)[-8:]
        name = 'vaccine_data_' + name

    return vaccine, name

# gets rki vaccine monitoring data for all states and extrapolates the values for counties according to their population
# Missing ratio values for the two different age groups are also estimated
def get_vaccine_data(read_data=dd.defaultDict['read_data'],
                     file_format=dd.defaultDict['file_format'],
                     out_folder=dd.defaultDict['out_folder'],
                     no_raw=dd.defaultDict['no_raw']):

    vaccine, name = download_vaccine_data()

    col_names = vaccine.columns
    vaccine = vaccine.replace('-', np.nan)
    vaccine[col_names[5]].fillna(vaccine[col_names[5]].mean(), inplace=True)
    
    # get population data for all countys
    population = getPopulationData.get_age_population_data(write_df=False)

    # create empty dataframe
    new_df = pandas.DataFrame([], columns=['ID_County', 'Administrated_Vaccines', 'First_Shot', 'Full_Vaccination',
                                           'Ratio_All', 'Ratio_Young', 'Ratio_Old'])

    # loop over all states
    for state, i in zip(vaccine['Unnamed: 0'], range(len(vaccine))):
        # get population data of current state
        pop_state = population[(population['ID_County'].values/1000).astype(int) == int(state)]

        # get county ids
        id_state = pop_state.ID_County

        # compute population ratio between countys and their state
        ratio = pop_state['Total'].values/sum(pop_state['Total'])

        # scale vaccine values to countys
        state_values = np.outer(ratio,vaccine.values[i, 1:4])

        # store scaled values in dataframe
        data = np.zeros((len(ratio), len(vaccine.values[i,:7])))
        data[:,0] = id_state
        data[:,1:4] = state_values
        data[:,4] =vaccine.values[i, 4]
        data[:,5] =vaccine.values[i, 5]
        if str(vaccine.values[i, 6]) == 'nan':
            pop_columns = ['<3 years', '3-5 years', '6-14 years', '15-17 years', '18-24 years',
              '25-29 years', '30-39 years', '40-49 years', '50-64 years',
              '65-74 years', '>74 years']
            data[:, 6] = np.divide(pop_state['Total'].values*vaccine.values[i, 4] \
                         - (np.sum([pop_state[x].values for x in pop_columns[:8]]) + pop_state[pop_columns[8]].values*2/3) * vaccine.values[i, 5],\
                         pop_state[pop_columns[8]].values*1/3 + np.sum(pop_state[x].values for x in pop_columns[9:]))
        else:
            data[:, 6] = vaccine.values[i,6]
        state_df = pandas.DataFrame(data, columns=new_df.columns)
        new_df = new_df.append(state_df)

    directory = out_folder
    directory = os.path.join(directory, 'Germany/')
    gd.check_dir(directory)

    gd.write_dataframe(new_df, directory, name, file_format)
    if file_format=='json_timeasstring':
        file_format = 'json'
    print('file written to:', directory + name + '.' + file_format)


def main():
    """! Main program entry."""

    arg_dict = gd.cli("vaccine")
    get_vaccine_data(**arg_dict)


if __name__ == "__main__":

   main()