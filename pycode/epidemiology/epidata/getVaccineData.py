import pandas
import requests
import io
import numpy as np
import os

from epidemiology.epidata  import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
import getPopulationData

def get_vaccine_data(read_data=dd.defaultDict['read_data'],
                       out_form=dd.defaultDict['out_form'],
                       out_folder=dd.defaultDict['out_folder']):

    # get rki vaccine data
    url = 'https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Daten/Impfquotenmonitoring.xlsx?__blob=publicationFile'
    header = {'User-Agent': 'Mozilla/5.0'}
    r = requests.get(url, headers=header)
    with io.BytesIO(r.content) as fh:
        df = pandas.io.excel.ExcelFile(fh, engine='openpyxl')
        sheet_names = df.sheet_names
        vaccine = pandas.read_excel(df, sheet_name=sheet_names[2], header=[1])[:16].drop(['Unnamed: 1', 'Unnamed: 10'], axis=1)

    # get population data for all countys
    population = getPopulationData.get_age_population_data(write_df=False)


    # create empty dataframe
    new_df = pandas.DataFrame([], columns=['Id_County', 'Alter1', 'Beruf1', 'Krank1', 'Pflegeheim1', 'Alter2', 'Beruf2', 'Krank2', 'Pflegeheim2'])

    # loop over all states
    for state, i in zip(vaccine['Unnamed: 0'], range(len(vaccine))):
        # get population data of current state
        pop_state = population[(population['ID_County'].values/1000).astype(int) == int(state)]

        # get county ids
        id_state = pop_state.ID_County

        # compute population ratio between countys and their state
        ratio = pop_state['Total'].values/sum(pop_state['Total'])

        # scale vaccine values to countys
        state_values = np.outer(ratio,vaccine.values[i, 1:9])

        # store scaled values in dataframe
        data = np.zeros((len(ratio), len(vaccine.values[i,:9])))
        data[:,0] = id_state
        data[:,1:] = state_values
        state_df = pandas.DataFrame(data, columns=new_df.columns)
        new_df = new_df.append(state_df)

    directory = out_folder
    directory = os.path.join(directory, 'Germany/')
    gd.check_dir(directory)

    name = sheet_names[2].replace('.', '_', 1).replace('.', '', 1)

    gd.write_dataframe(new_df, directory, name, out_form)
    print('file written to:', directory + name + '.' + out_form)



def main():
    """! Main program entry."""

    [read_data, out_form, out_folder] = gd.cli("population")
    get_vaccine_data(read_data, out_form, out_folder)


if __name__ == "__main__":

   main()