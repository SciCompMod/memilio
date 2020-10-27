import os
import sys
import pandas as pd
import numpy as np

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from epidemiology.epidata import getRKIData as grd
from epidemiology.epidata import getJHData as gjd


def get_rki_data_with_estimations(read_data=dd.defaultDict['read_data'],
                                  out_form=dd.defaultDict['out_form'],
                                  out_folder=dd.defaultDict['out_folder']):

    data_path = os.path.join(out_folder, 'Germany/')

    # get rki data
    grd.get_rki_data(read_data, out_form, out_folder, False)

    # get data from John Hopkins University
    gjd.get_jh_data(read_data, out_form, out_folder)

    # Now we now which data is generated and we can use it
    # read in jh data

    jh_data_file = os.path.join(data_path, 'whole_country_Germany_jh.json')
    df_jh = pd.read_json(jh_data_file)

    confirmed = dd.EngEng['confirmed']
    recovered = dd.EngEng['recovered']

    df_jh["Fraction"] =  df_jh[recovered] / df_jh[confirmed]

    # There might be sum NaN due to devision with 0
    df_jh["Fraction"] = df_jh["Fraction"].fillna(0)

    # get data from rki and make new data
    rki_files_to_change = ["all_germany_rki", "all_gender_rki", "all_age_rki",
                           "all_state_rki", "all_state_gender_rki", "all_state_age_rki",
                           "all_county_rki", "all_county_gender_rki", "all_county_age_rki"]

    rki_data_file = os.path.join(data_path, rki_files_to_change[0] + ".json")
    df_rki = pd.read_json(rki_data_file)

    #df.Date.dt.strftime('%Y-%m-%d')

    print(df_rki)

    # TODO: think about to add age dependent weight function

def main():
    [read_data, out_form, out_folder] = gd.cli("rkiest")
    get_rki_data_with_estimations(read_data, out_form, out_folder)


if __name__ == "__main__":
    main()


