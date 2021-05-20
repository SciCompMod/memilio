"""
@file getSimulationData.py

@brief Executes all data downloads which belong to the epidata package and downloads external data

The functions which are called are:
- getRKIData.get_rki_data
- getPopulationData.get_population_data
- getVacccineData.get_vaccine_data
- getDIVIData.get_divi_data
"""


from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from epidemiology.epidata import getVaccineData
from epidemiology.epidata import getPopulationData
from epidemiology.epidata import getRKIData
from epidemiology.epidata import getDIVIData


def get_simulation_data(read_data=dd.defaultDict['read_data'],
                 out_form=dd.defaultDict['out_form'],
                 out_folder=dd.defaultDict['out_folder'],
                 end_date=dd.defaultDict['end_date'],
                 fill_dates=dd.defaultDict['fill_dates'],
                 make_plot=dd.defaultDict['make_plot'],
                 moving_average=dd.defaultDict['moving_average'],
                 split_berlin=dd.defaultDict['split_berlin'],
                 start_date=dd.defaultDict['start_date'],
                 update_data=dd.defaultDict['update_data'],
                 ):
    """! Downloads all data from external sources

    The functions which are called are:
    - getRKIData.get_rki_data(read_data, out_form, out_folder, make_plot)
    - getPopulationData.get_population_data(read_data, out_form, out_folder)
    - getVaccineData.get_jh_data(read_data, out_form, out_folder)
    - getDIVIData.get_divi_data(read_data, out_form, out_folder, end_date, start_date, update_data)

    Keyword arguments:
    @param read_data False [Default] or True. Defines if data is read from file or downloaded.
    @param out_form File format which is used for writing the data. Default defined in defaultDict.
    @param out_folder Path to folder where data is written in folder out_folder/Germany.
    @param end_date [Optional] Date to stop to download data [Default = today].
    @param make_plot False [Default] or True. Defines if plots are generated with matplotlib.
    @param split_berlin True [Default] or False. Defines if Berlin counties is fused to just on county.
    @param start_date [Optional] Date to start to download data [Default = 2020-04-24].
    @param update_data "True" if existing data is updated or
    "False [Default]" if it is downloaded for all dates from start_date to end_date.
    """

    getRKIData.get_rki_data(read_data, out_form, out_folder, fill_dates, make_plot, moving_average, split_berlin)
    getPopulationData.get_population_data(read_data, out_form, out_folder)
    getPopulationData.get_age_population_data(read_data, out_form, out_folder)
    getDIVIData.get_divi_data(read_data, out_form, out_folder, end_date, start_date, update_data)
    getVaccineData.get_vaccine_data(read_data, out_form, out_folder)


def main():
    """! Main program entry."""

    [read_data, out_form, out_folder, end_date, fill_dates, make_plot, moving_average, split_berlin, start_date,
     update_data] = gd.cli("sim")
    get_simulation_data(read_data, out_form, out_folder, end_date, fill_dates, make_plot, moving_average, split_berlin,
                        start_date, update_data)


if __name__ == "__main__":
    main()
