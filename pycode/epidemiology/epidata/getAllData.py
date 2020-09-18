from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from epidemiology.epidata import getSpainData
from epidemiology.epidata import getJHData
from epidemiology.epidata import getPopulationData
from epidemiology.epidata import getRKIData
from epidemiology.epidata import getDIVIData


def get_all_data(read_data=dd.defaultDict['read_data'],
                 out_form=dd.defaultDict['out_form'],
                 out_folder=dd.defaultDict['out_folder'],
                 make_plot=dd.defaultDict['make_plot'],
                 end_date=dd.defaultDict['start_date'],
                 start_date=dd.defaultDict['end_date'],
                 update_data=dd.defaultDict['update_data'],
                 ):

    getRKIData.get_rki_data(read_data, out_form, out_folder, make_plot)
    getSpainData.get_spain_data(read_data, out_form, out_folder)
    getPopulationData.get_population_data(read_data, out_form, out_folder)
    getJHData.get_jh_data(read_data, out_form, out_folder)
    getDIVIData.get_divi_data(read_data, out_form, out_folder, end_date, start_date, update_data)


def main():
    [read_data, out_form, out_folder, make_plot, end_date, start_date, update_data] = gd.cli("all")
    get_all_data(read_data, out_form, out_folder, make_plot, end_date, start_date, update_data)


if __name__ == "__main__":

   main()