from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from epidemiology.epidata import getSpainData
from epidemiology.epidata import getJHData
from epidemiology.epidata import getPopulationData
from epidemiology.epidata import getRKIData


def get_all_data(read_data=dd.defaultDict['read_data'],
                 make_plot=dd.defaultDict['make_plot'],
                 out_form=dd.defaultDict['out_form'],
                 out_folder=dd.defaultDict['out_folder']):

    getRKIData.get_rki_data(read_data, make_plot, out_form, out_folder)
    getSpainData.get_spain_data(read_data, make_plot, out_form, out_folder)
    getPopulationData.get_population_data(read_data, make_plot, out_form, out_folder)
    getJHData.get_jh_data(read_data, make_plot, out_form, out_folder)


def main():
    [read_data, make_plot, out_form, out_folder] = gd.cli('Download all possible data')
    get_all_data(read_data, make_plot, out_form, out_folder)


if __name__ == "__main__":

   main()