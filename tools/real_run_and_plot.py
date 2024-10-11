
from plot_lct_realistic_scenario import *
import os


def run_simulation(R0s, simulationdays, save_dir=""):
    for R0 in range(len(R0s)):
        subprocess.call([f"./../build/bin/lct_realistic_scenario"])


def get_file_name(start_date, subcompartment, data_dir, boolagedistributed=False, boolsubcomp=False):
    filename = "real_" + start_date+"_" + f"{subcompartment}"+"_"
    if boolagedistributed:
        filename += "ageres"
    else:
        filename += "accumulated"
    if boolsubcomp:
        filename += "_subcompartments"
    return os.path.join(data_dir, filename)


def main():
    folder = "../data/simulation_lct_real/"
    start_date = '2020-9-1'
    num_subcomp = 10
    plot_new_infections_real([get_file_name(start_date, num_subcomp, folder)],
                             20000, legendplot=list(["LCT10"]),
                             filename_plot="real_new_infections_"+start_date)


if __name__ == "__main__":
    main()
