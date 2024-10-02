import subprocess
import os
import pandas as pd

from plot_lct_fictional_noage import plot_new_infections, get_maxpeak_incidence, plot_maxpeak_incidence, plot_subcompartments3D, plot_day_peak_incidence


def run_simulation(R0s, simulationdays, save_dir=""):
    for R0 in range(len(R0s)):
        subprocess.call([f"./../build/bin/lct_fictional_noage",
                        f"{R0s[R0]}", f"{simulationdays[R0]}", save_dir])


def plot_fixed_subcompartment(data_dir, R0s, subcompartment, filename):
    paths = []
    legendplot = list([])
    for R0 in R0s:
        paths.append(get_file_name(R0, subcompartment))
        legendplot.append(f"{int(R0)}"+".0")
    plot_new_infections(paths,
                        14e6, legendplot,
                        filename)


def get_file_name(R0, subcompartment, data_dir, boolsubcomp=False):
    filename = "fictional_lct_" + f"{int(R0)}"+".0_" + f"{subcompartment}"
    if boolsubcomp:
        filename += "_subcompartments"
    return os.path.join(data_dir, filename)


def main():
    data_dir = "../data/simulation_lct_noage/riseR0long/"
    # subcompartments = 5
    R0s = list([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    simulationdays = [140, 100, 90, 80, 60, 60, 60, 60, 60]
    # run_simulation(R0s, simulationdays, data_dir)
    plot_maxpeak_incidence(lambda R0, subcompartment: get_file_name(R0, subcompartment, data_dir=data_dir),
                           R0s, list([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), filename_plot="compare_peak_size")
    plot_day_peak_incidence(lambda R0, subcompartment: get_file_name(R0, subcompartment, data_dir=data_dir),
                            R0s, list([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), filename_plot="compare_peak_days")
    # plot_subcompartments3D(get_file_name(
    #     2, 50, "../data/simulation_lct_noage/riseR0short/", True), 50, 2, 4, filename_plot="a")


if __name__ == "__main__":

    main()
