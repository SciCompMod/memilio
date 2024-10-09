import subprocess
import os
import pandas as pd

from plot_lct_fictional_noage import *


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
    cases = [1, 2, 3]
    R0s = list([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])

    for case in cases:
        if case == 0:
            # Download data for long simulation with different R0s.
            data_dir = "../data/simulation_lct_noage/riseR0long/"
            simulationdays = [140, 100, 90, 80, 60, 60, 60, 60, 60]
            run_simulation(R0s, simulationdays, data_dir)
        elif case == 1:
            # Plots to compare time and size of epidemic peaks.
            data_dir = "../data/simulation_lct_noage/riseR0long/"
            plot_maxpeak_incidence(lambda R0, subcompartment: get_file_name(R0, subcompartment, data_dir=data_dir),
                                   R0s, list([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), filename_plot="compare_peak_size")
            plot_day_peak_incidence(lambda R0, subcompartment: get_file_name(R0, subcompartment, data_dir=data_dir),
                                    R0s, list([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), filename_plot="compare_peak_days")
        elif case == 2:
            # All 3d Plots for normal parameters.
            folder = "../data/simulation_lct_noage/riseR0short/"
            nums_subcomp = [10, 50]
            for num_subcomp in nums_subcomp:
                plot_subcompartments3D(get_file_name(
                    2, num_subcomp, folder, True), num_subcomp, 1, 1, filename_plot="subcompartments"+f"{num_subcomp}"+"_exposed")
                plot_subcompartments3D(get_file_name(
                    2, num_subcomp, folder, True), num_subcomp, 2, 1, filename_plot="subcompartments"+f"{num_subcomp}"+"_carrier")
                plot_subcompartments3D(get_file_name(
                    2, num_subcomp, folder, True), num_subcomp, 3, 1, filename_plot="subcompartments"+f"{num_subcomp}"+"_infected")
        elif case == 3:
            # All 3d Plots for swapped values TE and TC.
            folder = "../data/simulation_lct_noage/riseR0shortswappedTETC/"
            num_subcomp = 50
            plot_subcompartments3D(get_file_name(
                2, num_subcomp, folder, True), num_subcomp, 2, 1, filename_plot="subcompartments"+f"{num_subcomp}"+"_carrier_swappedTETC")
            plot_new_infections([get_file_name(2, 1, folder), get_file_name(2, 3, folder), get_file_name(2, 10, folder), get_file_name(2, 50, folder)],
                                20000, legendplot=list(["ODE", "LCT3", "LCT10", "LCT50"]),
                                filename_plot="new_infections_rise2.0_swappedTETC")


if __name__ == "__main__":

    main()
