import sys
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import h5py
from datetime import datetime

def plot_both_simulations_in_one_figure(path_memilio, path_panvXabmSim, colormap='Set1', xtick_step=150, show90=False):
    


def main():
    """ Main function for CLI usage. """
    parser = argparse.ArgumentParser(
        description="Plot infection state and location type results.")
    parser.add_argument("--name-of-simulation", type=str,
                        help="Name of the simulation")
    parser.add_argument("--path-to-memilio-sim",
                        help="Path to the Memilio simulation results file")
    parser.add_argument("--path-to-panvXabmSim",
                        help="Path to the panvXabmSim simulation results file")
    parser.add_argument("--colormap", type=str,
                        default='Set1', help="Matplotlib colormap")
    parser.add_argument("--xtick-step", type=int,
                        default=150, help="Step for x-axis ticks (usually hours)")
    parser.add_argument("--s90percentile", action="store_true",
                        help="If set, plot 90% percentile as well")
    args = parser.parse_args()

    if args.path_to_memilio_sim and args.path_to_panvXabmSim:
        print("Both Memilio and panvXabmSim paths provided. Plotting both results.")
        plot_both_simulations_in_one_figure(
            path_memilio=args.path_to_memilio_sim,
            path_panvXabmSim=args.path_to_panvXabmSim,
            colormap=args.colormap,
            xtick_step=args.xtick_step,
            show90=args.s90percentile
        )

    if not args.path_to_memilio_sim and not args.path_to_panvXabmSim:
        print("Please provide a path to Memilio or panvXabmSim results.")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("No arguments provided. Running in interactive mode.")
        main()
    else:
        print("Running in CLI mode with provided arguments.")
        main()
        sys.exit(0)
    sys.exit(1)
