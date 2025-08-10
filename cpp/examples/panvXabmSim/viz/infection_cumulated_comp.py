import sys
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import h5py
from datetime import datetime


def load_h5_results(base_path, percentile):
    """ Reads HDF5 results for a given group and percentile.

    @param[in] base_path Path to results directory.
    @param[in] percentile Subdirectory for percentile (e.g. 'p50').
    @return Dictionary with data arrays. Keys are dataset names from the HDF5 file 
            (e.g., 'Time', 'Total', age group names like 'Group1', 'Group2', etc.).
            Values are numpy arrays containing the corresponding time series data.
    """
    file_path = os.path.join(base_path, percentile, "Results.h5")
    with h5py.File(file_path, 'r') as f:
        data = {k: v[()] for k, v in f['0'].items()}
    return data


def calculate_cumulative_infections(data):
    """ Calculate cumulative infections from infection state data.

    @param[in] data Data array with infection states.
    @return Array with cumulative infections (all infected states summed and accumulated).
    """
    # Sum all infection states (excluding susceptible state 0)
    # States: 1=Exposed, 2=I_Asymp, 3=I_Symp, 4=I_Severe, 5=I_Critical, 7=Dead
    current_infections = data.sum(axis=1)
    # Calculate cumulative infections
    cumulative = current_infections

    return cumulative


def _format_x_axis_days(x, xtick_step):
    """ Helper to format x-axis as days. """
    plt.gca().set_xticks(x[::xtick_step])
    plt.gca().set_xticklabels([int(day) for day in x[::xtick_step]])


def plot_both_simulations_in_one_figure(path_memilio, path_panvXabmSim, output_path=None, colormap='Set1', xtick_step=150, show90=False):
    """ Plots cumulative infections for both simulations in one figure with percentiles.

    @param[in] path_memilio Path to Memilio simulation results.
    @param[in] path_panvXabmSim Path to panvXabmSim simulation results.
    @param[in] output_path Path where to save the comparison plots.
    @param[in] colormap Matplotlib colormap name.
    @param[in] xtick_step Step size for x-axis ticks.
    @param[in] show90 If True, plot 90% percentile bands in addition to 50% percentile.
    """

    # Load data for Memilio simulation
    print("Loading Memilio simulation data...")
    memilio_p50 = load_h5_results(path_memilio, "p50")
    memilio_p25 = load_h5_results(path_memilio, "p25")
    memilio_p75 = load_h5_results(path_memilio, "p75")

    time_memilio = memilio_p50['Time']
    memilio_total_50 = memilio_p50['Total']
    memilio_total_25 = memilio_p25['Total']
    memilio_total_75 = memilio_p75['Total']

    # Calculate cumulative infections for Memilio
    memilio_cum_50 = calculate_cumulative_infections(memilio_total_50)
    memilio_cum_25 = calculate_cumulative_infections(memilio_total_25)
    memilio_cum_75 = calculate_cumulative_infections(memilio_total_75)

    # Optional 90% percentiles for Memilio
    memilio_cum_05 = memilio_cum_95 = None
    if show90:
        memilio_p05 = load_h5_results(path_memilio, "p05")
        memilio_p95 = load_h5_results(path_memilio, "p95")
        memilio_total_05 = memilio_p05['Total']
        memilio_total_95 = memilio_p95['Total']
        memilio_cum_05 = calculate_cumulative_infections(memilio_total_05)
        memilio_cum_95 = calculate_cumulative_infections(memilio_total_95)

    # Load data for panvXabmSim simulation
    print("Loading panvXabmSim simulation data...")
    panv_p50 = load_h5_results(path_panvXabmSim, "p50")
    panv_p25 = load_h5_results(path_panvXabmSim, "p25")
    panv_p75 = load_h5_results(path_panvXabmSim, "p75")

    time_panv = panv_p50['Time']
    panv_total_50 = panv_p50['Total']
    panv_total_25 = panv_p25['Total']
    panv_total_75 = panv_p75['Total']

    # Calculate cumulative infections for panvXabmSim
    panv_cum_50 = calculate_cumulative_infections(panv_total_50)
    panv_cum_25 = calculate_cumulative_infections(panv_total_25)
    panv_cum_75 = calculate_cumulative_infections(panv_total_75)

    # Optional 90% percentiles for panvXabmSim
    panv_cum_05 = panv_cum_95 = None
    if show90:
        panv_p05 = load_h5_results(path_panvXabmSim, "p05")
        panv_p95 = load_h5_results(path_panvXabmSim, "p95")
        panv_total_05 = panv_p05['Total']
        panv_total_95 = panv_p95['Total']
        panv_cum_05 = calculate_cumulative_infections(panv_total_05)
        panv_cum_95 = calculate_cumulative_infections(panv_total_95)

    # Create the plot
    plt.figure('Cumulative_Infections_Comparison', figsize=(12, 8))

    title = 'Cumulative Infections: Memilio vs panvXabmSim'
    if show90:
        title += ' (with 90% percentiles)'
    plt.title(title, fontsize=14)

    # Use distinct colors for the two simulations
    color_memilio = 'blue'
    color_panv = 'red'

    # Plot Memilio simulation
    plt.plot(time_memilio, memilio_cum_50, color=color_memilio,
             linewidth=2.5, linestyle='-', label='Memilio (median)')
    plt.fill_between(time_memilio, memilio_cum_25, memilio_cum_75,
                     alpha=0.3, color=color_memilio, label='Memilio (25-75%)')

    # Optional 90% percentile for Memilio
    if show90 and memilio_cum_05 is not None and memilio_cum_95 is not None:
        plt.fill_between(time_memilio, memilio_cum_05, memilio_cum_95,
                         alpha=0.15, color=color_memilio, label='Memilio (5-95%)')

    # Plot panvXabmSim simulation
    plt.plot(time_panv, panv_cum_50, color=color_panv,
             linewidth=2.5, linestyle='-', label='panvXabmSim (median)')
    plt.fill_between(time_panv, panv_cum_25, panv_cum_75,
                     alpha=0.3, color=color_panv, label='panvXabmSim (25-75%)')

    # Optional 90% percentile for panvXabmSim
    if show90 and panv_cum_05 is not None and panv_cum_95 is not None:
        plt.fill_between(time_panv, panv_cum_05, panv_cum_95,
                         alpha=0.15, color=color_panv, label='panvXabmSim (5-95%)')

    # Format axes and labels
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    # Assuming similar time scales
    _format_x_axis_days(time_memilio, xtick_step)
    plt.xlabel('Days', fontsize=12)
    plt.ylabel('Cumulative number of infections', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save the plot with error handling
    saved_successfully = False

    # Try to save to specified output path first
    if output_path:
        try:
            if os.path.isdir(output_path):
                output_file = os.path.join(
                    output_path, "cumulative_infections_comparison.png")
            else:
                output_file = output_path

            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            plt.savefig(output_file, bbox_inches='tight', dpi=300)
            print(f"Comparison plot saved to {output_file}")
            saved_successfully = True
        except Exception as e:
            print(f"Error saving to specified output path: {e}")

    # Try to save to Memilio directory if output path failed or not specified
    if not saved_successfully:
        try:
            output_dir = os.path.dirname(path_memilio)
            output_file = os.path.join(
                output_dir, "cumulative_infections_comparison.png")
            plt.savefig(output_file, bbox_inches='tight', dpi=300)
            print(f"Comparison plot saved to {output_file}")
            saved_successfully = True
        except PermissionError:
            print(f"Permission denied: Cannot save to {output_dir}")
        except Exception as e:
            print(f"Error saving to Memilio directory: {e}")

    # Try to save to panvXabmSim directory
    if not saved_successfully:
        try:
            output_dir_panv = os.path.dirname(path_panvXabmSim)
            output_file_panv = os.path.join(
                output_dir_panv, "cumulative_infections_comparison.png")
            plt.savefig(output_file_panv, bbox_inches='tight', dpi=300)
            print(f"Comparison plot also saved to {output_file_panv}")
            saved_successfully = True
        except PermissionError:
            print(f"Permission denied: Cannot save to {output_dir_panv}")
        except Exception as e:
            print(f"Error saving to panvXabmSim directory: {e}")

    # Fallback: save to current working directory
    if not saved_successfully:
        try:
            fallback_file = os.path.join(
                os.getcwd(), "cumulative_infections_comparison.png")
            plt.savefig(fallback_file, bbox_inches='tight', dpi=300)
            print(
                f"Comparison plot saved to current directory: {fallback_file}")
            saved_successfully = True
        except Exception as e:
            print(f"Error saving to current directory: {e}")

    # Final fallback: save to user's home directory
    if not saved_successfully:
        try:
            home_dir = os.path.expanduser("~")
            home_file = os.path.join(
                home_dir, "cumulative_infections_comparison.png")
            plt.savefig(home_file, bbox_inches='tight', dpi=300)
            print(f"Comparison plot saved to home directory: {home_file}")
        except Exception as e:
            print(f"Error saving to home directory: {e}")
            print("Could not save plot to any location. Please check file permissions.")

    # Print summary statistics
    print("\nSummary Statistics:")
    print(
        f"Memilio - Final cumulative infections (median): {memilio_cum_50[-1]:.0f}")
    print(
        f"Memilio - Final cumulative infections (25-75%): {memilio_cum_25[-1]:.0f} - {memilio_cum_75[-1]:.0f}")
    print(
        f"panvXabmSim - Final cumulative infections (median): {panv_cum_50[-1]:.0f}")
    print(
        f"panvXabmSim - Final cumulative infections (25-75%): {panv_cum_25[-1]:.0f} - {panv_cum_75[-1]:.0f}")

    if show90:
        print(
            f"Memilio - Final cumulative infections (5-95%): {memilio_cum_05[-1]:.0f} - {memilio_cum_95[-1]:.0f}")
        print(
            f"panvXabmSim - Final cumulative infections (5-95%): {panv_cum_05[-1]:.0f} - {panv_cum_95[-1]:.0f}")


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
    parser.add_argument("--output-path", type=str,
                        help="Output path for comparison plots")
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
            output_path=args.output_path,
            colormap=args.colormap,
            xtick_step=args.xtick_step,
            show90=args.s90percentile
        )
    else:
        print("Please provide paths to both Memilio and panvXabmSim results.")
        print("Usage example:")
        print("python infection_cumulated_comp.py --path-to-memilio-sim /path/to/memilio/results --path-to-panvXabmSim /path/to/panv/results")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("No arguments provided. Running in interactive mode.")
        main()
    else:
        print("Running in CLI mode with provided arguments.")
        main()
        sys.exit(0)
    sys.exit(1)
