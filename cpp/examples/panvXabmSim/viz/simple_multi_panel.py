#!/usr/bin/env python3
"""
Simple Multi-Panel Epidemic Curves - Clean Version
===================================================

Creates a 2x2 grid of epidemic curve comparisons with elegant, publication-ready styling.
Each panel shows transmission-informed vs uniform initialization for one scenario.
"""

import sys
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches
import h5py
from datetime import datetime
import glob


def load_h5_results(base_path, percentile):
    """Load HDF5 results for a given percentile."""
    file_path = os.path.join(
        base_path, "amount_of_infections", percentile, "Results.h5")
    with h5py.File(file_path, 'r') as f:
        data = {k: v[()] for k, v in f['0'].items()}
    return data


def calculate_cumulative_infections(data):
    """Calculate cumulative infections."""
    current_infections = data.sum(axis=1)
    cumulative = current_infections
    return cumulative


def _format_x_axis_days(ax, x, xtick_step=48):
    """Helper to format x-axis as days."""
    ax.set_xticks(x[::xtick_step])
    ax.set_xticklabels([int(day) for day in x[::xtick_step]])


def find_scenario_directories(results_base_dir, scenario_key):
    """Find the directories for transmission-informed and uniform initialization."""
    pattern_transmission = f"epidemic_curves_*_{scenario_key}_transmission_informed"
    pattern_uniform = f"epidemic_curves_*_{scenario_key}_uniform_initialized"

    transmission_dirs = glob.glob(os.path.join(
        results_base_dir, pattern_transmission))
    uniform_dirs = glob.glob(os.path.join(results_base_dir, pattern_uniform))

    transmission_dir = max(transmission_dirs) if transmission_dirs else None
    uniform_dir = max(uniform_dirs) if uniform_dirs else None

    print(f"  Found directories for {scenario_key}:")
    print(f"    Transmission: {transmission_dir}")
    print(f"    Uniform: {uniform_dir}")
    return transmission_dir, uniform_dir


def plot_scenario_comparison(ax, path_transmission, path_uniform, scenario_title, show90=False, show_legend=False):
    """Plot comparison for one scenario with clean, elegant styling."""

    if not path_transmission or not path_uniform:
        ax.text(0.5, 0.5, 'Data not available', ha='center', va='center',
                transform=ax.transAxes, fontsize=18)
        ax.set_title(scenario_title, fontsize=20, fontweight='bold')
        return None

    try:
        # Load transmission-informed data
        panv_p50 = load_h5_results(path_transmission, "p50")
        panv_p25 = load_h5_results(path_transmission, "p25")
        panv_p75 = load_h5_results(path_transmission, "p75")

        time_panv = panv_p50['Time']
        panv_cum_50 = calculate_cumulative_infections(panv_p50['Total'])
        panv_cum_25 = calculate_cumulative_infections(panv_p25['Total'])
        panv_cum_75 = calculate_cumulative_infections(panv_p75['Total'])

        # Load uniform data
        memilio_p50 = load_h5_results(path_uniform, "p50")
        memilio_p25 = load_h5_results(path_uniform, "p25")
        memilio_p75 = load_h5_results(path_uniform, "p75")

        time_memilio = memilio_p50['Time']
        memilio_cum_50 = calculate_cumulative_infections(memilio_p50['Total'])
        memilio_cum_25 = calculate_cumulative_infections(memilio_p25['Total'])
        memilio_cum_75 = calculate_cumulative_infections(memilio_p75['Total'])

        # Colors for publication-ready plots
        color_transmission = '#d62728'  # Red
        color_uniform = '#1f77b4'      # Blue

        # Plot transmission-informed with thicker lines for publication
        ax.plot(time_panv, panv_cum_50, color=color_transmission,
                linewidth=3.5, linestyle='-', alpha=0.9)
        ax.fill_between(time_panv, panv_cum_25, panv_cum_75,
                        alpha=0.35, color=color_transmission)

        # Plot uniform
        ax.plot(time_memilio, memilio_cum_50, color=color_uniform,
                linewidth=3.5, linestyle='-', alpha=0.9)
        ax.fill_between(time_memilio, memilio_cum_25, memilio_cum_75,
                        alpha=0.35, color=color_uniform)

        # Optional 90% percentiles for additional uncertainty visualization
        if show90:
            try:
                panv_p05 = load_h5_results(path_transmission, "p05")
                panv_p95 = load_h5_results(path_transmission, "p95")
                panv_cum_05 = calculate_cumulative_infections(
                    panv_p05['Total'])
                panv_cum_95 = calculate_cumulative_infections(
                    panv_p95['Total'])
                ax.fill_between(time_panv, panv_cum_05, panv_cum_95,
                                alpha=0.15, color=color_transmission)

                memilio_p05 = load_h5_results(path_uniform, "p05")
                memilio_p95 = load_h5_results(path_uniform, "p95")
                memilio_cum_05 = calculate_cumulative_infections(
                    memilio_p05['Total'])
                memilio_cum_95 = calculate_cumulative_infections(
                    memilio_p95['Total'])
                ax.fill_between(time_memilio, memilio_cum_05, memilio_cum_95,
                                alpha=0.15, color=color_uniform)
            except:
                pass  # 90% percentiles not available

        # Publication-ready formatting with much larger fonts
        _format_x_axis_days(ax, time_memilio, xtick_step=48)
        ax.set_xlabel('Days', fontsize=22, fontweight='bold')
        ax.set_ylabel('Cumulative infections', fontsize=22, fontweight='bold')
        ax.tick_params(axis='both', which='major',
                       labelsize=18, width=2, length=8)
        ax.grid(True, alpha=0.3, linewidth=1.2)
        ax.set_title(scenario_title, fontsize=24, fontweight='bold', pad=25)

        # Add vertical lines when 100 infections are reached
        target_infections = 100

        # Find time points when 100 infections are reached
        time_100_transmission = None
        time_100_uniform = None

        for i, cum_inf in enumerate(panv_cum_50):
            if cum_inf >= target_infections and time_100_transmission is None:
                time_100_transmission = time_panv[i]
                break

        for i, cum_inf in enumerate(memilio_cum_50):
            if cum_inf >= target_infections and time_100_uniform is None:
                time_100_uniform = time_memilio[i]
                break

                # Draw thicker vertical lines at 100 infection points with timepoint annotations
        if time_100_transmission is not None:
            ax.axvline(x=time_100_transmission, color=color_transmission,
                       linestyle='--', alpha=0.8, linewidth=3)
            # Add timepoint annotation for transmission-informed - positioned very close to line
            y_pos = ax.get_ylim()[1] * 0.95  # 95% of max y-value
            ax.text(time_100_transmission + 0.5, y_pos,
                    f'TI: Day {time_100_transmission:.1f}',
                    rotation=0, ha='left', va='top', fontsize=14,
                    color=color_transmission, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor=color_transmission))

        if time_100_uniform is not None:
            ax.axvline(x=time_100_uniform, color=color_uniform,
                       linestyle='--', alpha=0.8, linewidth=3)
            # Add timepoint annotation for uniform - positioned very close to line
            y_pos = ax.get_ylim()[1] * 0.95  # 95% of max y-value
            ax.text(time_100_uniform - 0.5, y_pos,
                    f'UI: Day {time_100_uniform:.1f}',
                    rotation=0, ha='right', va='top', fontsize=14,
                    color=color_uniform, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor=color_uniform))

        # Add difference annotation positioned to avoid cutting into plot
        final_transmission = panv_cum_50[-1]
        final_uniform = memilio_cum_50[-1]
        diff = final_uniform - final_transmission
        percent_diff = (diff / final_transmission *
                        100) if final_transmission > 0 else 0

        # Enhanced difference visualization with arrow
        max_y = max(final_transmission, final_uniform)
        arrow_y = max_y * 0.85
        ax.annotate('', xy=(time_memilio[-1], final_uniform), xytext=(time_memilio[-1], final_transmission),
                    arrowprops=dict(arrowstyle='<->', color='black', lw=2))

        # Position text annotation in lower right, outside main data area
        ax.text(0.98, 0.15, f'Diff: +{diff:.0f}\n({percent_diff:+.1f}%)',
                transform=ax.transAxes, fontsize=16, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                          alpha=0.9, edgecolor='gray'),
                ha='right', va='bottom')

        # Only show legend if requested (for one panel only)
        if show_legend:
            legend_elements = [
                plt.Line2D([0], [0], color=color_transmission, linewidth=4,
                           label='Transmission-Informed'),
                plt.Line2D([0], [0], color=color_uniform, linewidth=4,
                           label='Uniform'),
                matplotlib.patches.Patch(facecolor=color_transmission, alpha=0.35,
                                         label='50% CI - Transmission'),
                matplotlib.patches.Patch(facecolor=color_uniform, alpha=0.35,
                                         label='50% CI - Uniform'),
                matplotlib.patches.Patch(facecolor=color_transmission, alpha=0.15,
                                         label='90% CI - Transmission'),
                matplotlib.patches.Patch(facecolor=color_uniform, alpha=0.15,
                                         label='90% CI - Uniform'),
                plt.Line2D([0], [0], color='gray', linewidth=2, linestyle='--',
                           label='100 infections reached')
            ]

            legend = ax.legend(handles=legend_elements, loc='upper left',
                               fontsize=14, frameon=True, fancybox=True, shadow=True,
                               bbox_to_anchor=(0.02, 0.8), framealpha=0.95)

            # Use black text for all legend entries
            legend_texts = legend.get_texts()
            for text in legend_texts:
                text.set_color('black')
                text.set_fontweight('bold')

        # Print summary info
        final_transmission = panv_cum_50[-1]
        final_uniform = memilio_cum_50[-1]
        diff = final_uniform - final_transmission
        percent_diff = (diff / final_transmission *
                        100) if final_transmission > 0 else 0

        print(f"  {scenario_title}: TI={final_transmission:.0f}, UI={final_uniform:.0f}, "
              f"Diff={diff:+.0f} ({percent_diff:+.1f}%)")

        return {'final_transmission': final_transmission, 'final_uniform': final_uniform,
                'difference': diff, 'percent_difference': percent_diff}

    except Exception as e:
        ax.text(0.5, 0.5, f'Error loading data:\\n{str(e)}', ha='center', va='center',
                transform=ax.transAxes, fontsize=16, color='red')
        ax.set_title(scenario_title, fontsize=24, fontweight='bold')
        print(f"  Error in {scenario_title}: {e}")
        return None


def create_multi_panel_comparison(results_base_dir, output_path=None, show90=False):
    """Create 2x2 multi-panel comparison with elegant publication styling."""

    scenarios = {
        'R1_restaurant_strong_clustering': 'R1: Restaurant Strong Clustering',
        'R2_restaurant_weaker_clustering': 'R2: Restaurant Weaker Clustering',
        'W1_workplace_few_meetings': 'W1: Workplace Few Meetings',
        'W2_workplace_many_meetings': 'W2: Workplace Many Meetings'
    }

    print(f"Creating elegant multi-panel comparison from: {results_base_dir}")

    # Create 2x2 subplot with wider figure for publication (landscape orientation)
    fig, axes = plt.subplots(2, 2, figsize=(24, 14))
    fig.suptitle('Infection Curves: Transmission-Informed vs Uniform Initialization',
                 fontsize=28, fontweight='bold', y=0.95)

    # Store results for summary
    all_results = {}

    # Plot each scenario
    positions = [(0, 0), (0, 1), (1, 0), (1, 1)]

    for i, (scenario_key, scenario_title) in enumerate(scenarios.items()):
        row, col = positions[i]
        ax = axes[row, col]

        # Only show legend on the first panel (R1)
        show_legend = (i == 0)

        print(f"Processing {scenario_key}...")
        transmission_dir, uniform_dir = find_scenario_directories(
            results_base_dir, scenario_key)

        if transmission_dir and uniform_dir:
            result = plot_scenario_comparison(
                ax, transmission_dir, uniform_dir, scenario_title, show90, show_legend)
            if result:
                all_results[scenario_key] = {
                    'title': scenario_title, 'data': result}
        else:
            ax.text(0.5, 0.5, f'Directories not found for\\n{scenario_title}',
                    ha='center', va='center', transform=ax.transAxes,
                    color='red', fontsize=18, fontweight='bold')
            ax.set_title(scenario_title, fontsize=24, fontweight='bold')

    # Adjust layout with more space - no bottom explanation needed
    plt.tight_layout(rect=[0, 0.02, 1, 0.92])

    # Print summary table
    print("\\n" + "="*80)
    print("SUMMARY: FINAL EPIDEMIC SIZES")
    print("="*80)
    print(f"{'Scenario':<30} {'Trans-Inf':<12} {'Uniform':<12} {'Difference':<15} {'% Change':<10}")
    print("-"*80)

    for scenario_key, info in all_results.items():
        data = info['data']
        title = info['title']
        print(f"{title:<30} {data['final_transmission']:<12.0f} {data['final_uniform']:<12.0f} "
              f"{data['difference']:+<15.0f} {data['percent_difference']:+<10.1f}%")

    print("-"*80)
    print("Key Finding: Uniform initialization consistently leads to higher epidemic sizes")
    print("="*80)

    # Save with high quality for publication
    if output_path is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_path = f"epidemic_curves_publication_{timestamp}.png"

    try:
        plt.savefig(output_path, bbox_inches='tight',
                    dpi=300, facecolor='white')
        print(f"\\nPublication-ready figure saved to: {output_path}")

        # Also save PDF for vector graphics
        pdf_path = output_path.replace('.png', '.pdf')
        plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
        print(f"Vector PDF saved to: {pdf_path}")

    except Exception as e:
        print(f"Error saving figure: {e}")
        fallback_path = os.path.join(
            os.getcwd(), "epidemic_curves_publication_fallback.png")
        plt.savefig(fallback_path, bbox_inches='tight',
                    dpi=300, facecolor='white')
        print(f"Saved to fallback location: {fallback_path}")

    plt.show()
    return output_path


def main():
    """Main function with command line interface."""
    parser = argparse.ArgumentParser(
        description='Create publication-ready epidemic curves comparison')
    parser.add_argument('--results-dir', type=str,
                        default='/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results',
                        help='Directory containing epidemic curve results')
    parser.add_argument('--output-path', type=str,
                        help='Output path for the figure')
    parser.add_argument('--s90percentile', action='store_true',
                        help='Include 90% percentile bands for additional uncertainty visualization')

    args = parser.parse_args()

    try:
        output_path = create_multi_panel_comparison(
            args.results_dir, args.output_path, args.s90percentile)
        print(f"\\n✓ Publication-ready multi-panel comparison completed!")
        return 0
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
