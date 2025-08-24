#!/usr/bin/env python3
"""
Comprehensive comparison script for all simulation scenarios.
Creates a 4-panel comparison visualization displaying all event types and simulation types.
"""

import argparse
import os
import sys
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import seaborn as sns

# Set style for better-looking plots
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")


class ScenarioComparator:
    def __init__(self, results_paths: List[str], labels: List[str], output_dir: str, timestamp: str):
        """
        Initialize the scenario comparator.

        Args:
            results_paths: List of paths to simulation results
            labels: List of labels for each simulation (format: simtype_eventtype)
            output_dir: Output directory for visualizations
            timestamp: Timestamp for file naming
        """
        self.results_paths = results_paths
        self.labels = labels
        self.output_dir = Path(output_dir)
        self.timestamp = timestamp

        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Parse labels to extract event types and simulation types
        self.scenarios = self._parse_labels()

        # Load all data
        self.data = self._load_all_data()

    def _parse_labels(self) -> Dict[str, Dict[str, str]]:
        """Parse labels to extract simulation and event type information."""
        scenarios = {}

        for i, label in enumerate(self.labels):
            if '_' in label:
                sim_type, event_type = label.split('_', 1)
                scenarios[label] = {
                    'sim_type': sim_type,
                    'event_type': event_type,
                    'path': self.results_paths[i],
                    'index': i
                }
            else:
                # Fallback for labels without underscore
                scenarios[label] = {
                    'sim_type': 'unknown',
                    'event_type': label,
                    'path': self.results_paths[i],
                    'index': i
                }

        return scenarios

    def _load_infection_data(self, results_path: str) -> Optional[pd.DataFrame]:
        """Load infection state data from a results directory."""
        infection_dir = Path(results_path) / "infection_state_per_age_group"

        if not infection_dir.exists():
            print(f"Warning: No infection state data found in {results_path}")
            return None

        # Try to find the main infection file
        possible_files = [
            "InfectionState.txt",
            "infection_state.txt",
            "infections.txt"
        ]

        data_file = None
        for filename in possible_files:
            file_path = infection_dir / filename
            if file_path.exists():
                data_file = file_path
                break

        if data_file is None:
            # Look for any .txt or .csv file
            txt_files = list(infection_dir.glob("*.txt"))
            csv_files = list(infection_dir.glob("*.csv"))

            if txt_files:
                data_file = txt_files[0]
            elif csv_files:
                data_file = csv_files[0]

        if data_file is None:
            print(f"Warning: No data files found in {infection_dir}")
            return None

        try:
            # Try to read the file (assuming it's tab or space separated)
            df = pd.read_csv(data_file, sep='\t', header=0)

            # If that doesn't work, try comma separated
            if df.shape[1] == 1:
                df = pd.read_csv(data_file, sep=',', header=0)

            # If still doesn't work, try space separated
            if df.shape[1] == 1:
                df = pd.read_csv(data_file, sep=' ', header=0)

            return df

        except Exception as e:
            print(f"Error reading data from {data_file}: {e}")
            return None

    def _load_all_data(self) -> Dict[str, pd.DataFrame]:
        """Load data for all scenarios."""
        data = {}

        for label, scenario in self.scenarios.items():
            df = self._load_infection_data(scenario['path'])
            if df is not None:
                data[label] = df
                print(f"Loaded data for {label}: {df.shape}")
            else:
                print(f"Failed to load data for {label}")

        return data

    def create_comprehensive_comparison(self, use_percentile: bool = True):
        """
        Create a comprehensive comparison visualization with all scenarios.

        Args:
            use_percentile: Whether to show 90th percentile bands
        """
        # Get unique event types and simulation types
        event_types = sorted(set(s['event_type']
                             for s in self.scenarios.values()))
        sim_types = sorted(set(s['sim_type'] for s in self.scenarios.values()))

        print(f"Event types: {event_types}")
        print(f"Simulation types: {sim_types}")

        # Create subplot layout - 2x3 grid to accommodate 5 event types
        n_events = len(event_types)
        n_cols = 3 if n_events > 4 else 2
        n_rows = (n_events + n_cols - 1) // n_cols  # Ceiling division

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
        if n_rows == 1:
            axes = axes.reshape(1, -1)
        elif n_cols == 1:
            axes = axes.reshape(-1, 1)

        # Flatten axes for easier indexing
        axes_flat = axes.flatten()

        # Colors for different simulation types
        colors = plt.cm.Set1(np.linspace(0, 1, len(sim_types)))
        sim_colors = dict(zip(sim_types, colors))

        for i, event_type in enumerate(event_types):
            ax = axes_flat[i]

            # Plot data for this event type
            for sim_type in sim_types:
                label_key = f"{sim_type}_{event_type}"

                if label_key in self.data:
                    df = self.data[label_key]

                    # Prepare data for plotting
                    if 'time' in df.columns or 'Time' in df.columns:
                        time_col = 'time' if 'time' in df.columns else 'Time'
                        time_data = df[time_col]
                    else:
                        # Create artificial time axis
                        time_data = range(len(df))

                    # Find columns that represent infection states (exclude time)
                    infection_cols = [col for col in df.columns
                                      if col.lower() not in ['time', 'day', 'date']
                                      and 'infect' in col.lower()]

                    if not infection_cols:
                        # If no specific infection columns, sum all numeric columns except time
                        numeric_cols = df.select_dtypes(
                            include=[np.number]).columns
                        infection_cols = [
                            col for col in numeric_cols if col != time_col]

                    if infection_cols:
                        # Calculate total infections (sum across age groups or infection states)
                        total_infections = df[infection_cols].sum(axis=1)

                        # Plot the main line
                        ax.plot(time_data, total_infections,
                                label=f"{sim_type}",
                                color=sim_colors[sim_type],
                                linewidth=2,
                                alpha=0.8)

                        # Add percentile bands if requested
                        if use_percentile and len(infection_cols) > 1:
                            # Calculate some spread measure (this is a simplified approach)
                            std_infections = df[infection_cols].std(axis=1)
                            ax.fill_between(time_data,
                                            total_infections - std_infections,
                                            total_infections + std_infections,
                                            color=sim_colors[sim_type],
                                            alpha=0.2)

            # Customize subplot
            ax.set_title(f"{event_type.replace('_', ' ').title()}",
                         fontsize=12, fontweight='bold')
            ax.set_xlabel("Time (days)")
            ax.set_ylabel("Total Infections")
            ax.grid(True, alpha=0.3)
            ax.legend()

            # Format x-axis if we have actual time data
            if hasattr(time_data, 'dtype') and 'datetime' in str(time_data.dtype):
                ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
                ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
                plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)

        # Hide extra subplots
        for j in range(len(event_types), len(axes_flat)):
            axes_flat[j].set_visible(False)

        # Adjust layout
        plt.tight_layout()

        # Save the comprehensive comparison
        output_file = self.output_dir / \
            f"comprehensive_comparison_{self.timestamp}.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Comprehensive comparison saved to: {output_file}")

        # Also save as PDF for publication quality
        pdf_file = self.output_dir / \
            f"comprehensive_comparison_{self.timestamp}.pdf"
        plt.savefig(pdf_file, dpi=300, bbox_inches='tight')
        print(f"PDF version saved to: {pdf_file}")

        plt.show()

    def create_summary_table(self):
        """Create a summary table with key metrics for each scenario."""
        summary_data = []

        for label, df in self.data.items():
            scenario = self.scenarios[label]

            # Calculate summary metrics
            try:
                # Find infection columns
                infection_cols = [col for col in df.columns
                                  if 'infect' in col.lower() and col.lower() != 'time']

                if infection_cols:
                    total_infections = df[infection_cols].sum(axis=1)

                    metrics = {
                        'Scenario': label,
                        'Simulation_Type': scenario['sim_type'],
                        'Event_Type': scenario['event_type'],
                        'Peak_Infections': total_infections.max(),
                        'Peak_Day': total_infections.argmax(),
                        'Final_Infections': total_infections.iloc[-1],
                        'Total_Infection_Days': (total_infections > 0).sum(),
                        'Growth_Rate': (total_infections.iloc[min(5, len(total_infections)-1)] /
                                        max(total_infections.iloc[0], 1)) if len(total_infections) > 5 else 1
                    }
                    summary_data.append(metrics)

            except Exception as e:
                print(f"Error calculating metrics for {label}: {e}")
                continue

        if summary_data:
            summary_df = pd.DataFrame(summary_data)

            # Save summary table
            summary_file = self.output_dir / \
                f"scenario_summary_{self.timestamp}.csv"
            summary_df.to_csv(summary_file, index=False)
            print(f"Summary table saved to: {summary_file}")

            # Display summary
            print("\n=== SCENARIO SUMMARY ===")
            print(summary_df.to_string(index=False))

            return summary_df

        return None

    def create_heatmap_comparison(self):
        """Create a heatmap comparing key metrics across scenarios."""
        summary_df = self.create_summary_table()

        if summary_df is not None and len(summary_df) > 1:
            # Pivot table for heatmap
            metrics_to_plot = ['Peak_Infections',
                               'Peak_Day', 'Final_Infections', 'Growth_Rate']

            # Create separate heatmaps for each metric
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            axes_flat = axes.flatten()

            for i, metric in enumerate(metrics_to_plot):
                if metric in summary_df.columns:
                    # Create pivot table
                    pivot_df = summary_df.pivot(index='Event_Type',
                                                columns='Simulation_Type',
                                                values=metric)

                    # Create heatmap
                    sns.heatmap(pivot_df,
                                annot=True,
                                fmt='.1f',
                                cmap='viridis',
                                ax=axes_flat[i])
                    axes_flat[i].set_xlabel("")
                    axes_flat[i].set_ylabel("")

            plt.tight_layout()

            # Save heatmap
            heatmap_file = self.output_dir / \
                f"metrics_heatmap_{self.timestamp}.png"
            plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
            print(f"Metrics heatmap saved to: {heatmap_file}")

            plt.show()


def main():
    """Main function to run the scenario comparison."""
    parser = argparse.ArgumentParser(
        description='Compare all simulation scenarios with comprehensive visualizations'
    )

    parser.add_argument(
        '--results-paths',
        nargs='+',
        required=True,
        help='List of paths to simulation results directories'
    )

    parser.add_argument(
        '--labels',
        nargs='+',
        required=True,
        help='List of labels for each simulation (format: simtype_eventtype)'
    )

    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for visualizations'
    )

    parser.add_argument(
        '--timestamp',
        required=True,
        help='Timestamp for file naming'
    )

    parser.add_argument(
        '--s90percentile',
        action='store_true',
        help='Show 90th percentile bands in plots'
    )

    args = parser.parse_args()

    # Validate inputs
    if len(args.results_paths) != len(args.labels):
        print("Error: Number of results paths must match number of labels")
        print(
            f"Got {len(args.results_paths)} paths and {len(args.labels)} labels")
        print("Paths:", args.results_paths)
        print("Labels:", args.labels)
        sys.exit(1)

    # Validate that result directories exist
    valid_paths = []
    valid_labels = []

    for path, label in zip(args.results_paths, args.labels):
        if os.path.exists(path):
            valid_paths.append(path)
            valid_labels.append(label)
            print(f"✓ Found results for {label}: {path}")
        else:
            print(
                f"⚠ Warning: Results path does not exist for {label}: {path}")

    if not valid_paths:
        print("Error: No valid results paths found")
        sys.exit(1)

    print(f"\nProcessing {len(valid_paths)} valid scenarios...")

    # Create comparator and run analyses
    try:
        comparator = ScenarioComparator(
            results_paths=valid_paths,
            labels=valid_labels,
            output_dir=args.output_dir,
            timestamp=args.timestamp
        )
    except Exception as e:
        print(f"Error creating comparator: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    if not comparator.data:
        print("Error: No data could be loaded from any scenarios")
        print("Available scenarios:", list(comparator.scenarios.keys()))
        sys.exit(1)

    print(f"Successfully loaded data for {len(comparator.data)} scenarios")
    print("Loaded scenarios:", list(comparator.data.keys()))

    try:
        # Create comprehensive comparison
        print("\n=== Creating Comprehensive Comparison ===")
        comparator.create_comprehensive_comparison(
            use_percentile=args.s90percentile)

        # Create summary table
        print("\n=== Creating Summary Table ===")
        summary_df = comparator.create_summary_table()

        # Create heatmap comparison (only if we have multiple scenarios)
        if summary_df is not None and len(summary_df) > 1:
            print("\n=== Creating Heatmap Comparison ===")
            comparator.create_heatmap_comparison()
        else:
            print("\n=== Skipping Heatmap (insufficient data) ===")

        print("\n=== Analysis Complete ===")
        print(f"All visualizations saved to: {args.output_dir}")

    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
