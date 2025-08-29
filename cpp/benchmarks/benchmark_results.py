import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogLocator, FuncFormatter
import json
import argparse
from pathlib import Path
from scipy import stats
import os

# Optional imports for enhanced functionality
try:
    import seaborn as sns
    sns.set_palette("husl")
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

# Set style for better-looking plots
if HAS_SEABORN:
    plt.style.use('seaborn-v0_8-whitegrid')
else:
    plt.style.use('default')
    plt.rcParams['grid.alpha'] = 0.3
    plt.rcParams['figure.facecolor'] = 'white'


class BenchmarkAnalyzer:
    """Enhanced benchmark analysis and visualization tool."""

    def __init__(self, fontsize=18):
        self.fontsize = fontsize
        self.colors = {'memilio': '#1f77b4', 'covasim': '#ff7f0e'}

        # Default data - can be overridden by loading from files
        self.population_sizes = [25, 50, 100,
                                 200, 400, 800, 1600]  # in thousands
        self.population_sizes = [p * 1000 for p in self.population_sizes]

        # Runtime data parallel (normalized per time step)
        self.memilio_times = np.array(
            [0.085, 0.16, 0.37, 0.75, 1.6, 3.1, 6.2]) * (1/120.0)
        self.covasim_times = np.array(
            [0.424, 0.757, 1.6, 3.4, 7.2, 16.2, 33.6]) * (1/120.0)

        # Runtime data single thread (normalized per time step)
        self.memilio_times = np.array(
            [0.1, 0.2, 0.37, 0.75, 1.6, 3.1, 6.2]) * (1/120.0)
        self.covasim_times = np.array(
            [0.424, 0.757, 1.6, 3.4, 7.2, 16.2, 33.6]) * (1/120.0)

    def load_data_from_file(self, filename):
        """Load benchmark data from JSON file."""
        if Path(filename).exists():
            with open(filename, 'r') as f:
                data = json.load(f)
                self.population_sizes = data.get(
                    'population_sizes', self.population_sizes)
                self.memilio_times = np.array(
                    data.get('memilio_times', self.memilio_times))
                self.covasim_times = np.array(
                    data.get('covasim_times', self.covasim_times))
                print(f"Loaded data from {filename}")
        else:
            print(f"File {filename} not found, using default data")

    def save_data_to_file(self, filename):
        """Save current benchmark data to JSON file."""
        data = {
            'population_sizes': self.population_sizes,
            'memilio_times': self.memilio_times.tolist(),
            'covasim_times': self.covasim_times.tolist()
        }
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"Data saved to {filename}")

    def calculate_metrics(self):
        """Calculate performance metrics."""
        # Time per 1000 agents
        memilio_per_1k = [1000 * t / p for t,
                          p in zip(self.memilio_times, self.population_sizes)]
        covasim_per_1k = [1000 * t / p for t,
                          p in zip(self.covasim_times, self.population_sizes)]

        # Speedup factor (Covasim time / MEmilio time)
        speedup = self.covasim_times / self.memilio_times

        # Scaling exponents (fit power law: time = a * pop_size^b)
        log_pop = np.log10(self.population_sizes)
        memilio_slope, memilio_intercept, memilio_r, _, _ = stats.linregress(
            log_pop, np.log10(self.memilio_times))
        covasim_slope, covasim_intercept, covasim_r, _, _ = stats.linregress(
            log_pop, np.log10(self.covasim_times))

        return {
            'memilio_per_1k': memilio_per_1k,
            'covasim_per_1k': covasim_per_1k,
            'speedup': speedup,
            'memilio_scaling': {'exponent': memilio_slope, 'r_squared': memilio_r**2},
            'covasim_scaling': {'exponent': covasim_slope, 'r_squared': covasim_r**2}
        }

    def plot_runtime_scaling(self, save_path=None):
        """Create the main runtime scaling plot."""
        fig, ax = plt.subplots(figsize=(12, 9))

        # Plot data with confidence intervals (assuming 10% uncertainty)
        memilio_err = self.memilio_times * 0.1
        covasim_err = self.covasim_times * 0.1

        ax.errorbar(self.population_sizes, self.memilio_times, yerr=memilio_err,
                    marker='o', linewidth=3, markersize=8, capsize=5,
                    color=self.colors['memilio'], label='MEmilio ABM')
        ax.errorbar(self.population_sizes, self.covasim_times, yerr=covasim_err,
                    marker='s', linewidth=3, markersize=8, capsize=5,
                    color=self.colors['covasim'], label='Covasim')

        # Add scaling lines
        metrics = self.calculate_metrics()
        log_pop_fine = np.linspace(np.log10(min(self.population_sizes)),
                                   np.log10(max(self.population_sizes)), 100)
        pop_fine = 10**log_pop_fine

        memilio_fit = 10**(metrics['memilio_scaling']['exponent'] * log_pop_fine +
                           np.log10(self.memilio_times[0]) -
                           metrics['memilio_scaling']['exponent'] * np.log10(self.population_sizes[0]))
        covasim_fit = 10**(metrics['covasim_scaling']['exponent'] * log_pop_fine +
                           np.log10(self.covasim_times[0]) -
                           metrics['covasim_scaling']['exponent'] * np.log10(self.population_sizes[0]))

        ax.plot(pop_fine, memilio_fit, '--', alpha=0.7,
                color=self.colors['memilio'])
        ax.plot(pop_fine, covasim_fit, '--', alpha=0.7,
                color=self.colors['covasim'])

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Population Size', fontsize=self.fontsize)
        ax.set_ylabel('Runtime per Time Step (seconds)',
                      fontsize=self.fontsize)
        ax.set_title('Epidemic Simulation Runtime Scaling Comparison',
                     fontsize=self.fontsize+4)

        # Format axes
        ax.tick_params(axis='both', which='major', labelsize=self.fontsize-2)
        ax.grid(True, alpha=0.3)

        # Custom legend with scaling information
        legend_text = [
            f'MEmilio ABM (scaling ∝ n^{metrics["memilio_scaling"]["exponent"]:.2f})',
            f'Covasim (scaling ∝ n^{metrics["covasim_scaling"]["exponent"]:.2f})'
        ]
        ax.legend(legend_text, fontsize=self.fontsize-2)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {save_path}")

        return fig, ax

    def plot_speedup_analysis(self, save_path=None):
        """Create speedup analysis plot."""
        metrics = self.calculate_metrics()

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

        # Top plot: Speedup factor
        ax1.semilogx(self.population_sizes, metrics['speedup'],
                     marker='o', linewidth=3, markersize=8, color='green')
        ax1.axhline(y=np.mean(metrics['speedup']), color='red', linestyle='--',
                    label=f'Average speedup: {np.mean(metrics["speedup"]):.1f}x')
        ax1.set_xlabel('Population Size', fontsize=self.fontsize)
        ax1.set_ylabel(
            'Speedup Factor\n(Covasim time / MEmilio time)', fontsize=self.fontsize)
        ax1.set_title('MEmilio Performance Advantage',
                      fontsize=self.fontsize+2)
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=self.fontsize-2)
        ax1.tick_params(axis='both', which='major', labelsize=self.fontsize-2)

        # Bottom plot: Time per 1000 agents
        ax2.semilogx(self.population_sizes, np.array(metrics['memilio_per_1k'])*1000,
                     marker='o', linewidth=3, markersize=8,
                     color=self.colors['memilio'], label='MEmilio ABM')
        ax2.semilogx(self.population_sizes, np.array(metrics['covasim_per_1k'])*1000,
                     marker='s', linewidth=3, markersize=8,
                     color=self.colors['covasim'], label='Covasim')
        ax2.set_xlabel('Population Size', fontsize=self.fontsize)
        ax2.set_ylabel(
            'Runtime per 1000 Agents\nper Time Step (milliseconds)', fontsize=self.fontsize)
        ax2.set_title('Computational Efficiency Comparison',
                      fontsize=self.fontsize+2)
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=self.fontsize-2)
        ax2.tick_params(axis='both', which='major', labelsize=self.fontsize-2)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Speedup analysis saved to {save_path}")

        return fig, (ax1, ax2)

    def print_summary(self):
        """Print summary statistics."""
        metrics = self.calculate_metrics()

        print("\n" + "="*60)
        print("BENCHMARK ANALYSIS SUMMARY")
        print("="*60)

        print(
            f"\nPopulation sizes tested: {[f'{p//1000}k' for p in self.population_sizes]}")

        print(f"\nScaling Analysis:")
        print(
            f"  MEmilio scaling exponent: {metrics['memilio_scaling']['exponent']:.3f} (R² = {metrics['memilio_scaling']['r_squared']:.3f})")
        print(
            f"  Covasim scaling exponent: {metrics['covasim_scaling']['exponent']:.3f} (R² = {metrics['covasim_scaling']['r_squared']:.3f})")

        print(f"\nPerformance Summary:")
        print(
            f"  Average speedup (MEmilio vs Covasim): {np.mean(metrics['speedup']):.1f}x faster")
        print(
            f"  MEmilio runtime range: {min(self.memilio_times)*1000:.1f} - {max(self.memilio_times)*1000:.1f} ms/timestep")
        print(
            f"  Covasim runtime range: {min(self.covasim_times)*1000:.1f} - {max(self.covasim_times)*1000:.1f} ms/timestep")

        print(f"\nEfficiency (ms per 1000 agents per timestep):")
        for i, pop in enumerate(self.population_sizes):
            memilio_eff = metrics['memilio_per_1k'][i] * 1000
            covasim_eff = metrics['covasim_per_1k'][i] * 1000
            print(
                f"  {pop//1000:4}k agents: MEmilio {memilio_eff:.2f}ms, Covasim {covasim_eff:.2f}ms")

        print("="*60)


def main():
    parser = argparse.ArgumentParser(
        description='Enhanced benchmark analysis tool')
    parser.add_argument('--save', type=str,
                        help='Save plots to specified directory')
    parser.add_argument('--data-file', type=str,
                        help='Save/load data to/from this JSON file')
    parser.add_argument('--no-show', action='store_true',
                        help='Don\'t display plots')

    args = parser.parse_args()

    # Create analyzer
    analyzer = BenchmarkAnalyzer(fontsize=18)

    # Create plots
    save_dir = Path(args.save) if args.save else None
    if save_dir:
        save_dir.mkdir(exist_ok=True)

    # Main scaling plot
    runtime_path = save_dir / 'runtime_scaling.png' if save_dir else None
    fig1, ax1 = analyzer.plot_runtime_scaling(runtime_path)

    # Speedup analysis
    speedup_path = save_dir / 'speedup_analysis.png' if save_dir else None
    fig2, (ax2, ax3) = analyzer.plot_speedup_analysis(speedup_path)

    # Print summary
    analyzer.print_summary()

    # Save data if requested
    if args.data_file:
        analyzer.save_data_to_file(args.data_file)

    # Show plots unless disabled
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
