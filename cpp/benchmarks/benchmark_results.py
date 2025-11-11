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


class rawData:
    """Class to hold raw benchmark data."""

    def __init__(self):
        # Population sizes in thousands
        self.population_sizes = [1000, 2000, 4000,
                                 8000, 16000, 32000, 64000, 128000, 256000]  # in thousands
        self.population_sizes = [p * 1000 for p in self.population_sizes]

        # Runtime memilio single thread (normalized per time step)
        self.memilio_times_single_core = np.array(
            [33594, 68616, 135322, 271997, 550301, 1130269, 2237273, 5551993, 11600000]) * (1/120.0) * (1/1000)

        # Runtime memilio four threads (normalized per time step)
        self.memilio_times_four_cores = np.array(
            [13812, 27802, 56257, 113088, 223888, 454182, 928952, 1944267, 4877925]) * (1/120.0)*(1/1000)

        self.memilio_times_sixteen_cores = np.array(
            [6905, 13902, 28128, 56388, 112888, 226091, 456476, 911234, 4877925]) * (1/120.0)*(1/1000)

        # covasim single thread (normalized per time step)
        self.covasim = np.array(
            [62415, 128200, 292536, 632381, 1298022, 2625049, 5385827]) * (1/120.0)*(1/1000)

        # now data for weak scaling

        # Runtime with 250.000 agents per core
        self.memilio_weak_scaling_250k = np.array(
            [8305, 9997, 13853, 22790, 40123, 77431]) * (1/120.0)*(1/1000)

        # Runtime with 500.000 agents per core
        self.memilio_weak_scaling_500k = np.array(
            [16685, 20424, 27934, 45817, 80567, 155877]) * (1/120.0)*(1/1000)

        # Runtime with 1.000.000 agents per core
        self.memilio_weak_scaling_1m = np.array(
            [34165, 41361, 56680, 91576, 159270, 315894]) * (1/120.0)*(1/1000)

        # Runtime with 2.000.000 agents per core
        self.memilio_weak_scaling_2m = np.array(
            [69787, 83565, 112698, 181267, 322380, 601845]) * (1/120.0)*(1/1000)

        # now data for strong scaling

        self.strong_scaling_cores = [1, 2, 4, 8, 16, 32, 64, 128]
        self.strong_scaling_nodes = [1, 2, 4, 8, 16, 32, 64, 128]

        # Runtime Strong scaling
        self.memilio_strong_scaling_128_runs_one_node = np.array(
            [1, 30345, 16234, 8956, 5234, 3120])
        self.memilio_strong_scaling_128_runs_multiple_nodes = np.array(
            [1, 15876, 8423, 4650, 2723, 1500])


class BenchmarkAnalyzer:
    """Enhanced benchmark analysis and visualization tool."""

    def __init__(self, fontsize=18):
        self.fontsize = fontsize
        self.colors = {'memilio': '#1f77b4', 'covasim': '#ff7f0e'}

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

    def plot_runtime_scaling_comparison(self, raw_data, save_path=None):
        """Create runtime scaling plot comparing single core, four cores, and Covasim with speedup."""
        fig, ax1 = plt.subplots(figsize=(12, 9))

        # Plot MEmilio single core (solid line)
        ax1.plot(raw_data.population_sizes, raw_data.memilio_times_single_core,
                 marker='o', linewidth=3, markersize=8,
                 color='#1f77b4', label='MEmilio ABM (1 core)', linestyle='-')

        # Plot MEmilio four cores (dotted line)
        ax1.plot(raw_data.population_sizes, raw_data.covasim_times_four_cores,
                 marker='s', linewidth=3, markersize=8,
                 color='#2ca02c', label='MEmilio ABM (4 cores)', linestyle=':')

        # Plot Covasim (solid line) - note: covasim has fewer data points
        covasim_pop_sizes = raw_data.population_sizes[:len(
            raw_data.covasim_times_parallel)]
        ax1.plot(covasim_pop_sizes, raw_data.covasim_times_parallel,
                 marker='^', linewidth=3, markersize=8,
                 color='#ff7f0e', label='Covasim', linestyle='-')

        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel('Population Size', fontsize=self.fontsize)
        ax1.set_ylabel('Runtime per Time Step (seconds)',
                       fontsize=self.fontsize)
        ax1.set_title('Runtime Scaling Comparison',
                      fontsize=self.fontsize+4)

        ax1.tick_params(axis='both', which='major', labelsize=self.fontsize-2)
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=self.fontsize-2, loc='upper left')

        # Create secondary y-axis for speedup
        ax2 = ax1.twinx()

        # Calculate speedup for 1 core (Covasim / MEmilio single core)
        speedup_1core = []
        for i, pop in enumerate(covasim_pop_sizes):
            speedup = raw_data.covasim_times_parallel[i] / \
                raw_data.memilio_times_single_core[i]
            speedup_1core.append(speedup)

        # Calculate speedup for 4 cores (Covasim / MEmilio four cores)
        speedup_4cores = []
        for i, pop in enumerate(covasim_pop_sizes):
            speedup = raw_data.covasim_times_parallel[i] / \
                raw_data.covasim_times_four_cores[i]
            speedup_4cores.append(speedup)

        # Plot both speedup lines on secondary y-axis
        ax2.plot(covasim_pop_sizes, speedup_1core,
                 marker='D', linewidth=2.5, markersize=7,
                 color='#d62728', label='Speedup (1 core vs Covasim)',
                 linestyle='--', alpha=0.8)

        ax2.plot(covasim_pop_sizes, speedup_4cores,
                 marker='v', linewidth=2.5, markersize=7,
                 color='#9467bd', label='Speedup (4 cores vs Covasim)',
                 linestyle='-.', alpha=0.8)

        ax2.set_ylabel('Speedup Factor (Covasim / MEmilio)',
                       fontsize=self.fontsize, color='#8b0000')
        ax2.tick_params(axis='y', labelcolor='#8b0000',
                        labelsize=self.fontsize-2)
        ax2.set_xscale('log')

        # Add legend for speedup lines
        ax2.legend(fontsize=self.fontsize-2, loc='center right')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Runtime scaling plot saved to {save_path}")

        return fig, ax1

    def plot_weak_scaling(self, raw_data, save_path=None):
        """Create weak scaling plot for different agent counts per core."""
        fig, ax = plt.subplots(figsize=(12, 9))

        # Number of cores for weak scaling
        num_cores = [1, 2, 4, 8, 16, 32]

        # Plot each configuration
        ax.plot(num_cores, raw_data.memilio_weak_scaling_250k,
                marker='o', linewidth=3, markersize=8,
                label='250k agents per core', linestyle='-')

        ax.plot(num_cores, raw_data.memilio_weak_scaling_500k,
                marker='s', linewidth=3, markersize=8,
                label='500k agents per core', linestyle='-')

        ax.plot(num_cores, raw_data.memilio_weak_scaling_1m,
                marker='^', linewidth=3, markersize=8,
                label='1M agents per core', linestyle='-')

        ax.plot(num_cores, raw_data.memilio_weak_scaling_2m,
                marker='D', linewidth=3, markersize=8,
                label='2M agents per core', linestyle='-')

        # Add ideal scaling line (constant runtime)
        ideal_runtime = raw_data.memilio_weak_scaling_250k[0]
        ax.plot(num_cores, [ideal_runtime] * len(num_cores),
                'k--', linewidth=2, alpha=0.5, label='Ideal (constant runtime)')

        ax.set_xlabel('Number of Cores', fontsize=self.fontsize)
        ax.set_ylabel('Runtime per Time Step (seconds)',
                      fontsize=self.fontsize)
        ax.set_title('Weak Scaling: MEmilio ABM', fontsize=self.fontsize+4)

        ax.tick_params(axis='both', which='major', labelsize=self.fontsize-2)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=self.fontsize-2)

        # Set x-axis ticks to actual core counts
        ax.set_xticks(num_cores)
        ax.set_xticklabels([str(c) for c in num_cores])

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Weak scaling plot saved to {save_path}")

        return fig, ax

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

    def print_weak_scaling_efficiency(self, raw_data):
        """Print a table showing weak scaling efficiency."""
        num_cores = [1, 2, 4, 8, 16, 32]

        # Calculate efficiency as percentage (ideal runtime / actual runtime * 100)
        # Ideal runtime is the runtime with 1 core

        def calculate_efficiency(runtimes):
            ideal = runtimes[0]
            return [(ideal / runtime * 100) for runtime in runtimes]

        eff_250k = calculate_efficiency(raw_data.memilio_weak_scaling_250k)
        eff_500k = calculate_efficiency(raw_data.memilio_weak_scaling_500k)
        eff_1m = calculate_efficiency(raw_data.memilio_weak_scaling_1m)
        eff_2m = calculate_efficiency(raw_data.memilio_weak_scaling_2m)

        print("\n" + "="*95)
        print("WEAK SCALING EFFICIENCY TABLE: MEmilio ABM")
        print("="*95)
        print(f"{'Cores':>8} | {'250k agents/core':>20} | {'500k agents/core':>20} | {'1M agents/core':>20} | {'2M agents/core':>20}")
        print(f"{'':>8} | {'Efficiency (%)':>20} | {'Efficiency (%)':>20} | {'Efficiency (%)':>20} | {'Efficiency (%)':>20}")
        print("-"*95)

        for i, cores in enumerate(num_cores):
            print(
                f"{cores:>8} | {eff_250k[i]:>19.1f}% | {eff_500k[i]:>19.1f}% | {eff_1m[i]:>19.1f}% | {eff_2m[i]:>19.1f}%")

        print("="*95)
        print("\nNote: Efficiency = (Runtime with 1 core / Runtime with N cores) × 100%")
        print(
            "      100% efficiency means perfect weak scaling (constant runtime per core)")
        print()


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

    # Create analyzer and raw data
    analyzer = BenchmarkAnalyzer(fontsize=18)
    raw_data = rawData()

    # Create plots
    save_dir = Path(args.save) if args.save else None
    if save_dir:
        save_dir.mkdir(exist_ok=True)

    # Runtime scaling plot (single core, four cores, Covasim)
    runtime_scaling_path = save_dir / 'runtime_scaling.png' if save_dir else None
    fig1, ax1 = analyzer.plot_runtime_scaling_comparison(
        raw_data, runtime_scaling_path)

    # Weak scaling plot
    weak_scaling_path = save_dir / 'weak_scaling.png' if save_dir else None
    fig2, ax2 = analyzer.plot_weak_scaling(raw_data, weak_scaling_path)

    # Print weak scaling efficiency table
    analyzer.print_weak_scaling_efficiency(raw_data)

    # Save data if requested
    if args.data_file:
        analyzer.save_data_to_file(args.data_file)

    # Show plots unless disabled
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
