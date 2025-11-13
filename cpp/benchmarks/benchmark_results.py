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
            [32947, 67147, 137760, 274461, 549385, 1119135, 2193598, 7112954, 17757387]) * (1/120.0) * (1/1000)

        # Runtime memilio four threads (normalized per time step)
        self.memilio_times_four_cores = np.array(
            [13812, 27802, 56257, 113088, 223888, 454182, 928952, 1944267, 4877925]) * (1/120.0)*(1/1000)

        self.memilio_times_sixteen_cores = np.array(
            [11571, 23337, 46865, 92908, 184811, 369252, 752734, 1609745, 3396075]) * (1/120.0)*(1/1000)

        # covasim single thread (normalized per time step)
        self.covasim = np.array(
            [37809, 75152, 198949, 409061, 863334, 1762960, 3618340]) * (1/120.0)*(1/1000)

        # opencovid single thread (normalized per time step)
        self.opencovid = np.array(
            [514, 861, 1584, 2959, 5918, 12816, 24275]) * (1/120.0)

        # now data for weak scaling
        self.weak_scaling_cores = [1, 2, 4, 8, 16, 32]

        # Runtime with 250.000 agents per core
        self.memilio_weak_scaling_250k = np.array(
            [8406, 10466, 15257, 25711, 46061, 97237]) * (1/120.0)*(1/1000)

        # Runtime with 500.000 agents per core
        self.memilio_weak_scaling_500k = np.array(
            [16674, 21146, 30163, 52217, 91833, 239517]) * (1/120.0)*(1/1000)

        # Runtime with 1.000.000 agents per core
        self.memilio_weak_scaling_1m = np.array(
            [33568, 42239, 61171, 103338, 184308, 484596]) * (1/120.0)*(1/1000)

        # Runtime with 2.000.000 agents per core
        self.memilio_weak_scaling_2m = np.array(
            [66931, 85021, 124221, 205287, 367955, 764577]) * (1/120.0)*(1/1000)

        # now data for strong scaling

        self.strong_scaling_cores = [1, 2, 4, 8, 16, 32, 64, 128]
        self.strong_scaling_nodes = [1, 2, 4, 8, 16, 32, 64, 128]

        # Runtime Strong scaling
        self.memilio_strong_scaling_128_runs_one_node = np.array(
            [2.646570e+04, 2.068018e+04,  1.078597e+04, 5.317400e+03, 2.740988e+03, 1.440401e+03, 8.953795e+02, 5.889476e+02])
        self.memilio_strong_scaling_128_runs_multiple_nodes = np.array(
            # only last data point available
            [6.66666e+04, 3.803566e+04, 1.878124e+04, 9.376841e+03, 4.762327e+03, 2.363056e+03, 1.186088e+03, 6.005680e+02])


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

    def plot_agent_scaling(self, raw_data, save_path=None):
        """Plot 1: Scaling with agents - memilio (1, 4, 16 cores), covasim, and opencovid."""
        fig, ax = plt.subplots(figsize=(12, 9))

        # Define colors for memilio lines (shades of blue)
        memilio_color_1 = '#0d47a1'   # dark blue
        memilio_color_4 = '#1976d2'   # medium blue
        memilio_color_16 = '#64b5f6'  # light blue
        covasim_color = '#ff7f0e'      # orange
        opencovid_color = '#2ca02c'    # green

        # Plot MEmilio with different core counts (all same color family)
        ax.plot(raw_data.population_sizes, raw_data.memilio_times_single_core,
                marker='o', linewidth=3, markersize=8,
                color=memilio_color_1, label='MEmilio ABM (1 core)', linestyle='-')

        ax.plot(raw_data.population_sizes, raw_data.memilio_times_four_cores,
                marker='o', linewidth=3, markersize=8,
                color=memilio_color_4, label='MEmilio ABM (4 cores)', linestyle='-')

        ax.plot(raw_data.population_sizes, raw_data.memilio_times_sixteen_cores,
                marker='o', linewidth=3, markersize=8,
                color=memilio_color_16, label='MEmilio ABM (16 cores)', linestyle='-')

        # Plot Covasim (fewer data points)
        covasim_pop_sizes = raw_data.population_sizes[:len(raw_data.covasim)]
        ax.plot(covasim_pop_sizes, raw_data.covasim,
                marker='s', linewidth=3, markersize=8,
                color=covasim_color, label='Covasim', linestyle='-')

        # Plot OpenCOVID (fewer data points)
        opencovid_pop_sizes = raw_data.population_sizes[:len(
            raw_data.opencovid)]
        ax.plot(opencovid_pop_sizes, raw_data.opencovid,
                marker='^', linewidth=3, markersize=8,
                color=opencovid_color, label='OpenCOVID', linestyle='-')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Number of Agents', fontsize=self.fontsize)
        ax.set_ylabel('Runtime per Time Step (seconds)',
                      fontsize=self.fontsize)
        ax.set_title('Agent-Based Model Runtime Scaling',
                     fontsize=self.fontsize+4)

        ax.tick_params(axis='both', which='major', labelsize=self.fontsize-2)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=self.fontsize-2, loc='upper left')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Agent scaling plot saved to {save_path}")

        return fig, ax

    def plot_one_node_strong_scaling(self, raw_data, save_path=None):
        """Plot 3: One node strong scaling."""
        fig, ax = plt.subplots(figsize=(12, 9))

        # Use only the data points that exist (skip the first element which is 1)
        cores = raw_data.strong_scaling_cores[0:len(
            raw_data.memilio_strong_scaling_128_runs_one_node)]
        runtimes = raw_data.memilio_strong_scaling_128_runs_one_node[0:]

        # Plot runtime
        ax.plot(cores, runtimes,
                marker='o', linewidth=3, markersize=8,
                color='#1f77b4', label='MEmilio ABM (One Node)', linestyle='-')

        # Calculate and plot ideal scaling
        if len(runtimes) > 0:
            ideal_scaling = [runtimes[0] / c * cores[0] for c in cores]
            ax.plot(cores, ideal_scaling,
                    'k--', linewidth=2, alpha=0.5, label='Ideal Scaling')

        ax.set_xlabel('Number of Cores', fontsize=self.fontsize)
        ax.set_ylabel('Runtime per Time Step (seconds)',
                      fontsize=self.fontsize)
        ax.set_title('Strong Scaling: One Node - MEmilio ABM',
                     fontsize=self.fontsize+4)

        ax.tick_params(axis='both', which='major', labelsize=self.fontsize-2)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=self.fontsize-2)

        # Set x-axis to log scale if useful, or linear
        ax.set_xscale('log', base=2)
        ax.set_yscale('log')
        ax.set_xticks(cores)
        ax.set_xticklabels([str(c) for c in cores])

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"One node strong scaling plot saved to {save_path}")

        return fig, ax

    def plot_multi_node_strong_scaling(self, raw_data, save_path=None):
        """Plot 4: Multi-node strong scaling."""
        fig, ax = plt.subplots(figsize=(12, 9))

        # Use only the data points that exist (skip the first element which is 1)
        nodes = raw_data.strong_scaling_nodes[0:len(
            raw_data.memilio_strong_scaling_128_runs_multiple_nodes)]
        runtimes = raw_data.memilio_strong_scaling_128_runs_multiple_nodes[0:]

        # Plot runtime
        ax.plot(nodes, runtimes,
                marker='s', linewidth=3, markersize=8,
                color='#2ca02c', label='MEmilio ABM (Multi-Node)', linestyle='-')

        # Calculate and plot ideal scaling
        if len(runtimes) > 0:
            ideal_scaling = [runtimes[0] / n * nodes[0] for n in nodes]
            ax.plot(nodes, ideal_scaling,
                    'k--', linewidth=2, alpha=0.5, label='Ideal Scaling')

        ax.set_xlabel('Number of Nodes', fontsize=self.fontsize)
        ax.set_ylabel('Runtime per Time Step (seconds)',
                      fontsize=self.fontsize)
        ax.set_title('Strong Scaling: Multi-Node - MEmilio ABM',
                     fontsize=self.fontsize+4)

        ax.tick_params(axis='both', which='major', labelsize=self.fontsize-2)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=self.fontsize-2)

        # Set x-axis to log scale if useful, or linear
        ax.set_xscale('log', base=2)
        ax.set_yscale('log')
        ax.set_xticks(nodes)
        ax.set_xticklabels([str(n) for n in nodes])

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Multi-node strong scaling plot saved to {save_path}")

        return fig, ax

    def plot_weak_scaling(self, raw_data, save_path=None):
        """Plot 2: Weak scaling plot for different agent counts per core."""
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

        ax.set_yscale('log')

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

    def print_strong_scaling_efficiency(self, raw_data):
        """Print a table showing strong scaling efficiency."""
        # One node strong scaling
        cores_one_node = raw_data.strong_scaling_cores[1:len(
            raw_data.memilio_strong_scaling_128_runs_one_node)]
        runtimes_one_node = raw_data.memilio_strong_scaling_128_runs_one_node[1:]

        # Multi-node strong scaling
        nodes_multi = raw_data.strong_scaling_nodes[1:len(
            raw_data.memilio_strong_scaling_128_runs_multiple_nodes)]
        runtimes_multi = raw_data.memilio_strong_scaling_128_runs_multiple_nodes[1:]

        def calculate_strong_efficiency(baseline, runtimes, units):
            """Calculate strong scaling efficiency: (baseline / (runtime * units)) * 100"""
            return [(baseline / (runtime * unit) * 100) for runtime, unit in zip(runtimes, units)]

        print("\n" + "="*70)
        print("STRONG SCALING EFFICIENCY TABLE: MEmilio ABM")
        print("="*70)

        # One node table
        if len(runtimes_one_node) > 0:
            print("\nOne Node Strong Scaling:")
            print(
                f"{'Cores':>10} | {'Runtime (s)':>15} | {'Speedup':>10} | {'Efficiency (%)':>15}")
            print("-"*70)
            baseline_one = runtimes_one_node[0] * cores_one_node[0]
            for i, cores in enumerate(cores_one_node):
                speedup = baseline_one / runtimes_one_node[i]
                efficiency = (speedup / cores) * 100
                print(
                    f"{cores:>10} | {runtimes_one_node[i]:>15.2f} | {speedup:>10.2f} | {efficiency:>15.1f}%")

        # Multi-node table
        if len(runtimes_multi) > 0:
            print("\nMulti-Node Strong Scaling:")
            print(
                f"{'Nodes':>10} | {'Runtime (s)':>15} | {'Speedup':>10} | {'Efficiency (%)':>15}")
            print("-"*70)
            baseline_multi = runtimes_multi[0] * nodes_multi[0]
            for i, nodes in enumerate(nodes_multi):
                speedup = baseline_multi / runtimes_multi[i]
                efficiency = (speedup / nodes) * 100
                print(
                    f"{nodes:>10} | {runtimes_multi[i]:>15.2f} | {speedup:>10.2f} | {efficiency:>15.1f}%")

        print("="*70)
        print("\nNote: Efficiency = (Speedup / Number of cores or nodes) × 100%")
        print("      100% efficiency means perfect linear scaling")
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

    args.save = "/Users/saschakorf/Nosynch/Arbeit/memilio/example_results"

    # Create analyzer and raw data
    analyzer = BenchmarkAnalyzer(fontsize=18)
    raw_data = rawData()

    # Create plots
    save_dir = Path(args.save) if args.save else None
    if save_dir:
        save_dir.mkdir(exist_ok=True)

    # Plot 1: Agent scaling with memilio (1, 4, 16 cores), covasim, and opencovid
    agent_scaling_path = save_dir / 'agent_scaling.png' if save_dir else None
    fig1, ax1 = analyzer.plot_agent_scaling(raw_data, agent_scaling_path)

    # Plot 2: Weak scaling
    weak_scaling_path = save_dir / 'weak_scaling.png' if save_dir else None
    fig2, ax2 = analyzer.plot_weak_scaling(raw_data, weak_scaling_path)

    # Plot 3: One node strong scaling
    one_node_strong_path = save_dir / \
        'strong_scaling_one_node.png' if save_dir else None
    fig3, ax3 = analyzer.plot_one_node_strong_scaling(
        raw_data, one_node_strong_path)

    # Plot 4: Multi-node strong scaling
    multi_node_strong_path = save_dir / \
        'strong_scaling_multi_node.png' if save_dir else None
    fig4, ax4 = analyzer.plot_multi_node_strong_scaling(
        raw_data, multi_node_strong_path)

    # Print efficiency tables
    analyzer.print_weak_scaling_efficiency(raw_data)
    analyzer.print_strong_scaling_efficiency(raw_data)

    # Save data if requested
    if args.data_file:
        analyzer.save_data_to_file(args.data_file)

    # Show plots unless disabled
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
