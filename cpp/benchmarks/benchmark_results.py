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

# Set professional plot style
plt.style.use('default')
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['figure.facecolor'] = 'white'

# DPI for high-quality output
DPI = 300

# Base fontsize settings
BASE_FONTSIZE = 17
TICK_FONTSIZE = int(0.8 * BASE_FONTSIZE)
LEGEND_FONTSIZE = int(0.8 * BASE_FONTSIZE)


def set_fontsize(base_fontsize=17):
    """Set font sizes for all plot elements."""
    fontsize = base_fontsize
    plt.rcParams.update({
        'font.size': fontsize,
        'axes.titlesize': fontsize * 1,
        'axes.labelsize': fontsize ,
        'xtick.labelsize': fontsize * 0.8,
        'ytick.labelsize': fontsize * 0.8,
        'legend.fontsize': fontsize * 0.8,
        'font.family': "Arial"
    })


class rawData:
    """Class to hold raw benchmark data."""

    def __init__(self):
        # Population sizes in thousands
        self.population_sizes = [1000, 2000, 4000,
                                 8000, 16000, 32000, 64000, 128000, 256000]  # in thousands
        self.population_sizes = [p * 1000 for p in self.population_sizes]

        # Runtime memilio single thread (normalized per time step)
        self.memilio_times_single_core = np.array(
            [31023, 62863, 124198, 248658, 498499, 1003958, 1992911, 6342611, 15568413]) * (1/120.0) * (1/1000)

        # Runtime memilio four threads (normalized per time step)
        self.memilio_times_four_cores = np.array(
            [14744, 30224, 60749, 119645, 240825, 496879, 1031504, 2561311, 5752128]) * (1/120.0)*(1/1000)

        self.memilio_times_sixteen_cores = np.array(
            [11818, 23459, 47549, 94222, 187279, 378181, 776280, 1656997, 3477961]) * (1/120.0)*(1/1000)

        # covasim single thread (normalized per time step)
        self.covasim = np.array(
            [37809, 75152, 198949, 409061, 863334, 1762960, 3618340]) * (1/120.0)*(1/1000)

        # opencovid single thread (normalized per time step) (in seconds!!!!)
        self.opencovid = np.array(
            [380, 622, 1069, 2009, 3947, 8569, 16142]) * (1/120.0)

        # now data for weak scaling
        self.weak_scaling_cores = [1, 2, 4, 8, 16, 32]

        # Runtime with 250.000 agents per core
        self.memilio_weak_scaling_250k = np.array(
            [7625, 10058, 14891, 25895, 46850, 124340]) * (1/120.0)*(1/1000)

        # Runtime with 500.000 agents per core
        self.memilio_weak_scaling_500k = np.array(
            [15343, 20838, 30058, 52553, 93834, 237016]) * (1/120.0)*(1/1000)

        # Runtime with 1.000.000 agents per core
        self.memilio_weak_scaling_1m = np.array(
            [31135, 41238, 60003, 103803, 188685, 483887]) * (1/120.0)*(1/1000)

        # Runtime with 2.000.000 agents per core
        self.memilio_weak_scaling_2m = np.array(
            [61156, 82754, 119354, 208809, 377700, 891970]) * (1/120.0)*(1/1000)

        # now data for strong scaling

        self.strong_scaling_cores = [1, 2, 4, 8, 16, 32, 64, 128]
        self.strong_scaling_nodes = [1, 2, 4, 8, 16, 32, 64, 128]
        # multiply each by 128 for total cores
        self.strong_scaling_nodes = [
            n * 128 for n in self.strong_scaling_nodes]

        # Runtime Strong scaling
        self.memilio_strong_scaling_128_runs_one_node = np.array(
            [2.512869e+04, 1.911620e+04,  9.766570e+03, 4.919744e+03,  2.492793e+03, 1.299748e+03,  8.281876e+02,  5.801488e+02])
        self.memilio_strong_scaling_128_runs_multiple_nodes = np.array(
            # only last data point available
            [7.413365e+04, 3.721441e+04, 1.865978e+04, 9.374108e+03, 4.684981e+03, 2.357355e+03, 1.187477e+03, 6.005680e+02])
        

        # "Graph-ODE w/ mobility: [19566.8, 9718.44, 4847.9, 2417.92, 1252.95, 615.575, 329.08, 190.454],
        # "LCT": [746.859, 373.144, 185.975, 92.8182, 46.6567, 24.8617, 12.8418, 7.05302],
        # "IDE": [5358.52, 2672.95, 1341, 666.594, 332.608, 176.283, 91.0512, 46.2337]
        self.memilio_graph_ode_strong_scaling_128_runs_one_node = np.array(
            [19566.8, 9718.44, 4847.9, 2417.92, 1252.95, 615.575, 329.08, 190.454])
        self.memilio_lct_strong_scaling_128_runs_one_node = np.array(
            [746.859, 373.144, 185.975, 92.8182, 46.6567, 24.8617, 12.8418, 7.05302])
        self.memilio_ide_strong_scaling_128_runs_one_node = np.array(
            [5358.52, 2672.95, 1341, 666.594, 332.608, 176.283, 91.0512, 46.2337])
        


class BenchmarkAnalyzer:
    """Enhanced benchmark analysis and visualization tool."""

    def __init__(self, fontsize=None):
        self.fontsize = fontsize if fontsize is not None else BASE_FONTSIZE
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

    def _annotate_speedups_minmax(self, ax, x_data, y_data_fast, y_data_slow, color, fontsize=None):
        """Draw triangle markers and speedup labels for min and max speedup only."""
        if fontsize is None:
            fontsize = max(int(TICK_FONTSIZE * 0.85), 6)

        try:
            import matplotlib.patheffects as mpatheffects
        except Exception:
            mpatheffects = None

        # Get axis limits for label positioning
        try:
            x_limits = ax.get_xlim()
            is_log_x = ax.get_xscale() == 'log'
        except Exception:
            x_limits = None
            is_log_x = False

        def _shift_label_x(x_value):
            if x_limits is None:
                return x_value
            xmin, xmax = x_limits
            if not (np.isfinite(xmin) and np.isfinite(xmax)):
                return x_value
            if is_log_x and x_value > 0:
                candidate = x_value * 1.4  # Moderate shift to the right
                max_allowed = xmax / 1.01 if xmax > 0 else xmax
                if x_value >= max_allowed:
                    return x_value
                return candidate if candidate < max_allowed else max_allowed
            span = xmax - xmin
            if span <= 0:
                return x_value
            # Moderate shift to the right
            candidate = x_value + 0.10 * span
            max_allowed = xmax - 0.01 * span
            if x_value >= max_allowed:
                return x_value
            return candidate if candidate < max_allowed else max_allowed

        # First pass: calculate all speedups to find min and max
        min_len = min(len(x_data), len(y_data_fast), len(y_data_slow))
        speedups = []
        valid_indices = []

        for i in range(min_len):
            y_fast = y_data_fast[i]
            y_slow = y_data_slow[i]

            if np.isfinite(y_fast) and np.isfinite(y_slow) and y_fast > 0 and y_slow > 0:
                speedup = y_slow / y_fast
                if np.isfinite(speedup):
                    speedups.append(speedup)
                    valid_indices.append(i)

        if not speedups:
            return

        # Find indices of min and max speedup
        min_speedup_idx = valid_indices[np.argmin(speedups)]
        max_speedup_idx = valid_indices[np.argmax(speedups)]

        # Draw annotations only for min and max
        for i in [min_speedup_idx, max_speedup_idx]:
            x = x_data[i]
            y_fast = y_data_fast[i]
            y_slow = y_data_slow[i]

            lower, upper = (y_fast, y_slow) if y_fast < y_slow else (
                y_slow, y_fast)

            # Draw triangle markers along the vertical line
            n_points = 12
            y_points = np.logspace(np.log10(lower), np.log10(upper), n_points)
            x_points = np.full(n_points, x)

            ax.plot(
                x_points,
                y_points,
                color=color,
                linewidth=0,
                linestyle='None',
                marker='^',
                markersize=4,
                markerfacecolor=color,
                markeredgecolor=color,
                alpha=0.45,
                zorder=1.5,
            )

            # Calculate and display speedup
            speedup = y_slow / y_fast

            y_mid = np.sqrt(y_fast * y_slow)
            # Position slightly above center
            shift_factor = 1.4
            offset_upper = upper * 0.99
            offset_lower = lower * 1.02
            if y_fast < y_slow:
                # Position slightly above center
                y_pos = min(y_mid * shift_factor, offset_upper)
            else:
                y_pos = max(y_mid * shift_factor, offset_lower)

            label = f"{speedup:.0f}×"
            # Increase fontsize - was using TICK_FONTSIZE * 0.85, now use BASE_FONTSIZE * 0.9
            label_fontsize = int(BASE_FONTSIZE * 0.9)
            text_obj = ax.text(
                _shift_label_x(x),
                y_pos,
                label,
                ha="center",
                va="center",
                color=color,
                fontsize=label_fontsize,
                zorder=4,
                clip_on=False,
            )
            if mpatheffects is not None:
                text_obj.set_path_effects(
                    [mpatheffects.withStroke(linewidth=2.0, foreground="white")])

    def plot_agent_scaling(self, raw_data, save_path=None):
        """Plot 1: Scaling with agents - memilio (1, 4, 16 cores), covasim, and opencovid."""
        set_fontsize()
        figsize = (8, 5)
        panel = (0.2, 0.2, 0.78, 0.75)
        fig = plt.figure(figsize=figsize, dpi=DPI)
        ax = fig.add_axes(panel)

        # Define colors for memilio lines (shades of blue)
        memilio_color_1 = '#741194'   # dark blue
        memilio_color_4 = '#741194'   # medium blue
        memilio_color_16 = '#741194'  # light blue
        covasim_color = '#5D8A2B'      # orange
        opencovid_color = '#E89A63'    # green

        # Plot MEmilio with different core counts (all same color family)
        ax.plot(raw_data.population_sizes, raw_data.memilio_times_single_core,
                marker='o',  
                color=memilio_color_1, label='MEmilio ABM (1 core)', linestyle='-')

        ax.plot(raw_data.population_sizes, raw_data.memilio_times_four_cores,
                marker='o',  
                color=memilio_color_4, label='MEmilio ABM (4 cores)', linestyle='--')
        ax.plot(raw_data.population_sizes, raw_data.memilio_times_sixteen_cores,
                marker='o',  
                color=memilio_color_16, label='MEmilio ABM (16 cores)', linestyle='dotted')

        # Plot Covasim (fewer data points)
        covasim_pop_sizes = raw_data.population_sizes[:len(raw_data.covasim)]
        ax.plot(covasim_pop_sizes, raw_data.covasim,
                marker='o',  
                color=covasim_color, label='Covasim', linestyle='-')

        # Plot OpenCOVID (fewer data points)
        opencovid_pop_sizes = raw_data.population_sizes[:len(
            raw_data.opencovid)]
        ax.plot(opencovid_pop_sizes, raw_data.opencovid,
                marker='o',  
                color=opencovid_color, label='OpenCOVID', linestyle='-')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Agents [#]')
        ax.set_ylabel('Runtime [s]')

        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=LEGEND_FONTSIZE, loc='lower right',
                  framealpha=0.95, edgecolor='gray', fancybox=False)

        # Add speedup annotations between different implementations (min and max only)
        # MEmilio 16 cores vs Covasim (use covasim color)
        self._annotate_speedups_minmax(ax, covasim_pop_sizes,
                                       raw_data.memilio_times_sixteen_cores[:len(
                                           raw_data.covasim)],
                                       raw_data.covasim, covasim_color)

        # MEmilio 16 cores vs OpenCOVID (use opencovid color)
        self._annotate_speedups_minmax(ax, opencovid_pop_sizes,
                                       raw_data.memilio_times_sixteen_cores[:len(
                                           raw_data.opencovid)],
                                       raw_data.opencovid, opencovid_color)

        if save_path:
            fig.savefig(save_path, dpi=DPI)
            # Also save PDF version
            pdf_path = str(save_path).replace('.png', '.pdf')
            fig.savefig(pdf_path, dpi=DPI)
            print(
                f"Agent scaling plots saved to:\n  {save_path}\n  {pdf_path}")

        return fig, ax

    def plot_one_node_strong_scaling(self, raw_data, save_path=None):
        """Plot 3: One node strong scaling."""
        set_fontsize()
        figsize = (8, 5)
        panel = (0.2, 0.2, 0.78, 0.75)
        fig = plt.figure(figsize=figsize, dpi=DPI)
        ax = fig.add_axes(panel)

        # Set tick label sizes
        ax.tick_params(axis='both', which='both', labelsize=TICK_FONTSIZE)

        # Use only the data points that exist (skip the first element which is 1)
        cores = raw_data.strong_scaling_cores[0:len(
            raw_data.memilio_strong_scaling_128_runs_one_node)]
        runtimes = raw_data.memilio_strong_scaling_128_runs_one_node[0:]
        runtimes_graph_ode = raw_data.memilio_graph_ode_strong_scaling_128_runs_one_node[0:]
        runtimes_lct = raw_data.memilio_lct_strong_scaling_128_runs_one_node[0:]
        runtimes_ide = raw_data.memilio_ide_strong_scaling_128_runs_one_node[0:]

        # Plot runtime
        ax.plot(cores, runtimes,
                marker='o',  
                color='#741194', label='ABM', linestyle='-')
        ax.plot(cores, runtimes_graph_ode,
                marker='o',  
                color="#0D47A1", label='Graph-ODE\nw/ mobility', linestyle='-')
        ax.plot(cores, runtimes_lct,
                marker='o',  
                color='#5D8A2B', label='LCT', linestyle='-')
        ax.plot(cores, runtimes_ide,
                marker='o',  
                color='#E89A63', label='IDE', linestyle='-')
        

        # Calculate and plot ideal scaling
        if len(runtimes) > 0:
            ideal_scaling = [((runtimes[1] / c * cores[1])*0.26)for c in cores]
            ax.plot(cores, ideal_scaling,
                    'k--', linewidth=2, alpha=0.5, label='Ideal scaling')

        ax.set_xlabel('Cores [#]')
        ax.set_ylabel('Runtime [s]')

        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=LEGEND_FONTSIZE,
                  framealpha=0.95, edgecolor='gray', fancybox=False)

        # Set x-axis to log scale if useful, or linear
        ax.set_xscale('log', base=2)
        ax.set_yscale('log')
        ax.set_xticks(cores)
        ax.set_xticklabels([str(c) for c in cores])

        if save_path:
            fig.savefig(save_path, dpi=DPI)
            pdf_path = str(save_path).replace('.png', '.pdf')
            fig.savefig(pdf_path, dpi=DPI)
            print(
                f"One node strong scaling plots saved to:\n  {save_path}\n  {pdf_path}")

        return fig, ax

    def plot_multi_node_strong_scaling(self, raw_data, save_path=None):
        """Plot 4: Multi-node strong scaling."""
        set_fontsize()
        figsize = (8, 5)
        panel = (0.2, 0.2, 0.78, 0.75)
        fig = plt.figure(figsize=figsize, dpi=DPI)
        ax = fig.add_axes(panel)

        # Set tick label sizes
        ax.tick_params(axis='both', which='both', labelsize=TICK_FONTSIZE)

        # Use only the data points that exist (skip the first element which is 1)
        nodes = raw_data.strong_scaling_nodes[0:len(
            raw_data.memilio_strong_scaling_128_runs_multiple_nodes)]
        runtimes = raw_data.memilio_strong_scaling_128_runs_multiple_nodes[0:]

        # Plot runtime
        ax.plot(nodes, runtimes,
                marker='o',  
                color='#741194', label='ABM', linestyle='-')

        # Calculate and plot ideal scaling
        if len(runtimes) > 0:
            ideal_scaling = [runtimes[1] / n * nodes[1] for n in nodes]
            ax.plot(nodes, ideal_scaling,
                    'k--', linewidth=2, alpha=0.5, label='Ideal scaling')

        ax.set_xlabel('Cores [#]')
        ax.set_ylabel('Runtime [s]')

        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=LEGEND_FONTSIZE,
                  framealpha=0.95, edgecolor='gray', fancybox=False)

        # Set x-axis to log scale if useful, or linear
        ax.set_xscale('log', base=2)
        ax.set_yscale('log')
        ax.set_xticks(nodes)
        ax.set_xticklabels([str(n) for n in nodes])

        if save_path:
            fig.savefig(save_path, dpi=DPI)
            pdf_path = str(save_path).replace('.png', '.pdf')
            fig.savefig(pdf_path, dpi=DPI)
            print(
                f"Multi-node strong scaling plots saved to:\n  {save_path}\n  {pdf_path}")

        return fig, ax

    def plot_weak_scaling(self, raw_data, save_path=None):
        """Plot 2: Weak scaling plot for different agent counts per core."""
        set_fontsize()
        figsize = (8, 5)
        panel = (0.2, 0.2, 0.78, 0.75)
        fig = plt.figure(figsize=figsize, dpi=DPI)
        ax = fig.add_axes(panel)

        # Number of cores for weak scaling
        num_cores = [1, 2, 4, 8, 16, 32]

        # Plot each configuration
        ax.plot(num_cores, raw_data.memilio_weak_scaling_250k,
                marker='o',  
                label='250k agents per core', linestyle='-')

        ax.plot(num_cores, raw_data.memilio_weak_scaling_500k,
                marker='s',  
                label='500k agents per core', linestyle='-')

        ax.plot(num_cores, raw_data.memilio_weak_scaling_1m,
                marker='^',  
                label='1M agents per core', linestyle='-')

        ax.plot(num_cores, raw_data.memilio_weak_scaling_2m,
                marker='D',  
                label='2M agents per core', linestyle='-')

        # Add ideal scaling line (constant runtime)
        # ideal_runtime = raw_data.memilio_weak_scaling_250k[0]
        # ax.plot(num_cores, [ideal_runtime] * len(num_cores),
        #         'k--', linewidth=2, alpha=0.5, label='Ideal (constant runtime)')

        ax.set_xlabel('Cores [#]')
        ax.set_ylabel('Runtime [s]')

        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=LEGEND_FONTSIZE,
                  framealpha=0.95, edgecolor='gray', fancybox=False)

        # Set x-axis ticks to actual core counts
        ax.set_xticks(num_cores)
        ax.set_xticklabels([str(c) for c in num_cores])

        ax.set_yscale('log')

        if save_path:
            fig.savefig(save_path, dpi=DPI)
            pdf_path = str(save_path).replace('.png', '.pdf')
            fig.savefig(pdf_path, dpi=DPI)
            print(f"Weak scaling plots saved to:\n  {save_path}\n  {pdf_path}")

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
        print("WEAK SCALING EFFICIENCY TABLE: ABM")
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
        cores_one_node = raw_data.strong_scaling_cores[0:len(
            raw_data.memilio_strong_scaling_128_runs_one_node)]
        runtimes_one_node = raw_data.memilio_strong_scaling_128_runs_one_node[0:]

        # Multi-node strong scaling
        nodes_multi = raw_data.strong_scaling_nodes[0:len(
            raw_data.memilio_strong_scaling_128_runs_multiple_nodes)]
        runtimes_multi = raw_data.memilio_strong_scaling_128_runs_multiple_nodes[0:]

        def calculate_strong_efficiency(baseline, runtimes, units):
            """Calculate strong scaling efficiency: (baseline / (runtime * units)) * 100"""
            return [(baseline / (runtime * unit) * 100) for runtime, unit in zip(runtimes, units)]

        print("\n" + "="*70)
        print("STRONG SCALING EFFICIENCY TABLE: ABM")
        print("="*70)

        # One node tablez
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
    analyzer = BenchmarkAnalyzer(fontsize=20)
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
