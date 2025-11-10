# import os
# # Choose sensible values for your machine; start with physical core count
# os.environ["NUMBA_NUM_THREADS"] = "8"       # Numba worker pool size
# os.environ["NUMBA_THREADING_LAYER"] = "omp"  # or 'tbb'/'omp' if available

# # BLAS/LAPACK backends (pick what your NumPy uses)
# os.environ["OMP_NUM_THREADS"] = "8"         # OpenMP
# os.environ["OPENBLAS_NUM_THREADS"] = "8"    # OpenBLAS
# os.environ["MKL_NUM_THREADS"] = "8"         # Intel MKL

import covasim as cv
import time
import sys


def covasim_benchmark(pop_size, n_days=120):
    """
    Benchmark covasim simulation with specified population size.

    Args:
        pop_size: Number of agents in the simulation
        n_days: Number of days to simulate (default 120 to match ABM 5 days * 24)

    Returns:
        Runtime in seconds
    """
    # Calculate initial infections as 0.05% of population (matching ABM benchmark)
    pop_infected = max(1, int(pop_size * 0.0005))

    pars = dict(
        pop_size=pop_size,
        pop_infected=pop_infected,
        start_day='2025-01-01',
        n_days=n_days,
        pop_type='hybrid',
        location='Germany',
        rescale=False
    )

    sim = cv.Sim(pars)
    cv.options.set(numba_parallel='full')
    sim.initialize()

    t_start = time.process_time()
    sim.run(verbose=0)
    t_end = time.process_time()

    runtime = t_end - t_start
    print(f"Population size: {pop_size:,}, Runtime: {runtime:.3f}s")
    return runtime


if __name__ == "__main__":
    # Population sizes matching the ABM benchmark
    pop_sizes = [1000000, 2000000, 4000000, 8000000,
                 16000000, 32000000, 64000000]
    # pop_sizes = [256000000, 512000000, 1024000000]

    # If command line argument provided, use that specific size
    if len(sys.argv) > 1:
        try:
            pop_size = int(sys.argv[1])
            covasim_benchmark(pop_size)
        except ValueError:
            print("Usage: python covasim_benchmark.py [population_size]")
            sys.exit(1)
    else:
        # Run all benchmark sizes
        print("Running Covasim benchmarks with population sizes matching ABM benchmark:")
        print("=" * 60)

        total_time = 0
        for pop_size in pop_sizes:
            runtime = covasim_benchmark(pop_size)
            total_time += runtime

        print("=" * 60)
        print(f"Total runtime for all benchmarks: {total_time:.3f}s")
