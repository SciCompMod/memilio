"""Helpers for running many independent ABM forward passes in parallel.

- `run_batch` uses a thread pool. forward_pass() releases the GIL for the duration of the
  C++ call, so real OS threads can run it concurrently, sharing one ABMPopulation/TestingBudget
  in memory.

- `run_batch_processes` uses a process pool instead. Each worker process has its own
  independent heap/allocator and builds its own copy of the population once
"""
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

from memilio.simulation import set_log_level
from memilio.simulation.abm import ABMPopulation, forward_pass


def run_batch(params_list, population, design, max_workers=None, show_progress=True):
    """Run forward_pass(population, params, design) for each params dict, across worker threads.

    Args:
        params_list: sequence of parameter dicts (e.g. {"beta": 1.0, "t_exposed": 5.0}); see
            forward_pass() for the recognised keys.
        population: a single, shared ABMPopulation.
        design: a single, shared TestingBudget.
        max_workers: worker thread count; defaults to os.cpu_count() + 4 (ThreadPoolExecutor's
            default)
        show_progress: show a tqdm progress bar.

    Returns:
        List of results (forward_pass dicts), in the same order as params_list.
    """
    results = [None] * len(params_list)
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(forward_pass, population, params, design): i
            for i, params in enumerate(params_list)
        }
        iterator = as_completed(futures)
        if show_progress:
            from tqdm.auto import tqdm
            iterator = tqdm(iterator, total=len(futures), unit="sim")
        for future in iterator:
            results[futures[future]] = future.result()
    return results


def run_batch_designs(params_list, design_list, population, max_workers=None, show_progress=True):
    """Like run_batch(), but each simulation uses its own TestingBudget (design_list[i]).

    Used to amortize inference over the testing design: each training example is generated at a
    different budget drawn from a design prior, so a single network learns p(theta | X, design).

    Args:
        params_list: sequence of parameter dicts.
        design_list: sequence of TestingBudget, one per params entry (same length/order).
        population: a single, shared ABMPopulation.
        max_workers, show_progress: as in run_batch().

    Returns:
        List of results (forward_pass dicts), in the same order as params_list.
    """
    assert len(params_list) == len(design_list), "params_list and design_list must be the same length"
    results = [None] * len(params_list)
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(forward_pass, population, params, design): i
            for i, (params, design) in enumerate(zip(params_list, design_list))
        }
        iterator = as_completed(futures)
        if show_progress:
            from tqdm.auto import tqdm
            iterator = tqdm(iterator, total=len(futures), unit="sim")
        for future in iterator:
            results[futures[future]] = future.result()
    return results


_pool_population = None
_pool_design = None


def _init_process_worker(total_population, design, log_level=None):
    """ProcessPoolExecutor initializer: builds one ABMPopulation per worker process.
    """
    global _pool_population, _pool_design
    _pool_population = ABMPopulation(total_population=total_population)
    _pool_design = design
    if log_level is not None:
        forward_pass(_pool_population, {"beta": 1.0, "t_exposed": 5.0}, _pool_design)
        set_log_level(log_level)


def _run_one_process(params):
    return forward_pass(_pool_population, params, _pool_design)


def run_batch_processes(params_list, total_population, design, max_workers=None,
                        show_progress=True, log_level=None):
    """Like run_batch(), but spreads work across worker *processes* instead of threads.

    Args:
        params_list: sequence of parameter dicts (e.g. {"beta": 1.0, "t_exposed": 5.0}); see
            forward_pass() for the recognised keys.
        total_population: passed to ABMPopulation(total_population=...) in each worker.
        design: a TestingBudget; pickled and sent to each worker as-is.
        max_workers: worker process count; defaults to os.cpu_count().
        show_progress: show a tqdm progress bar.
        log_level: if given (a memilio.simulation.LogLevel).

    Returns:
        List of results, in the same order as params_list.
    """
    results = [None] * len(params_list)
    with ProcessPoolExecutor(
        max_workers=max_workers, initializer=_init_process_worker,
        initargs=(total_population, design, log_level),
    ) as pool:
        futures = {pool.submit(_run_one_process, p): i for i, p in enumerate(params_list)}
        iterator = as_completed(futures)
        if show_progress:
            from tqdm.auto import tqdm
            iterator = tqdm(iterator, total=len(futures), unit="sim")
        for future in iterator:
            results[futures[future]] = future.result()
    return results
