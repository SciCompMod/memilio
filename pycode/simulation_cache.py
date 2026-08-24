import os

if "KERAS_BACKEND" not in os.environ:
    os.environ["KERAS_BACKEND"] = "jax"

import json
from datetime import datetime
from pathlib import Path

import numpy as np

from memilio.simulation.abm import ABMPopulation
from abm_batch import run_batch_designs
from abm_inference import (
    TOTAL_POPULATION, PARAM_KEYS, OBSERVED_DAYS, N_CT_BINS,
    SOURCES_SURVEY, SOURCES_DIAGNOSTIC, SOURCES_POOLED,
    prior_theta, make_design, total_budget_for_rate, _assemble_dual_stream,
)

# ─────────────────────────────────────────────────────────────────────────────
CACHE_DIR = Path(__file__).parent / "simulation_cache_data"
SHARD_SIZE = 10_000         
TARGET_N = 200_000    

TOTAL_POPULATION_ = TOTAL_POPULATION
DESIGN_RATE_LO, DESIGN_RATE_HI = 50, 2000   # tests / 100k people / day
PERIODS = [1, 2, 3]
# ─────────────────────────────────────────────────────────────────────────────


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def _manifest_path(cache_dir):
    return Path(cache_dir) / "manifest.json"


def _this_config():
    return {
        "total_population": TOTAL_POPULATION_,
        "design_rate_lo": DESIGN_RATE_LO,
        "design_rate_hi": DESIGN_RATE_HI,
        "periods": list(PERIODS),
        "observed_days": OBSERVED_DAYS,
        "n_bins": N_CT_BINS,
        "shard_size": SHARD_SIZE,
    }


def _existing_shards(cache_dir):
    """Sorted list of shard file Paths already on disk."""
    return sorted(Path(cache_dir).glob("shard_*.npz"))


def read_manifest_config(cache_dir):
    """The `config` dict a cache was generated under."""
    return json.loads(_manifest_path(cache_dir).read_text())["config"]


def cache_size(cache_dir):
    """Total number of examples currently available at `cache_dir` (0 if it doesn't exist yet, or
    is still empty/mid-first-shard)."""
    cache_dir = Path(cache_dir)
    if not cache_dir.exists():
        return 0
    total = 0
    for path in _existing_shards(cache_dir):
        with np.load(path) as d:
            total += len(d["total_budget"])
    return total


def _load_or_init_manifest(cache_dir):
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = _manifest_path(cache_dir)
    this_config = _this_config()
    if path.exists():
        manifest = json.loads(path.read_text())
        if manifest["config"] != this_config:
            raise ValueError(
                f"{path} was generated with a different config than the one currently set:\n"
                f"  existing: {manifest['config']}\n  current:  {this_config}\n"
                "Point CACHE_DIR at a new directory instead of mixing configs into one cache.")
        return manifest
    manifest = {"config": this_config, "created": datetime.now().isoformat()}
    path.write_text(json.dumps(manifest, indent=2))
    return manifest


def generate_shard(shard_idx, cache_dir, population, rng):
    """Simulate SHARD_SIZE fresh examples and write them as shard_{idx:05d}.npz. rng is shared
    across shards so examples are never repeated."""
    draws = [prior_theta(rng) for _ in range(SHARD_SIZE)]
    totals = np.exp(rng.uniform(np.log(total_budget_for_rate(DESIGN_RATE_LO)),
                                np.log(total_budget_for_rate(DESIGN_RATE_HI)), SHARD_SIZE))
    periods = rng.choice(PERIODS, SHARD_SIZE)
    designs = [make_design(t, p) for t, p in zip(totals, periods)]
    outs = run_batch_designs(draws, designs, population, show_progress=True, max_workers=14)
    data = _assemble_dual_stream(draws, totals, periods, outs)

    path = Path(cache_dir) / f"shard_{shard_idx:05d}.npz"
  
    tmp_path = path.with_name(f".tmp_{path.name}")
    np.savez_compressed(tmp_path, **data)  
                                           
    tmp_path.rename(path) 
                           
    return path


def main():
    _load_or_init_manifest(CACHE_DIR)  # create/validate; raises if CACHE_DIR's config doesn't match
    shards = _existing_shards(CACHE_DIR)
    n_existing = len(shards) * SHARD_SIZE
    log(f"cache: {CACHE_DIR}  |  {len(shards)} shard(s), ~{n_existing:,} examples so far "
       f"(target {TARGET_N:,})")

    if n_existing >= TARGET_N:
        log("target already reached; nothing to do.")
        return

    rng = np.random.default_rng() 
                                   
    population = ABMPopulation(total_population=TOTAL_POPULATION_)

    next_idx = len(shards)
    n_total = n_existing
    while n_total < TARGET_N:
        log(f"generating shard {next_idx} ({SHARD_SIZE:,} examples)...")
        path = generate_shard(next_idx, CACHE_DIR, population, rng)
        n_total += SHARD_SIZE
        log(f"  wrote {path.name}  |  ~{n_total:,}/{TARGET_N:,} examples total")
        next_idx += 1

    log("target reached.")


def load_cached_dataset(cache_dir, n=None, sources=SOURCES_POOLED, scale_by_tests=False, rng=None):
    """Load (a random subset of, if `n` given) the cache at `cache_dir`.
    Loads shard-by-shard (shuffled) rather than the whole cache at once when `n` is well under the
    cache's total size, so sampling a modest subset from a huge cache doesn't require reading all
    of it into memory first."""
    cache_dir = Path(cache_dir)
    manifest = json.loads(_manifest_path(cache_dir).read_text())
    if manifest["config"]["observed_days"] != OBSERVED_DAYS or manifest["config"]["n_bins"] != N_CT_BINS:
        raise ValueError(
            f"{cache_dir} was written with observed_days/n_bins="
            f"{manifest['config']['observed_days']}/{manifest['config']['n_bins']}, but this "
            f"module currently uses {OBSERVED_DAYS}/{N_CT_BINS} -- the cache predates a schema "
            "change (e.g. OBSERVED_DAYS's margin) and can't be safely reconstructed as-is.")
    shards = _existing_shards(cache_dir)
    if not shards:
        raise FileNotFoundError(f"no shards found in {cache_dir}")

    rng = rng or np.random.default_rng()

    order = rng.permutation(len(shards)) if n is not None else np.arange(len(shards))

    collected = []
    n_collected = 0
    for shard_i in order:
        with np.load(shards[shard_i]) as d:
            city = np.zeros_like(d["histogram_ct_survey"])
            tot = np.zeros_like(d["daily_test_total_survey"])
            if SOURCES_SURVEY.issubset(sources):
                city = city + d["histogram_ct_survey"]
                tot = tot + d["daily_test_total_survey"]
            if SOURCES_DIAGNOSTIC.issubset(sources):
                city = city + d["histogram_ct_diagnostic"]
                tot = tot + d["daily_test_total_diagnostic"]

            if scale_by_tests:
                city = np.divide(city, tot[..., None], out=np.zeros_like(city), where=tot[..., None] > 0)

            shard_data = {k: d[k] for k in PARAM_KEYS}
            shard_data["total_budget"] = d["total_budget"]
            shard_data["period"] = d["period"]
            shard_data["n_ever_infected"] = d["n_ever_infected"]
            shard_data["histogram_ct"] = city

        shard_n = len(shard_data["total_budget"])
        if n is not None and n_collected + shard_n > n:
            keep = rng.choice(shard_n, size=n - n_collected, replace=False)
            shard_data = {k: v[keep] for k, v in shard_data.items()}
            shard_n = n - n_collected

        collected.append(shard_data)
        n_collected += shard_n
        if n is not None and n_collected >= n:
            break

    merged = {k: np.concatenate([c[k] for c in collected], axis=0) for k in collected[0]}
    data = {k: merged[k].reshape(-1, 1) for k in PARAM_KEYS}
    data.update({
        "log_budget": np.log(merged["total_budget"]).reshape(-1, 1),
        "frequency": (1.0 / merged["period"]).reshape(-1, 1),
        "histogram_ct": merged["histogram_ct"],
        "histogram_bin": merged["histogram_ct"].sum(axis=-1),
        "n_ever_infected": merged["n_ever_infected"],
    })
    return data


if __name__ == "__main__":
    main()
