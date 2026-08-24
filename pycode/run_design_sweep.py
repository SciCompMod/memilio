"""Design-amortized CT-vs-binary information-gain sweep over the surveillance design.
"""
import os

if "KERAS_BACKEND" not in os.environ:
    os.environ["KERAS_BACKEND"] = "jax"

import csv
import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import keras
import bayesflow as bf

from memilio.simulation.abm import ABMPopulation
from abm_batch import run_batch_designs
from abm_inference import (
    PARAM_KEYS, PAR_NAMES, TOTAL_POPULATION, EPIDEMIC_PARAM_KEYS,
    SOURCES_SURVEY, SOURCES_DIAGNOSTIC, SOURCES_POOLED,
    prior_theta, make_design, make_val_at_design, build_workflow, _assemble,
    estimate_mutual_information, sample_posterior, local_jacobian_svd,
    shared_bounds, apply_recovery_bounds, apply_pairgrid_bounds, plot_posterior_comparison_grid,
    total_budget_for_rate, rate_per_100k_per_day,
)
import simulation_cache

# Set None to not use Cache
CACHE_DIR = simulation_cache.CACHE_DIR

# Fixed parent for every run's results folder 
OUTPUT_DIR = Path(__file__).parent / "design_sweep_results"

SOURCES_NAMES = {SOURCES_SURVEY: "survey", SOURCES_DIAGNOSTIC: "diagnostic", SOURCES_POOLED: "pooled"}


# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────
@dataclass
class SweepConfig:
    """Class to hold main knobs for a single sweep run."""
    n_train: int = 100000
    n_val_per_design: int = 1000

    total_population: int = TOTAL_POPULATION
    early_stop_patience: int = 15
    early_stop_start_from_epoch: int = 50
    train_epochs: int = 120

    sources: frozenset = SOURCES_SURVEY
    scale_by_tests: bool = True

    design_rate_lo: float = 50       # tests / 100k people / day
    design_rate_hi: float = 2000     # tests / 100k people / day
    periods: list = field(default_factory=lambda: [1, 2, 3])  # cadence: sample every N days

    # Grid at which to report Delta.
    rate_grid: list = field(default_factory=lambda: [50, 100, 200, 400, 600, 800, 1000, 1300, 1600, 2000])
    period_grid: list = field(default_factory=lambda: [1, 2, 3])

    @property
    def design_total_lo(self):
        return total_budget_for_rate(self.design_rate_lo)

    @property
    def design_total_hi(self):
        return total_budget_for_rate(self.design_rate_hi)

    @property
    def total_budget_grid(self):
        return [total_budget_for_rate(r) for r in self.rate_grid]

    def tag(self):
        """Short string identifying this config, for use in results folder names."""
        parts = [
            SOURCES_NAMES.get(self.sources, "sources" + "".join(str(s) for s in sorted(self.sources))),
            "rate" if self.scale_by_tests else "count",
            f"train{self.n_train // 1000}k",
            f"ep{self.train_epochs}",
            f"r{self.design_rate_lo}-{self.design_rate_hi}",
        ]
        if self.total_population != TOTAL_POPULATION:
            parts.append(f"pop{self.total_population}")
        if self.early_stop_patience != 15:
            parts.append(f"pat{self.early_stop_patience}")
        return "_".join(parts)


def make_training_dataset(n, population, rng, config, show_progress=True):
    """Mixed-design training set: theta AND (total_budget, period) drawn per simulation."""
    draws = [prior_theta(rng) for _ in range(n)]
    totals = np.exp(rng.uniform(np.log(config.design_total_lo), np.log(config.design_total_hi), n))
    periods = rng.choice(config.periods, n)
    # make_val_at_design() assumes one fixed design.
    designs = [make_design(t, p) for t, p in zip(totals, periods)]
    outs = run_batch_designs(draws, designs, population, show_progress=show_progress, max_workers=os.cpu_count())
    return _assemble(draws, totals, periods, outs, sources=config.sources, scale_by_tests=config.scale_by_tests)


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def _training_data_from_cache_or_simulate(config, population, rng):
    """Training set only. From the on-disk cache if it currently has enough examples, else
    simulated fresh. The early-stopping validation set is always simulated fresh."""
    available = simulation_cache.cache_size(CACHE_DIR) if CACHE_DIR is not None else 0
    if available >= config.n_train:
        manifest = simulation_cache.read_manifest_config(CACHE_DIR)
        if manifest["total_population"] != config.total_population:
            log(f"cache at {CACHE_DIR} was generated with population={manifest['total_population']}, "
               f"but this config wants {config.total_population} -- simulating fresh instead.")
        else:
            log(f"loading training data from cache ({CACHE_DIR}, {available:,} examples available, "
               f"cache design range {manifest['design_rate_lo']}-{manifest['design_rate_hi']}/100k/day "
               f"periods={manifest['periods']} -- check these cover what this config needs)...")
            return simulation_cache.load_cached_dataset(CACHE_DIR, n=config.n_train, sources=config.sources,
                                                         scale_by_tests=config.scale_by_tests, rng=rng)
    elif CACHE_DIR is not None:
        log(f"cache at {CACHE_DIR} has {available:,}/{config.n_train:,} examples needed -- simulating fresh instead.")

    log(f"Simulating mixed-design training data (n={config.n_train}, population={config.total_population})...")
    return make_training_dataset(config.n_train, population, rng, config)


def main(config=None):
    """Run one full sweep for `config` (a SweepConfig; defaults to SweepConfig()).
    Returns the output_dir where results were written."""
    config = config or SweepConfig()
    output_dir = OUTPUT_DIR / f"{config.tag()}_{datetime.now():%Y%m%d_%H%M%S}"
    output_dir.mkdir(parents=True, exist_ok=True)
    log(f"Writing results to {output_dir}")

    rng = np.random.default_rng(20260813)
    population = ABMPopulation(total_population=config.total_population)

    raw_train = _training_data_from_cache_or_simulate(config, population, rng)
    log("Simulating a mixed-design validation set (for early stopping)...")
    raw_val = make_training_dataset(config.n_val_per_design * 2, population, rng, config)

    early_stop = keras.callbacks.EarlyStopping(monitor="val_loss", patience=config.early_stop_patience,
                                               restore_best_weights=True,
                                               start_from_epoch=config.early_stop_start_from_epoch)

    _prior_draws = [prior_theta(rng) for _ in range(2000)]
    prior_draws = {k: np.array([d[k] for d in _prior_draws]) for k in PARAM_KEYS}

    KEY_LABEL = [("histogram_ct", "ct"), ("histogram_bin", "bin")]

    workflows, posts = {}, {}
    for key, label in KEY_LABEL:
        log(f"Training {label} workflow (design-amortized)...")
        wf = build_workflow(key, checkpoint_dir=output_dir)
        hist = wf.fit_offline(data=raw_train, epochs=config.train_epochs, batch_size=64,
                              validation_data=raw_val, callbacks=[early_stop])
        f = bf.diagnostics.plots.loss(hist)
        f.savefig(output_dir / f"loss_{label}.png", dpi=150, bbox_inches="tight")
        plt.close(f)

        log(f"Sampling {label} posterior on the held-out mixed-design validation set...")
        workflows[label] = wf
        posts[label] = sample_posterior(wf, raw_val, key)

    bounds = shared_bounds(raw_val, posts)

    example_indices = [0, 1, 2]
    for key, label in KEY_LABEL:
        wf, post = workflows[label], posts[label]
        log(f"Posterior diagnostics ({label})...")

        f = bf.diagnostics.plots.recovery(
            post, raw_val, variable_keys=PARAM_KEYS, variable_names=PAR_NAMES)
        apply_recovery_bounds(f, bounds)
        f.suptitle(f"Recovery -- {label} (design-amortized)", y=1.02)
        f.savefig(output_dir / f"recovery_{label}.png", dpi=150, bbox_inches="tight")
        plt.close(f)

        f = bf.diagnostics.plots.calibration_ecdf(
            post, raw_val, variable_keys=PARAM_KEYS,
            variable_names=PAR_NAMES, difference=True)
        f.suptitle(f"ECDF calibration -- {label} (design-amortized)", y=1.02)
        f.savefig(output_dir / f"calibration_ecdf_{label}.png", dpi=150, bbox_inches="tight")
        plt.close(f)

        for example_idx in [0, 1]:
            g = bf.diagnostics.plots.pairs_posterior(
                estimates=post, targets=raw_val, priors=prior_draws, dataset_id=example_idx,
                variable_keys=PARAM_KEYS, variable_names=PAR_NAMES)
            apply_pairgrid_bounds(g, bounds)
            g.figure.suptitle(f"{label} -- prior vs posterior (example #{example_idx})", y=1.02)
            g.figure.savefig(output_dir / f"pairs_posterior_{label}_{example_idx}.png",
                             dpi=150, bbox_inches="tight")
            plt.close(g.figure)

        log(f"Local Jacobian identifiability ({label})...")
        for idx in example_indices:
            S, U, _ = local_jacobian_svd(wf, raw_val, idx, key)
            ratio = S[0] / S[-1] if S[-1] > 0 else np.inf
            sloppy_dir = U[:, np.argmax(S)]
            log(f"    instance #{idx}: local FRACTIONAL posterior std={S.round(3)} (ratio={ratio:.1f}x), "
                f"sloppiest direction in log-theta ({','.join(PARAM_KEYS)})={sloppy_dir.round(3)}")

        workflows[label] = wf


    log("Posterior comparison grids (CT vs binary)...")
    for idx in example_indices:
        f = plot_posterior_comparison_grid(posts, raw_val, idx, bounds=bounds)
        f.savefig(output_dir / f"posterior_comparison_{idx}.png", dpi=150, bbox_inches="tight")
        plt.close(f)

    total_budget_grid = config.total_budget_grid
    rows = []
    for total_budget in total_budget_grid:
        rate = rate_per_100k_per_day(total_budget)
        for period in config.period_grid:
            log(f"[rate={rate:.0f}/100k/day period={period}] simulating validation set and estimating Delta...")
            val = make_val_at_design(config.n_val_per_design, population, total_budget, period, rng,
                                     sources=config.sources, scale_by_tests=config.scale_by_tests)
            kl_ct,  frac_invalid_ct  = estimate_mutual_information(workflows["ct"],  val, "histogram_ct",
                                                                    keep_keys=EPIDEMIC_PARAM_KEYS)
            kl_bin, frac_invalid_bin = estimate_mutual_information(workflows["bin"], val, "histogram_bin",
                                                                    keep_keys=EPIDEMIC_PARAM_KEYS)
            I_ct, I_bin = kl_ct.mean(), kl_bin.mean()
            delta = I_ct - I_bin
            se_delta = (kl_ct - kl_bin).std(ddof=1) / np.sqrt(len(kl_ct))
            rows.append({
                "total_budget": total_budget, "rate_per_100k_per_day": rate,
                "period": period, "frequency": 1.0 / period,
                "I_ct": I_ct, "I_bin": I_bin, "delta": delta, "se_delta": se_delta,
                "n_ever_infected": float(val["n_ever_infected"].mean()),
                "frac_invalid_ct": frac_invalid_ct, "frac_invalid_bin": frac_invalid_bin,
            })
            log(f"    I_ct={I_ct:.3f} I_bin={I_bin:.3f} delta={delta:.3f} +/- {se_delta:.3f} "
                f"| n_ever_infected~{val['n_ever_infected'].mean():.0f} "
                f"| out-of-support: ct={frac_invalid_ct:.2%} bin={frac_invalid_bin:.2%}")

    with open(output_dir / "design_sweep.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    with open(output_dir / "design_sweep.json", "w") as fh:
        json.dump(rows, fh, indent=2)

    def grid(field_name):
        g = np.full((len(total_budget_grid), len(config.period_grid)), np.nan)
        for r in rows:
            g[total_budget_grid.index(r["total_budget"]), config.period_grid.index(r["period"])] = r[field_name]
        return g

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, (field_name, title) in zip(axes, [("delta", "Delta(d) [nats] (epidemic params only)"),
                                               ("I_ct", "I_ct [nats] (epidemic params only)"),
                                               ("n_ever_infected", "epidemic size (ever infected)")]):
        im = ax.imshow(grid(field_name), origin="lower", aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(config.period_grid))); ax.set_xticklabels(config.period_grid)
        ax.set_yticks(range(len(config.rate_grid))); ax.set_yticklabels(config.rate_grid)
        ax.set_xlabel("period (days between surveillance events)")
        ax.set_ylabel("tests / 100k people / day")
        ax.set_title(title)
        fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(output_dir / "design_sweep.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"Sweep complete. Delta surface -> {output_dir / 'design_sweep.png'}")
    return output_dir


if __name__ == "__main__":
    main()
