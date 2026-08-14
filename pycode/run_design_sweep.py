"""Design-amortized CT-vs-binary information-gain sweep over the 2-D surveillance design
(total testing budget x testing frequency).

budget_fraction is the fraction sampled per surveillance event. A design is (total_budget,
period): total tests-per-person over the run, and the cadence (test every `period` days). Both
`total_budget` and `frequency = 1/period` are direct conditions in the neural network,
so one CT workflow and one binary workflow learn p(theta | X, d) amortized over designs.

Usage: python run_design_sweep.py
"""
import os

if "KERAS_BACKEND" not in os.environ:
    os.environ["KERAS_BACKEND"] = "jax"

import csv
import json
import warnings
from collections.abc import Collection
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import keras
import bayesflow as bf
from scipy import stats

from memilio.simulation.abm import ABMPopulation, TestingBudget
from abm_batch import run_batch_designs

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────
TOTAL_POPULATION = 10000
QUARANTINE_COMPLIANCE = 0.6     # city-level P(isolate | positive), fixed per training run
N_TRAIN = 16000
N_VAL_PER_DESIGN = 600

SIM_DAYS = 30
SURVEY, DIAGNOSTIC = 0, 1        # source tags from forward_pass
N_CT_BINS = 41                   # CT 0..40

# Design priors (amortised over): total_budget ~ log-uniform, period from a discrete set.
DESIGN_TOTAL_LO, DESIGN_TOTAL_HI = 0.1, 1.5   # total survey tests per person over the run
PERIODS = [1, 2, 3, 5, 7]                      # cadence: sample every `period` days

# Grid at which to report Delta (the surface).
TOTAL_BUDGET_GRID = [0.15, 0.3, 0.6, 1.0, 1.5]
PERIOD_GRID = [1, 2, 3, 5, 7]

# theta priors (full-support lognormal, matching joint_inference).
PRIOR_BETA_S,  PRIOR_BETA_SCALE  = 0.5, 1.0
PRIOR_TEXP_S, PRIOR_TEXP_SCALE = 0.4, 5.0

NUM_POSTERIOR_SAMPLES = 500
OUTPUT_DIR = Path(__file__).parent / "design_sweep_results"


# ─────────────────────────────────────────────────────────────────────────────
# Data-generating process
# ─────────────────────────────────────────────────────────────────────────────
def n_events(period, sim_days=SIM_DAYS):
    """Number of surveillance events (testing days 0, period, 2*period, ... < sim_days)."""
    return (sim_days - 1) // period + 1


def records_to_histogram(out, sources: Collection[int] = frozenset({SURVEY}), sim_days=SIM_DAYS, n_bins=N_CT_BINS):
    """forward_pass() dict -> (sim_days, n_bins) city-wide daily CT histogram of positives.
    `sources` is a set of stream tags to pool: {0}=survey, {1}=diagnostic, {0, 1}=both.
    Positives are out['positives'] = [day, person, age, loc, ct, source]."""
    city = np.zeros((sim_days, n_bins))
    pos = out["positives"]
    if pos.shape[0]:
        m = np.isin(pos[:, 5], list(sources))
        days = np.round(pos[m, 0]).astype(int)
        cts = np.clip(np.round(pos[m, 4]).astype(int), 0, n_bins - 1)
        valid = (days >= 0) & (days < sim_days)
        np.add.at(city, (days[valid], cts[valid]), 1)
    return city


def prior_theta(rng):
    return {
        "beta":  rng.lognormal(np.log(PRIOR_BETA_SCALE),  PRIOR_BETA_S),
        "t_exposed": rng.lognormal(np.log(PRIOR_TEXP_SCALE), PRIOR_TEXP_S),
    }


def prior_log_prob(beta, t_exposed):
    return (stats.lognorm.logpdf(beta,  s=PRIOR_BETA_S,  scale=PRIOR_BETA_SCALE)
            + stats.lognorm.logpdf(t_exposed, s=PRIOR_TEXP_S, scale=PRIOR_TEXP_SCALE))


def make_design(total_budget, period):
    """(total_budget, period) -> a TestingBudget whose per-event fraction sums to total_budget."""
    per_event = total_budget / n_events(period)
    return TestingBudget(budget_fraction=float(per_event), test_period_days=int(period))


def _assemble(draws, totals, periods, outs):
    city = np.stack([records_to_histogram(o, {SURVEY}) for o in outs])
    return {
        "beta":  np.array([d["beta"]  for d in draws], dtype=np.float64).reshape(-1, 1),
        "t_exposed": np.array([d["t_exposed"] for d in draws], dtype=np.float64).reshape(-1, 1),
        # design conditions, already in network-friendly form:
        "log_budget": np.log(np.asarray(totals, dtype=np.float64)).reshape(-1, 1),
        "frequency":  (1.0 / np.asarray(periods, dtype=np.float64)).reshape(-1, 1),
        "histogram_ct":  city,
        "histogram_bin": city.sum(axis=-1),
        "n_ever_infected": np.array([o["n_ever_infected"][0, 0] for o in outs], dtype=np.float64),
    }


def make_training_dataset(n, population, rng, show_progress=True):
    """Mixed-design training set: theta AND (total_budget, period) drawn per simulation."""
    draws = [prior_theta(rng) for _ in range(n)]
    totals = np.exp(rng.uniform(np.log(DESIGN_TOTAL_LO), np.log(DESIGN_TOTAL_HI), n))
    periods = rng.choice(PERIODS, n)
    designs = [make_design(t, p) for t, p in zip(totals, periods)]
    outs = run_batch_designs(draws, designs, population, show_progress=show_progress)
    return _assemble(draws, totals, periods, outs)


def make_val_at_design(n, population, total_budget, period, rng):
    """Validation set at a single fixed design."""
    draws = [prior_theta(rng) for _ in range(n)]
    designs = [make_design(total_budget, period)] * n
    outs = run_batch_designs(draws, designs, population, show_progress=False)
    return _assemble(draws, np.full(n, total_budget), np.full(n, period), outs)


# ─────────────────────────────────────────────────────────────────────────────
# Model
# ─────────────────────────────────────────────────────────────────────────────
class GRU(bf.networks.SummaryNetwork):
    def __init__(self, hidden_dim=128, summary_dim=16, dropout=0.2, recurrent_dropout=0.1, **kwargs):
        super().__init__(**kwargs)
        self.gru1 = keras.layers.GRU(hidden_dim, return_sequences=True,
                                     dropout=dropout, recurrent_dropout=recurrent_dropout)
        self.gru2 = keras.layers.GRU(hidden_dim, dropout=dropout, recurrent_dropout=recurrent_dropout)
        self.summary_stats = keras.layers.Dense(summary_dim)

    def call(self, time_series, **kwargs):
        training = kwargs.get("stage") == "training"
        x = self.gru1(time_series, training=training)
        x = self.gru2(x, training=training)
        return self.summary_stats(x)


def build_workflow(summary_key):
    """Infers (beta, t_exposed) from a histogram summary AND the design (log_budget, frequency) as
    direct conditions -- so one trained network handles any (total_budget, period)."""
    adapter = (
        bf.adapters.Adapter()
        .convert_dtype("float64", "float32")
        .as_time_series(summary_key)
        .concatenate(["beta", "t_exposed"], into="inference_variables")
        .concatenate(["log_budget", "frequency"], into="inference_conditions")
        .rename(summary_key, "summary_variables")
        .log(["inference_variables", "summary_variables"], p1=True)
    )
    return bf.BasicWorkflow(
        adapter=adapter,
        inference_network=bf.networks.CouplingFlow(),
        summary_network=GRU(),
        standardize=["inference_variables", "inference_conditions"],
    )


def estimate_mutual_information(workflow, val, x_key, num_posterior_samples=NUM_POSTERIOR_SAMPLES, val_batch=32):
    """I(theta;X|d) for a validation set at a fixed design; the design (log_budget, frequency,
    carried in `val`) is passed as a condition to sample() and log_prob()."""
    n_val = len(val["beta"])
    kls = np.empty(n_val)
    n_invalid, n_total = 0, 0
    cond_keys = ("log_budget", "frequency")
    for start in range(0, n_val, val_batch):
        idx = slice(start, min(start + val_batch, n_val))
        cond = {x_key: val[x_key][idx], **{k: val[k][idx] for k in cond_keys}}
        post = workflow.sample(conditions=cond, num_samples=num_posterior_samples)
        beta_s, t_exposed_s = post["beta"].reshape(-1, 1), post["t_exposed"].reshape(-1, 1)
        b = val["beta"][idx].shape[0]
        tiled = {x_key: np.repeat(val[x_key][idx], num_posterior_samples, axis=0)}
        for k in cond_keys:
            tiled[k] = np.repeat(val[k][idx], num_posterior_samples, axis=0)
        log_q = workflow.log_prob(data={**tiled, "beta": beta_s, "t_exposed": t_exposed_s})
        log_q = log_q.reshape(b, num_posterior_samples)
        log_p = prior_log_prob(beta_s, t_exposed_s).reshape(b, num_posterior_samples)
        diff = np.where(np.isfinite(log_q) & np.isfinite(log_p), log_q - log_p, np.nan)
        n_invalid += np.isnan(diff).sum()
        n_total += diff.size
        kls[idx] = np.nanmean(diff, axis=1)
    return kls, n_invalid / n_total


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    global OUTPUT_DIR
    OUTPUT_DIR = OUTPUT_DIR.parent / f"{OUTPUT_DIR.name}_{datetime.now():%Y%m%d_%H%M%S}"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    log(f"Writing results to {OUTPUT_DIR}")

    rng = np.random.default_rng(20260813)
    population = ABMPopulation(total_population=TOTAL_POPULATION, quarantine_compliance=QUARANTINE_COMPLIANCE)

    log(f"Simulating mixed-design training data (n={N_TRAIN})...")
    raw_train = make_training_dataset(N_TRAIN, population, rng)
    log("Simulating a mixed-design validation set (for early stopping)...")
    raw_val = make_training_dataset(N_VAL_PER_DESIGN * 2, population, rng)

    early_stop = keras.callbacks.EarlyStopping(monitor="val_loss", patience=15,
                                               restore_best_weights=True, start_from_epoch=50)
    workflows = {}
    for key, label in [("histogram_ct", "ct"), ("histogram_bin", "bin")]:
        log(f"Training {label} workflow (design-amortized)...")
        wf = build_workflow(key)
        hist = wf.fit_offline(data=raw_train, epochs=300, batch_size=64,
                              validation_data=raw_val, callbacks=[early_stop])
        f = bf.diagnostics.plots.loss(hist)
        f.savefig(OUTPUT_DIR / f"loss_{label}.png", dpi=150, bbox_inches="tight")
        plt.close(f)
        workflows[label] = wf

    rows = []
    for total_budget in TOTAL_BUDGET_GRID:
        for period in PERIOD_GRID:
            log(f"[budget={total_budget} period={period}] simulating validation set and estimating Delta...")
            val = make_val_at_design(N_VAL_PER_DESIGN, population, total_budget, period, rng)
            kl_ct,  _ = estimate_mutual_information(workflows["ct"],  val, "histogram_ct")
            kl_bin, _ = estimate_mutual_information(workflows["bin"], val, "histogram_bin")
            I_ct, I_bin = kl_ct.mean(), kl_bin.mean()
            delta = I_ct - I_bin
            se_delta = (kl_ct - kl_bin).std(ddof=1) / np.sqrt(len(kl_ct))
            rows.append({
                "total_budget": total_budget, "period": period, "frequency": 1.0 / period,
                "I_ct": I_ct, "I_bin": I_bin, "delta": delta, "se_delta": se_delta,
                "n_ever_infected": float(val["n_ever_infected"].mean()),
            })
            log(f"    I_ct={I_ct:.3f} I_bin={I_bin:.3f} delta={delta:.3f} +/- {se_delta:.3f} "
                f"| n_ever_infected~{val['n_ever_infected'].mean():.0f}")

    # ── Save + plot Delta(total_budget, period) surface ────────────────────────
    with open(OUTPUT_DIR / "design_sweep.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    with open(OUTPUT_DIR / "design_sweep.json", "w") as fh:
        json.dump(rows, fh, indent=2)

    def grid(field):
        g = np.full((len(TOTAL_BUDGET_GRID), len(PERIOD_GRID)), np.nan)
        for r in rows:
            g[TOTAL_BUDGET_GRID.index(r["total_budget"]), PERIOD_GRID.index(r["period"])] = r[field]
        return g

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, (field, title) in zip(axes, [("delta", "Delta(d) [nats]"),
                                         ("I_ct", "I_ct [nats]"),
                                         ("n_ever_infected", "epidemic size (ever infected)")]):
        im = ax.imshow(grid(field), origin="lower", aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(PERIOD_GRID))); ax.set_xticklabels(PERIOD_GRID)
        ax.set_yticks(range(len(TOTAL_BUDGET_GRID))); ax.set_yticklabels(TOTAL_BUDGET_GRID)
        ax.set_xlabel("period (days between surveillance events)")
        ax.set_ylabel("total_budget (tests/person)")
        ax.set_title(title)
        fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "design_sweep.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"Sweep complete. Delta surface -> {OUTPUT_DIR / 'design_sweep.png'}")


if __name__ == "__main__":
    main()
