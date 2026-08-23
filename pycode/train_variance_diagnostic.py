import os

if "KERAS_BACKEND" not in os.environ:
    os.environ["KERAS_BACKEND"] = "jax"

from datetime import datetime

import numpy as np
import keras

from memilio.simulation.abm import ABMPopulation
from abm_inference import (
    build_workflow, estimate_mutual_information, make_val_at_design, total_budget_for_rate,
)
from run_design_sweep import SweepConfig, make_training_dataset

# ─────────────────────────────────────────────────────────────────────────────
CONFIG = SweepConfig(early_stop_patience=25, train_epochs=160)  
N_REPEATS = 5  
EVAL_POINTS = [(50, 1), (500, 1), (2000, 1)]  
# ─────────────────────────────────────────────────────────────────────────────


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    rng = np.random.default_rng(20260813)
    population = ABMPopulation(total_population=CONFIG.total_population)

    log(f"Simulating training data ONCE (n={CONFIG.n_train})...")
    raw_train = make_training_dataset(CONFIG.n_train, population, rng, CONFIG)
    log("Simulating early-stopping validation data ONCE...")
    raw_val = make_training_dataset(CONFIG.n_val_per_design * 2, population, rng, CONFIG)

    log(f"Simulating {len(EVAL_POINTS)} evaluation point(s) ONCE...")
    eval_sets = {}
    for rate, period in EVAL_POINTS:
        total_budget = total_budget_for_rate(rate)
        eval_sets[(rate, period)] = make_val_at_design(
            CONFIG.n_val_per_design, population, total_budget, period, rng,
            sources=CONFIG.sources, scale_by_tests=CONFIG.scale_by_tests)

    results = {point: {"ct": [], "bin": []} for point in EVAL_POINTS}
    for rep in range(N_REPEATS):
        log(f"=== training replicate {rep + 1}/{N_REPEATS} (fresh init, SAME data) ===")

        wf_ct = build_workflow("histogram_ct")
        early_stop = keras.callbacks.EarlyStopping(monitor="val_loss", patience=CONFIG.early_stop_patience,
                                                   restore_best_weights=True,
                                                   start_from_epoch=CONFIG.early_stop_start_from_epoch)
        wf_ct.fit_offline(data=raw_train, epochs=CONFIG.train_epochs, batch_size=64,
                          validation_data=raw_val, callbacks=[early_stop], verbose=0)

        wf_bin = build_workflow("histogram_bin")
        early_stop = keras.callbacks.EarlyStopping(monitor="val_loss", patience=CONFIG.early_stop_patience,
                                                   restore_best_weights=True,
                                                   start_from_epoch=CONFIG.early_stop_start_from_epoch)
        wf_bin.fit_offline(data=raw_train, epochs=CONFIG.train_epochs, batch_size=64,
                           validation_data=raw_val, callbacks=[early_stop], verbose=0)

        for point in EVAL_POINTS:
            val = eval_sets[point]
            kl_ct,  _ = estimate_mutual_information(wf_ct,  val, "histogram_ct")
            kl_bin, _ = estimate_mutual_information(wf_bin, val, "histogram_bin")
            results[point]["ct"].append(kl_ct.mean())
            results[point]["bin"].append(kl_bin.mean())
            log(f"    rate={point[0]:>5} period={point[1]}: I_ct={kl_ct.mean():.3f} "
                f"I_bin={kl_bin.mean():.3f} delta={kl_ct.mean() - kl_bin.mean():.3f}")

    log("=== summary: Delta variance from training stochasticity ALONE (same data every time) ===")
    for point in EVAL_POINTS:
        ct = np.array(results[point]["ct"])
        bin_ = np.array(results[point]["bin"])
        deltas = ct - bin_
        log(f"  rate={point[0]:>5} period={point[1]}: delta = {deltas.mean():.4f} +/- "
            f"{deltas.std(ddof=1):.4f}  (range [{deltas.min():.4f}, {deltas.max():.4f}], n={N_REPEATS})")


if __name__ == "__main__":
    main()
