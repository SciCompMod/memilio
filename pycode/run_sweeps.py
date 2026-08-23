"""Run run_design_sweep.main() for a batch of SweepConfigs back to back, in one process.

Edit CONFIGS below to list the designs (experiment configurations) you want to compare, then:
    python run_sweeps.py

Each config gets its own tagged, timestamped results folder under design_sweep_results/ (see
SweepConfig.tag() / run_design_sweep.main()) -- nothing here needs to track that itself.

Runs sequentially, on purpose: run_batch_designs() already saturates every CPU core for a single
sweep's simulation phase, so running configs concurrently would only contend for the same cores,
not speed anything up.

A failure in one config is logged and does not stop the rest of the batch -- see the summary
printed at the end for which configs actually finished.
"""
import traceback
from datetime import datetime

import run_design_sweep as rds
from run_design_sweep import SweepConfig
from abm_inference import SOURCES_SURVEY, SOURCES_DIAGNOSTIC, SOURCES_POOLED

CONFIGS = [
    SweepConfig(sources=SOURCES_SURVEY, scale_by_tests=False),
    SweepConfig(sources=SOURCES_SURVEY, scale_by_tests=False),
    SweepConfig(sources=SOURCES_SURVEY, scale_by_tests=True),
    SweepConfig(sources=SOURCES_SURVEY, scale_by_tests=True)
    # SweepConfig(sources=SOURCES_SURVEY, scale_by_tests=False),          # rate vs raw counts
    # SweepConfig(design_rate_lo=100, design_rate_hi=1000),               # a narrower budget range
    # SweepConfig(n_train=60000, train_epochs=150),                       # more training data
]


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    log(f"Running {len(CONFIGS)} config(s)...")
    results = []
    for i, config in enumerate(CONFIGS):
        log(f"=== config {i + 1}/{len(CONFIGS)}: {config.tag()} ===")
        try:
            output_dir = rds.main(config)
            results.append((config.tag(), output_dir, None))
            log(f"    done -> {output_dir}")
        except Exception as exc:
            results.append((config.tag(), None, exc))
            log(f"    FAILED: {exc!r}")
            traceback.print_exc()

    log("=== batch summary ===")
    for tag, output_dir, exc in results:
        status = f"-> {output_dir}" if exc is None else f"FAILED: {exc!r}"
        log(f"  {tag}: {status}")


if __name__ == "__main__":
    main()
