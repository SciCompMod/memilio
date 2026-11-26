#!/usr/bin/env python3
"""Fit the contact weights of the location split simulation to an observed distribution of infection sites.

The contact rate of each location type is a heuristic factor times a weight (--contact_weights). This
script adjusts the weights of the fitted location types and the scale of the viral shed k until

  1. the share of infections per fitted location type, among infections of persons aged 15+ in those
     location types, matches the target shares, and
  2. the epidemic grows as fast as a reference: by default the original heuristic model (HEURISTIC
     weights, the given k), so that rebalancing the settings does not change the overall speed.

Default targets: the French ComCor study (Galmiche et al., BMC Public Health 2024, 584,846 adults,
October 2020 to August 2022). Known or suspected sources were household members 30.8%, extended family
15.6%, work colleagues 15.0%, friends 11.0% (other/multiple 27.6% left out). Mapped to the ABM and
renormalized: Home 42.5%, SocialEvent 36.7% (extended family + friends), Work 20.7%.

Identifiability: three shares (two independent) plus the growth rate pin down three numbers. The script
fits the home and social event weights and k; work, school and shop keep their weights as the
reference. ComCor has no data for schools and shops, so their weights stay relative to work.

Each evaluation runs short simulations (default 21 days) and measures in the exponential phase, from
day 7 on, so that the random placement of the initial infections has washed out. Shares and growth use
whole weeks because of the weekly rhythm of the mobility (school on weekdays, events on weekends). All
evaluations use the same seed, so differences between iterations come from the parameters, not noise.

Household transmission is close to saturation in this model, so the shares respond nonlinearly to the
weights. The update is a damped step in log space whose sensitivity (d log share / d log weight) is
estimated from the previous iterations, a secant method. The sensitivity is bounded below and the
steps are limited, so that a weight the shares barely respond to (the home weight, because households
saturate) is not pushed far away from its start.

Usage:
    python fit_contact_weights.py --binary ../../../../build/bin/location_split_simulation
    python fit_contact_weights.py --binary ... --target home=0.5,social_event=0.3,work=0.2
"""

import argparse
import glob
import json
import math
import os
import subprocess

import numpy as np
import pandas as pd

LOCATION_TYPE = {"home": "Home", "work": "Work", "social_event": "SocialEvent", "school": "School",
                 "shop": "BasicsShop"}
ALL_WEIGHTS = ["home", "school", "work", "social_event", "shop"]
# Original heuristic factors of the panvXabmSim branch, expressed as weights on top of 24 h / hours per
# visit day (ContactWeights::heuristic() in parameter_setter.h). The fit starts from these.
HEURISTIC = {"home": 16.8875, "school": 1.6, "work": 0.625, "social_event": 2.4, "shop": 0.1584}
COMCOR_TARGET = {"home": 30.8, "social_event": 15.6 + 11.0, "work": 15.0}
ADULT_AGE_GROUPS = [2, 3, 4, 5]  # 15-34, 35-59, 60-79, 80+


def parse_mapping(text):
    mapping = {}
    for entry in text.split(","):
        key, value = entry.split("=")
        if key not in LOCATION_TYPE:
            raise ValueError(f"Unknown location type '{key}', expected one of {list(LOCATION_TYPE)}.")
        mapping[key] = float(value)
    return mapping


def normalized(mapping):
    total = sum(mapping.values())
    return {key: value / total for key, value in mapping.items()}


def weights_string(weights):
    return ",".join(f"{key}={weights[key]:.6g}" for key in ALL_WEIGHTS)


def simulate(args, weights, k, directory):
    """Run the simulation and return (shares of the target types, growth rate per day, all-age shares)."""
    if not os.path.exists(os.path.join(directory, "summary.txt")):
        command = [args.binary, "--runs", str(args.runs), "--days", str(args.days),
                   "--n_persons", str(args.n_persons), "--initial_infections", str(args.initial_infections),
                   "--infection_k", f"{k:.6g}", "--seed", str(args.seed),
                   "--contact_weights", weights_string(weights), "--output_dir", directory]
        if args.threads:
            command += ["--threads", str(args.threads)]
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL)

    events = pd.concat([pd.read_csv(path) for path in glob.glob(os.path.join(directory, "runs", "run_*_infection_events.csv"))])
    events = events[events["Location_Type"] != "Initial"]
    day = (events["Timestep"] - 1) // 24

    window = events[day >= args.burn_in]
    adults = window[window["Age_Group"].isin(ADULT_AGE_GROUPS)]
    counts = {key: int((adults["Location_Type"] == LOCATION_TYPE[key]).sum()) for key in args.target}
    total = sum(counts.values())
    if total == 0:
        raise RuntimeError("No adult infections in the target settings; increase --days or --initial_infections.")
    shares = {key: count / total for key, count in counts.items()}

    # Growth from the last two whole weeks, pooled over runs.
    last_week = ((day >= args.days - 7) & (day < args.days)).sum()
    previous_week = ((day >= args.days - 14) & (day < args.days - 7)).sum()
    growth = math.log(max(last_week, 1) / max(previous_week, 1)) / 7.0

    all_shares = window["Location_Type"].value_counts(normalize=True).to_dict()
    return shares, growth, all_shares, total


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary", required=True, help="Path to location_split_simulation.")
    parser.add_argument("--target", type=parse_mapping, default=normalized(COMCOR_TARGET),
                        help="Target shares, e.g. home=0.425,social_event=0.367,work=0.207 (normalized).")
    parser.add_argument("--fit", default="home,social_event",
                        help="Weights to fit; the others stay fixed as the reference.")
    parser.add_argument("--start_weights", type=parse_mapping, default={},
                        help="Starting weights; missing ones are the heuristic values.")
    parser.add_argument("--k", type=float, default=2.4, help="Starting scale of the viral shed.")
    parser.add_argument("--target_growth", type=float, default=None,
                        help="Target growth rate per day. Default: that of the heuristic weights with --k.")
    parser.add_argument("--n_persons", type=int, default=10000)
    parser.add_argument("--runs", type=int, default=20)
    parser.add_argument("--days", type=int, default=21, help="Simulated days; at least --burn_in + 14.")
    parser.add_argument("--burn_in", type=int, default=7, help="Days left out of the share measurement.")
    parser.add_argument("--initial_infections", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--max_iterations", type=int, default=10)
    parser.add_argument("--share_tolerance", type=float, default=0.02, help="Absolute, per share.")
    parser.add_argument("--growth_tolerance", type=float, default=0.005, help="Per day.")
    parser.add_argument("--generation_time", type=float, default=5.0,
                        help="Days; first guess of how growth responds to k (r changes by log(k ratio) / T).")
    parser.add_argument("--output", default="calibration_output")
    args = parser.parse_args()

    args.target = normalized(args.target)
    fitted = [key for key in args.fit.split(",") if key]
    if args.days < args.burn_in + 14:
        parser.error("--days must be at least --burn_in + 14.")
    for key in fitted:
        if key not in args.target:
            parser.error(f"Fitted weight '{key}' has no target share.")
    os.makedirs(args.output, exist_ok=True)

    weights = dict(HEURISTIC)
    weights.update(args.start_weights)
    k = args.k

    if args.target_growth is None:
        _, args.target_growth, _, _ = simulate(args, dict(HEURISTIC), args.k,
                                               os.path.join(args.output, "reference"))
    print(f"Target shares (adults): {', '.join(f'{key} {value:.3f}' for key, value in args.target.items())}")
    print(f"Target growth rate: {args.target_growth:.4f} per day\n")

    # Sensitivities, start with the linear guess.
    elasticity = {key: 1.0 for key in fitted}
    growth_per_log_k = 1.0 / args.generation_time
    history = []
    previous = None

    for iteration in range(args.max_iterations):
        shares, growth, all_shares, n_events = simulate(args, weights, k, os.path.join(args.output, f"iteration_{iteration}"))
        error = max(abs(shares[key] - args.target[key]) for key in args.target)
        history.append({"iteration": iteration, "weights": dict(weights), "k": k, "shares": shares,
                        "growth": growth, "all_age_shares": all_shares, "adult_events": n_events,
                        "max_share_error": error})
        print(f"iteration {iteration}: k={k:.4g}, {weights_string(weights)}")
        print(f"    adult shares {', '.join(f'{key} {shares[key]:.3f}' for key in args.target)} "
              f"(max error {error:.3f}, {n_events} events), growth {growth:.4f}/day")

        if error <= args.share_tolerance and abs(growth - args.target_growth) <= args.growth_tolerance:
            print("Converged.")
            break

        # Secant update of the sensitivities from the last step, kept in a sane range.
        if previous is not None:
            for key in fitted:
                step = math.log(weights[key] / previous["weights"][key])
                if abs(step) > 1e-6:
                    change = math.log(shares[key] / previous["shares"][key])
                    elasticity[key] = float(np.clip(change / step, 0.3, 2.0))
            log_k_step = math.log(k / previous["k"])
            if abs(log_k_step) > 1e-6:
                growth_per_log_k = float(np.clip((growth - previous["growth"]) / log_k_step, 0.05, 0.5))
        previous = {"weights": dict(weights), "shares": shares, "growth": growth, "k": k}

        for key in fitted:
            log_step = math.log(args.target[key] / shares[key]) / elasticity[key]
            weights[key] *= math.exp(float(np.clip(log_step, -0.5, 0.5)))
        log_k = (args.target_growth - growth) / growth_per_log_k
        k *= math.exp(float(np.clip(log_k, -0.3, 0.3)))
    else:
        print("Not converged within --max_iterations; the last evaluated values are written.")
    weights, k = history[-1]["weights"], history[-1]["k"]

    result = {"weights": weights, "infection_k": k, "contact_weights_option": weights_string(weights),
              "target_shares": args.target, "target_growth": args.target_growth, "fitted": fitted,
              "settings": {key: getattr(args, key) for key in
                           ["n_persons", "runs", "days", "burn_in", "initial_infections", "seed"]},
              "history": history}
    path = os.path.join(args.output, "fit_result.json")
    with open(path, "w") as file:
        json.dump(result, file, indent=2)
    print(f"\nWrote {path}")
    print(f"Use: --infection_k {k:.4g} --contact_weights {weights_string(weights)}")


if __name__ == "__main__":
    main()
