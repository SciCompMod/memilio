#!/usr/bin/env python3
#############################################################################
# Copyright (C) 2020-2026 MEmilio
#
# Authors: Sascha Korf
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
"""Neural posterior estimation for the Halle ABM.

Trains a neural posterior estimator on the ensemble that abm_halle_pais wrote in "ensemble" mode, then
evaluates the posterior at the observed deaths. This replaces the grid search of the ABM paper: instead
of one point estimate on a 6-per-dimension grid, it returns a full posterior over all fitted parameters
from the same simulation budget.

Only the channels that have real data enter the fit. The PAIS channels are carried through the ensemble
as a forward projection and are deliberately excluded here: with no PAIS data, the posterior would be
flat in any PAIS parameter, so those are not fitted at all (see fit_parameters() on the C++ side).

Requires torch and sbi:
    pip install sbi

Usage:
    # 1. write the prior bounds and generate the ensemble (the expensive step, run under MPI)
    abm_halle_pais --mode '"priors"'   --output_dir '"out"'
    mpirun -n 64 abm_halle_pais --mode '"ensemble"' --num_runs 3000 --num_days 90 ... --output_dir '"out"'
    # 2. train the estimator and sample the posterior
    python3 fit_npe.py --ensemble-dir out --deaths data/Halle/halle_deaths.csv \
        --start-date 2022-07-01 --num-days 90 --out out
"""

import argparse
import csv
import datetime
import glob
import os


def read_priors(path):
    """Read the prior bounds written by the "priors" mode."""
    names, lower, upper = [], [], []
    with open(path) as handle:
        for row in csv.DictReader(handle):
            names.append(row['name'])
            lower.append(float(row['lower']))
            upper.append(float(row['upper']))
    if not names:
        raise SystemExit(f'{path} holds no parameters.')
    return names, lower, upper


def read_ensemble(ensemble_dir, parameter_names, fit_channel_prefix):
    """Read every ensemble shard and split it into parameters and the fitted part of the observable.

    Returns (theta, x, channel_columns) where theta has one row of parameters per run, x the matching
    values of the channels whose name starts with fit_channel_prefix, and channel_columns their names.
    """
    shards = sorted(glob.glob(os.path.join(ensemble_dir, 'ensemble_rank*.csv')))
    if not shards:
        raise SystemExit(f'No ensemble_rank*.csv in {ensemble_dir}.')

    theta, x, channel_columns = [], [], None
    for shard in shards:
        with open(shard) as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                continue
            columns = [name for name in reader.fieldnames
                       if name.startswith(fit_channel_prefix)]
            if channel_columns is None:
                channel_columns = columns
            elif columns != channel_columns:
                raise SystemExit(f'{shard} has different channel columns than the earlier shards.')
            for row in reader:
                theta.append([float(row[name]) for name in parameter_names])
                x.append([float(row[name]) for name in columns])

    if not theta:
        raise SystemExit(f'The shards in {ensemble_dir} hold no runs.')
    if not channel_columns:
        raise SystemExit(f'No columns starting with "{fit_channel_prefix}" in the ensemble.')
    print(f'read {len(theta)} runs from {len(shards)} shard(s), '
          f'{len(parameter_names)} parameters, {len(channel_columns)} observable values')
    return theta, x, channel_columns


def read_observed_deaths(path, start_date, num_days):
    """Read the real cumulative deaths and rebase them onto the start of the fit window.

    The simulation reports deaths that occur after its start, so the data has to be rebased the same way,
    by subtracting the cumulative count on the start date.
    """
    cumulative = {}
    with open(path) as handle:
        for row in csv.DictReader(handle):
            cumulative[datetime.date.fromisoformat(row['date'])] = float(row['deaths'])

    if start_date not in cumulative:
        raise SystemExit(f'{path} has no entry for the start date {start_date}. It covers '
                         f'{min(cumulative)} to {max(cumulative)}.')
    baseline = cumulative[start_date]

    observed = []
    for day in range(num_days + 1):
        date = start_date + datetime.timedelta(days=day)
        if date not in cumulative:
            raise SystemExit(f'{path} has no entry for {date}, which is day {day} of the fit window. '
                             f'It covers {min(cumulative)} to {max(cumulative)}.')
        observed.append(cumulative[date] - baseline)
    return observed


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--ensemble-dir', required=True,
                        help='Directory holding the ensemble shards and priors.csv.')
    parser.add_argument('--deaths', required=True,
                        help='CSV with the real cumulative deaths, from prepare_data.py.')
    parser.add_argument('--start-date', required=True,
                        help='First day of the fit window, matching the ensemble.')
    parser.add_argument('--num-days', type=int, required=True,
                        help='Length of the fit window in days, matching the ensemble.')
    parser.add_argument('--out', required=True, help='Directory to write the posterior to.')
    parser.add_argument('--num-posterior-samples', type=int, default=10000,
                        help='Number of samples drawn from the posterior.')
    parser.add_argument('--fit-channel-prefix', default='deaths_',
                        help='Only observable columns with this prefix enter the fit.')
    args = parser.parse_args()

    try:
        import torch
        from sbi.inference import NPE
        from sbi.utils import BoxUniform
    except ImportError as error:
        raise SystemExit(f'This script needs torch and sbi ({error}). Install them with: pip install sbi')

    names, lower, upper = read_priors(os.path.join(args.ensemble_dir, 'priors.csv'))
    theta, x, channel_columns = read_ensemble(args.ensemble_dir, names, args.fit_channel_prefix)
    observed = read_observed_deaths(
        args.deaths, datetime.date.fromisoformat(args.start_date), args.num_days)

    if len(observed) != len(channel_columns):
        raise SystemExit(f'The ensemble has {len(channel_columns)} fitted values per run but the data '
                         f'gives {len(observed)}. Check that --num-days matches the ensemble.')

    theta_tensor = torch.tensor(theta, dtype=torch.float32)
    x_tensor = torch.tensor(x, dtype=torch.float32)
    observed_tensor = torch.tensor([observed], dtype=torch.float32)

    prior = BoxUniform(low=torch.tensor(lower, dtype=torch.float32),
                       high=torch.tensor(upper, dtype=torch.float32))

    # A single round of NPE. The ensemble is drawn from the prior, so it is exactly the training set the
    # first round of NPE expects. For a sequential fit, sample the next round's parameters from the
    # posterior below, run the ensemble mode on them, and append the shards.
    inference = NPE(prior=prior)
    inference.append_simulations(theta_tensor, x_tensor)
    print('training the density estimator ...')
    density_estimator = inference.train()
    posterior = inference.build_posterior(density_estimator)

    samples = posterior.sample((args.num_posterior_samples,), x=observed_tensor)

    os.makedirs(args.out, exist_ok=True)
    samples_path = os.path.join(args.out, 'posterior_samples.csv')
    with open(samples_path, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        writer.writerows(samples.tolist())
    print(f'wrote {samples_path}')

    print('\nposterior marginals:')
    print(f'{"parameter":<28}{"median":>10}{"2.5%":>10}{"97.5%":>10}{"prior":>18}')
    for index, name in enumerate(names):
        column = samples[:, index]
        quantiles = torch.quantile(column, torch.tensor([0.025, 0.5, 0.975]))
        print(f'{name:<28}{quantiles[1]:>10.4f}{quantiles[0]:>10.4f}{quantiles[2]:>10.4f}'
              f'{f"[{lower[index]}, {upper[index]}]":>18}')

    # A marginal that still fills its prior range is the signature of a parameter the data cannot
    # constrain, which is worth knowing before reading anything into the fit.
    for index, name in enumerate(names):
        column = samples[:, index]
        width = (torch.quantile(column, torch.tensor(0.975))
                 - torch.quantile(column, torch.tensor(0.025))).item()
        if width > 0.9 * (upper[index] - lower[index]):
            print(f'\nWARNING: the posterior of {name} covers almost its whole prior range. '
                  f'The deaths data does not constrain it.')


if __name__ == '__main__':
    main()
