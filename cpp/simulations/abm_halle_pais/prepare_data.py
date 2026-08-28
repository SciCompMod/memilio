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
"""Turn the RKI data downloaded by memilio-epidata into the three CSVs the Halle simulation reads.

The simulation needs daily values per age group, while the RKI files hold cumulative counts, so the
cumulative series are differenced here rather than in C++.

Produces, in the output directory:
  halle_cases.csv         date,age_group,new_cases   - reported cases, scaled by dark_figure in C++
  halle_vaccinations.csv  date,age_group,new_doses   - first doses, used unscaled
  halle_deaths.csv        date,deaths                - the fit target, cumulative

Age groups are the indices 0-5 of the six RKI groups (0-4, 5-14, 15-34, 35-59, 60-79, 80+).

Usage:
    python3 prepare_data.py --data-dir data/pydata/Germany --out-dir data/Halle
"""

import argparse
import csv
import json
import os

#: County id of Halle (Saale), Stadt.
HALLE_COUNTY_ID = 15002
#: State id of Sachsen-Anhalt, used to distribute county vaccination totals over age groups.
SACHSEN_ANHALT_STATE_ID = 15

#: Age group labels of the case files, in the order of the age group indices used by the simulation.
CASE_AGE_GROUPS = ['A00-A04', 'A05-A14', 'A15-A34', 'A35-A59', 'A60-A79', 'A80+']
#: Age group labels of the vaccination files, same order. Note that the oldest group is labelled
#: differently from the case files, where it is "A80+".
VACC_AGE_GROUPS = ['0-4', '5-14', '15-34', '35-59', '60-79', '80-99']


def load_json(path):
    with open(path) as handle:
        return json.load(handle)


def to_daily(cumulative_by_date):
    """Difference a {date: cumulative value} mapping into daily new values, clipped at zero.

    The RKI series are revised downwards occasionally, which would otherwise yield negative counts.
    """
    daily = {}
    previous = None
    for date in sorted(cumulative_by_date):
        value = cumulative_by_date[date]
        daily[date] = max(0.0, value - previous) if previous is not None else 0.0
        previous = value
    return daily


def write_csv(path, header, rows):
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    print(f'wrote {path} ({len(rows)} rows)')


def prepare_cases(data_dir, out_dir):
    """Write the daily new cases per age group and the daily deaths of Halle."""
    records = load_json(os.path.join(data_dir, 'cases_all_county_age_ma7.json'))
    records = [r for r in records if r['ID_County'] == HALLE_COUNTY_ID]
    if not records:
        raise SystemExit(f'No rows for county {HALLE_COUNTY_ID} in the case data.')

    case_rows = []
    deaths_by_date = {}
    for index, label in enumerate(CASE_AGE_GROUPS):
        by_date = {r['Date']: r['Confirmed'] for r in records if r['Age_RKI'] == label}
        deaths = {r['Date']: r['Deaths'] for r in records if r['Age_RKI'] == label}
        if not by_date:
            raise SystemExit(f'No rows for age group {label} in the case data.')
        for date, value in sorted(to_daily(by_date).items()):
            case_rows.append([date, index, round(value, 6)])
        # Deaths stay cumulative: the simulation also reports deaths cumulatively, and both sides are
        # rebased onto the start of the fit window, which is more stable than daily counts for one city.
        for date, value in deaths.items():
            deaths_by_date[date] = deaths_by_date.get(date, 0.0) + value

    write_csv(os.path.join(out_dir, 'halle_cases.csv'),
              ['date', 'age_group', 'new_cases'], case_rows)
    write_csv(os.path.join(out_dir, 'halle_deaths.csv'), ['date', 'deaths'],
              [[date, round(deaths_by_date[date], 6)] for date in sorted(deaths_by_date)])


def prepare_vaccinations(data_dir, out_dir):
    """Write the daily first vaccine doses per age group of Halle.

    Prefers the county-and-age file that memilio-epidata can produce. If only the county totals and the
    state-level age breakdown are available, the county total of each day is split over the age groups
    using that day's state age distribution.
    """
    county_age_path = os.path.join(data_dir, 'vacc_county_ageinf_ma7.json')
    rows = []

    if os.path.exists(county_age_path):
        records = [r for r in load_json(county_age_path)
                   if r['ID_County'] == HALLE_COUNTY_ID]
        for index, label in enumerate(VACC_AGE_GROUPS):
            by_date = {r['Date']: r['Vacc_partially']
                       for r in records if r['Age_RKI'] == label}
            for date, value in sorted(to_daily(by_date).items()):
                rows.append([date, index, round(value, 6)])
    else:
        print(f'{county_age_path} not found, splitting county totals by the state age distribution.')
        county = load_json(os.path.join(data_dir, 'vacc_county_ma7.json'))
        county_daily = to_daily({r['Date']: r['Vacc_partially'] for r in county
                                 if r['ID_County'] == HALLE_COUNTY_ID})
        state = [r for r in load_json(os.path.join(data_dir, 'vacc_states_ageinf_ma7.json'))
                 if r['ID_State'] == SACHSEN_ANHALT_STATE_ID]
        available = {r['Age_RKI'] for r in state}
        missing = [label for label in VACC_AGE_GROUPS if label not in available]
        if missing:
            raise SystemExit(f'Age groups {missing} are not in the state vaccination data, which has '
                             f'{sorted(available)}. They would silently receive zero vaccinations.')
        state_daily = {label: to_daily({r['Date']: r['Vacc_partially'] for r in state
                                        if r['Age_RKI'] == label})
                       for label in VACC_AGE_GROUPS}
        for date, total in sorted(county_daily.items()):
            shares = [state_daily[label].get(date, 0.0)
                      for label in VACC_AGE_GROUPS]
            total_share = sum(shares)
            if total_share <= 0.0:
                continue
            for index, share in enumerate(shares):
                rows.append([date, index, round(total * share / total_share, 6)])

    if not rows:
        raise SystemExit('No vaccination rows produced.')
    write_csv(os.path.join(out_dir, 'halle_vaccinations.csv'),
              ['date', 'age_group', 'new_doses'], rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data-dir', default=os.path.join('data', 'pydata', 'Germany'),
                        help='Directory holding the JSON files downloaded by memilio-epidata.')
    parser.add_argument('--out-dir', default=os.path.join('data', 'Halle'),
                        help='Directory to write the CSVs to.')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    prepare_cases(args.data_dir, args.out_dir)
    prepare_vaccinations(args.data_dir, args.out_dir)


if __name__ == '__main__':
    main()
