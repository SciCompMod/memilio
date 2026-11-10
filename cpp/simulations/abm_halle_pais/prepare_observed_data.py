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
"""Build the observed data the Halle ABM is fitted to.

Writes one row per day of the fit window with the real counts for Halle (Saale),
county 15002, in the same channels as the model output (halle_model.h):

  deaths      cumulative deaths since the start of the window. RKI dates a death
              by the reference date of its case, so the series is shifted by
              DEATH_DELAY_DAYS, as in paper_abm_bs_testing.cpp.
  admissions  new hospital admissions. Only published per state, so the
              Sachsen-Anhalt series is scaled by the population share of Halle.
              RKI dates an admission by the reporting date of its case.
  icu         COVID-19 patients in intensive care (DIVI), Halle.
  reported    newly reported cases by reporting date, Halle.

All values are daily and unsmoothed; the fit compares 7-day sums of the flows.
A JSON file with the populations needed to scale the model output is written
next to the CSV.

Input files are the outputs of the memilio-epidata console scripts, run with
-o data/pydata (they write into its Germany/ subfolder):
    getcasedata --rep-date            -> cases_all_county_repdate.json
    getcasedata                       -> cases_all_county.json
    getdividata                       -> county_divi.json
    gethospitalizationdata            -> hospit_state.json
    getpopuldata                      -> county_current_population.json
The *_ma1 variants of the case files are accepted as well.
"""
import argparse
import csv
import datetime
import json
import os

HALLE_COUNTY_ID = 15002
SACHSEN_ANHALT_STATE_ID = 15
#: Days from the reference date of a case to the death, see paper_abm_bs_testing.cpp.
DEATH_DELAY_DAYS = 18


def load(data_dir, *names):
    """Load the first of the given JSON files that exists in data_dir."""
    for name in names:
        path = os.path.join(data_dir, name)
        if os.path.exists(path):
            with open(path) as handle:
                return name, json.load(handle)
    raise SystemExit(f'None of {", ".join(names)} found in {data_dir}. '
                     'Generate them with memilio-epidata, see the docstring of this script.')


def series(records, key, value, filter_key, filter_value, source):
    """Extract a {date: value} mapping for one county or state."""
    out = {}
    for record in records:
        if record.get(filter_key) == filter_value:
            out[datetime.date.fromisoformat(record[key][:10])] = float(record[value])
    if not out:
        raise SystemExit(f'{source} has no entries for {filter_key} {filter_value}.')
    return out


def at(values, date, source):
    if date not in values:
        raise SystemExit(f'{source} has no entry for {date}; it covers {min(values)} to {max(values)}. '
                         'Download more recent data.')
    return values[date]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--data-dir', default='data/pydata/Germany',
                        help='Directory holding the memilio-epidata output.')
    parser.add_argument('--start-date', default='2022-07-01')
    parser.add_argument('--num-days', type=int, default=90)
    parser.add_argument('--out', default='data/Germany/halle/halle_observed.csv')
    args = parser.parse_args()

    start = datetime.date.fromisoformat(args.start_date)
    dates = [start + datetime.timedelta(days=d) for d in range(args.num_days)]

    name, records = load(args.data_dir, 'county_current_population.json')
    population = {r['ID_County']: r['Population'] for r in records}
    halle_population = population[HALLE_COUNTY_ID]
    state_population = sum(p for county, p in population.items()
                           if county // 1000 == SACHSEN_ANHALT_STATE_ID)

    name, records = load(args.data_dir, 'cases_all_county_repdate.json', 'cases_all_county_repdate_ma1.json')
    confirmed = series(records, 'Date', 'Confirmed', 'ID_County', HALLE_COUNTY_ID, name)
    reported = [at(confirmed, d, name) - at(confirmed, d - datetime.timedelta(days=1), name) for d in dates]

    name, records = load(args.data_dir, 'cases_all_county.json', 'cases_all_county_ma1.json')
    deaths_by_refdate = series(records, 'Date', 'Deaths', 'ID_County', HALLE_COUNTY_ID, name)
    delay = datetime.timedelta(days=DEATH_DELAY_DAYS)
    deaths_at_start = at(deaths_by_refdate, start - delay - datetime.timedelta(days=1), name)
    deaths = [at(deaths_by_refdate, d - delay, name) - deaths_at_start for d in dates]

    name, records = load(args.data_dir, 'county_divi.json')
    icu_by_date = series(records, 'Date', 'ICU', 'ID_County', HALLE_COUNTY_ID, name)
    icu = [at(icu_by_date, d, name) for d in dates]

    name, records = load(args.data_dir, 'hospit_state.json')
    state_admissions = series(records, 'Date', 'Hospitalized', 'ID_State', SACHSEN_ANHALT_STATE_ID, name)
    share = halle_population / state_population
    admissions = [share * at(state_admissions, d, name) for d in dates]

    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    with open(args.out, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['date', 'day', 'deaths', 'admissions', 'icu', 'reported'])
        for day, row in enumerate(zip(dates, deaths, admissions, icu, reported)):
            writer.writerow([row[0].isoformat(), day] + [round(v, 4) for v in row[1:]])

    meta = {'start_date': args.start_date, 'num_days': args.num_days,
            'halle_population': halle_population, 'sachsen_anhalt_population': state_population,
            'death_delay_days': DEATH_DELAY_DAYS}
    meta_path = os.path.splitext(args.out)[0] + '_meta.json'
    with open(meta_path, 'w') as handle:
        json.dump(meta, handle, indent=2)

    print(f'Wrote {args.out} and {meta_path}: {args.num_days} days from {start}, '
          f'{deaths[-1]:.0f} deaths, {sum(admissions):.0f} admissions, '
          f'ICU {min(icu):.0f}-{max(icu):.0f}, {sum(reported):.0f} reported cases.')


if __name__ == '__main__':
    main()
