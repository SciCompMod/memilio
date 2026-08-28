# Halle ABM with PAIS, fitted by simulation based inference

The Halle agent based model with post-acute infection syndrome (PAIS), and the simulator side of a
neural posterior estimation (NPE) fit against real data.

The split is deliberate: the C++ binary is the simulator and generates the ensemble under MPI, which is
where essentially all of the compute sits. The density estimator lives in Python, because that is where
the mature normalizing flow implementations are.

```
prepare_data.py   RKI JSON            -> the three CSVs the simulation reads
abm_halle_pais    population + data   -> observable time series, or an ensemble of (parameters, observable)
fit_npe.py        ensemble + deaths   -> posterior over the fitted parameters
```

## What is fitted

Four parameters, with the prior bounds of the grid search in the ABM paper
(`cpp/simulations/paper_abm_bs_testing.cpp` on branch `abm_paper_test_bs`):

| Parameter | Prior | Paper optimum | Acts on |
| --- | --- | --- | --- |
| `viral_shedding_rate` | 1.4 – 2.0 | 1.52 | `InfectionRateFromViralShed` |
| `dark_figure` | 2.5 – 5.5 | 4.3 | scale of the seeded infection history |
| `testing_probability_sympt` | 0.01 – 0.045 | 0.024 | testing scheme for symptomatic persons |
| `ratio_asympt_to_sympt` | 2 – 15 | 4.6 | testing scheme for asymptomatic persons |

The paper's fifth parameter, `contact_red_lockdown`, is deliberately absent: testing is the only measure
in this setup, so there is no contact reduction to fit.

`fit_parameters()` in `halle_model.h` is the single place where the fitted parameters and their priors are
defined. Adding one there changes the dimension of the fit with no other code change, on either side.

### What is not fitted, and why

The PAIS parameters are held at the values of `cpp/examples/abm_aims_halle.cpp`. They act only on the PAIS
output channels, so without PAIS data the posterior would be flat in every one of them and the fit would
return the prior while looking like a result. PAIS is therefore a forward projection here: it is simulated
and reported, but not inferred. `fit_npe.py` prints a warning for any parameter whose posterior still
covers its prior range, which is the empirical version of this check.

## Data

| File | Columns | Source |
| --- | --- | --- |
| `halle_population_data.csv` | `age,home_id,school_id,event_id,shopping_id,work_id` | provided; `age` is an index into an 11-group scheme |
| `halle_cases.csv` | `date,age_group,new_cases` | `prepare_data.py` |
| `halle_vaccinations.csv` | `date,age_group,new_doses` | `prepare_data.py` |
| `halle_deaths.csv` | `date,deaths` (cumulative) | `prepare_data.py` |

`age_group` is the index 0–5 of the six RKI groups (0-4, 5-14, 15-34, 35-59, 60-79, 80+). The population
file's `age` column uses a different, finer scheme, which `age_group_from_input()` aggregates onto these
six; every parameter value in the simulation is given for the six.

Download the RKI inputs with `memilio-epidata` first. The county-and-age vaccination file
(`vacc_county_ageinf_ma7.json`) is preferred; if only the county totals and the state age breakdown are
available, `prepare_data.py` splits the county total by the state age distribution and says so.

## Running

The memilio CLI parses values as JSON, so string arguments need quotes inside quotes. Use `--help` for the
full list of options, and `--read_from_json` to keep a run's configuration in a file.

Write the prior bounds:

```bash
abm_halle_pais --mode '"priors"' --output_dir '"out"'
```

A single run:

```bash
abm_halle_pais --mode '"simulate"' \
  --start_date '"2022-07-01"' --num_days 90 \
  --cases_file '"data/Halle/halle_cases.csv"' \
  --vaccinations_file '"data/Halle/halle_vaccinations.csv"' \
  --output_dir '"out"'
```

The ensemble, which is the expensive step. Runs are handed out by index modulo rank, so no communication
is needed and each rank writes its own shard:

```bash
mpirun -n 64 abm_halle_pais --mode '"ensemble"' --num_runs 3000 \
  --start_date '"2022-07-01"' --num_days 90 \
  --cases_file '"data/Halle/halle_cases.csv"' \
  --vaccinations_file '"data/Halle/halle_vaccinations.csv"' \
  --output_dir '"out"'
```

Then the fit:

```bash
pip install sbi
python3 fit_npe.py --ensemble-dir out --deaths data/Halle/halle_deaths.csv \
  --start-date 2022-07-01 --num-days 90 --out out
```

For a sequential fit, sample the next round's parameters from the posterior, run `ensemble` mode on them
and append the shards.

## Notes and caveats

- **Simulation length is a plain option.** `--num_days` and `--start_date` set the window; nothing about
  the length is baked in. A 90 day window costs roughly 90 s per run at 286k agents on one core, so an
  ensemble of 3000 is on the order of 75 core-hours. A full run to the end of 2025 (about 1280 days) is
  roughly 20 minutes, which is fine for posterior predictive draws but not for an ensemble.
- **Deaths are rebased onto the start of the window** on both sides. The seeded infection history kills
  part of the population before `t0`, so the raw cumulative count is dominated by deaths that already
  happened. `fit_npe.py` rebases the real data the same way; the two must stay consistent.
- **Deaths alone are a weak signal for a single city, and a weak signal for two of the four parameters.**
  `testing_probability_sympt` and `ratio_asympt_to_sympt` act on detection, which affects deaths only
  indirectly through quarantine. If their posteriors come back flat, add a detected-cases channel to
  `observable_channels()` and the `LogChannels` logger, and point `--fit-channel-prefix` at it.
- **A missing input file is an error, not a fallback.** `allow_missing_history` exists only for smoke
  tests: without the histories the population has no immunity and no pre-existing PAIS at `t0`, which
  silently changes the epidemic being fitted.
- **The PAIS transition rates are not the colleague's matrix.** `PAISTransitionMatrix` is consumed as
  rates in 1/day, while the matrix in `abm_aims_halle.cpp` is documented as per-day probabilities, and its
  medium row (`Medium->Severe` 0.99998/day) empties the medium state within a day under either reading.
  `set_pais_parameters()` uses plausible literature durations instead and says so; the semantics of the
  original numbers need to be settled before they are used.
- **`PAISProbabilitySeverityFactor` is below 1** (0.267 / 0.247 / 0.217) while `PAIS::init_or_refresh`
  applies it as an increase for severe acute infections. Either the values or the direction is wrong.
