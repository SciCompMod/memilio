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

**No intervention of any kind is modelled at present** - no contact reduction, and no testing either. The
paper's `contact_red_lockdown` is therefore absent, and so are its two testing parameters
(`testing_probability_sympt`, `ratio_asympt_to_sympt`): with no testing scheme in the model they would act
on nothing and the fit would treat them as pure noise dimensions. `add_testing_strategy()` and the two
parameters are recoverable from git history when testing is reinstated; `make_model` still takes `tmax`
for exactly that reason.

Measured while testing was still in: sweeping both testing parameters across their full priors moved
deaths by only about 10% (2183 to 1967 over 90 days), so they were weakly identified by deaths even then.

`MaximumContacts` is not set on any location. It rescales each `ContactRates` row to sum to at most its
value, so setting it would partly undo the contact matrices of `set_contact_rates()`;
`set_local_parameters_ger()` on `abmXpanvadere` likewise sets contact rates and never caps them. Note that
the previous cap of 5 on `SocialEvent` was strongly binding - removing it raised 90-day deaths from about
2100 to 3516.

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
| `priors.csv` | `name,lower,upper` | `abm_halle_pais --mode '"priors"'` |

The RKI inputs in `data/pydata/Germany` were refreshed to cover 2021-01-01 – 2022-12-01 (vaccinations from
2020-12-01, which is the start of the campaign), so the 2022-07-01 fit window and its history lookback are
covered. `vacc_county_ageinf_ma7.json` is now present, so the state-level split fallback is not used.

`halle_cases.csv` is an **input, not a fit target**: `dark_figure` times the reported cases is what seeds
the pre-`t0` infection history, and therefore the immunity and the pre-existing PAIS at `t0`. It is needed
even though no case channel enters the objective.

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

## Open issue: the prior box cannot reach the observed deaths

**This is unresolved and blocks a meaningful fit.** With the four parameters above, no point in the prior
box reproduces the deaths observed in the target window. The smallest epidemic the box allows, at the
corner `1.4, 5.5, 0.045, 2.0`, still gives **1584 deaths over 2022-07-01 to 2022-12-01 against 67
observed**; the centre of the prior gives 2182. Removing the contact caps and the testing schemes moved
the centre further up, and no longer splitting oversized locations moved it further still: 4062 deaths
over the first 90 days against roughly 40 observed. An NPE fit run
in this state returns a posterior pressed against the edge of the prior, which looks like a result but is
not one.

The likely cause is that the severity cascade is Alpha-era: `SeverePerInfectedSymptoms`,
`CriticalPerInfectedSevere` and `DeathsPerInfectedCritical` are the values the ABM paper calibrated on
Braunschweig in spring 2021, when Alpha met a largely non-immune population and the infection fatality
rate was on the order of a percent. Halle in the second half of 2022 saw BA.5 meet a population with
hybrid immunity, at roughly a tenth of that.

A scaling factor on the entry into the severe branch brackets the target cleanly, which is evidence for
this diagnosis: 0.02 gives 43 deaths and 0.05 gives 105, against 67 observed, so the implied factor is
about 0.03. Whether to express this as a fitted parameter, as a hardcoded Omicron cascade, or some other
way is a design question and is deliberately left open here.

## Contact exposure is not normalised by location size

`total_exposure_by_contacts()` in `abm/model_functions.cpp` sums the raw viral shed of everyone present and
multiplies by `ContactRates`, with **no division by the number of Person%s at the Location**. A Location's
transmission therefore scales linearly with its absolute headcount, so large Locations are super-spreader
engines by construction rather than by parameterisation.

That matters here because the population file's Locations are very unevenly sized:

| Location | ids | mean | median | p90 | p99 | max |
| --- | --- | --- | --- | --- | --- | --- |
| home | 174568 | 1.6 | 1 | 3 | 5 | 6 |
| school | 235 | 123 | 104 | 250 | 567 | 700 |
| work | 5250 | 31 | 11 | 57 | 308 | 4464 |
| event | 313 | 914 | 470 | 1789 | 6743 | 16340 |
| shop | 15 | 19068 | 9164 | 54975 | 58519 | 58519 |

Schools and workplaces are plausible as institutions, so nothing is done to them. The 15 shops at 19000
Person%s each and the events at up to 16000 are the outliers, and with no per-capita normalisation they
dominate transmission. `MaximumContacts` is the mechanism the library provides to compensate, which is
worth bearing in mind now that it is switched off: the cap of 5 on `SocialEvent` was throttling exactly
this effect.

## Notes and caveats

- **Simulation length is a plain option.** `--num_days` and `--start_date` set the window; nothing about
  the length is baked in. A 90 day window costs roughly 90 s per run at 286k agents on one core, so an
  ensemble of 3000 is on the order of 75 core-hours. A full run to the end of 2025 (about 1280 days) is
  roughly 20 minutes, which is fine for posterior predictive draws but not for an ensemble.
- **The history lookback must stay short enough not to saturate.** A Person is seeded at most once, while
  the reported cases they are drawn from include reinfections. Over 365 days around the Omicron wave,
  reported Halle cases times a dark figure of 4 come to 97% of the population, which seeds nearly everyone
  as recovered and leaves no susceptible for the epidemic to run in: deaths then flatline after 60 days.
  The default is therefore 90 days, which stays unsaturated across the whole `dark_figure` prior, and
  `seed_history()` now warns loudly whenever a request exceeds an age group's population.
- **Locations are never split.** One model Location per input id, for every type. An earlier version
  capped schools at 600 and workplaces at 100 and split anything larger, on the assumption that the input
  assigned implausibly many people per location. Measurement did not support that: it split exactly two of
  235 schools and broke up the workplace tail, both of which are realistic institution sizes. Removing the
  splitting raised 90-day deaths from 3516 to 4062.
- **Contact rates are read from `data/Germany/contacts`.** `ContactRates` defaults to an all-zero baseline,
  so a model that does not set it has no contact transmission whatsoever. The per-location factors come
  from `set_local_parameters_ger()` on the `abmXpanvadere` branch. Both that file and
  `cpp/examples/abm_aims_halle.cpp` hardcode the matrices and share a transcription error in the home
  matrix (0.0504 where `baseline_home.txt` has 0.4504), which reading the file avoids.
- **Deaths are rebased onto the start of the window** on both sides. The seeded infection history kills
  part of the population before `t0`, so the raw cumulative count is dominated by deaths that already
  happened. `fit_npe.py` rebases the real data the same way; the two must stay consistent.
- **Deaths alone are a weak signal for a single city.** Only 67 deaths fall in the whole 2022-07-01 to
  2022-12-01 window. If a posterior comes back flat, add a detected-cases or ICU channel to
  `observable_channels()` and the `LogChannels` logger and point `--fit-channel-prefix` at it.
- **`QuarantineEffectiveness` defaults to 0.0**, i.e. a quarantined Person sheds virus at full rate. Only
  the mobility block in `Model::perform_mobility` restrains them. This is inherited from the library
  default and matters as soon as testing or isolation is reinstated.
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
