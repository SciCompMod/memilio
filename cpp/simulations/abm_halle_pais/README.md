# Halle ABM with PAIS

The Halle agent-based model with post-acute infection syndrome (PAIS). It is used from Python through
`memilio.simulation.halle_abm.run_once`, which `pycode/examples/simulation/fitting_bayesflow.py` calls as
the simulator of a BayesFlow fit.

```python
days, values = halle_abm.run_once(person_file, contact_dir, history_file, theta, start_date, num_days)
# values[day] = [deaths since start_date, new admissions, ICU occupancy, new positive tests]
```

`theta` holds, in the order of `halle_abm.parameter_names()`: `viral_shedding_rate`
(`InfectionRateFromViralShed`), `dark_figure` (>= 1), `testing_probability_symptomatic` and
`ratio_symptomatic_to_asymptomatic` (>= 1). The values are for the model population (286,022 agents, which
include persons outside Halle); the fit scales them by the population of Halle over the number of agents.

Fitting workflow (from the repository root):

1. Download the RKI/DIVI data with memilio-epidata (see the docstring of `prepare_observed_data.py`).
2. `python3 cpp/simulations/abm_halle_pais/prepare_observed_data.py` writes
   `data/Germany/halle/halle_observed.csv` and `halle_observed_meta.json` for 2022-07-01 + 90 days.
3. `python3 pycode/examples/simulation/fitting_bayesflow.py`.

`analyze_population.py` plots the population file to check its plausibility.

## Files

| File | Content |
| --- | --- |
| `halle_model.h` / `halle_model.cpp` | `make_model` (builds the model) and `run_once` (builds, simulates, logs the daily outputs) |
| `prepare_observed_data.py` | Real deaths, admissions, ICU and reported cases for Halle in the model's channels |
| `analyze_population.py` | Summary and plots of the population file |
| `pycode/.../bindings/models/abm_pais.cpp` | pybind11 module `_simulation_halle_abm`, compiles `halle_model.cpp` directly |

## Input data

| Argument | File | Columns |
| --- | --- | --- |
| `person_file` | `halle_population_data.csv` | `age,home_id,school_id,event_id,shopping_id,work_id` (extra columns are ignored); `age` is an index into an 11-group scheme, aggregated onto the six RKI groups |
| `history_file` | `260304_infections_vaccines_halle.csv` | `age_group,ih_infection_{1..4}_date,ih_vaccine_{1..3}_date`, dates as `YYYY-MM-DD` or `NA` |
| `contact_dir` | `data/Germany/contacts` | `baseline_{home,school_pf_aver,work,other}.txt`, 6x6 |

The population and cohort files are not in git (`*.csv` is ignored); locally they are in `data_additional/`.

## What the model does

- One agent per row of the population file, one Location per input id, never split. Sex is drawn from
  the age-dependent female share of Halle, since the PAIS parameters are stratified by sex.
- Every agent is assigned a Hospital and an ICU so that severe and critical cases are moved there.
- Contact rates are the German baseline matrices times per-location factors from
  `set_local_parameters_ger()` on the `abmXpanvadere` branch. No `MaximumContacts`, no aerosol
  transmission.
- Movement as on `abmXpanvadere`: work starts between 6:00 and 8:00, school at 8:00, no shopping for
  ages 0-4, shopping takes 2 hours, weekday events start from 20:00 and last at least 5 hours, and
  symptomatic agents isolate at home 2 days after symptom onset. These use the ABM parameters
  `BasicShoppingDuration`, `SocialEventDuration`, `SocialEventWeekdayStart` and
  `SymptomaticIsolationDelay`, whose defaults keep the previous behaviour for other models. Unlike on
  `abmXpanvadere`, an event stay that runs past midnight ends after its minimum duration instead of
  lasting until 20:00 the next day.
- Testing as in `paper_abm_bs_testing.cpp` (branch `abm_paper_test_bs`): on entering a school, workplace,
  shop or social event, symptomatic agents are tested with `testing_probability_symptomatic`, all others
  with that probability divided by `ratio_symptomatic_to_asymptomatic`. Generic test (sensitivity 0.71,
  specificity 0.996), results valid for 3 days; a positive test denies entry and starts quarantine.
- Fixed infection parameters (transition times, age-dependent severity) and PAIS parameters follow
  `cpp/examples/abm_aims_halle.cpp`, except for the PAIS transition rates (see below).
- The infection and vaccination history before the start date is replayed from the cohort: every agent
  gets the full history of one cohort member of its age, in date order, with courses fatal before `t0`
  redrawn and PAIS advanced daily up to `t0`. Drawing a member per agent is what scales the 5,000-person
  cohort up to the population. The cohort has no one under 18, so age groups 0-4 and 5-14 start fully
  susceptible (to be discussed).
- The initial dark figure scales the cohort infections of the 30 days before `t0`: agents whose donor has
  no recent infection get an additional one with the probability that makes the recent share of the age
  bin `dark_figure` times the cohort's, dated like the recent cohort infections. Older history is unscaled.
- Outputs per day: deaths since `t0` and ICU (InfectedCritical) occupancy at the end of the day, new
  InfectedSevere (admissions) and new positive tests during the day.
- `DeathsPerInfectedSevere` is not set, so its default of 0.1 applies on top of the deaths through ICU.

## Open points

- `PAISTransitionMatrix` is consumed as rates in 1/day. The matrix in `abm_aims_halle.cpp` is documented
  as per-day probabilities and empties the medium state within a day, so plausible literature durations
  are used instead until the semantics are settled.
- `PAISProbabilitySeverityFactor` is below 1 (0.267 / 0.247 / 0.217) while `PAIS::init_or_refresh`
  applies it as an increase for severe acute infections. Either the values or the direction is wrong.
- Transmission at a Location scales with its absolute headcount, so the 15 shops (median about 9000
  agents) and the large events dominate transmission.
- `run_once` uses a freshly seeded random number generator, so runs are not reproducible from a seed.
