# Halle ABM with PAIS

The Halle agent-based model with post-acute infection syndrome (PAIS). It is used from Python through
`memilio.simulation.halle_abm.run_once`, which `pycode/examples/simulation/fitting_bayesflow.py` calls as
the simulator of a BayesFlow fit.

```python
days, values = halle_abm.run_once(person_file, contact_dir, history_file, theta, start_date, num_days)
# values[day] = [deaths since start_date, InfectedSevere, InfectedCritical]
```

`theta[0]` is `InfectionRateFromViralShed` (the viral shedding rate). Any further entries are ignored.

## Files

| File | Content |
| --- | --- |
| `halle_model.h` / `halle_model.cpp` | `make_model` (builds the model) and `run_once` (builds, simulates, logs the daily outputs) |
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
  transmission, no testing and no other intervention.
- Fixed infection parameters (transition times, age-dependent severity) and PAIS parameters follow
  `cpp/examples/abm_aims_halle.cpp`, except for the PAIS transition rates (see below).
- The infection and vaccination history before the start date is replayed from the cohort: every agent
  gets the full history of one cohort member of its age, in date order, with fatal courses redrawn and
  PAIS advanced daily up to `t0`. The cohort has no one under 18, so age groups 0-4 and 5-14 start fully
  susceptible.
- Deaths are rebased onto the start date, since the seeded history already contains deaths.

## Open points

- `PAISTransitionMatrix` is consumed as rates in 1/day. The matrix in `abm_aims_halle.cpp` is documented
  as per-day probabilities and empties the medium state within a day, so plausible literature durations
  are used instead until the semantics are settled.
- `PAISProbabilitySeverityFactor` is below 1 (0.267 / 0.247 / 0.217) while `PAIS::init_or_refresh`
  applies it as an increase for severe acute infections. Either the values or the direction is wrong.
- Transmission at a Location scales with its absolute headcount, so the 15 shops (median about 9000
  agents) and the large events dominate transmission.
- `run_once` uses a freshly seeded random number generator, so runs are not reproducible from a seed.
