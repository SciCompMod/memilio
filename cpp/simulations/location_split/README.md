# Location split simulation

A multi-run agent based simulation on a synthetic city. The city is generated from the demographic
and infrastructure statistics of a region (Germany, France or the United States) and the epidemic is
started from a number of randomly drawn infected persons.

This is the port of `cpp/examples/panvXabmSim` from the `abmXpanvadere` branch onto the ABM of `main`.

## What is different from panvXabmSim

* The ABM was renamed and reworked between the two branches. Everything here is written against
  `mio::abm::Model`; no file under `cpp/models/abm` is modified.
* The Panvadere coupling is gone. On the old branch the initial infections came either from an
  external microsimulation of a restaurant or work meeting, or from a small dedicated event ABM.
  Here they are drawn uniformly at random from the population (`--initial_infections`).
* The three sets of region parameters used to live in three headers that all reused the namespace
  `CityParameters`, so only one could be compiled at a time. They are now `RegionParameters` values
  selected with `--region`.
* The branch logged the site of every infection with `Location::get_infected_persons()`, which only
  exists there. Here the simulation advances hour by hour and compares every person with its state
  of the previous hour: a person newly infected in step *k* was infected at the location it occupied
  in step *k-1*, because the ABM lets people interact before it moves them.
* `mio::save_results` is not used, because with percentiles enabled it also serializes the model
  parameters as a graph, which `abm::Model` does not support. The h5 files are written directly.

## Calibrating the contact weights

The contact rate of each location type is its contact matrix times 24 h / the hours a person spends
there on a day it goes there (`HoursPerVisitDay`, measured from the mobility rules: home 19.3, school 8,
work 10, social event 4, shop 1.2), times a weight, set with
`--contact_weights home=...,school=...,work=...,social_event=...,shop=...` (defaults: the calibrated values below).
`calibration/fit_contact_weights.py` fits the weights to an observed distribution of infection
sites:

```bash
python calibration/fit_contact_weights.py --binary build/bin/location_split_simulation --output calibration_output
```

It runs short simulations (20 runs, 10,000 persons, 21 days, measured from day 7), compares the shares
of infections of persons aged 15+ at Home, SocialEvent and Work with the targets, and adjusts the home
and social event weights and `k` until the shares match and the epidemic grows as fast as with the
start weights (by default the current defaults). Work, school and shop keep their weights.

The default targets come from the French ComCor study (Galmiche et al., BMC Public Health 2024): known
or suspected sources of 584,846 infected adults were household 30.8%, extended family 15.6%, work
colleagues 15.0% and friends 11.0%. Mapped to the ABM (extended family and friends to SocialEvent) and
renormalized: Home 42.5%, SocialEvent 36.7%, Work 20.7%.

With `--fix_k`, `k` stays at `--k` and only the weights are fitted; the growth rate is then not
matched. The defaults were fitted that way with `k = 1.596` (seed 42), converged within 0.5 percentage
points per share in three iterations; the early growth rate is 0.131 per day (doubling time 5.3 days).
These values are the hardcoded defaults (`ContactWeights` in `parameter_setter.h`, `DEFAULT_INFECTION_K`
in `defaults.h`):

```
python calibration/fit_contact_weights.py --binary build/bin/location_split_simulation --k 1.596 --fix_k
--infection_k 1.596 --contact_weights home=20.3425,school=1.6,work=0.625,social_event=6.08233,shop=0.1584
```

With these defaults, 10,000 persons and 10 initial infections, about 90% of the population is infected
by day 90 (seed 42, 20 runs). The previous defaults (`k = 2.11`, home 21.616, social event 4.8) gave
about 92%.

Caveats:

* The home weight is poorly determined. Households are close to saturation, so the home share barely
  responds to it: home weights from 16.4 to 21.6 give home shares of 40.7% to 41.3% (at `k = 2.11`). The social
  event weight and `k` are well determined.
* School and shop have no target data. Their weights stay relative to work, so the school share of all
  infections (15% after the fit, 20% before) is an assumption, not a fitted value. Studies whose
  result depends on school transmission should vary the school weight.
* ComCor covers adults under the NPIs of 2020-2022 in France; it is a first target, not German data.

## Splitting locations into rooms

`--school_room_size`, `--work_room_size` and `--event_room_size` split schools, workplaces and social
event locations in the same way; each type has its own shuffle stream. The description below uses
schools as the example.

With `--school_room_size n` every school (elementary and secondary alike) is split into
`ceil(pupils / n)` rooms of near equal size. Each room is a School location of its own, not a cell,
because the mobility rules of the ABM always place persons in cell 0. The pupils of a school are
shuffled at random and every pupil keeps its room for the whole simulation. Every room gets the same
school contact rates as an unsplit school, so a pupil has as many contacts per hour as before, only
with fewer persons.

The shuffle uses a stream with the counter of the run but its own key, so it does not draw from the
Model's generator. Arms that differ only in `--school_room_size` and share `--seed` therefore have the
same city, the same initially infected persons and the same person streams in every run; the rooms
are the only difference. The additional rooms are added after all other locations, so all other
location ids are the same as well.

```bash
for n in 0 50 25 10 5; do
  ./build/bin/location_split_simulation --runs 50 --n_persons 10000 --seed 42 \
      --school_room_size $n --output_dir results/rooms_$n
done
```

## Building

The simulation needs HDF5. It is not part of the default build. Independent runs are executed in
parallel with OpenMP if it is found; otherwise (e.g. Apple Clang without `libomp`) one after another.
With Apple Clang, point CMake to Homebrew's `libomp` (`brew install libomp`):

```bash
cmake -S cpp -B build -DMEMILIO_BUILD_SIMULATIONS=ON -DOpenMP_ROOT=$(brew --prefix libomp)
cmake --build build --target location_split_simulation
```

Keep `MEMILIO_ENABLE_OPENMP` off (the default). If it is on, the ABM library parallelises every
model itself, and its `omp single` blocks would bind to the loop over runs, so the runs are then
executed one after another with a parallel model each. That is much slower for many small runs.

## Running

```bash
./build/bin/location_split_simulation --runs 50 --n_persons 10000
```

| Option | Meaning |
| --- | --- |
| `--region <name>` | `germany` (default), `france` or `usa` |
| `--runs <n>` | number of runs, default 5 |
| `--days <n>` | simulated days, default 90 |
| `--n_persons <n>` | total population, default 1000 |
| `--initial_infections <n>` | randomly drawn initially infected persons, default 10 |
| `--infection_k <x>` | scale of the viral shed, default 1.596 (calibrated together with the contact weights) |
| `--school_room_size <n>` | split every school into rooms of at most `n` pupils, default 0 (no split) |
| `--work_room_size <n>` | split every workplace into rooms of at most `n` employees, default 0 |
| `--event_room_size <n>` | split every social event location into rooms of at most `n` persons, default 0 |
| `--contact_weights <w>` | relative transmission weight per location type; defaults are the calibrated values `home=20.3425,school=1.6,work=0.625,social_event=6.08233,shop=0.1584` |
| `--threads <n>` | runs executed in parallel, default: all cores |
| `--output_dir <path>` | output directory, default `results/run_<timestamp>` |
| `--seed <n>` | seed for reproducibility |
| `--write_person_locations` | write where every person was at every hour; needed for the transmission tree plots |
| `--redraw_school_rooms` | with `--school_room_size`: assign the pupils of every school to its rooms anew every midnight (same room sizes, own random stream, so day 0 is identical to fixed rooms) |
| `--write_school_rooms` | write, per school room and hour, the persons N, susceptible s and infectious i present and the new infections, aggregated over rooms with the same day, N, s and i; used by `study/theory_check.py` |
| `--print_city` | print the derived city infrastructure and exit |

Two invocations with the same `--seed` and the same options produce identical results, independent
of `--threads`. Every run builds its own city; run *i* of two invocations with the same seed and city
options has the same city.

## Output

```
<output_dir>/
  summary.txt                          run configuration
  scenario.json                        label and school room size of this arm, read by viz/
  city_config.csv                      derived city statistics
  location_id_and_type.csv             location -> type (the same in every run)
  runs/run_<i>_persons.csv             age group, home, workplace and school (room) of every person in run i
  runs/run_<i>_infection_events.csv    every infection with the location it happened at
  runs/run_<i>_person_locations.csv    location of every person at every hour (only with --write_person_locations)
  runs/run_<i>_school_rooms.csv        day, N, s, i, room hours and new infections in school rooms (only with --write_school_rooms)
  infection_state_per_age_group/       h5, hourly, per run and p05/p25/p50/p75/p95, age * InfectionState::Count + state
  infection_per_location_type_per_age_group/  h5, new infections per day, age * LocationType::Count + type
  total_infections/                    h5, cumulative number of non-susceptible persons
```

## Memory

Infection events are detected while the simulation runs, so a run only keeps the state of the
previous hour (a few bytes per person) and the hourly aggregated time series. A run of 10,000 persons
over 90 days needs about 30 MB. Parallel runs each need their own memory. `--write_person_locations`
streams to disk and grows with persons times hours.

## Model limitations

* Households are filled by drawing age groups in proportion to the persons left, so the age mix is
  the same in every household size. There is no family structure: households with children are only
  guaranteed one person of age 15+, and about a fifth of them have no adult of age 15-59.
* The contact weights are fitted to the ComCor infection sites for Home, SocialEvent and Work only;
  the school and shop weights are not calibrated.
* Nobody of age 60+ works (see `calc_num_workplaces_and_worker`), the whole 5-14 age group attends
  school, there is no daycare, and every city has one hospital and one ICU.
* Every person has one fixed social event location, so social events are fixed groups of about 15.
* There is no testing or isolation; symptomatic persons keep going to school and work.

## Layout

| File | Contents |
| --- | --- |
| `include/city_parameters.h`, `src/city_parameters.cpp` | demographics and infrastructure ratios per region |
| `include/city_builder.h`, `src/city_builder.cpp` | builds the `abm::Model` and assigns people to locations |
| `include/parameter_setter.h`, `src/parameter_setter.cpp` | infection parameters and contact matrices |
| `include/multi_run_simulator.h`, `src/multi_run_simulator.cpp` | runs the ensemble and writes the results |
| `include/file_utils.h`, `src/file_utils.cpp` | output directories and timestamps |
| `main.cpp` | command line |
