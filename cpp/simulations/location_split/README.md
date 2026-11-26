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
  exists there. `LogPersonStates` records the location and the infection state of every person at
  every time step instead, and `derive_infection_events` reconstructs the infection sites from that:
  a person newly infected in step *k* was infected at the location it occupied in step *k-1*, because
  the ABM lets people interact before it moves them.
* `mio::save_results` is not used, because with percentiles enabled it also serializes the model
  parameters as a graph, which `abm::Model` does not support. The h5 files are written directly.

## Splitting schools into rooms

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
  ./build/bin/location_split_simulation --runs 50 --days 60 --n_persons 10000 --infection_k 1.25 \
      --seed 42 --school_room_size $n --output_dir results/rooms_$n
done
```

## Building

The simulation needs HDF5. It is not part of the default build:

```bash
cmake -S cpp -B build -DMEMILIO_BUILD_SIMULATIONS=ON
cmake --build build --target location_split_simulation
```

## Running

```bash
./build/bin/location_split_simulation --runs 50 --days 14 --n_persons 10000
```

| Option | Meaning |
| --- | --- |
| `--region <name>` | `germany` (default), `france` or `usa` |
| `--runs <n>` | number of runs, default 5 |
| `--days <n>` | simulated days, default 7 |
| `--n_persons <n>` | total population, default 1000 |
| `--initial_infections <n>` | randomly drawn initially infected persons, default 10 |
| `--infection_k <x>` | scale of the viral shed, default 22.6 |
| `--school_room_size <n>` | split every school into rooms of at most `n` pupils, default 0 (no split) |
| `--output_dir <path>` | output directory, default `results/run_<timestamp>` |
| `--seed <n>` | seed for reproducibility |
| `--no_infection_events` | do not write the per run infection event csv |
| `--write_contacts` | write the pairwise contact hours per run; quadratic in the location size |
| `--write_person_locations` | write where every person was at every hour; needed for the transmission tree plots |
| `--print_city` | print the derived city infrastructure and exit |

Two invocations with the same `--seed` and the same options produce identical results.

## Output

```
<output_dir>/
  summary.txt                          run configuration
  scenario.json                        label and school room size of this arm, read by viz/
  city_config.csv                      derived city statistics
  household_id.csv                     person -> home
  work_id.csv                          person -> workplace, empty for non-workers
  location_id_and_type.csv             location -> type
  initial_infection_households.csv     household size of the initially infected, per run
  runs/run_<i>_infection_events.csv    every infection with the location it happened at
  runs/run_<i>_contact_hours.csv       pairwise contact hours (only with --write_contacts)
  infection_state_per_age_group/       h5, per run and p05/p25/p50/p75/p95, age * InfectionState::Count + state
  infection_per_location_type_per_age_group/  h5, same layout, age * LocationType::Count + type
  total_infections/                    h5, cumulative number of non-susceptible persons
```

## Memory

`LogPersonStates` keeps one record per person and simulated hour, about 24 bytes each. A run with
10^6 persons over 14 days therefore needs several GB. `--write_contacts` is additionally quadratic in
the number of people per location and should only be used on small populations.

## Layout

| File | Contents |
| --- | --- |
| `include/city_parameters.h`, `src/city_parameters.cpp` | demographics and infrastructure ratios per region |
| `include/city_builder.h`, `src/city_builder.cpp` | builds the `abm::Model` and assigns people to locations |
| `include/parameter_setter.h`, `src/parameter_setter.cpp` | infection parameters and contact matrices |
| `include/custom_loggers.h`, `src/custom_loggers.cpp` | loggers |
| `include/multi_run_simulator.h`, `src/multi_run_simulator.cpp` | runs the ensemble and writes the results |
| `include/file_utils.h`, `src/file_utils.cpp` | output directories and timestamps |
| `main.cpp` | command line |
