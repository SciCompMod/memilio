# Location split study

Scripts that run the school, workplace and social event splits with the calibrated defaults and
evaluate them. Run from the repository root after building `location_split_simulation`:

```bash
STUDY_DIR=study_output CALIBRATION_RESULT=calibration_output/fit_result.json \
    bash cpp/simulations/location_split/study/pipeline.sh
```

| Script | What it does |
| --- | --- |
| `pipeline.sh` | runs all arms one after another (finished arms are skipped), then the analysis and figures |
| `theory_check.py` | compares the school room counts (persons, susceptible and infectious present, new infections per room and hour) with the one-step formula of the theory note, and fixed with daily redrawn rooms |
| `analyze.py` | paired comparison of every arm with the no-split arm of its group: infected by day 90, delay to 20% infected, peak of the 7-day mean incidence, infection shares per setting; writes `paired_results_day90.csv` |
| `make_figures.py` | the figures of the presentation |

Environment variables: `STUDY_DIR` (simulation output, default `study_output`), `STUDY_OUT` (tables and
figures, default `$STUDY_DIR/results`), `SIM_BINARY` (default `build/bin/location_split_simulation`),
`CALIBRATION_RESULT` (output of `calibration/fit_contact_weights.py`, for the calibration figure; skipped if
missing), `THREADS` (runs in parallel, default: all cores), `PYTHON` (interpreter with numpy, pandas, h5py
and matplotlib, default `python3`).

A higher school share (e.g. 15% or 25% of all infections) cannot be reached in this model: only about
10% of the persons go to school and almost all of them are infected anyway, so schools can cause at most
about 12% of the infections at the default attack rate.

Arms, all with seed 42, 10,000 persons, 90 days and 10 initial infections:

| Arms | Runs | Compared with |
| --- | --- | --- |
| `sw1_no_split`, `sw1_rooms_50/25/10/5`: fixed school rooms | 50 | `sw1_no_split` |
| `sw1_redrawn_rooms_50/25/10/5`: school rooms redrawn every day | 50 | `sw1_no_split` |
| `work_rooms_5/3`, `event_rooms_5/3`: workplace and social event rooms, only with `WORK_EVENT_SPLITS=1` | 30 | `sw1_no_split` |
