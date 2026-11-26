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
| `scan_school.py` | fits the school weight and `k` so that schools cause 15% or 25% of all infections at the growth rate of the calibrated model |
| `analyze.py` | paired comparison of every arm with the no-split arm of its group: infected by day 90, delay to 20% infected, peak of the 7-day mean incidence, infection shares per setting; writes `paired_results_day90.csv` |
| `make_figures.py` | the figures of the presentation |

Environment variables: `STUDY_DIR` (simulation output, default `study_output`), `STUDY_OUT` (tables and
figures, default `$STUDY_DIR/results`), `SIM_BINARY` (default `build/bin/location_split_simulation`),
`CALIBRATION_RESULT` (output of `calibration/fit_contact_weights.py`, for the calibration figure).
