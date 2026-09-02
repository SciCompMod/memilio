# Comparison plots for the location split simulation

Every script here compares an arbitrary number of **arms** of the location split study. An arm is one
output directory of `location_split_simulation`. Arms are ordinal, they differ in how much of the
split location is handled outside the ABM:

```
no split (baseline)  ->  quarter split  ->  half split  ->  full split
```

The plots order the arms along that axis, colour them along a sequential ramp, and draw the baseline
in black as the reference everything else is measured against. Two arms, four arms or seven arms all
work; nothing is hardcoded to a pair.

## Quick start

```bash
python -m pip install -r requirements.txt

python compare_scenarios.py \
    --scenario 'no split=results/no_split' \
    --scenario 'quarter split=results/quarter_split' \
    --scenario 'half split=results/half_split' \
    --scenario 'full split=results/full_split' \
    --output figures
```

That writes every figure plus `summary_statistics.csv` and prints a comparison table.

## How an arm learns its split level

In order of precedence:

1. **`scenario.json` in the results directory** — the intended mechanism once the split is
   implemented:

   ```json
   {"label": "half split", "split_fraction": 0.5, "split_location_type": "SocialEvent"}
   ```

   `split_fraction` is the share of the location handled outside the ABM: `0.0` is the baseline,
   `1.0` is a full split.

2. **The label** on the command line, or the directory name if no label is given. Understood:
   `no`, `quarter`, `third`, `half`, `three_quarter`, `full`, or a number (`split_0.75`, `75%`).

3. **The order on the command line**, if neither of the above resolves.

`--baseline 'LABEL'` overrides which arm is the reference; by default it is the arm with the
smallest split fraction.

## The scripts

| Script | Figure | Replaces on the `abmXpanvadere` branch |
| --- | --- | --- |
| `compare_scenarios.py` | everything below + `summary_statistics.csv` | `compare_all_scenarios.py`, `combined_scenario_comparison.py` |
| `plot_epidemic_curves.py` | cumulative infections with IQR band, daily incidence, difference to baseline | `infection_cumulated_comp.py`, `simple_multi_panel.py`, `multi_seed_comparison.py` |
| `plot_location_breakdown.py` | composition and absolute infections per location type, difference to baseline; `--pies` for the original pie layout | `location_infection_pie_chart.py` |
| `plot_location_timeline.py` | location type x time heatmap per arm, `--differences` against baseline, plus a stacked area view | `temporal_infection_heatmap.py`, `time_loc_type.py` |
| `plot_household_heatmap.py` | infection burden per household over time, arms as rows; `--use-workplaces` | `comparative_temporal_heatmap.py` |
| `plot_transmission_tree.py` | transmission chains over time and the generation / reproduction structure | `infection_timeline.py`, `comparative_infection_trees.py`, `contact_network.py` |
| `scenario_io.py` | shared loading, ordering and colours | — |

Common options: `--scenario`, `--baseline`, `--output`, `--dpi`, `--format {png,pdf,svg}`.

## Data each script needs

| Script | Needs |
| --- | --- |
| `plot_epidemic_curves.py` | `total_infections/` (always written) |
| `plot_location_breakdown.py`, `plot_location_timeline.py` | `runs/run_*_infection_events.csv` (default on, disabled by `--no_infection_events`) |
| `plot_household_heatmap.py` | the above plus `household_id.csv` / `work_id.csv` |
| `plot_transmission_tree.py` | the above plus `runs/run_*_person_locations.csv`, so the simulation must be run with `--write_person_locations` |

`compare_scenarios.py --skip-trees` leaves out the part that needs `--write_person_locations`.

## How the transmission tree is reconstructed

The simulation records, for every infection, the exact hour and Location it happened at, and with
`--write_person_locations` also where every Person was at every hour. The infector of somebody
infected during hour *t* at Location *L* must have been at *L* during that hour and already
infectious, which is usually a handful of candidates; among those the one furthest into its own
infectious window is picked. On the test data every single infection resolved to a candidate.

This is a reconstruction, not ground truth: the ABM computes an infection risk per person rather
than naming an infector, so the assignment within a location is a plausible attribution, not a
record. The generation depth and the offspring distribution are the parts that compare meaningfully
across arms; individual edges are not.

## Notes

- Arms may have different city realizations. Person and location ids are therefore **not**
  comparable between arms, which is why every plot aggregates before comparing.
- The household heatmap sorts households by size and separates the size blocks with thin lines. The
  branch encoded size as a per-size colour; sorting keeps the same information and frees the colour
  channel for the infection rate.
- Pie charts stop being readable past two arms, so `plot_location_breakdown.py` defaults to stacked
  and grouped bars. `--pies` still produces the original layout.
