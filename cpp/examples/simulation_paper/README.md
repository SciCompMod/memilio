# Example intervention strategies

Run the example graph_germany_nuts3_ode.py with one of the four scenarios (0-4).
From build dir:

```bash
./bin/graph_germany_nuts3_ode -TestCase 0
```

Results will automatically be written into the results folder. 
After running all 4 test cases you can create plots by running the files in the plots_paper folder starting with "nuts3_"
I would suggest:
- nuts3_single_plot_all.py
- nuts3_triple_plot_simple.py
- nuts3_triple_plot_no_open.py
- nuts3_dynamic_regional_plot.py

And if you need legends for these plots run:
- nuts3_legend.py


## State Epidemiological parameters and states
Additionally, you will get terminal output of current Epidemiological parameters and states. You can save them by running

```bash
./bin/graph_germany_nuts3_ode -TestCase 0 > ../examples/simulation_paper/results/Output_open.txt
```
Change name of file according to your test case (0: "open", 1:"same", 2:"lockdown", 3"dynamic").
After that you can process the terminal output via the script nuts3_process_contacts.py.
