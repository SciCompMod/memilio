from plotting_settings import plotting_dir, set_fontsize, colors, dpi
import matplotlib.pyplot as plt
import pandas as pd
import os

save_dir = os.path.join(plotting_dir, "plots")

models = {"osecir": "ODE", "lsecir": "LCT", "isecir": "IDE"}
model_colors = {"osecir": colors['Purple'], "lsecir": colors['Teal'], "isecir": colors['Orange']}

data = {
    "regions": [10, 100, 1000],
    "osecir": [],
    "lsecir": [],
    "isecir": []
}

for model in models.keys():
    folder_path  = os.path.join(plotting_dir, "..", f"simulation_paper_{models[model].lower()}", "results_runtime")
    if not os.path.isdir(folder_path):
        continue

    for region in data['regions']:
        df = pd.read_csv(os.path.join(folder_path, f"{model}_{str(region)}regions_sim_time.csv"))
        data[model].append(df.C1.mean())

set_fontsize()        
figsize = (5, 3.5)
panel = (0.2, 0.2, 0.78, 0.75)
fig = plt.figure(figsize=figsize)
ax = fig.add_axes(panel)

for model in models.keys():
    ax.plot(data['regions'], data[model], marker = 'o', color = model_colors[model], label = models[model])

ax.set_xlabel('Regions [#]')
ax.set_ylabel('Runtime [s]')
ax.legend()
fig.savefig(os.path.join(save_dir, 'mean_runtimes.png'), dpi=dpi)

ax.set_yscale('log')
ax.set_xscale('log')
ax.legend(ncol = 1)
fig.savefig(os.path.join(save_dir, 'mean_runtimes_log.png'), dpi=dpi)