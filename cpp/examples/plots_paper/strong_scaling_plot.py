from plotting_settings import plotting_dir, set_fontsize, colors, dpi
import matplotlib.pyplot as plt
import pandas as pd
import os

save_dir = os.path.join(plotting_dir, "plots")

models = {"osecir": "ODE", "lsecir": "LCT", "isecir": "IDE"}
model_colors = {"osecir": colors['Purple'], "lsecir": colors['Teal'], "isecir": colors['Orange']}

data = {
    "processors": [1, 2, 4, 8, 16, 32, 64, 128],
    "osecir": [19566.8, 9718.44, 4847.9, 2417.92, 1252.95, 615.575, 329.08, 190.454],
    "lsecir": [746.859, 373.144, 185.975, 92.8182, 46.6567, 24.8617, 12.8418, 7.05302],
    "isecir": [5358.52, 2672.95, 1341, 666.594, 332.608, 176.283, 91.0512, 46.2337]
}

set_fontsize()        
figsize = (8, 5)
panel = (0.2, 0.2, 0.78, 0.75)
fig = plt.figure(figsize=figsize)
ax = fig.add_axes(panel)

for model in models.keys():
    ax.plot(data['processors'], data[model], marker = 'o', color = model_colors[model], label = models[model])

ax.set_xlabel('Cores [#]')
ax.set_xticks(data['processors'])
ax.set_ylabel('Runtime [s]')
ax.legend()
fig.savefig(os.path.join(save_dir, 'Scaling.png'), dpi=dpi)

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xticks(data['processors'])
ax.set_xticklabels(data['processors'])
ax.legend(ncol = 1)
fig.savefig(os.path.join(save_dir, 'Scaling_log.png'), dpi=dpi)
