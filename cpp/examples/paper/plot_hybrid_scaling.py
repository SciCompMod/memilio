from plotting_settings import *
import os

dir = 'V:/bick_ju/MemilioPaper/ScalingHybrid'
save_dir = "H:/Documents/MEmilioPaper/ScalingHybrid/"

models = {"2_hybrid": "Hybrid (2)", "5_hybrid": "Hybrid (5)", "10_hybrid": "Hybrid (10)", "dabm": "dABM", "osecir": "ODE"}
model_colors = {"dabm": colors['Teal'], "osecir": colors['Purple'], "2_hybrid": colors['Red'], "5_hybrid": colors['Orange'], "10_hybrid": colors['Yellow']}

data = {
    "pop": [],
    "dabm": [],
    "osecir": [],
    "2_hybrid": [],
    "5_hybrid": [],
    "10_hybrid": []
}

l = sorted(os.listdir(dir))
l_int = [int(v) for v in l]
l_int = sorted(l_int)
for folder in l_int:
    folder_path = os.path.join(dir, str(folder))
    if not os.path.isdir(folder_path):
        continue
    
    data['pop'].append(folder)
    for model in models.keys():
        df = pd.read_csv(os.path.join(folder_path, f"pop_{str(folder)}_{model}_sim_time.csv"))
        data[model].append(df.C1.mean())

set_fontsize()        
figsize = (5, 3.5)
panel = (0.2, 0.2, 0.78, 0.75)
fig = plt.figure(figsize=figsize)
ax = fig.add_axes(panel)

for model in models.keys():
    ax.plot(data['pop'], data[model], marker = 'o', color = model_colors[model], label = models[model])
    
ax.set_xlabel('Population size [#]')
ax.set_ylabel('Runtime [s]')
ax.legend()
fig.savefig(save_dir + 'mean_runtimes.png', dpi=dpi)

ax.set_yscale('log')
ax.set_xscale('log')
ax.legend(ncol = 1)
fig.savefig(save_dir + 'mean_runtimes_log.png', dpi=dpi)


