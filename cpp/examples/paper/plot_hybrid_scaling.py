from plotting_settings import *
import os
import numpy as np

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
figsize = (8, 5)
panel = (0.2, 0.2, 0.78, 0.75)
fig = plt.figure(figsize=figsize)
ax = fig.add_axes(panel)

for model in models.keys():
    ax.plot(data['pop'], data[model], marker = 'o', color = model_colors[model], label = models[model])

x_last = data['pop'][-1]
y1_last = data['dabm'][-1]
y2_last = data['10_hybrid'][-1]

x_mid = x_last
y_mid = np.sqrt(y1_last * y2_last)
N = 15
transform = ax.transData
inv_transform = ax.transData.inverted()
p1 = transform.transform((x_last, y1_last))
p2 = transform.transform((x_last, y2_last))
ys_disp = np.linspace(p1[1], p2[1], N)
ys = [inv_transform.transform((p1[0], yd))[1] for yd in ys_disp]

speedup = y1_last / y2_last
label = f"{speedup:.0f}x"

ax.scatter(
    [x_last] * N,   # all at same x
    ys,
    marker="^",
    color=colors['Dark grey'],
    s=20
)

ax.text(
    x_mid*0.96, y_mid,
    label,
    ha='center', va='bottom',
    color=colors['Dark grey']
)
    
ax.set_xlabel('Population size [#]')
ax.set_ylabel('Runtime [s]')
ax.legend()
ax.grid(True, alpha=0.3) 
fig.savefig(save_dir + 'mean_runtimes.png', dpi=dpi)

ax.set_yscale('log')
ax.set_xscale('log')
#ax.legend(ncol = 1)
fig.savefig(save_dir + 'mean_runtimes_log.png', dpi=dpi)


