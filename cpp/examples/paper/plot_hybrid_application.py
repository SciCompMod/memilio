import pandas as pd
import matplotlib.patches as mpatches
import numpy as np

from plotting_settings import *

def plot_ts_output(axis, filename, percentiles, color, compartments, model):
    while len(percentiles) > 0:
        if len(percentiles) == 1:
            df = pd.read_csv(filename + f"_p{percentiles[0]}.csv")
            axis.plot(df["Time"], df[compartments].sum(axis=1), color=color)
            del percentiles[0]
        else:
            df_low = pd.read_csv(filename + f"_p{percentiles[0]}.csv")
            df_high = pd.read_csv(filename + f"_p{percentiles[-1]}.csv")
            axis.plot(df_low["Time"], df_low[compartments].sum(axis=1), color=color, alpha=0.2)
            axis.plot(df_high["Time"], df_high[compartments].sum(axis=1), color=color, alpha=0.2)
            if(model == "dabm"):
                axis.fill_between(df_low["Time"], df_low[compartments].sum(axis=1), df_high[compartments].sum(axis=1), alpha=0.3, facecolor='none', edgecolor=color, hatch='////', linewidth=0.0)
            elif(model == "hybrid"):
                axis.fill_between(df_low["Time"], df_low[compartments].sum(axis=1), df_high[compartments].sum(axis=1), alpha=0.3, facecolor='none', edgecolor=color, hatch='\\\\', linewidth=0.0)
            else:
                axis.fill_between(df_low["Time"], df_low[compartments].sum(axis=1), df_high[compartments].sum(axis=1), alpha=0.3, color=color)
            del percentiles[0]
            del percentiles[-1]
    return axis

def plot_extinction(axis, filename, percentiles, color, compartments, model):
    while len(percentiles) > 0:
        if len(percentiles) == 1:
            linestyle = "--" if model=="dabm" else "-"
            df = pd.read_csv(filename + f"_extinction_p{percentiles[0]}.csv")
            axis.plot(df["Time"], df[compartments].sum(axis=1), color=color, linestyle=linestyle)
            del percentiles[0]
        else:
            df_low = pd.read_csv(filename + f"_extinction_p{percentiles[0]}.csv")
            df_high = pd.read_csv(filename + f"_extinction_p{percentiles[-1]}.csv")
            axis.plot(df_low["Time"], df_low[compartments].sum(axis=1), color=color, alpha=0.2)
            axis.plot(df_high["Time"], df_high[compartments].sum(axis=1), color=color, alpha=0.2)
            if(model == "dabm"):
                axis.fill_between(df_low["Time"], df_low[compartments].sum(axis=1), df_high[compartments].sum(axis=1), alpha=0.3, facecolor='none', edgecolor=color, hatch='////', linewidth=0.0)
            elif(model == "hybrid"):
                axis.fill_between(df_low["Time"], df_low[compartments].sum(axis=1), df_high[compartments].sum(axis=1), alpha=0.3, facecolor='none', edgecolor=color, hatch='\\\\', linewidth=0.0)
            else:
                axis.fill_between(df_low["Time"], df_low[compartments].sum(axis=1), df_high[compartments].sum(axis=1), alpha=0.3, color=color)
            del percentiles[0]
            del percentiles[-1]
    return axis

result_dir = "V:/bick_ju/MemilioPaper/HybridApplication/"
save_dir = "H:/Documents/MEmilioPaper/HybridApplication/"
model_colors = {"dabm": colors["Teal"], "hybrid": colors["Red"], "osecir": colors["Purple"], "init": colors["Orange"], "sim": colors["Blue"]}
set_fontsize()

# Plot time series
panel = [0.18, 0.18, 0.8, 0.8]
figsize = (5, 3.5)
fig = plt.figure(figsize=figsize)
ax = fig.add_axes(panel)
ax = plot_ts_output(ax, result_dir + "dabm", ["05", "50", "95"], model_colors["dabm"], ["E", "Ins", "InsC", "Isy", "IsyC", "Isev", "Icri"], "dabm")
ax = plot_ts_output(ax, result_dir + "osecir", ["05", "50", "95"], model_colors["osecir"], ["E", "Ins", "InsC", "Isy", "IsyC", "Isev", "Icri"], "osecir")
ax = plot_ts_output(ax, result_dir + "hybrid", ["05", "50", "95"], model_colors["hybrid"], ["E", "Ins", "InsC", "Isy", "IsyC", "Isev", "Icri"], "hybrid")
ax.set_ylabel("Infected [#]")
#ax.set_xlabel("Time [days]")
fig.savefig(save_dir + "hybrid_application_sim_results.png", dpi=dpi)

fig = plt.figure(figsize=figsize)
fig.legend(handles=[mpatches.Patch(color=model_colors["dabm"], label="ABM"), mpatches.Patch(color=model_colors["osecir"], label="PBM"), mpatches.Patch(color=model_colors["hybrid"], label="Hybrid")], loc='center', ncol=1)
fig.savefig(save_dir + "hybrid_application_legend.png", dpi=dpi)

# Plot extinction for ABM and Hybrid
fig = plt.figure(figsize=figsize)
ax = fig.add_axes(panel)
ax = plot_extinction(ax, result_dir + "hybrid", ["05", "50", "95"], model_colors["hybrid"], ["E", "Ins", "InsC", "Isy", "IsyC", "Isev", "Icri"], "hybrid")
ax = plot_extinction(ax, result_dir + "dabm", ["05", "50", "95"], model_colors["dabm"], ["E", "Ins", "InsC", "Isy", "IsyC", "Isev", "Icri"], "dabm")
ax.set_ylabel("Infected [#]")
ax.set_xlabel("Time [days]")
fig.savefig(save_dir + "hybrid_application_sim_results_ext_surv.png", dpi=dpi)

# Runtimes (mean)
sim_time = []
init_time = []

df_sim_all = pd.DataFrame()
df_init_all = pd.DataFrame()

for key in model_colors.keys():
    if key in ["init", "sim"]:
        continue
    df_init = pd.read_csv(result_dir + key + "_init_time.csv")
    df_sim = pd.read_csv(result_dir + key + "_sim_time.csv")
    df_sim_all[key] = df_sim["C1"]
    df_init_all[key] = df_init["C1"]
    print("Median:")
    print(f"{key} init: {df_init['C1'].median()} s, sim: {df_sim['C1'].median()} s")
    print("Mean:")
    print(f"{key} init: {df_init['C1'].mean()} s, sim: {df_sim['C1'].mean()} s")
    sim_time.append(df_sim["C1"].mean())
    init_time.append(df_init["C1"].mean())
    
labels = ["ABM", "Hybrid", "PBM"]
# Pie chart for ABM and Hybrid
fig = plt.figure(figsize=(figsize[0]/2., figsize[1]/2.))
ax = fig.add_axes([0.1, 0.35, 0.8, 0.7])
ax.pie(sim_time[:2], startangle=90, colors=[model_colors["dabm"], model_colors["hybrid"]])
ax.legend(handles=[mpatches.Patch(color=model_colors["dabm"], label=f"ABM ({100*sim_time[0]/sum(sim_time[:2]):.2f}%)"), mpatches.Patch(color=model_colors["hybrid"], label=f"Hybrid ({100*sim_time[1]/sum(sim_time[:2]):.2f}%)")], loc="lower center",bbox_to_anchor=(0.5, -0.5), ncol=1)
fig.savefig(save_dir + "hybrid_application_sim_time_pie.png", dpi=dpi)
plt.close(fig)


# Boxplots for all three models
fig = plt.figure(figsize=figsize)
ax = fig.add_axes(panel)
positions1 = np.arange(len(df_sim_all.columns))
positions2 = positions1 + 0.3
ax.boxplot([df_sim_all[col] for col in df_sim_all.columns], positions=positions1, widths=0.25, patch_artist=True,
           boxprops=dict(facecolor=model_colors["sim"]))
ax.boxplot([df_init_all[col] for col in df_init_all.columns], positions=positions2, widths=0.25, patch_artist=True,
           boxprops=dict(facecolor=model_colors["init"]))

ax.set_xticks(positions1 + 0.15)
ax.set_xticklabels(df_sim_all.columns)
ax.set_ylabel('Runtime [s]')
ax.set_yscale('log')
fig.savefig(save_dir + "hybrid_application_runtimes_boxplot.png", dpi=dpi)
    
fig = plt.figure(figsize=(figsize[0], figsize[1]))
ax = fig.add_axes(panel)
x = np.arange(len(labels))
width = 0.35

ax.bar(x - width/2, init_time, width, label='Initialization', color = model_colors["init"])
ax.bar(x + width/2, sim_time, width, label='Simulation', color = model_colors["sim"])

ax.set_ylabel('Runtime [s]')
ax.set_yscale('log')
ax.set_xticks(x)
ax.set_xticklabels(labels)
fig.savefig(save_dir + "hybrid_application_runtimes_mean.png", dpi=dpi)

fig = plt.figure(figsize=figsize)
fig.legend(handles=[mpatches.Patch(color=model_colors["init"], label="Initialization"), mpatches.Patch(color=model_colors["sim"], label="Simulation")], loc='center', ncol=2)
fig.savefig(save_dir + "hybrid_application_runtimes_legend.png", dpi=dpi)
