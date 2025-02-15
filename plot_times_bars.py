import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

cwd = os.getcwd()
save_path = os.path.join(cwd, "saves")
os.makedirs(save_path, exist_ok=True)

df = pd.read_csv("times.csv")

df_long = pd.melt(
    df,
    id_vars=["Num_Pred", "Days", "Dampings"],
    value_vars=["Mean_Time_Per_Run", "Mean_time_ODE"],
    var_name="Method",
    value_name="Time"
)
df_long["DaysStr"] = df_long["Days"].astype(str)

# get declarations for the legend


def combine_dampings_method(row):
    method_label = "GNN Surrogate" if row["Method"] == "Mean_Time_Per_Run" else "Original"
    damping_label = "no contact change" if row[
        "Dampings"] == 0 else f"{row['Dampings']} contact changes"
    return f"{method_label} ({damping_label})"


df_long["Dampings_Method"] = df_long.apply(combine_dampings_method, axis=1)

sns.set_style("whitegrid")

num_preds = df_long["Num_Pred"].unique()
for num_pred in num_preds:
    df_subset = df_long[df_long["Num_Pred"] == num_pred]
    g = sns.catplot(
        data=df_subset,
        x="Time",
        y="DaysStr",
        hue="Dampings_Method",
        kind="bar",
        palette="husl",
        sharex=True
    )
    g.set_axis_labels("Time (s)", "Days")
    g.set(xscale="log")
    g.set_titles(f"Num_Pred = {num_pred}")

    for ax in g.axes.flat:
        for p in ax.patches:
            width = p.get_width()
            ax.text(width, p.get_y() + p.get_height() / 2,
                    f'{width:.3f}', ha='left', va='center')

    # delete legend
    g._legend.remove()
    g.fig.tight_layout()
    g.fig.savefig(os.path.join(save_path, f'plots_Num_Pred_{num_pred}.png'))

# Separate legend plot
# Create a dummy plot to extract the legend
fig, ax = plt.subplots()
sns.barplot(
    data=df_long,
    x="Time",
    y="DaysStr",
    hue="Dampings_Method",
    palette="Set2",
    ax=ax
)
handles, labels = ax.get_legend_handles_labels()
plt.close(fig)  # Close the dummy plot

# Create a new figure for the legend
fig, ax = plt.subplots()
fig.legend(handles, labels, loc='center', ncol=4)
plt.axis('off')
fig.savefig(os.path.join(save_path, 'legend.png'), bbox_inches='tight')
plt.close(fig)
