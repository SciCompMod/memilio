import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

cwd = os.getcwd()
save_path = os.path.join(cwd, "saves")
os.makedirs(save_path, exist_ok=True)

# timings GNN vs Graph-ODE
df = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "times.csv"))

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
    g.fig.savefig(os.path.join(save_path, f'plots_Num_Pred_{num_pred}.png'), dpi=500)

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
fig.savefig(os.path.join(save_path, 'legend.png'), bbox_inches='tight', dpi=500)
plt.close(fig)


# MAPEs LSTM
df = pd.DataFrame(data=[[30, 0, 0.65, 7.34], [60, 0, 0.58, 5.62], [90, 0, 1.24, 9.75], [30, 1, 0.82, 8.87], [
                  30, 2, 0.73, 8.02], [30, 3, 0.99, 10.53]], columns=['Days', 'Dampings', 'Test  MAPE (log-scale)', 'Test MAPE (orig. scale)'])

columns = ['Days', 'Dampings'] 
column_labels = ['Days', 'Contact changes'] 
dfs_plot = [df[df['Dampings']==0], df[df['Days']==30]] # inverted order
for i in range(len(dfs_plot)):
    g = sns.catplot(
        data=dfs_plot[i],
        y="Test MAPE (orig. scale)",
        x=columns[i],
        # hue="Dampings_Method",
        kind="bar",
        palette="husl",
        sharex=True,
        legend=False
    )
    g.set_axis_labels(column_labels[i], "Test MAPE (orig. scale, in %)")
    g.tight_layout()
    g.savefig(os.path.join(save_path, 'results_lstm_Ibased_withagegroups_testMAPE_'+str(i)+'.png'), bbox_inches='tight', dpi=500)                  


# MAPEs GNN
df = pd.DataFrame(data=[[30, 0, 5.83, 15.05], [60, 0, 4.14, 10.33], [90, 0, 4.29, 10.69], [30, 1, 6.08, 17.89], [
                  30, 2, 6.74, 16.48], [30, 3, 8.75, 27.11]], columns=['Days', 'Dampings', 'Test  MAPE (log-scale)', 'Test MAPE (orig. scale)'])


columns = ['Days', 'Dampings'] 
column_labels = ['Days', 'Contact changes'] 
dfs_plot = [df[df['Dampings']==0], df[df['Days']==30]] # inverted order
for i in range(len(dfs_plot)):
    g = sns.catplot(
        data=dfs_plot[i],
        y="Test MAPE (orig. scale)",
        x=columns[i],
        # hue="Dampings_Method",
        kind="bar",
        palette="husl",
        sharex=True,
        legend=False
    )
    g.set_axis_labels(column_labels[i], "Test MAPE (orig. scale, in %)")
    g.tight_layout()
    g.savefig(os.path.join(save_path, 'results_gnn_Ibased_nodeswithvariance_testMAPE_'+str(i)+'.png'), bbox_inches='tight', dpi=500)                  
