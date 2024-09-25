import geopandas
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from matplotlib.colors import LogNorm

path_cwd = os.getcwd()
data = {
    "t$^{max}_{tr}$": [1, 0.1, 0.01, 0.001, 0.0001, 1e-05, 1e-06, 1e-07, 1e-08, 1e-09, 1e-10, 9.999999999999999e-12],
    "Error Total": [110006.855862158, 7809.678580792167, 758.4331693110972, 75.5635651796531, 7.553568871060584,
                    0.7553290186968746, 0.07553262348223019, 0.007553257381912364, 0.0007553240172079833, 7.553329773282593e-05,
                    7.553509614659763e-06, 7.555947320214238e-07],
    "Error Node 1": [42652.74831245717, 3116.170021125493, 291.6104624439494, 29.05456909470188, 2.904395888854338,
                     0.2904289812222902, 0.02904279257331477, 0.002904275449273624, 0.0002904251050594491,
                     2.904295436473338e-05, 2.903962008063293e-06, 2.905622041721283e-07],
    "Error Node 2": [11842.27890009439, 868.8034841553293, 91.80981365107719, 9.145218934132371, 0.9141655800415017,
                     0.09141299606120287, 0.009141264136454205, 0.0009141262820766135, 9.141279536149776e-05,
                     9.141163362729119e-06, 9.138012270565464e-07, 9.114645517884858e-08],
    "Error Node 3": [14805.40108249544, 611.9039145924595, 104.8323831287708, 10.44402371852397, 1.044011616163256,
                     0.1043972558069759, 0.01043968692584691, 0.001043968956469398, 0.0001043973541520541, 1.044021767281218e-05,
                     1.044495057324539e-06, 1.048376146005816e-07],
    "Error Node 4": [24586.05619540355, 1731.906475605919, 169.2367481621293, 16.86220954590431, 1.685608444948962,
                     0.1685547207953988, 0.01685540911503605, 0.001685539680499141, 0.0001685535981285378,
                     1.685497043444722e-05, 1.685289124599266e-06, 1.678938101250984e-07],
    "Error Node 5": [16129.53819924653, 1483.263487910087, 101.1521940589988, 10.07828946191161, 1.007460927209842,
                     0.1007424134903008, 0.01007420540869014, 0.001007420326301663, 0.0001007422525061959,
                     1.00742758822658e-05, 1.007916879253823e-06, 1.013182525500787e-07]
}

df = pd.DataFrame(data)
df = df.applymap(lambda x: f"{x:.2e}" if x >= 0.0001 else f"{x:.2e}")

# Exportieren als LaTeX
latex_code = df.to_latex(index=False, column_format='|c' *
                         (len(df.columns) + 1) + '|', escape=True)

# Speichern in einer Datei
with open('table_err.tex', 'w') as f:
    f.write(latex_code)

# plot a table
fig, ax = plt.subplots(figsize=(12, 4))
ax.axis('off')
ax.axis('tight')
table = ax.table(cellText=df.values, colLabels=df.columns,
                 cellLoc='center', loc='center', colLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.2)
for key, cell in table.get_celld().items():
    if key[0] == 0:
        cell.set_text_props(weight='bold')
plt.savefig(path_cwd + '/table_err.png')
plt.show()

# # # plot a graph
# plt.figure(figsize=(12, 4))
# plt.plot(df["t$^{max}_{tr}$"], df["Error Total"], label='Error Total')
# plt.plot(df["t$^{max}_{tr}$"], df["Error Node 1"], label='Error Node 1')
# plt.plot(df["t$^{max}_{tr}$"], df["Error Node 2"], label='Error Node 2')
# plt.plot(df["t$^{max}_{tr}$"], df["Error Node 3"], label='Error Node 3')
# plt.plot(df["t$^{max}_{tr}$"], df["Error Node 4"], label='Error Node 4')
# plt.plot(df["t$^{max}_{tr}$"], df["Error Node 5"], label='Error Node 5')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('t$^{max}_{tr}$')
# plt.ylabel('Error')
# plt.legend()
# plt.grid()
# plt.savefig(path_cwd + '/plot_err.png')
# plt.show()


# also a simple map plot
# map_data = geopandas.read_file(
#     os.path.join(
#         os.getcwd(),
#         'tools/vg2500_12-31.utm32s.shape/vg2500/VG2500_KRS.shp'))

# map_data = map_data[['ARS', 'GEN', 'NUTS', 'geometry']]

# county_names = ['Köln', 'Leverkusen',
#                 'Rheinisch-Bergischer Kreis', 'Rhein-Sieg-Kreis', 'Bonn']
# county_ids = ['05315', '05316', '05378', '05382', '05314']
# map_data = map_data[map_data['ARS'].isin(county_ids)]

# # for each df["t$^{max}_{tr}$"] i have the error for each county
# for i, t_max_tr in enumerate(df["t$^{max}_{tr}$"]):
#     plt.figure(figsize=(10, 10))
#     plt.title(f"Errors for t$^{{max}}_{{tr}} = {t_max_tr}$")

#     # Error values for this t$^{max}_{tr}$
#     errors = [
#         df["Error Node 1"][i],
#         df["Error Node 2"][i],
#         df["Error Node 3"][i],
#         df["Error Node 4"][i],
#         df["Error Node 5"][i],
#     ]

#     # Create a mapping of counties to their respective errors
#     error_mapping = dict(zip(county_names, errors))

#     # Add the error data to the GeoDataFrame
#     map_data['Error'] = map_data['GEN'].map(error_mapping)

#     # Plotting the map with the error data
#     ax = map_data.plot(column='Error', cmap='viridis', legend=True)
#     ax.set_axis_off()  # Turn off the axis
#     # Save the plot if needed
#     plt.savefig(path_cwd + f'/map_plot_tmax_{t_max_tr}.png')

#     # plt.show()

# # Calculate global min and max error values for the colorbar
# global_min = df.loc[:, ["Error Node 1", "Error Node 2",
#                         "Error Node 3", "Error Node 4", "Error Node 5"]].min().min()
# global_max = df.loc[:, ["Error Node 1", "Error Node 2",
#                         "Error Node 3", "Error Node 4", "Error Node 5"]].max().max()

# # Number of subplots
# n_subplots = len(df["t$^{max}_{tr}$"])
# n_cols = 3  # Number of columns in the grid
# # Calculate the required number of rows
# n_rows = (n_subplots + n_cols - 1) // n_cols

# fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols,
#                          figsize=(15, 10), constrained_layout=True)

# # Flatten the axes array to easily iterate over it
# axes = axes.flatten()

# # Plot each t$^{max}_{tr}$ in a subplot
# for i, (t_max_tr, ax) in enumerate(zip(df["t$^{max}_{tr}$"], axes)):
#     ax.set_title(f"t$^{{max}}_{{tr}}$ = {t_max_tr:.1e}")

#     # Error values for this t$^{max}_{tr}$
#     errors = [
#         df["Error Node 1"][i],
#         df["Error Node 2"][i],
#         df["Error Node 3"][i],
#         df["Error Node 4"][i],
#         df["Error Node 5"][i],
#     ]

#     # Create a mapping of counties to their respective errors
#     error_mapping = dict(zip(county_names, errors))

#     # Add the error data to the GeoDataFrame
#     map_data['Error'] = map_data['GEN'].map(error_mapping)

#     # Plotting the map with the error data using global min and max for color scaling (logarithmic)
#     map_data.plot(column='Error', cmap='viridis', legend=False, ax=ax,
#                   norm=LogNorm(vmin=global_min, vmax=global_max))
#     map_data.boundary.plot(ax=ax, linewidth=0.8, linestyle='--', color='black')

#     ax.set_axis_off()  # Turn off the axis

# for j in range(i + 1, n_rows * n_cols):
#     fig.delaxes(axes[j])

# # Create a single logarithmic colorbar for all subplots
# cbar = fig.colorbar(plt.cm.ScalarMappable(norm=LogNorm(vmin=global_min, vmax=global_max), cmap='viridis'),
#                     ax=axes, orientation='vertical', fraction=0.02, pad=0.04)
# cbar.set_label('Error (log scale)')

# # Save the plot grid if needed
# plt.savefig(path_cwd + '/grid_map_plots_log.png')
