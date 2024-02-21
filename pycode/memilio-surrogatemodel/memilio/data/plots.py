import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt

### get all neccessary plots
df_nodamp_MLP_part1 = pd.read_csv('pycode/memilio-surrogatemodel/memilio/secir_groups_nodamp_100days/datfarame_secirgroups_100days_nodamp')
df_nodamp_MLP_part2 = pd.read_csv('pycode/memilio-surrogatemodel/memilio/secir_groups_nodamp_100days/datfarame_secirgroups_100days_nodamp_MLP')
df_nodamp_MLP_part2['model']=['Dense', 'Dense', 'Dense']
df_nodamp_CNN_LSTM = pd.read_csv('pycode/memilio-surrogatemodel/memilio/secir_groups_nodamp_100days/datfarame_secirgroups_100days_nodamp_CNN_LSTM')
df_nodamp_100days = pd.concat([df_nodamp_MLP_part1 , df_nodamp_MLP_part2, df_nodamp_CNN_LSTM], axis=0)


df_damps_100days = pd.read_csv('pycode/memilio-surrogatemodel/memilio/secir_groups_2345damp_test/dataframe_secirgroups_2345damp_test')
df_damps_100days['dampings']=[2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5]

df_nodamp_100days.rename(columns={'model': 'Model'}, inplace=True)
df_nodamp_100days['dampings']=np.full(len(df_nodamp_100days),0)

# load data with one damping
df_dc_LSTM = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/secir_groups_onedamp_day_and_matrix/datarame_secirgroups_onedamp_noinfo_day_and_matrix_LSTM')
df_dc_CNN = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/secir_groups_onedamp_day_and_matrix/datarame_secirgroups_onedamp_noinfo_day_and_matrix_CNN')
df_dc_MLP = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/secir_groups_onedamp_day_and_matrix/datarame_secirgroups_onedamp_noinfo_day_and_matrix_MLP')

means_day_and_contact = []
for i in [df_dc_MLP, df_dc_CNN, df_dc_LSTM]:
            means_day_and_contact.append(i['MAPE'].mean().round(4))

df_one_damp = pd.DataFrame(data = {'Model':['Dense', 'CNN', 'LSTM'],'dampings':np.full(len(means_day_and_contact),1), 'MAPE':means_day_and_contact})

# cancat all three dfs together
df_all = pd.concat([df_damps_100days,df_nodamp_100days, df_one_damp])







# before implementing dampings we wanted a baseline value: Test MAPE for 100 days prediction without dampings
def barplot(df):
    df_mean= df.groupby('model').mean()['MAPE']
    fig, ax = plt.subplots()
    ax.bar(df_mean.index.values, df_mean.values)
    ax.set_ylabel('Test MAPE')
    ax.set_xlabel('Model')
    ax.set_title('Test MAPE per model for 100 days prediction ')
    ax.bar_label(ax.containers[0], label_type='edge')
    # pad the spacing between the number and the edge of the figure
    ax.margins(y=0.1)

    plt.show()
    plt.savefig("secir_groups_bar_100days.png")


def plot_multipledampings(df):
    df_grouped = pd.DataFrame(data=df[['MAPE','Model', 'dampings']].groupby(['Model', 'dampings']).mean())
    df_grouped.reset_index(inplace=True)
    df_plot = df_grouped.pivot(index='Model', columns='dampings', values='MAPE')

    MAPE = {
        'MLP': df_plot.loc['Dense'].values.round(2).tolist(), 
        'CNN':df_plot.loc['CNN'].values.round(2).tolist(), 
        'LSTM': df_plot.loc['LSTM'].values.round(2).tolist(), 
    }

    
    #df_bar=df_opt[['optimizer',  'kfold_test']]


    x = df_grouped['dampings'].unique() # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained')

    for attribute, measurement in MAPE.items():
        offset = width * multiplier
        rects = ax.bar(x +offset, measurement, width, label=attribute)
        ax.bar_label(rects, padding=3, fontsize = 8)
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.

    #ax.set_yticks(x+width, layers)
    ax.legend(loc='upper left', ncols=3, fontsize="8",)

    ax.set_ylabel('test MAPE')
    ax.set_xlabel('Number of dampings')
    ax.set_title('Mean Test MAPE for different number of dampings')
    plt.savefig("groups_multipledamp.png")


barplot(df_nodamp_100days)
plot_multipledampings(df_all)