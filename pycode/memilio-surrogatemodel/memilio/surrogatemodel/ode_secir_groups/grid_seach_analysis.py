import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib as mpl
import matplotlib.colors as colors


df = pd.read_csv("/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/secir_groups_grid_search/dataframes_concatenated_groups")

df_hyper = pd.read_csv("/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/secir_groups_LSTM_hyperparamter/dataframe_0_512_LSTM_opt")

def heatmap(df):

    # create a colormap ranging from green to red --> green for low MAPE, red gor high MAPE 


    # This dictionary defines the colormap
    cdict = {'red':  ((0.0, 0.0, 0.0),   # no red at 0
                    (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
                    (1.0, 0.8, 0.8)),  # set to 0.8 so its not too bright at 1

            'green': ((0.0, 0.8, 0.8),   # set to 0.8 so its not too bright at 0
                    (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
                    (1.0, 0.0, 0.0)),  # no green at 1

            'blue':  ((0.0, 0.0, 0.0),   # no blue at 0
                    (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
                    (1.0, 0.0, 0.0))   # no blue at 1
        }

    # Create the colormap using the dictionary
    GnRd = colors.LinearSegmentedColormap('GnRd', cdict)


    plt.figure().clf() 
    df_heatmap1 = pd.DataFrame(data =  df.loc[(df['model'] == 'Dense')][['number_of_hidden_layers', 'number_of_neurons', 'kfold_test']])
    df_heatmap1= df_heatmap1.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='kfold_test')

    df_heatmap2 = pd.DataFrame(data =  df.loc[(df['model'] == 'CNN')][['number_of_hidden_layers', 'number_of_neurons', 'kfold_test']])
    df_heatmap2= df_heatmap2.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='kfold_test')

    df_heatmap3 = pd.DataFrame(data =  df.loc[(df['model'] == 'LSTM')][['number_of_hidden_layers', 'number_of_neurons', 'kfold_test']])
    df_heatmap3= df_heatmap3.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='kfold_test')

    fig, axs = plt.subplots(nrows = 2, ncols = 2, sharex=False, figsize = (20,20), constrained_layout = True)

    for ax, df_heatmap, name  in zip(axs.flat, [df_heatmap1 ,df_heatmap2, df_heatmap3], ['MLP', 'CNN', 'LSTM']):
        
        im = ax.imshow(df_heatmap.values, cmap ='Greens_r' )
        plt.rcParams.update({'font.size': 30})
        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(len(df_heatmap.columns)), labels=df_heatmap.columns, fontsize = 25)
        ax.set_yticks(np.arange(len(df_heatmap.index)), labels=df_heatmap.index, fontsize = 25)

        ax.set_ylabel('number of hidden layers', fontsize = 25)
        ax.set_xlabel('number of neurons per layer', fontsize = 25)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                rotation_mode="anchor" )
        

        # Loop over data dimensions and create text annotations.
        for i in range(len(df_heatmap.index)):
            for j in range(len(df_heatmap.columns)):
                text = ax.text(j, i, np.around(df_heatmap.values, decimals=2)[i, j],
                            ha="center", va="center", color="k", fontsize = 25)
                
        ax.set_title('Model = '+name, fontsize = 30) 


        # cbar_kw = {}        
        # cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
        # cbar.ax.set_ylabel('MAPE', rotation=-90, va="bottom")
               

    #cax,kw = mpl.colorbar.make_axes([ax for ax in axs.flat])
    #plt.colorbar(im, cax=cax, **kw)
    #fig.subplots_adjust(right=0.8)
    #fig.colorbar(im, ax=axs.ravel().tolist(), location = 'right')  
    fig.colorbar(im, ax = axs, shrink=0.75, label = 'Test MAPE')
    #fig.tight_layout()
    fig.delaxes(axs[1][1])
    #plt.subplots_adjust(wspace=0.1, hspace=0.1)
    
    plt.show()
    plt.savefig("heatmap_layers_neurons_all_secir_groups_test.png")




# mean MAPE for number of layers


def barplot(df):
    plt.figure().clf() 
   

    df_1 = pd.DataFrame(data =  df.loc[(df['model'] == 'Dense')][['number_of_hidden_layers', 'number_of_neurons', 'kfold_test']])
    df_1= df_1.pivot(index='number_of_neurons', columns='number_of_hidden_layers', values='kfold_test')
    
    df_2 = pd.DataFrame(data =  df.loc[(df['model'] == 'CNN')][['number_of_hidden_layers', 'number_of_neurons', 'kfold_test']])
    df_2= df_2.pivot(index='number_of_neurons', columns='number_of_hidden_layers', values='kfold_test')
    
    df_3 = pd.DataFrame(data =  df.loc[(df['model'] == 'LSTM')][['number_of_hidden_layers', 'number_of_neurons', 'kfold_test']])
    df_3= df_3.pivot(index='number_of_neurons', columns='number_of_hidden_layers', values='kfold_test')
    

    MAPE = {
        'MLP': df_1.mean().values.round(2).tolist(), 
        'CNN':df_2.mean().values.round(2).tolist(), 
        'LSTM':df_3.mean().values.round(2).tolist(), 
    }

    layers = df_1.mean().index.values
    #df_bar=df_opt[['optimizer',  'kfold_test']]


    x = np.arange(len(layers))  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained')

    for attribute, measurement in MAPE.items():
        offset = width * multiplier
        rects = ax.bar(x +offset, measurement, width, label=attribute)
        ax.bar_label(rects, padding=3)
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.

    #ax.set_yticks(x+width, layers)
    ax.legend(loc='upper right', ncols=3)

    ax.set_ylabel('test MAPE')
    ax.set_xlabel('Number of layers')
    ax.set_title('Mean Test MAPE for different number of layers')
    plt.savefig("bar_groups_layer_correctTestValues.png")




# mean MAPE for number of neurons per layer
def lineplot(df):
       
    plt.figure().clf() 
    df_1 = pd.DataFrame(data =  df.loc[(df['model'] == 'Dense')][['number_of_hidden_layers', 'number_of_neurons', 'kfold_test']])
    df_1= df_1.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='kfold_test')
    
    df_2 = pd.DataFrame(data =  df.loc[(df['model'] == 'CNN')][['number_of_hidden_layers', 'number_of_neurons', 'kfold_test']])
    df_2= df_2.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='kfold_test')
    
    df_3 = pd.DataFrame(data =  df.loc[(df['model'] == 'LSTM')][['number_of_hidden_layers', 'number_of_neurons', 'kfold_test']])
    df_3= df_3.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='kfold_test')
    

    MAPE = {
        'MLP': df_1.mean().values.round(2).tolist(), 
        'CNN':df_2.mean().values.round(2).tolist(), 
        'LSTM':df_3.mean().values.round(2).tolist(), 
    }

    fig, ax = plt.subplots(figsize=(8, 5), layout='constrained')
    x = df_1.mean().index.values
    
    ax.plot(x,df_1.mean().values.round(2).tolist(), label='MLP') 
    ax.plot(x,df_2.mean().values.round(2).tolist(), label='CNN') 
    ax.plot(x,df_3.mean().values.round(2).tolist(), label='LSTM') 


    ax.set_xticks(x, labels=x, fontsize = 12)
    

    ax.set_ylabel('Test MAPE', fontsize = 15)
    ax.set_xlabel('number of neurons per layer', fontsize = 15)


    ax.legend()  # Add a legend.

    ax.set_ylabel('test MAPE')
    ax.set_xlabel('Number of Neurons')
    ax.set_title('Mean Test MAPE for different number of neurons')
    plt.savefig("line_groups_neurons_correctTestValues.png")





def simple_barplot(df):
    plt.figure().clf() 
    df_bar=df[['optimizer',  'kfold_test']]
    df_bar.sort_values(by='kfold_test', inplace = True, ascending=False)
    

    #df_bar=df_opt[['optimizer',  'kfold_test']]

    fig, ax = plt.subplots()

    rects = ax.barh(df_bar['optimizer'], df_bar['kfold_test'].round(4))

    ax.set_ylabel('Test MAPE')
    ax.set_xlabel('Opitmizer')
    ax.set_title('Optimizer for LSTM ')


    large_bars = [p if p > 2 else '' for p in df_bar['kfold_test'].round(4)]
    small_bars = [p if p <= 2 else '' for p in df_bar['kfold_test'].round(4)]


    ax.bar_label(rects, small_bars,
                  padding=5, color='black')
    ax.bar_label(rects, large_bars,
                  padding=-40, color='white')

    #ax.bar_label(rects[:4], ax.containers[:4], label_type='edge')
    #ax.bar_label(rects, ax.containers[4:], label_type='edge', padding=-32, color='white', fontweight='bold')
    # pad the spacing between the number and the edge of the figure
    ax.margins(y=0.1)

    plt.show()
    plt.savefig("barh_optimizer_LSTM.png")
    ax.set_ylabel('Optimizer')
    ax.set_xlabel('Test MAPe')
    ax.set_title('Test MAPE for different optimizer')
    plt.savefig("secirgroups_optimizer_LSTM_barplot.png")




############# other plots ##############
# plot for one damp

def plot_one_damp_secir_groups():
    # load data from models with damping but no information about damping provided 
    # we performed five runs and will utiliz the average value 
    # one random damping, 30 days prediction, baseline matrix is correct, minimum matrix is 0 

    df_noinfo_CNN = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/secir_groups_onedamp_noinfo/datfarame_secirgroups_onedamp_noinfo_CNN')
    df_noinfo_LSTM = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/secir_groups_onedamp_noinfo/datfarame_secirgroups_onedamp_noinfo')
    df_noinfo_MLP = pd.read_csv('pycode/memilio-surrogatemodel/memilio/secir_groups_onedamp_noinfo/datfarame_secirgroups_onedamp_noinfo_MLP')

    df_dc_LSTM = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/secir_groups_onedamp_day_and_matrix/datarame_secirgroups_onedamp_noinfo_day_and_matrix_LSTM')
    df_dc_CNN = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/secir_groups_onedamp_day_and_matrix/datarame_secirgroups_onedamp_noinfo_day_and_matrix_CNN')
    df_dc_MLP = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/secir_groups_onedamp_day_and_matrix/datarame_secirgroups_onedamp_noinfo_day_and_matrix_MLP')

    means_noinfo= []
    for i in [df_noinfo_MLP, df_noinfo_CNN, df_noinfo_LSTM]:
        means_noinfo.append(i['MAPE'].mean().round(4))

    means_day_and_contact = []
    for i in [df_dc_MLP, df_dc_CNN, df_dc_LSTM]:
            means_day_and_contact.append(i['MAPE'].mean().round(4))

    models = ("MLP", "CNN", "LSTM")
    penguin_means = {
        'Model 1': means_noinfo,
        'Model 2': means_day_and_contact}

    x = np.arange(len(models))  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained')

    for attribute, measurement in penguin_means.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute)
        ax.bar_label(rects, padding=3)
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Test MAPE')
    ax.set_title('Model comparison for prediction with variable damping')
    ax.set_xticks(x + width, models)
    ax.legend(loc='upper right', ncols=3)
    plt.savefig("secirgroups_onedamp_modelcomparison.png")





heatmap(df)
barplot(df)
simple_barplot(df_hyper)
