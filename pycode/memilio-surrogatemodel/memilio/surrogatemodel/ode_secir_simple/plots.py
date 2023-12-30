import os 
import tensorflow as tf 
import pickle 
from  memilio.surrogatemodel.ode_secir_simple.model import split_data, get_test_statistic
import numpy as np 
from memilio.simulation.secir import InfectionState
import matplotlib.pyplot as plt
import pandas as pd 
#load data 

path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')
    
filename = "data_secir_simple.pickle"

if not os.path.isfile(os.path.join(path_data, 'data_secir_simple.pickle')):
        ValueError("no dataset found in path: " + path_data)

file = open(os.path.join(path_data,filename), 'rb')

data = pickle.load(file)
data_splitted = split_data(data['inputs'], data['labels'])

test_inputs = data_splitted['test_inputs']
test_labels = data_splitted['test_labels']


df_gridsearch = pd.read_csv('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/secir_simple_grid_searchdataframe')
# load trained model 
#new_model = tf.keras.models.load_model('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple')
new_model = tf.keras.models.load_model('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_150days')

pred = new_model.predict(test_inputs)

# transform arrays, so we assign the data to the right compartment 
# for our plot we only need one sample, not all 1000 test samples
sample_id = 0
pred_transformed =  pred[sample_id].transpose()
test_labels_transformed =  np.array(test_labels)[sample_id].transpose()

# we have to reverse the log normalization if we want to show real numbers in our plot 
pred_reversed = np.expm1(pred_transformed)
labels_reversed = np.expm1(test_labels_transformed)


def lineplots_pred_label(pred_reversed, labels_reversed):
        pred = pred_reversed
        labels = labels_reversed

        compartment_array = []
        for compartment in InfectionState.values():
                compartment_array.append(compartment) 
        index=[str(compartment).split('.')[1]
               for compartment in compartment_array]
        
        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(10,13))
        
        #fig, axes = plt.subplots(nrows=2, ncols=4, sharey=False)
        axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
        for ax, c, p, l in zip(axes, index, pred, labels):
                ax.plot(p)
                ax.plot(l)
                ax.set_title(c, fontsize = 10)
                
                ax.legend(loc='upper right', ncols=3)
        
        ax7.set_xlabel('Time')
        ax8.set_xlabel('Time')





        fig.suptitle('Predicted values and labels for compartments', fontsize=16)
        fig.legend(loc='upper right', ncols=3)
                
        plt.savefig("secir_simple_compartment_lines.png")

        
def lineplot_number_of_days():
        model30 = tf.keras.models.load_model('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple')
        model60 = tf.keras.models.load_model('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_60days')
        model90 = tf.keras.models.load_model('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_90days')
        model120 = tf.keras.models.load_model('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_120')
        model150 = tf.keras.models.load_model('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_150days')

        models = [model30, model60, model90, model120, model150]
        filenames =["data_secir_simple.pickle", "data_secir_simple_60days.pickle","data_secir_simple_90days.pickle","data_secir_simple_120days.pickle","data_secir_simple_150days.pickle"]
        days = [30,60,90,120,150]
        MAPE = []
        for model, file in zip(models, filenames): 
                model = model 
                path = os.path.dirname(os.path.realpath(__file__))
                path_data = os.path.join(os.path.dirname(os.path.realpath(
                        os.path.dirname(os.path.realpath(path)))), 'data')
                
                filename = file

                if not os.path.isfile(os.path.join(path_data, 'data_secir_simple.pickle')):
                        ValueError("no dataset found in path: " + path_data)

                file = open(os.path.join(path_data,filename), 'rb')

                data = pickle.load(file)
                data_splitted = split_data(data['inputs'], data['labels'])

                test_inputs = data_splitted['test_inputs']
                test_labels = data_splitted['test_labels']

                mape = get_test_statistic(test_inputs, test_labels, model)
                mean_mape = mape.mean()
                MAPE.append(mean_mape)


                
      

                
        plt.figure().clf()
        fig, ax = plt.subplots()
        ax.plot(days, MAPE,  marker = 'o' )
        ax.set_xticks(days)
        ax.set_xlabel('Number of days')
        ax.set_ylabel('MAPE loss')
        ax.set_title('MAPE loss for number of days to be predicted')
        plt.savefig("plot_days_secirsimple.png")


def heatmap(df_gridsearch):
    df = df_gridsearch

    plt.figure().clf() 
    df_heatmap1 = pd.DataFrame(data =  df.loc[(df['model'] == 'Dense')][['number_of_hidden_layers', 'number_of_neurons', 'kfold_test']])
    df_heatmap1= df_heatmap1.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='kfold_test')

    df_heatmap2 = pd.DataFrame(data =  df.loc[(df['model'] == 'CNN')][['number_of_hidden_layers', 'number_of_neurons', 'kfold_test']])
    df_heatmap2= df_heatmap2.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='kfold_test')

    df_heatmap3 = pd.DataFrame(data =  df.loc[(df['model'] == 'LSTM')][['number_of_hidden_layers', 'number_of_neurons', 'kfold_test']])
    df_heatmap3= df_heatmap3.pivot(index='number_of_hidden_layers', columns='number_of_neurons', values='kfold_test')

    fig, axs = plt.subplots(nrows = 2, ncols = 2, sharex=False, figsize = (20,20), constrained_layout = True)

    for ax, df_heatmap, name  in zip(axs.flat, [df_heatmap1 ,df_heatmap2, df_heatmap3], ['MLP', 'CNN', 'LSTM']):
        
        im = ax.imshow(df_heatmap.values, cmap ='Blues_r' )
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
    plt.savefig("heatmap_layers_neurons_all_secir_simple.png")
