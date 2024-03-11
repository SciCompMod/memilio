import os 
import tensorflow as tf 
import pickle 
from  memilio.surrogatemodel.ode_secir_simple.model import split_data, get_test_statistic
from memilio.surrogatemodel.ode_secir_groups.data_generation_nodamp import get_population
import numpy as np 
from memilio.simulation.secir import InfectionState
import matplotlib.pyplot as plt
import pandas as pd 
#load data 
path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')
    
filename = "data_secir_simple_150days.pickle"

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
new_model = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_simple_bestLSTM_2024_150days')

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
                ax.plot(l, label ='label')
                ax.plot(p, label='pred', linestyle = '--')

                ax.set_title(c, fontsize = 10)
                
                #ax.legend(loc='upper right', ncols=3)
        
        ax7.set_xlabel('Time')
        ax8.set_xlabel('Time')


        lines = [] 
        line_labels = [] 
        for ax in fig.axes: 
                Line, Label = ax.get_legend_handles_labels() 
                # print(Label) 
                lines.extend(Line) 
                line_labels.extend(Label) 

        fig.legend(lines[:2], line_labels[:2], loc='lower center') 
        fig.suptitle('Predicted values and labels for compartments', fontsize=16)
        #fig.legend(loc='upper right', ncols=3)
        
                
        plt.savefig("secir_simple_compartment_lines.png")


def closeup_susceptible(pred_reversed, labels_reversed):
        pred = pred_reversed
        labels = labels_reversed

        compartment_array = []
        for compartment in InfectionState.values():
                compartment_array.append(compartment) 
        index=[str(compartment).split('.')[1]
               for compartment in compartment_array]
        index_s = [0]
        pred_s = pred[0][:20] #only 20 first days of susceptible

        days = np.arange(1,21)
        plt.figure().clf()
        fig, ax = plt.subplots()
        ax.plot(days, labels[0][:20], label = 'labels')
        ax.plot(days, pred_s, label = 'pred', linestyle = '--')
        
        ax.set_xticks(days)
        ax.set_xlabel('Number of days')
        ax.set_ylabel('MAPE loss')
        ax.legend(loc='upper right')
        ax.set_title('Closeup of first 20 days of Susceptible compartment')
        plt.savefig("susceptible_closeup_secir_simple_dotted.png")

        
def lineplot_number_of_days():
       #model30 = tf.keras.models.load_model('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple')
        model60 = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_60days')
        model90 = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_90days')
        model120 = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_120days')
        model150 = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_150days')

        models = [model60, model90, model120, model150]
        filenames =[ "data_secir_simple_60days.pickle","data_secir_simple_90days.pickle","data_secir_simple_120days.pickle","data_secir_simple_150days.pickle"]
        days = [30,60,90,120,150]
        MAPE = []
        MAPE.append(0.1162)
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
                mean_mape = mape.mean()[0]
                MAPE.append(mean_mape)
     
                
        plt.figure().clf()
        fig, ax = plt.subplots(figsize =(12,7))
        ax.plot(days, MAPE,  marker = 'o' )
        ax.set_xticks(days)
        ax.set_xlabel('Number of days')
        ax.set_ylabel('MAPE')
        ax.set_title('MAPE for number of days to be predicted')
        plt.savefig("plot_days_secirsimple_long.png")


def lineplot_number_of_days_secir_groups():
        model30 = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_groups_best_LSTM')
        model60 = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_groups_best_LSTM_60days')
        model90 = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_groups_best_LSTM_90days')
        model120 = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_groups_best_LSTM_120days')
        model150 = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_groups_best_LSTM_150days')

        models = [model30, model60, model90, model120, model150]
        filenames =[ 'data_secir_groups_30days_nodamp.pickle',"data_secir_groups_60days_nodamp.pickle","data_secir_groups_90days_nodamp.pickle",
                    "data_secir_groups_120days_nodamp.pickle","data_secir_groups_150days_nodamp.pickle"]
        days = [30,60,90,120,150]
        MAPE = []

        for model, file in zip(models, filenames): 
                model = model 
                path = os.path.dirname(os.path.realpath(__file__))
                path_data = os.path.join(os.path.dirname(os.path.realpath(
                        os.path.dirname(os.path.realpath(path)))), 'data')
                
                filename = file

                file = open(os.path.join(path_data,filename), 'rb')

                data = pickle.load(file)
                data_splitted = split_data(data['inputs'], data['labels'])

                test_inputs = data_splitted['test_inputs']
                test_labels = data_splitted['test_labels']

                mape = get_test_statistic(test_inputs, test_labels, model)
                mean_mape = mape.mean()[0]
                MAPE.append(mean_mape)
     
                
        plt.figure().clf()
        fig, ax = plt.subplots(figsize =(12,7))
        ax.plot(days, MAPE,  marker = 'o' )
        ax.set_xticks(days)
        ax.set_xlabel('Number of days')
        ax.set_ylabel('MAPE')
        ax.set_title('MAPE for number of days to be predicted')
        plt.savefig("plot_days_secirgroups_long.png")



def lineplot_days_simple_and_groups():
        path = os.path.dirname(os.path.realpath(__file__))
        path_data = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'data')
        
        filenames_simple=['data_secir_simple.pickle', "data_secir_simple_60days.pickle","data_secir_simple_90days.pickle","data_secir_simple_120days.pickle",
                          "data_secir_simple_150days.pickle"]
        filenames_groups = ['data_secir_groups_30days_nodamp.pickle',"data_secir_groups_60days_nodamp.pickle","data_secir_groups_90days_nodamp.pickle",
                    "data_secir_groups_120days_nodamp.pickle","data_secir_groups_150days_nodamp.pickle"]
        
        path = os.path.dirname(os.path.realpath(__file__))
        path_models = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'saved_models')
        
        modelnames_simple = ['saved_models_secir_simple_bestLSTM_2024','saved_models_secir_simple_bestLSTM_2024_60days', 
                             'saved_models_secir_simple_bestLSTM_2024_90days', 'saved_models_secir_simple_bestLSTM_2024_120days',
                             'saved_models_secir_simple_bestLSTM_2024_150days']
        modelnames_groups = ['saved_models_secir_groups_best_LSTM', 'saved_models_secir_groups_best_LSTM_60days', 'saved_models_secir_groups_best_LSTM_90days',
                             'saved_models_secir_groups_best_LSTM_120days', 'saved_models_secir_groups_best_LSTM_150days']
        
        days = [30,60,90,120,150]

        MAPE_simple = []
        MAPE_groups = []
        for filenames, modelnames, MAPE in zip([filenames_simple, filenames_groups], [modelnames_simple, modelnames_groups],[MAPE_simple, MAPE_groups]):
                
                for modelname, filename in zip(modelnames, filenames): 

                        
                        file = open(os.path.join(path_data,filename), 'rb')
                        model = tf.keras.models.load_model(os.path.join(path_models, modelname))
                        

                        data = pickle.load(file)
                        data_splitted = split_data(data['inputs'], data['labels'])

                        test_inputs = data_splitted['test_inputs']
                        test_labels = data_splitted['test_labels']

                        mape = get_test_statistic(test_inputs, test_labels, model)
                        mean_mape = mape.mean()[0]
                        MAPE.append(mean_mape)

        # version 1 
        # plt.figure().clf()
        # fig, ax = plt.subplots(figsize =(8,4))
        # ax.plot(days, MAPE_simple,  marker = 'o', label='simple' )
        # ax.plot(days, MAPE_groups,  marker = 'x', label='groups' )
        # ax.set_xticks(days)
        # ax.set_xlabel('Number of days')
        # ax.set_ylabel('MAPE')
        # ax.set_title('MAPE for number of days to be predicted')
        # plt.savefig("plot_days_simple_and_groups.png")

        # version two with two axes 
                               
        plt.figure().clf()
        fig, ax1 = plt.subplots(figsize =(8,4))

        color = 'tab:red'
        ax1.plot(days, MAPE_simple, label ='simple', color = color, marker = 'o' )
        ax1.set_xticks(days)
        ax1.set_xlabel('Number of days')
        ax1.set_ylabel('MAPE', color = color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.set_title('MAPE for number of days to be predicted')
        


        ax2 = ax1.twinx()
        color = 'tab:blue'
        ax2.set_ylabel('MAPE', color = color)
        ax2.plot(days, MAPE_groups, color = color, label = 'groups', marker = 'x')
        ax2.tick_params(axis='y', labelcolor=color)
        fig.legend(loc= 'lower left', fontsize = 9)
        fig.tight_layout()
        

        plt.savefig("plot_days_simple_and_groups_v2.png")





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


def plot_losses(): # for secir groups
        path = os.path.dirname(os.path.realpath(__file__))
        path_data = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'data')
        
        #filenames = ["data_secir_simple_150days.pickle", 'data_secir_groups_150days_nodamp.pickle']
        filename = 'data_secir_groups_150days_nodamp.pickle'

      
        path_models = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'saved_models')
        
        #modelnames = ['saved_models_secir_simple_bestLSTM_2024_150days', 'saved_models_secir_groups_best_LSTM_150days' ]
        modelname = 'saved_models_secir_groups_best_LSTM_150days'
        
        #for filename, modelname in zip(filenames, modelnames):

        if not os.path.isfile(os.path.join(path_data, 'data_secir_simple.pickle')):
                        ValueError("no dataset found in path: " + path_data)

        file = open(os.path.join(path_data,filename), 'rb')

        data = pickle.load(file)
        data_splitted = split_data(data['inputs'], data['labels'])

        test_inputs = data_splitted['test_inputs']
        test_labels = data_splitted['test_labels']

                # load trained model 
                
        new_model = tf.keras.models.load_model(os.path.join(path_models, modelname))

        pred = new_model.predict(test_inputs)

        # we have to reverse the log normalization if we want to show real numbers in our plot 
        pred_reversed = np.expm1(pred)
        labels_reversed = np.expm1(test_labels)

        MAE_150 = []
        for abse in (abs(labels_reversed - pred_reversed).reshape([150,1000,48])):
                        MAE_150.append(np.mean(abse))

        from tensorflow.python.ops.numpy_ops import np_config
        np_config.enable_numpy_behavior()

        MAPE_150 = []
        for abse in ((abs((test_labels - pred)/test_labels)).reshape([150,1000,48])):
                        MAPE_150.append(100*np.mean(abse))
                
        days = np.arange(0,150)
        plt.figure().clf()
        fig, ax1 = plt.subplots(figsize =(10,5))

        color = 'tab:red'
        ax1.plot(days, MAE_150, label ='MSE', color = color )
        ax1.set_xticks(days)
        ax1.set_xlabel('Number of days')
        ax1.set_ylabel('MAE', color = color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.set_xticks(np.arange(0,150,5))
        ax1.set_title('MAE and MAPE for each day')

        ax2 = ax1.twinx()
        color = 'tab:blue'
        ax2.set_ylabel('MAPE', color = color)
        ax2.plot(days, MAPE_150, color = color)
        ax2.tick_params(axis='y', labelcolor=color)
        fig.tight_layout()

        plt.savefig("MAE_MAPE_150_secir_groups.png")

        
# take closer look on day 56 which has highest MSE and MAPE 
pred_56 = pred.reshape([150,1000,48])[56]
labels_56 = test_labels.reshape([150,1000,48])[56]     

MAPE_56 = 100*np.mean(abs((labels_56 - pred_56)/labels_56))



# compartment error plots 
# for secir simple
def compartment_error_simple():
                #load data 
        path = os.path.dirname(os.path.realpath(__file__))
        path_data = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'data')
        
        filename = "data_secir_simple_150days.pickle"

        file = open(os.path.join(path_data,filename), 'rb')

        data = pickle.load(file)
        data_splitted = split_data(data['inputs'], data['labels'])

        test_inputs = data_splitted['test_inputs']
        test_labels = data_splitted['test_labels']

        new_model = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_simple_bestLSTM_2024_150days')

        pred = new_model.predict(test_inputs)

        pred_reversed = np.expm1(pred)
        labels_reversed = np.expm1(test_labels)

        compartment_array = []
        for compartment in InfectionState.values():
                compartment_array.append(compartment) 
        index=[str(compartment).split('.')[1]
               for compartment in compartment_array]
        
        mae =  np.mean((abs(labels_reversed - pred_reversed)), axis = 0).transpose()
        mape = 100*np.mean(abs((test_labels - pred)/test_labels), axis = 0).transpose()
                
        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(10,13))
        
        #fig, axes = plt.subplots(nrows=2, ncols=4, sharey=False)
        axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
        for ax, c, ms, ma in zip(axes, index, mae, mape):
                
                color = 'tab:blue'
                ax.plot(ms, label ='MAE', color = color)
                ax.set_xlabel('Number of days')
                ax.set_ylabel('MAE', color = color)
                ax.tick_params(axis='y', labelcolor=color)
                ax.set_title(c, fontsize = 10)

                ax2 = ax.twinx()
                color = 'tab:green'
                ax2.set_ylabel('MAPE', color = color)
                ax2.plot(ma, color = color, label = 'MAPE' , linestyle = '--')
                ax2.tick_params(axis='y', labelcolor=color)
                fig.tight_layout()


                #ax.plot(ms, label ='MSE')
                #ax.plot(ma, label='MAPE', linestyle = '--')#

                #ax.set_title(c, fontsize = 10)
                
                #ax.legend(loc='upper right', ncols=3)
        
        ax7.set_xlabel('Days')
        ax8.set_xlabel('Days')


        #lines = [] 
        #line_labels = [] 
        #for ax in fig.axes: 
        #        Line, Label = ax.get_legend_handles_labels() 
        #        # print(Label) 
        #        lines.extend(Line) 
        #        line_labels.extend(Label) 

        #fig.legend(lines[:2], line_labels[:2], loc='lower center') 
        #fig.('MSE and MAPE for compartments', fontsize=16)                
        plt.savefig("150_simple_MSE_and_MAPE.png")


def compartment_error_groups():
        path = os.path.dirname(os.path.realpath(__file__))
        path_data = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'data')
        
        filename = "data_secir_groups_150days_nodamp.pickle"

        file = open(os.path.join(path_data,filename), 'rb')

        data = pickle.load(file)
        data_splitted = split_data(data['inputs'], data['labels'])

        test_inputs = data_splitted['test_inputs']
        test_labels = data_splitted['test_labels']

        new_model = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_groups_best_LSTM_150days')

        pred = new_model.predict(test_inputs)

        pred_reversed = np.expm1(pred)
        labels_reversed = np.expm1(test_labels)

        compartment_array = []
        for compartment in InfectionState.values():
                compartment_array.append(compartment) 
        index=[str(compartment).split('.')[1]
               for compartment in compartment_array]




        # calcuate average of all age gorups 
        #average_labels_reversed = []
        #for sample in labels_reversed = 
        
        mae =  np.mean((abs(labels_reversed - pred_reversed)), axis = 0).transpose()
        mape = 100*np.mean(abs((test_labels - pred)/test_labels), axis = 0).transpose()
                
        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(10,13))
        
        #fig, axes = plt.subplots(nrows=2, ncols=4, sharey=False)
        axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
        for ax, c, ms, ma in zip(axes, index, mae, mape):
                
                color = 'tab:blue'
                ax.plot(ms, label ='MAE', color = color)
                ax.set_xlabel('Number of days')
                ax.set_ylabel('MAE', color = color)
                ax.tick_params(axis='y', labelcolor=color)
                ax.set_title(c, fontsize = 10)

                ax2 = ax.twinx()
                color = 'tab:orange'
                ax2.set_ylabel('MAPE', color = color)
                ax2.plot(ma, color = color, label = 'MAPE' , linestyle = '--')
                ax2.tick_params(axis='y', labelcolor=color)
                fig.tight_layout()


                #ax.plot(ms, label ='MSE')
                #ax.plot(ma, label='MAPE', linestyle = '--')#

                #ax.set_title(c, fontsize = 10)
                
                #ax.legend(loc='upper right', ncols=3)
        
        ax7.set_xlabel('Days')
        ax8.set_xlabel('Days')


        lines = [] 
        line_labels = [] 
        for ax in fig.axes: 
                Line, Label = ax.get_legend_handles_labels() 
                # print(Label) 
                lines.extend(Line) 
                line_labels.extend(Label) 

        fig.legend(lines[:2], line_labels[:2], loc='lower center') 
        #fig.('MSE and MAPE for compartments', fontsize=16)                
        plt.savefig("150_simple_MSE_and_MAPE.png")

def plot_30days_differentmodels():
        path = os.path.dirname(os.path.realpath(__file__))
        path_data = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'data')
        
        filenames_simple=['data_secir_simple.pickle', "data_secir_simple_60days.pickle","data_secir_simple_90days.pickle","data_secir_simple_120days.pickle",
                          "data_secir_simple_150days.pickle"]

        path = os.path.dirname(os.path.realpath(__file__))
        path_models = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'saved_models')
        
        modelnames_simple = ['saved_models_secir_simple_bestLSTM_2024','saved_models_secir_simple_bestLSTM_2024_60days', 
                             'saved_models_secir_simple_bestLSTM_2024_90days', 'saved_models_secir_simple_bestLSTM_2024_120days',
                             'saved_models_secir_simple_bestLSTM_2024_150days']
        
        days = ['30days','60days','90days','120days','150days']


        compartment_array = []
        for compartment in InfectionState.values():
                        compartment_array.append(compartment) 
        index=[str(compartment).split('.')[1]
        for compartment in compartment_array]
       

        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(10,13))

        
        #fig, axes = plt.subplots(nrows=2, ncols=4, sharey=False)
        axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]

                
        for modelname, filename , name_plot in zip( modelnames_simple, filenames_simple, days):
                file = open(os.path.join(path_data,filename), 'rb')

                data = pickle.load(file)
                data_splitted = split_data(data['inputs'], data['labels'])

                test_inputs = data_splitted['test_inputs']
                test_labels = data_splitted['test_labels']

                new_model = tf.keras.models.load_model(os.path.join(path_models, modelname))

                pred = new_model.predict(test_inputs)

                pred_reversed = np.expm1(pred)
                labels_reversed = np.expm1(test_labels)

                mape = 100*np.mean(abs((test_labels - pred)/test_labels), axis = 0).transpose()
                for ax, c, ma in zip(axes, index, mape):
                        
                        ax.set_ylabel('Test MAPE')
                        ax.set_xlabel('Days')
                        ax.plot(ma[:30], label = name_plot)
                        ax.set_title(c, fontsize = 10)
                                
                        fig.tight_layout()

                
                ax7.set_xlabel('Days')
                ax8.set_xlabel('Days')

        lines = [] 
        line_labels = [] 
        for ax in fig.axes: 
                Line, Label = ax.get_legend_handles_labels() 
                # print(Label) 
                lines.extend(Line) 
                line_labels.extend(Label) 

        fig.legend(lines[:(len(days))], line_labels[:len(days)], loc='upper center',  bbox_to_anchor=(0.5, -0.05), shadow=True, ncol=3)
        #fig.('MSE and MAPE for compartments', fontsize=16)                
        plt.savefig("MAPE_30_days_allmodels_secirsimple.png")



def plot_inputs():
        path = os.path.dirname(os.path.realpath(__file__))
        path_data = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'data')
        filename = "data_secir_simple_90days.pickle"

        file = open(os.path.join(path_data,filename), 'rb')

        data = pickle.load(file)
        data_splitted = split_data(data['inputs'], data['labels'])
        train_inputs = np.transpose(data_splitted['train_inputs']) # reshape to [n_compartments, n_days, n_samples]

        p50 = []
        p25 = []        
        p75 = []
        p05 = []
        p95 = []

        percentiles = [50,25,75,5,95]
        arrays = [p50, p25, p75, p05, p95]
        labels = ['percentile p50', 'percentile p25', 'percentile p75', 'percentile p05', 'percentile p95']

        inputs_reversed = np.expm1(train_inputs)
        
        for comp in inputs_reversed:
                for day, a  in zip(comp, arrays):
                        array_ = []
                        for p in zip( percentiles): 
                                array_.append(np.percentile(day, p))
                        a.append(array_)

        compartment_array = []
        for compartment in InfectionState.values():
                        compartment_array.append(compartment) 
        index=[str(compartment).split('.')[1]
        for compartment in compartment_array]
       

        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(8,10))
        axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
        
        for array, label in zip(arrays, labels):
                array = np.squeeze(array)

                for ax, c, values in zip(axes, index, array):
                       
       
                #for values, label in zip(array, labels):                       
                        ax.set_ylabel('Number of individuals')
                        ax.set_xlabel('Days')
                        ax.plot(values, label = label)
                        ax.set_title(c, fontsize = 10)
                                                
                        fig.tight_layout()

                        
                        ax7.set_xlabel('Days')
                        ax8.set_xlabel('Days')

                # lines = [] 
                # line_labels = [] 
                # for ax in fig.axes: 
                #         Line, Label = ax.get_legend_handles_labels() 
                #         # print(Label) 
                #         lines.extend(Line) 
                #         line_labels.extend(Label) 

        #fig.legend(lines[:(len(days))], line_labels[:len(days)], loc='upper center',  bbox_to_anchor=(0.5, -0.05), shadow=True, ncol=3)
        fig.legend()
        plt.savefig("inputs_secir_simple.png")



def hist_plot_populations():
        path_population = os.path.abspath(
                r"data//pydata//Germany//county_population.json")
        populations = get_population(path_population )
        populations = np.asarray(populations).sum(axis = 1)

       
        plt.figure().clf()
        
        fig, axs= plt.subplots(1,2 ,figsize=(8,5))
        #axs[0].ticklabel_format(style='plain')
        axs[0].hist(populations, bins = 100)    
        axs[0].set_xlabel('Population size')
        axs[0].set_title('Histogram')
        axs[0].set_xticks(np.arange(0,3500000,500000))

        #axs[1].ticklabel_format(style='plain')
        axs[1].boxplot(populations, vert  = False)
        axs[1].set_title('Boxplot')
        axs[1].set_xticks(np.arange(0,3500000,500000))
        plt.savefig("population_hist_and_boxplot.png")

