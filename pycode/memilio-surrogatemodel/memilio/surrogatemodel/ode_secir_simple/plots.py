import os 
import tensorflow as tf 
import pickle 
from  memilio.surrogatemodel.ode_secir_simple.model import split_data, get_test_statistic
#from memilio.surrogatemodel.ode_secir_groups.data_generation_nodamp import get_population
import numpy as np 
from memilio.simulation.secir import InfectionState
import matplotlib.pyplot as plt
import pandas as pd 
import seaborn as sns
#load data 

path = os.path.dirname(os.path.realpath(__file__))
path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')
    
filename = "data_secir_simple_150days_w2.pickle"

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
new_model = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_150days_w')

pred = new_model.predict(test_inputs)

# transform arrays, so we assign the data to the right compartment 
# for our plot we only need one sample, not all 1000 test samples
sample_id = 894
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
        infectionstates = ['Susceptible','Exposed', 'InfectedNoSymptoms', 'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']
        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(10,13), constrained_layout = True)
        
        #fig, axes = plt.subplots(nrows=2, ncols=4, sharey=False)
        axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
        for ax, c, p, l in zip(axes, infectionstates, pred, labels):
                ax.plot(l, label ='label')
                ax.plot(p, label='pred', linestyle = '--')
                ax.set_xlabel('Number of days')
                ax.set_ylabel('Number of individuals')
                
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
        
                
        plt.savefig("secir_simple_compartment_lines_w.png")


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
        plt.savefig("susceptible_closeup_secir_simple_dotted_w.png")

        
def lineplot_number_of_days():
       #model30 = tf.keras.models.load_model('/home/schm_a45/Documents/code3/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple')
        model60 = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_simple_60days')
        model90 = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_simple_90days')
        model120 = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_simple_120days')
        model150 = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_simple_150days')

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
        
        filenames_simple=['data_secir_simple_30days_widerinput_2.pickle', "data_secir_simple_60days_w2.pickle","data_secir_simple_90days_w2.pickle","data_secir_simple_120days_w2.pickle",
                          "data_secir_simple_150days_w2.pickle"]
        filenames_groups = ['data_secir_groups_30days_nodamp_w.pickle',"data_secir_groups_60days_nodamp_w.pickle","data_secir_groups_90days_nodamp_w.pickle",
                    "data_secir_groups_120days_nodamp_w.pickle","data_secir_groups_150days_nodamp_w.pickle"]
        
        # filenames_simple=['data_secir_simple.pickle', "data_secir_simple_60days.pickle","data_secir_simple_90days.pickle","data_secir_simple_120days.pickle",
        #                   "data_secir_simple_150days.pickle"]
        # filenames_groups = ['data_secir_groups_30days_nodamp.pickle',"data_secir_groups_60days_nodamp.pickle","data_secir_groups_90days_nodamp.pickle",
        #             "data_secir_groups_120days_nodamp.pickle","data_secir_groups_150days_nodamp.pickle"]
        
        path = os.path.dirname(os.path.realpath(__file__))
        path_models = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'saved_models')
        
        # modelnames_simple = ['saved_models_secir_simple_bestLSTM_2024','saved_models_secir_simple_bestLSTM_2024_60days', 
        #                      'saved_models_secir_simple_bestLSTM_2024_90days', 'saved_models_secir_simple_bestLSTM_2024_120days',
        #                      'saved_models_secir_simple_bestLSTM_2024_150days']
        # modelnames_groups = ['saved_models_secir_groups_best_LSTM', 'saved_models_secir_groups_best_LSTM_60days', 'saved_models_secir_groups_best_LSTM_90days',
        #                      'saved_models_secir_groups_best_LSTM_120days', 'saved_models_secir_groups_best_LSTM_150days']

        modelnames_simple = ['saved_models_secir_simple_30days_w','saved_models_secir_simple_60days_w', 
                             'saved_models_secir_simple_90days_w', 'saved_models_secir_simple_120days_w',
                             'saved_models_secir_simple_150days_w']
        modelnames_groups = ['saved_models_groups_nodamp_30days_w', 'saved_models_groups_nodamp_60days_w', 'saved_models_groups_nodamp_90days_w',
                             'saved_models_groups_nodamp_120days_w', 'saved_models_groups_nodamp_150days_w']
        
        days = [30,60,90,120,150]

        #MAPE_simple = [0.0922, 0.1485, 0.2049, 0.1680, 0.3642]
        #MAPE_groups = [0.4123, 0.3225, 0.4032, 0.4176, 0.7222]
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
        

        plt.savefig("plot_days_simple_and_groups_v2_w.png")





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
        filename = 'data_secir_groups_150days_nodamp_w.pickle'

      
        path_models = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'saved_models')
        
        #modelnames = ['saved_models_secir_simple_bestLSTM_2024_150days', 'saved_models_secir_groups_best_LSTM_150days' ]
        modelname = 'saved_models_groups_nodamp_150days_w'
        
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

        plt.savefig("MAE_MAPE_150_secir_groups_w.png")

        
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
        
        filename = "data_secir_simple_150days_w2.pickle"

        file = open(os.path.join(path_data,filename), 'rb')

        data = pickle.load(file)
        data_splitted = split_data(data['inputs'], data['labels'])

        test_inputs = data_splitted['test_inputs']
        test_labels = data_splitted['test_labels']

        new_model = tf.keras.models.load_model('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_simple_150days_w')

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
        infectionstates = ['Susceptible','Exposed', 'InfectedNoSymptoms', 'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']
        for ax, c, ms, ma in zip(axes, infectionstates, mae, mape):
                
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
        plt.savefig("150_simple_MAE_and_MAPE_w.png")


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



        infectionstates = ['Susceptible','Exposed', 'InfectedNoSymptoms', 'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']
        # calcuate average of all age gorups 
        #average_labels_reversed = []
        #for sample in labels_reversed = 
        
        mae =  np.mean((abs(labels_reversed - pred_reversed)), axis = 0).transpose()
        mape = 100*np.mean(abs((test_labels - pred)/test_labels), axis = 0).transpose()
                
        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(10,13))
        
        #fig, axes = plt.subplots(nrows=2, ncols=4, sharey=False)
        axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
        for ax, c, ms, ma in zip(axes, infectionstates, mae, mape):
                
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
        
        # filenames_simple=['data_secir_simple.pickle', "data_secir_simple_60days.pickle","data_secir_simple_90days.pickle","data_secir_simple_120days.pickle",
        #                   "data_secir_simple_150days.pickle"]
        filenames_simple=['data_secir_simple_30days_widerinput_2.pickle', "data_secir_simple_60days_w2.pickle","data_secir_simple_90days_w2.pickle","data_secir_simple_120days_w2.pickle",
                          "data_secir_simple_150days_w2.pickle"]

        path = os.path.dirname(os.path.realpath(__file__))
        path_models = os.path.join(os.path.dirname(os.path.realpath(
                #os.path.dirname(os.path.realpath(path)))), 'saved_models')
                os.path.dirname(os.path.realpath(path)))))
        path_models = '/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio'
        modelnames_simple = ['saved_models_secir_simple_30days_w','saved_models_secir_simple_60days_w', 
                             'saved_models_secir_simple_90days_w', 'saved_models_secir_simple_120days_w',
                             'saved_models_secir_simple_150days_w']
        
        # modelnames_simple = ['saved_models_secir_simple_bestLSTM_2024','saved_models_secir_simple_bestLSTM_2024_60days', 
        #                      'saved_models_secir_simple_bestLSTM_2024_90days', 'saved_models_secir_simple_bestLSTM_2024_120days',
        #                      'saved_models_secir_simple_bestLSTM_2024_150days']
        
        days = ['30days','60days','90days','120days','150days']


        # compartment_array = []
        # for compartment in InfectionState.values():
        #                 compartment_array.append(compartment) 
        # index=[str(compartment).split('.')[1]
        # for compartment in compartment_array]
       

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
                infectionstates = ['Susceptible','Exposed', 'InfectedNoSymptoms', 'InfectedSymptoms', 'InfectedSevere', 'InfectedCritical', 'Recovered', 'Dead']
                mape = 100*np.mean(abs((test_labels - pred)/test_labels), axis = 0).transpose()
                for ax, c, ma in zip(axes, infectionstates, mape):
                        
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
        plt.savefig("MAPE_30_days_allmodels_secirsimple_w.png")



def plot_inputs():
        path = os.path.dirname(os.path.realpath(__file__))
        path_data = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'data')
        filename = "data_secir_simple_5daysinput_90.pickle"

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
                for day in comp:
                        for a, p in zip(arrays, percentiles):
                                #array_ = []
                                #for p in zip( percentiles): 
                                 
                                a.append(np.percentile(day, p))
                                #a.append(array_)

        compartment_array = []
        for compartment in InfectionState.values():
                        compartment_array.append(compartment) 
        index=[str(compartment).split('.')[1]
        for compartment in compartment_array]

        #linestyle = [':', '--', '-', '--', ':']
        #colors = 
       

        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(8,10))
        axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
        
        for array, label, in zip(arrays, labels):
                #array = np.squeeze(array)
                n_compartments = 8 
                n_days = 5
                array = np.asarray(array).reshape([n_compartments,n_days])
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
        #fig.legend()
        plt.savefig("inputs_secir_simple_matplotlib.png")





def plot_inputs_withseaborn():
        path = os.path.dirname(os.path.realpath(__file__))
        path_data = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'data')
        filename = "data_secir_simple_30days_widerinput_3.pickle"

        file = open(os.path.join(path_data,filename), 'rb')

        data = pickle.load(file)
        data_splitted = split_data(data['inputs'], data['labels'])
        train_inputs = np.transpose(data_splitted['train_inputs']) # reshape to [n_compartments, n_days, n_samples]

        #labels = ['percentile p50', 'percentile p25', 'percentile p75', 'percentile p05', 'percentile p95']
        labels = ['percentile p50', 'percentile p25/75', 'percentile p05/95']

        inputs_reversed = np.expm1(train_inputs)        

        compartment_array = []
        for compartment in InfectionState.values():
                        compartment_array.append(compartment) 
        index=[str(compartment).split('.')[1]
        for compartment in compartment_array]


        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(8,10))
        axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]

        
        intervals = [75,95]
        colors = ['tab:red', 'tab:green', 'tab:blue', 'tab:red', 'tab:green', 'tab:blue', 'tab:red', 'tab:green']   

        
        for ax, c, values , color in zip(axes, index, inputs_reversed, colors): 
                     
       
                #for values, label in zip(array, labels): 
                df_inputs = pd.DataFrame(data =values.transpose())
                df_inputs.columns = ['1', '2', '3', '4', '5']
                
                df2 = pd.melt(df_inputs,  
                  var_name="Day")                      
                
                ax.set_ylabel('Number of individuals')
                ax.set_xlabel('Days')
                for interval in intervals:
                        sns.lineplot(ax = ax, data =df2, x = 'Day', y = 'value',  estimator="median", errorbar=("pi", interval), 
                                      color=color, legend = 'auto')


                ax.set_title(c, fontsize = 10)
                #handles, labels = ax.get_legend_handles_labels()
                #ax.legend(labels=labels, loc='upper right')
                                                
                fig.tight_layout()
                        
                ax7.set_xlabel('Days')
                ax8.set_xlabel('Days')

        plt.savefig("inputs_secir_simple_w3.png")





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



def plot_input_days():
        df = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/secir_simple_grid_search_input_width/dataframe_simple_30days')
        df= df[['input_days', 'kfold_train', 'kfold_test']]

        df_groups = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/secir_simple_grid_search_input_width/dataframe_groups30_test')
        df_groups= df_groups[['input_days', 'kfold_train', 'kfold_test']]


        input_days = [1,2,3,4,5]
        penguin_means = {
        'secir simple': df['kfold_test'],
        'secir groups': df_groups['kfold_test'][:5]}

        x = np.arange(len(input_days))  # the label locations
        width = 0.25  # the width of the bars
        multiplier = 0

        fig, ax = plt.subplots(layout='constrained')

        for attribute, measurement in penguin_means.items():
                offset = width * multiplier
                rects = ax.bar(x + offset, measurement.round(4), width, label=attribute)
                ax.bar_label(rects, padding=3)
                multiplier += 1

                # Add some text for labels, title and custom x-axis tick labels, etc.
                ax.set_ylabel('Test MAPE')
                ax.set_xlabel('Number of input days')
                #ax.set_title('')
                ax.set_xticks(x + width, input_days)
                ax.legend(loc='upper left', ncols=3)
                #ax.set_ylim(0, 250)

        plt.savefig('input_days_simple_groups.png')
  

def boxplot_inputs():
        # path = os.path.dirname(os.path.realpath(__file__))
        # path_data = os.path.join(os.path.dirname(os.path.realpath(
        #         os.path.dirname(os.path.realpath(path)))), 'data')
        # filename = "data_secir_simple_30days_widerinput.pickle"

        file = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/data/data_secir_simple.pickle', 'rb')
        file_w = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/data/data_secir_simple_30days_widerinput.pickle', 'rb')
        file_w2 = open('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/data/data_secir_simple_30days_widerinput_3.pickle', 'rb') 
       
        data = pickle.load(file)
        data = np.expm1(data['inputs'])
        data_w = pickle.load(file_w)
        data_w = np.expm1(data_w['inputs'])
        data_w2 = pickle.load(file_w2)
        data_w2 = np.expm1(data_w2['inputs'])

        compartment_array = []
        for compartment in InfectionState.values():
                compartment_array.append(compartment) 
        index = [str(compartment).split('.')[1] for compartment in compartment_array]

        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2, sharey=False, figsize=(8, 10))
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

        for ax, compartment, d, dw1, dw2 in zip(axes, index, np.asarray(data).transpose(), np.asarray(data_w).transpose(), np.asarray(data_w2).transpose()):
      
                d_df = pd.DataFrame(data = d)
                d_df = pd.melt(d_df.transpose(),var_name="Day") 
                d_df['type'] = 'b'  

                dw1_df = pd.DataFrame(data = dw1)
                dw1_df = pd.melt(dw1_df.transpose(), var_name="Day") 
                dw1_df['type'] = 'w1'  
        
                dw2_df = pd.DataFrame(data = dw2)
                dw2_df = pd.melt(dw2_df.transpose(), var_name="Day") 
                dw2_df['type'] = 'w2' 
                df_all = pd.concat([d_df, dw1_df, dw2_df], ignore_index=True)

                sns.boxplot(ax = ax, x='Day', y='value', data = df_all, hue = 'type', palette = 'Set1', width = 0.8, legend = 'auto')
                ax.set_title(compartment, fontsize = 10)
                ax.legend().set_visible(False)


                handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper right', ncol=3, bbox_to_anchor=(0.7, 1), frameon=False)
        plt.tight_layout()
                 
        plt.savefig("boxplot_input_compartments_b_w1_w3.png")


def plot_results_on_b_w1_w3():
        df_b = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/secir_simple_dataframessecir_simple_baselineLSTM_30')
        df_w1 = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/secir_simple_dataframessecir_simple_w1LSTM_30_w1')
        #df_w2 = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/secir_simple_dataframessecir_simple_w2LSTM_30_w2')
        df_w3 = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/secir_simple_dataframessecir_simple_w3LSTM_30_w3')

        df_all = pd.concat([df_b, df_w1, df_w3])


        #df = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/secir_simple_grid_search_input_width/dataframe')
        #df= df[['input_days', 'kfold_train', 'kfold_test']]

        #df_groups = pd.read_csv('/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/secir_simple_grid_search_input_width/dataframe_90days')
        #df_groups= df_groups[['input_days', 'kfold_train', 'kfold_test']]


        input_days = ['b', 'w1', 'w2']
        penguin_means = {
        'train': df_all['kfold_train'],
        'test': df_all['kfold_test']}

        x = np.arange(len(input_days))  # the label locations
        width = 0.25  # the width of the bars
        multiplier = 0

        fig, ax = plt.subplots(layout='constrained')

        for attribute, measurement in penguin_means.items():
                offset = width * multiplier
                rects = ax.bar(x + offset, measurement.round(4), width, label=attribute)
                ax.bar_label(rects, padding=3)
                multiplier += 1

                # Add some text for labels, title and custom x-axis tick labels, etc.
                ax.set_ylabel('Test MAPE')
                #ax.set_xlabel('')
                #ax.set_title('')
                ax.set_xticks(x + width, input_days)
                ax.legend(loc='upper left', ncols=3)
                #ax.set_ylim(0, 250)

        plt.savefig('train_test_secir_simple_v_w1_w2.png')




def compartment_error_simple_barplot_comarison_MAEmodel():
        # this plot compares our classic model 30 day secir simple model with a model that used MAE as a loss metric

                #load data 
        #path = os.path.dirname(os.path.realpath(__file__))
        #path_data = os.path.join(os.path.dirname(os.path.realpath(
        #        os.path.dirname(os.path.realpath(path)))), 'data')
        
        #filename = "data_secir_simple_150days.pickle"

        #file = open(os.path.join(path_data,filename), 'rb')
        filename = '/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/data/data_secir_simple.pickle'
        modelnames = ['/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_simple_bestLSTM_2024',
                      '/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models/saved_models_secir_simple_bestLSTM_2024_30days_mae']

        all_mape =[]
        for modelname in modelnames:
        
                file = open(filename, 'rb')
        

                data = pickle.load(file)
                data_splitted = split_data(data['inputs'], data['labels'])

                test_inputs = data_splitted['test_inputs']
                test_labels = data_splitted['test_labels']

                new_model = tf.keras.models.load_model(modelname)

                pred = new_model.predict(test_inputs)

                #mae =  np.mean((abs(labels_reversed - pred_reversed)), axis = 0).transpose()
                mape = 100*np.mean(abs((test_labels - pred)/test_labels), axis = 0).transpose()
                all_mape.append(np.mean(mape, axis = 1))
                
        compartment_array = []
        for compartment in InfectionState.values():
                compartment_array.append(compartment) 
        index = [3,5]
        compartments_cleaned= np.delete(compartment_array, index)

        compartmentnames=[str(compartment).split('.')[1]
               for compartment in compartments_cleaned]

        compartment_errors = {
        'MAPEmodel': all_mape[0],
        'MAEmodel': all_mape[1]}


        x = np.arange(-1,len(compartmentnames)-1)  # the label locations
        width = 0.25  # the width of the bars
        multiplier = 0

        fig, ax = plt.subplots(layout='constrained')

        for attribute, measurement in compartment_errors.items():

                offset = width * multiplier
                rects = ax.barh(x + offset, measurement.round(4), width, label=attribute)

                large_bars = [p if p > 0.28 else '' for p in measurement.round(4)]
                small_bars = [p if p <= 0.28 else '' for p in measurement.round(4)]

                ax.bar_label(rects, small_bars,
                  padding=3, color='black', fontsize = 8 )
                ax.bar_label(rects, large_bars,
                                        padding=-40, color='white', fontsize = 8)


                #ax.bar_label(rects, padding=3, fontsize = 8)
                multiplier += 1

                # Add some text for labels, title and custom x-axis tick labels, etc.
                ax.set_xlabel('Test MAPE')
                #ax.set_ylabel('')
                #ax.set_title('')
                ax.set_yticks(x + width, compartmentnames)
                ax.legend(loc='lower right', ncols=3)
                #ax.set_ylim(0, 250)

        plt.savefig("secir_simple_compartments_distribution_MAPE_MAEmodels_.png")


def compartment_error_simple_barplot():
                #load data 
        #path = os.path.dirname(os.path.realpath(__file__))
        #path_data = os.path.join(os.path.dirname(os.path.realpath(
        #        os.path.dirname(os.path.realpath(path)))), 'data')
        
        #filename = "data_secir_simple_150days.pickle"

        #file = open(os.path.join(path_data,filename), 'rb')
        filenames = ['/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/data/data_secir_simple_30days_widerinput_2.pickle', 
                     '/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/data/data_secir_simple_90days_w2.pickle',
                     '/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/data/data_secir_simple_150days_w2.pickle']
        modelnames = ['/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_30days_w',
                      '/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_90days_w',
                      '/home/schm_a45/Documents/Code/memilio/memilio/pycode/memilio-surrogatemodel/memilio/saved_models_secir_simple_150days_w']

        all_mape =[]
        for filename, modelname in zip(filenames, modelnames):
        
                file = open(filename, 'rb')
        

                data = pickle.load(file)
                data_splitted = split_data(data['inputs'], data['labels'])

                test_inputs = data_splitted['test_inputs']
                test_labels = data_splitted['test_labels']

                new_model = tf.keras.models.load_model(modelname)

                pred = new_model.predict(test_inputs)
                pred_reversed = np.expm1(pred)
                labels_reversed = np.expm1(test_labels)

                #mae =  np.mean((abs(labels_reversed - pred_reversed)), axis = 0).transpose()
                mape = 100*np.mean(abs((test_labels - pred)/test_labels), axis = 0).transpose()
                all_mape.append(np.mean(mape, axis = 1))
                #all_mape.append(np.mean(mae, axis = 1))
                
        compartment_array = []
        for compartment in InfectionState.values():
                compartment_array.append(compartment) 
        index = [3,5]
        compartments_cleaned= np.delete(compartment_array, index)

        compartmentnames=[str(compartment).split('.')[1]
               for compartment in compartments_cleaned]

        compartment_errors = {
        '30d': all_mape[0],
        '90d': all_mape[1], 
        '150d': all_mape[2]}


        x = np.arange(-1,len(compartmentnames)-1)  # the label locations
        width = 0.25  # the width of the bars
        multiplier = 0

        fig, ax = plt.subplots(layout='constrained')

        for attribute, measurement in compartment_errors.items():

                offset = width * multiplier
                rects = ax.barh(x + offset, measurement.round(4), width, label=attribute)

                large_bars = [p if p > 0.28 else '' for p in measurement.round(4)] # for MAPE
                small_bars = [p if p <= 0.28 else '' for p in measurement.round(4)] # for MAPE
                
                # large_bars = [p if p > 80 else '' for p in measurement.round(2)] # for MAE
                # small_bars = [p if p <= 80 else '' for p in measurement.round(2)] # for MAE



                ax.bar_label(rects, small_bars,
                  padding=3, color='black', fontsize = 8 )
                ax.bar_label(rects, large_bars,
                                        padding=-40, color='white', fontsize = 8)


                #ax.bar_label(rects, padding=3, fontsize = 8)
                multiplier += 1

                # Add some text for labels, title and custom x-axis tick labels, etc.
                ax.set_xlabel('Test MAPE')
                #ax.set_xlabel('Test MAE')
                #ax.set_ylabel('')
                #ax.set_title('')
                ax.set_yticks(x + width, compartmentnames)
                ax.legend(loc='lower right', ncols=3)
                #ax.set_ylim(0, 250)

        plt.savefig("secir_simple_compartments_distribution_MAPE_w.png")


def lineplot_days_simple_and_groups_MAE():
        path = os.path.dirname(os.path.realpath(__file__))
        path_data = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'data')
        
        # filenames_simple=['data_secir_simple.pickle', "data_secir_simple_60days.pickle","data_secir_simple_90days.pickle","data_secir_simple_120days.pickle",
        #                   "data_secir_simple_150days.pickle"]
        
        filenames_simple=['data_secir_simple_30days_widerinput_2.pickle', "data_secir_simple_60days_w2.pickle","data_secir_simple_90days_w2.pickle","data_secir_simple_120days_w2.pickle",
                          "data_secir_simple_150days_w2.pickle"]
        
        filenames_groups = ['data_secir_groups_30days_nodamp_w.pickle',"data_secir_groups_60days_nodamp_w.pickle","data_secir_groups_90days_nodamp_w.pickle",
                    "data_secir_groups_120days_nodamp_w.pickle","data_secir_groups_150days_nodamp_w.pickle"]

        # filenames_groups = ['data_secir_groups_30days_nodamp.pickle',"data_secir_groups_60days_nodamp.pickle","data_secir_groups_90days_nodamp.pickle",
        #             "data_secir_groups_120days_nodamp.pickle","data_secir_groups_150days_nodamp.pickle"]
        
        path = os.path.dirname(os.path.realpath(__file__))
        path_models = os.path.join(os.path.dirname(os.path.realpath(
                os.path.dirname(os.path.realpath(path)))), 'saved_models')
        
        # modelnames_simple = ['saved_models_secir_simple_bestLSTM_2024','saved_models_secir_simple_bestLSTM_2024_60days', 
        #                      'saved_models_secir_simple_bestLSTM_2024_90days', 'saved_models_secir_simple_bestLSTM_2024_120days',
        #                      'saved_models_secir_simple_bestLSTM_2024_150days']
        
        modelnames_simple = ['saved_models_secir_simple_30days_w','saved_models_secir_simple_60days_w', 
                             'saved_models_secir_simple_90days_w', 'saved_models_secir_simple_120days_w',
                             'saved_models_secir_simple_150days_w']
        
        modelnames_groups = ['saved_models_groups_nodamp_30days_w', 'saved_models_groups_nodamp_60days_w', 'saved_models_groups_nodamp_90days_w',
                             'saved_models_groups_nodamp_120days_w', 'saved_models_groups_nodamp_150days_w']
        
        
        # modelnames_simple = ['saved_models_secir_simple_30days_w','saved_models_secir_simple_60days_w', 
        #                      'saved_models_secir_simple_90days_w', 'saved_models_secir_simple_120days_w',
        #                      'saved_models_secir_simple_150days_w']
        
        days = [30,60,90,120,150]

        MAE_simple = []
        MAE_groups = []
        for filenames, modelnames, MAE in zip([filenames_simple, filenames_groups], [modelnames_simple, modelnames_groups],[MAE_simple, MAE_groups]):
                
                for modelname, filename in zip(modelnames, filenames): 

                        
                        file = open(os.path.join(path_data,filename), 'rb')
                        model = tf.keras.models.load_model(os.path.join(path_models, modelname))
                        

                        data = pickle.load(file)
                        data_splitted = split_data(data['inputs'], data['labels'])

                        test_inputs = data_splitted['test_inputs']
                        test_labels = data_splitted['test_labels']
                        pred = model.predict(test_inputs)
                        pred_reversed = np.expm1(pred)
                        labels_reversed = np.expm1(test_labels)

                        mae =  np.mean((abs(labels_reversed - pred_reversed)), axis = 0).transpose()

                        mean_mae = mae.mean()
                        MAE.append(mean_mae)

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
        ax1.plot(days, MAE_simple, label ='simple', color = color, marker = 'o' )
        ax1.set_xticks(days)
        ax1.set_xlabel('Number of days')
        ax1.set_ylabel('MAE', color = color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.set_title('MAE for number of days to be predicted')
        


        ax2 = ax1.twinx()
        color = 'tab:blue'
        ax2.set_ylabel('MAE', color = color)
        ax2.plot(days, MAE_groups, color = color, label = 'groups', marker = 'x')
        ax2.tick_params(axis='y', labelcolor=color)
        fig.legend(loc= 'lower left', fontsize = 9)
        fig.tight_layout()
        

        plt.savefig("plot_days_simple_and_groups_MAE_w.png")
