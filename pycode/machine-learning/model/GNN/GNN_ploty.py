import os 
import tensorflow as tf 
import pickle 
import matplotlib.pyplot as plt
import pandas as pd

def lineplot_number_of_days():

        filenames =[ "GNNtype1_ARMA_60days.csv", "GNNtype1_ARMA_90days.csv", "GNNtype1_ARMA_120days.csv"]
        days = [30,60,90,120]
        MAPE = []
        MAPE.append(2.0147)  #  for 30 days experiment
        for file in filenames: 
                
                path = os.path.dirname(os.path.realpath(__file__))
                path_data = os.path.join(os.path.dirname(os.path.realpath(
                        os.path.dirname(os.path.realpath(path)))), 'dataframe_gridsearch_2024')
                
                filename = os.path.join(path_data, file)

                df = pd.read_csv(filename)

                MAPE.append(df['kfold_test'][0])
     
                
        plt.figure().clf()
        fig, ax = plt.subplots(figsize =(8,5))
        ax.plot(days, MAPE,  marker = 'o' )
        ax.set_xticks(days)
        ax.set_xlabel('Number of days')
        ax.set_ylabel('MAPE')
        ax.set_title('MAPE for number of days to be predicted')
        plt.savefig("plot_days_GNN.png")

