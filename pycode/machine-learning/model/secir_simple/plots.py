import pandas  as pd
import matplotlib.pyplot as plt 
from matplotlib.lines import Line2D


data = [[0.1645, 0.2988, 0.3428], [0.1186, 0.1907, 0.1686], [0.2042, 0.1498, 0.2169]]
df = pd.DataFrame(data = data, columns = ['2','3','4'], index=["MLP", "LSTM", "CNN"])


def plot_layer(df):
    plt.figure().clf() 
       
    linestyles = ['--', '-',  ':']

    markers = []
    for m in Line2D.markers:
        try:
            if len(m) == 1 and m != " ":
                markers.append(m)
        except TypeError:
            pass

    for values,label, ls, m in zip(df.values, df.index, linestyles, markers): 
        plt.plot(values, label = label,  linestyle=ls, marker = m )


    plt.title('Network architecture: influence of number of hidden layers')
    plt.legend()
    plt.xlabel('Number of layers')
    #plt.xticks([2,3,4])
    plt.ylabel('Test MAPE error')
    plt.savefig("layers_secir_simple.png")

plot_layer(df)