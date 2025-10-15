import matplotlib
import matplotlib.pyplot as plt
import pandas as pd


def set_fontsize(base_fontsize=17):
    fontsize = base_fontsize
    plt.rcParams.update({
        'font.size': fontsize,
        'axes.titlesize': fontsize * 1,
        'axes.labelsize': fontsize,
        'xtick.labelsize': fontsize * 0.8,
        'ytick.labelsize': fontsize * 0.8,
        'legend.fontsize': fontsize * 0.8,
        'font.family': "Arial"
    })


plt.style.use('default')

dpi = 300

colors = {"Light blue": "#48ADFF", "Blue": "#155489", "Red": "#BA421E", "Purple": "#8200AE", "Light green": "#8EBE4B", "Green": "#5D8A2B",
          "Teal": "#20A398", "Yellow": "#FFD45E", "Orange": "#E5945A", "Rose": "#CF688F", "Grey": "#C0BFBF", "Dark grey": "#616060", "Light grey": "#E2E2E2"}

if __name__ == '__main__':
    set_fontsize()
