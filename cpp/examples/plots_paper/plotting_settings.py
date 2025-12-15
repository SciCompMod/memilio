import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import os


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

colors = {"Blue": "#155489",
          "Medium blue": "#64A7DD",
          "Light blue": "#B4DCF6",
          "Lilac blue": "#AECCFF",
          "Turquoise": "#76DCEC",
          "Light green": "#B6E6B1",
          "Medium green": "#54B48C",
          "Green": "#5D8A2B",
          "Teal": "#20A398",
          "Yellow": "#FBD263",
          "Orange": "#E89A63",
          "Rose": "#CF7768",
          "Red": "#A34427",
          "Purple": "#741194",
          "Grey": "#C0BFBF",
          "Dark grey": "#616060",
          "Light grey": "#F1F1F1"}

plotting_dir = os.path.dirname(__file__)

if __name__ == '__main__':
    set_fontsize()