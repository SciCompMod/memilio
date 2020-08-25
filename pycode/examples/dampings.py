import numpy as np
import pandas as pd


def create_dampings(path, days = 120):
    comix_all = pd.read_csv(path + 'images/Comix.csv', index_col=0).values
    comix_home = pd.read_csv(path + 'images/Comix_Home.csv', index_col=0).values
    comix_work = pd.read_csv(path + 'images/Comix_Work.csv', index_col=0).values
    comix_other = pd.read_csv(path + 'images/Comix_Other.csv', index_col=0).values
    comix_school = pd.read_csv(path + 'images/Comix_School.csv', index_col=0).values
    polymod_all = pd.read_csv(path + 'images/Polymod.csv', index_col=0).values
    polymod_home = pd.read_csv(path + 'images/Polymod_Home.csv', index_col=0).values
    polymod_work = pd.read_csv(path + 'images/Polymod_Work.csv', index_col=0).values
    polymod_other = pd.read_csv(path + 'images/Polymod_Other.csv', index_col=0).values
    polymod_school = pd.read_csv(path + 'images/Polymod_School.csv', index_col=0).values

    home_diff = polymod_home - comix_home
    work_diff = polymod_work - comix_work
    other_diff = polymod_other - comix_other
    school_diff = polymod_school - comix_school

    other_range = 13
    work_range = 14
    home_range = 14
    school_range = 12

    other_start = 49
    work_start = 30
    home_start = 30
    school_start = 43

    other_weights = np.zeros(other_range)
    other_weights[0] = 23
    other_weights[5] = 36
    other_weights[12] = 34

    school_weights = np.zeros(school_range)
    school_weights[6] = 2
    school_weights[-1] = 1

    other_weights_reverse = np.ones(50)

    damp_dates = np.repeat(polymod_all[:, :, np.newaxis], days, axis=2)
    for i in range(days - work_start):
        if i < work_range:
            damp_dates[:,:, work_start + i] -= work_diff * (i + 1)/work_range
        else:
            damp_dates[:, :, work_start + i] -= work_diff
    for i in range(days - home_start):
        if i < home_range:
            damp_dates[:, :, home_start + i] -= home_diff * (i + 1) / home_range
        else:
            damp_dates[:, :, home_start + i] -= home_diff
    for i in range(days - other_start):
        if i < other_range:
            damp_dates[:,:,other_start + i] -= other_diff*sum(other_weights[:i+1])/(sum(other_weights))
        else:
            damp_dates[:, :,other_start + i] -= other_diff
    for i in range(days - school_start):
        if i < school_range:
            damp_dates[:,:,school_start +i] -= school_diff*sum(school_weights[:i+1])/sum(school_weights)
        else:
            damp_dates[:, :, school_start + i] -= school_diff
    '''for i in range(70):
        if i < 50:
            damp_dates[:,:,50+i] += other_diff*sum(other_weights_reverse[:i+1])/(sum(other_weights_reverse)*2)
        else:
            damp_dates[:, :, 50+i] += other_diff/2'''




    ratios = np.zeros(damp_dates.shape)
    for i in range(8):
        for j in range(8):
            for d in range(damp_dates.shape[2]):
                ratios[i,j,d] = damp_dates[i,j,d]/polymod_all[i,j]

    dfs = []
    for i in range(ratios.shape[2]):
        dfs.append(pd.DataFrame(data=ratios[:,:,i], columns=['0-4', '5-17', '18-29', '30-39', '40-49', '50-59', '60-69', '70+']))
    df = pd.concat(dfs, axis=1, keys=np.arange(ratios.shape[2]))
    df.to_csv(path + 'dampings.csv')


