import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as pclrs
import pandas as pd
import datetime as dt
import numpy as np

def fix_eisenach(map_data : gpd.GeoDataFrame):
    #merge geometries for eisenach with wartburgkreis
    wartburg = map_data.ARS == '16063'
    eisenach = map_data.ARS == '16056'
    map_data.at[wartburg, 'geometry'] = [map_data[wartburg].geometry.values[0].union(map_data[eisenach].geometry.values[0])]
    # remove eisenach 
    return map_data.drop(map_data[eisenach].index.values[0])

def plot(
    date : dt.date,
    filter_state : int = None, # 1-16: id of state, None: all states
    filter_agegroup : int = None, # 0: 0-14, 1: 15-59, 2: 60+, None: all groups
    filter_vacc_state : int = 1, # 0: partial, 1: full, 2: booster 
):
    date_str = date.strftime('%Y-%m-%d')

    #load and filter RKI data
    rki_data = pd.read_json('data/all_county_ageinf_vacc_all_dates_old.json')
    rki_filter = rki_data.Date == date_str
    if not filter_agegroup is None:
        if filter_agegroup == 0:
            age_strs = ['0-4', '5-14']
        elif filter_agegroup == 1:
            age_strs = ['15-34', '35-59']
        else:
            age_strs = ['60-79', '80-99']
        rki_filter = rki_filter & (rki_data.Age_RKI.isin(age_strs))
    rki_data = rki_data[rki_filter]
    if filter_vacc_state == 0:
        col_name = 'Vacc_partially'
    elif filter_vacc_state == 1:
        col_name = 'Vacc_completed'
    else:
        col_name = 'Vacc_refreshed'
    rki_data['Count'] = rki_data[col_name]

    #load and filter sanitized data
    new_data = pd.read_json('data/all_county_ageinf_vacc_all_dates_sanitized.json')
    new_filter = new_data.Date == date_str
    if not filter_agegroup is None:
        if filter_agegroup == 0:
            age_strs = ['0-4', '5-14']
        elif filter_agegroup == 1:
            age_strs = ['15-34', '35-59']
        else:
            age_strs = ['60-79', '80-99']
        new_filter = new_filter & (new_data.Age_RKI.isin(age_strs))
    new_data = new_data[new_filter]
    if filter_vacc_state == 0:
        col_name = 'Vacc_partially'
    elif filter_vacc_state == 1:
        col_name = 'Vacc_completed'
    else:
        col_name = 'Vacc_refreshed'
    new_data['Count'] = new_data[col_name]

    #read and filter pop data
    pop_data = pd.read_json('data/county_current_population.json')
    pop_filter = pd.Series([True for i in range(len(pop_data))])
    if not filter_agegroup is None:
        if filter_agegroup == 0:
            age_strs = {'<3 years' : 1, '3-5 years' : 1, '6-14 years' : 1}
        elif filter_agegroup == 1:
            age_strs = {'15-17 years' : 1, '18-24 years' : 1, '25-29 years' : 1, '30-39 years' : 1, '40-49 years' : 1, '50-64 years' : 2/3 }
        else:
            age_strs = {'50-64 years' : 1/3, '65-74 years' : 1, '>74 years' : 1}
        pop_data['Count'] = pop_data.apply(lambda e: sum([f * e[g] for (g, f) in age_strs.items()]), axis = 1)
    else:
        pop_data['Count'] = pop_data['Total']
    if not pop_filter is None:
        pop_data = pop_data[pop_filter]
    
    #read and filter map data
    map_data = gpd.read_file("shapes/vg2500_01-01.utm32s.shape/vg2500/vg2500_krs.shp")
    map_data = fix_eisenach(map_data)
    map_filter = pd.Series([True for i in range(len(map_data))], index = map_data.index)
    if not filter_state is None:
        map_filter = map_filter & (map_data.ARS.str.startswith('{:02}'.format(filter_state)))
    map_data = map_data[map_filter]

    # print('rki data:')
    # print(rki_data)
    # print('new data:')
    # print(new_data)
    # print('pop data:')
    # print(pop_data)
    # print('map_data:')
    # print(map_data)

    #add vaccination rates to map
    def get_county(df, id):
        df = df[df.ID_County == id]
        if not df.empty:
            return df.Count.sum()
        else:
            print('Value for County {} not found'.format(id))
            return None
    def get_vacc_rates(e):
        pop = get_county(pop_data, int(e.ARS))
        abs_new = get_county(new_data, int(e.ARS))
        abs_rki = get_county(rki_data, int(e.ARS))
        rel_new = abs_new / pop if (abs_new is not None and pop is not None) else None
        rel_rki = abs_rki / pop if (abs_rki is not None and pop is not None) else None
        return (rel_rki, rel_new)
    map_data[['VaccRateRKI', 'VaccRateNew']] = map_data.apply(get_vacc_rates, axis=1, result_type='expand')

    #calc min and max and map to non-uniform scale
    map_data_vacc_columns = map_data[['VaccRateRKI', 'VaccRateNew']]
    min_value = min(map_data_vacc_columns.min().values)
    max_value = max(map_data_vacc_columns.max().values)
    s = 0.7
    
    def normalize_lin(x, src, dest):
        return (x - src[0]) / (src[1] - src[0]) * (dest[1] - dest[0]) + dest[0]
    def normalize_sine(x):
        if (x < s):
            i = [min_value, s]
            j = [-np.pi / 2, 0]
        else:
            i = [s, max_value]
            j = [0, np.pi / 2]
        map_x = normalize_lin(x, i, j)
        sin_map_x = np.sin(map_x)
        y = normalize_lin(sin_map_x, [-1, 1], [0, 1])
        return y
    def normalize_sine_inv(y):
        map_y = normalize_lin(y, [0, 1], [-1, 1])
        asin_map_y = np.arcsin(map_y)
        if (asin_map_y < 0):
            i = [min_value, s]
            j = [-np.pi / 2, 0]
        else:
            i = [s, max_value]
            j = [0, np.pi / 2]
        x = normalize_lin(asin_map_y, j, i)
        return x
    e = 2
    def normalize_skewed_sine(x):
        x = normalize_lin(x, (min_value, max_value), (-np.pi/2, np.pi/2))
        x = np.sin(x)
        def p(x):
            return (-1)**(e-1) * (0.5**e) * (x - 1)**e + 1
        x = p(x)
        return x
    def normalize_skewed_sine_inv(y):
        def p_inv(y):
            return -((-1)**e * (1-y) / (0.5**e))**(1/e) + 1
        y = p_inv(y)
        y = np.arcsin(y)
        y = normalize_lin(y, (-np.pi/2, np.pi/2), (min_value, max_value))
        return y

    color_norm = pclrs.FuncNorm((np.vectorize(normalize_sine), np.vectorize(normalize_sine_inv)), vmin = min_value, vmax = max_value)

    #file name extensions for filters
    part_filename_age = 'allages'
    if not filter_agegroup is None:
        part_filename_age = 'age{0}'.format(filter_agegroup)    
    part_filename_state = 'allstates'
    if not filter_state is None:
        part_filename_state = 'state{0}'.format(filter_state)
    part_filename_vacc = 'vacc{0}'.format(filter_vacc_state)
    part_filename_date = date.strftime('%Y%m%d')
    part_filename = '{}_{}_{}_{}'.format(part_filename_date, part_filename_vacc, part_filename_age, part_filename_state)

    #interactive html files
    def save_interactive(col, filename):
        map_data.explore(col, legend = True, vmin = min_value, vmax = max_value).save(filename) # TODO: norm

    save_interactive('VaccRateRKI', 'vaccrate_{}_rki.html'.format(part_filename))
    save_interactive('VaccRateNew', 'vaccrate_{}_new.html'.format(part_filename))

    from matplotlib.gridspec import GridSpec
    fig = plt.figure(constrained_layout=True)
    gs = GridSpec(3, 3, figure=fig, wspace = 0.1, width_ratios = [0.05, 1, 1], height_ratios = [0.15, 1, 0.15])
    cax = fig.add_subplot(gs[1,0])
    ax1 = fig.add_subplot(gs[:,1])
    ax2 = fig.add_subplot(gs[:,2])
    map_data.plot('VaccRateRKI', ax = ax1, cax = cax, legend = True, vmin = min_value, vmax = max_value, norm = color_norm)
    ax1.set_axis_off()
    map_data.plot('VaccRateNew', ax = ax2, legend = False, vmin = min_value, vmax = max_value, norm = color_norm)
    ax2.set_axis_off()

    plt.savefig('vaccrate_' + part_filename + "_compare.svg")

    # plt.show()

if __name__ == '__main__':
    for vacc in [1]:
        for state in [None, 9]:
            for age in [None, 1]:
                plot(dt.date(2021, 10, 1), filter_vacc_state=vacc, filter_state=state, filter_agegroup=age)