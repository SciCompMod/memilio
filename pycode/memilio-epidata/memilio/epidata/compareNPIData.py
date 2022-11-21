import os
import csv
import pandas as pd
import numpy as np

# directory = '/home/wend_aa/Documents/PSS/NPIs'
directory = 'c:\\work\\projets\\epidemiology\\code\\memilio\\data/pydata\\Germany/'

#numberofcities = 2

df_npis_old = pd.read_csv(
    os.path.join(directory, 'kr_massnahmen_unterkategorien.csv'),
    sep=',')  # , nrows=numberofcities*1248
print(df_npis_old)
numberofcities = int(len(df_npis_old) / 1248)
print('Number of cities', numberofcities)

codelist = ['m01a', 'm01b', 'm02a', 'm02b', 'm03', 'm04', 'm05', 'm06', 'm07', 'm08', 'm09',
        'm10', 'm11', 'm12', 'm13', 'm14', 'm15', 'm16', 'm17', 'm18', 'm19', 'm20', 'm21']

#list = ['m01a', 'm01b']

# df_npis_new = pd.read_csv(
#  os.path.join(directory, 'new', 'kr_massn_unterkat_m01a.csv'),
# sep=',', nrows=4000)

# print(df_npis_old)

# print(pd.testing.assert_frame_equal(df_npis_old, df_npis_new))
# print(pd.testing.assert_frame_equal(df_npis_old, df_npis_new))


# number of days that we have data for in new data set
numberofdays = 883

# create data frame that contains structure of old data
df = df_npis_old.iloc[:, :6]

start_county = 1

df_local = [pd.DataFrame() for i in range(401)]
counter_col = 0
for code in codelist:
    # print(code)
    df_npis_new = pd.read_csv(
        os.path.join(directory, 'new',
                        'kr_massn_unterkat_{}.csv'.format(code)),
        sep=',')  # , skiprows=1248
    # print(df_npis_new.iloc[0:numberofdays])
    #df_npis_new = df_npis_new.iloc[1152:]
    #print('df_new', df_npis_new)
    counties = np.sort(df_npis_new.ags5.unique())
    if len(df_npis_new) / len(counties) != numberofdays:
        print('error')
    if len(counties) != 401:
        print('error')

    # extract dates from new data
    dates = df_npis_new.iloc[:numberofdays, 5]
    # rename dates so that they match dates from old npi dataframe
    dates_new = ['d' + date.replace('-', '') for date in dates]

    for i in range(0,401):
        if counter_col == 0:
            df_local[i] = pd.DataFrame(columns=list(df_npis_new.columns[0:5]) + ['code'] + dates_new)
        else:
            print('todo')

        df_npis_new[df_npis_new.ags5 == counties[i]].iloc[:, 6:].T.reset_index()

    if counter_col == 0:
        df = df_npis_new.copy()
    else:
        # extend data frame with columns for each date that is included in new data set
        df = pd.concat([df, pd.DataFrame(columns=list(df_npis_new.columns[0:5])+ dates_new)])

    for i in range(start_county, numberofcities):
        print('County number', i)
        counter_col = 0



    # get number of codes for current subcategory, we substract 6 because we have 6 columns for bundesland, kreis, etc.
    numberofcodes = len(df_npis_new.columns) - 6
    #print('numberofcodes', numberofcodes)

    # insert values from new dataset into data frame
    df_local[0].iloc[:, 6:] = df_npis_new.iloc[0:numberofdays, 6:].T

    #print(df.iloc[i*1248 + counter: i*1248 + counter + numberofcodes])

    counter_col += numberofcodes
    counter_counties += 1
print(df)

df_dropped = df
df_npis_old_dropped = df_npis_old
df_lost_subcat = df_npis_old


for i in range(numberofcities):
    # drop entries for subcategories 22, 23, and 24 because they are not identified for new data set
    #print((i+1)*counter, (i+1)*1248)
    # df_dropped = df_dropped.drop(
    #   df_dropped.index[range((i+1)*counter, i*counter + 1248)])
    df_dropped = df_dropped.drop(
        df_dropped.index[range((i+1)*counter_col, i*counter_col + 1248)])
    #print('df_dropped in loop', df_dropped.iloc[:(i+1)*counter + 3])
    # df_npis_old_dropped = df_npis_old_dropped.drop(df_npis_old_dropped.index[range(
    #   (i+1)*counter, i*counter + 1248)])
    df_npis_old_dropped = df_npis_old_dropped.drop(df_npis_old_dropped.index[range(
        (i+1)*counter_col, i*counter_col + 1248)])

    # extract all dropped rows from df_npis_old
   # print('indices', (i)*1248, (i+1)*counter)
    df_lost_subcat = df_lost_subcat.drop(
        df_lost_subcat.index[range(i*(1248-counter_col), i*(1248-counter_col) + counter_col)])
    # print('df_dropped in loop lost subcat',
    #     df_lost_subcat.iloc[:(i)*(1248-counter) + 3])

#print('df', df.iloc[:1248-counter])
df_dropped = df_dropped.iloc[:, :723]
print('df_dropped', df_dropped)
print('df_npis_old_dropped', df_npis_old_dropped)


#df = df.iloc[:counter, :723]

#df_npis_old = df_npis_old[:counter]
#print(df_npis_old_dropped.iloc[2*counter:], df_dropped.iloc[2*counter:])
# print('Differences between old and new data:', pd.testing.assert_frame_equal(
#   df_npis_old_dropped.iloc[2*counter:, :], df_dropped.iloc[2*counter:, :], check_dtype=False))
print('df_npis_old_dropped columns', df_npis_old_dropped.columns)
print('df_dropped columns', df_dropped.columns)

# Check differences columnwise
nodiffcolumns = []
for i in range(df_npis_old_dropped.shape[1]):
    try:
        if pd.testing.assert_frame_equal(pd.DataFrame(df_npis_old_dropped.iloc[:, i]), pd.DataFrame(df_dropped.iloc[:, i]), check_dtype=False) == None:
            nodiffcolumns.append(df_npis_old_dropped.columns[i])
            #print('No difference in column', df_npis_old_dropped.columns[i])

    except AssertionError as e:
        print(e, "\n")


# Check differences rowwise
nodiffrows = []
for i in range(df_npis_old_dropped.shape[0]):
    try:
        if pd.testing.assert_frame_equal(pd.DataFrame(df_npis_old_dropped.iloc[i]), pd.DataFrame(df_dropped.iloc[i]), check_dtype=False) == None:
            nodiffrows.append(
                df_npis_old_dropped.iloc[i, :6])
            #print('No difference in column', df_npis_old_dropped.columns[i])

    except AssertionError as e:
        print(e, "\n")

# Check differences citywise
nodiffcities = []
for i in range(numberofcities):
    try:
        if pd.testing.assert_frame_equal(pd.DataFrame(df_npis_old_dropped.iloc[i*counter_col: (i+1)*counter_col, :]), pd.DataFrame(df_dropped.iloc[i*counter_col: (i+1)*counter_col, :]), check_dtype=False) == None:
            nodiffcities.append(df_npis_old_dropped.iloc[i*counter_col]['kreis'])
            #print('No difference in column', df_npis_old_dropped.columns[i])

    except AssertionError as e:
        print(e, "\n")


# check if all values in dropped subcategories are = -99
print('df_lostsubcat', df_lost_subcat.iloc[:, :])
print(np.sum(np.where(df_lost_subcat.iloc[:, 6:] != -99)))
columnwisecheck = (df_lost_subcat.iloc[:, 6:] == -99).all()
checkallcolumns = (columnwisecheck == True).all()
print('Dropped subcategories have never been active:', checkallcolumns)


# print(nodiffcolumns)
print('Number of no diff columns', len(nodiffcolumns))
print('Total number of columns', df_npis_old_dropped.shape[1])


#print('df_npis_old_dropped columns', df_npis_old_dropped.columns)
#print('df_dropped columns', df_dropped.columns)
print('Column names are equal',
      (df_npis_old_dropped.columns == df_dropped.columns).all())


# print(nodiffrows)
print('Number of no diff rows', len(nodiffrows))
print('Total number of rows', df_npis_old_dropped.shape[0])


print(nodiffcities)
print('Number of no diff cities', len(nodiffcities))
print('Total number of cities', numberofcities)


# save results in csv file
pd.DataFrame(nodiffcolumns).to_csv(os.path.join(
    directory, 'comparedata_columns.csv'))  # , nodiffrows, nodiffcities
pd.DataFrame(nodiffrows).to_csv(
    os.path.join(directory, 'comparedata_rows.csv'))
pd.DataFrame(nodiffcities).to_csv(
    os.path.join(directory, 'comparedata_cities.csv'))
