import datetime
import numpy as np
import pandas as pd
import unittest


from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd
from epidemiology.epidata import getRKIData


class TestGetRKIData(unittest.TestCase):

    def test_moving_average(self):
        Date = ['2020-01-02', '2020-01-03','2020-01-20', '2020-01-25', '2020-01-30']
        Date = np.array([np.datetime64(x) for x in Date])
        ID_State = np.ones(len(Date))
        State = np.repeat('Schleswig-Holstein', len(Date))
        Confirmed = np.arange(1,len(Date)+1)

        df = pd.DataFrame(np.array([ID_State, State, Date, Confirmed]).T, columns = ['ID_State', 'State', 'Date', 'Confirmed'])

        Date = ['2020-01-03', '2020-01-04','2020-01-21', '2020-01-26', '2020-01-31']
        Date = np.array([np.datetime64(x) for x in Date])
        ID_State = np.ones(len(Date))*3
        State = np.repeat('Niedersachsen', len(Date))
        Confirmed = np.arange(len(Date)) + 5

        df = df.append(pd.DataFrame(np.array([ID_State, State, Date, Confirmed]).T, columns = ['ID_State', 'State', 'Date', 'Confirmed']))

        Date = ['2020-01-' + str(x).zfill(2) for x in range(2,32)]
        Date = np.array([np.datetime64(x) for x in Date])
        ID_State = np.ones(len(Date))
        State = np.repeat('Schleswig-Holstein', len(Date))
        Confirmed = np.ones(len(Date))*2
        for i in range(7):
            Confirmed[0 + i] = (2*i +1)/(1+i)
        for i in range(5):
            Confirmed[18 + i] = (15 + i)/7
        Confirmed[23] = (21/7)
        for i in range(1,5):
            Confirmed[23 + i] = (22 + i)/7

        Confirmed[28] = (28/7)
        Confirmed[29] = (30/7)

        result_df = pd.DataFrame(np.array([ID_State, State, Date, Confirmed]).T, columns = ['ID_State', 'State', 'Date', 'Confirmed'])

        Date = ['2020-01-' + str(x).zfill(2) for x in range(2,32)]
        Date = np.array([np.datetime64(x) for x in Date])
        ID_State = np.ones(len(Date))*3
        State = np.repeat('Niedersachsen', len(Date))
        Confirmed = np.ones(len(Date)) * 6
        Confirmed[0] = 0.0
        Confirmed[1] = 5.0
        for i in range(1,7):
            Confirmed[0 + i] = (5*i +i -1)/(1+i)
        Confirmed[7] = 41 / 7
        for i in range(5):
            Confirmed[19 + i] = (43 + i) / 7
        Confirmed[24] = 49/7
        for i in range(1,5):
            Confirmed[24 + i] = (50 + i) / 7

        Confirmed[29] = (56 / 7)

        result_df = result_df.append(pd.DataFrame(np.array([ID_State, State, Date, Confirmed]).T, columns = ['ID_State', 'State', 'Date', 'Confirmed']))
        result_df = result_df.astype({'ID_State': 'float64', 'Confirmed': 'float64', 'Date': 'datetime64[ns]'})
        result_df.index = (range(len(result_df)))

        test_df = getRKIData.fill_df(df, ['ID_State', 'State'], ['Confirmed'], True)

        pd.testing.assert_frame_equal(test_df, result_df)

    def test_fill_df(self):
        Date = ['2020-01-02', '2020-01-03','2020-01-20', '2020-01-25', '2020-01-30']
        Date = np.array([np.datetime64(x) for x in Date])
        ID_State = np.ones(len(Date))
        State = np.repeat('Schleswig-Holstein', len(Date))
        Confirmed = np.arange(1,len(Date)+1)

        df = pd.DataFrame(np.array([ID_State, State, Date, Confirmed]).T, columns = ['ID_State', 'State', 'Date', 'Confirmed'])

        Date = ['2020-01-03', '2020-01-04','2020-01-21', '2020-01-26', '2020-01-31']
        Date = np.array([np.datetime64(x) for x in Date])
        ID_State = np.ones(len(Date))*3
        State = np.repeat('Niedersachsen', len(Date))
        Confirmed = np.arange(len(Date)) + 5

        df = df.append(pd.DataFrame(np.array([ID_State, State, Date, Confirmed]).T, columns = ['ID_State', 'State', 'Date', 'Confirmed']))

        Date = ['2020-01-' + str(x).zfill(2) for x in range(2,32)]
        Date = np.array([np.datetime64(x) for x in Date])
        ID_State = np.ones(len(Date))
        State = np.repeat('Schleswig-Holstein', len(Date))
        Confirmed = np.ones(len(Date))*2
        Confirmed[0] = 1
        Confirmed[18:23] = 3
        Confirmed[23:28] = 4
        Confirmed[28:30] = 5

        result_df = pd.DataFrame(np.array([ID_State, State, Date, Confirmed]).T, columns = ['ID_State', 'State', 'Date', 'Confirmed'])



        Date = ['2020-01-' + str(x).zfill(2) for x in range(2,32)]
        Date = np.array([np.datetime64(x) for x in Date])
        ID_State = np.ones(len(Date))*3
        State = np.repeat('Niedersachsen', len(Date))
        Confirmed = np.ones(len(Date)) * 6
        Confirmed[0] = 0
        Confirmed[1] = 5
        Confirmed[19:24] = 7
        Confirmed[24:29] = 8
        Confirmed[29:30] = 9


        result_df = result_df.append(pd.DataFrame(np.array([ID_State, State, Date, Confirmed]).T, columns = ['ID_State', 'State', 'Date', 'Confirmed']))
        result_df = result_df.astype({'ID_State': 'float64', 'Date': 'datetime64[ns]'})
        result_df.index = (range(len(result_df)))

        test_df = getRKIData.fill_df(df, ['ID_State', 'State'], ['Confirmed'], False)

        pd.testing.assert_frame_equal(test_df, result_df)

    def test_fuse_berlin(self):
        Date = ['2020-01-02', '2020-01-03','2020-01-20', '2020-01-25', '2020-01-30']
        Date = np.array([np.datetime64(x) for x in Date])
        ID_County = np.ones(len(Date))*11001
        County = np.repeat('Berlin-Mitten', len(Date))
        ID_State = np.ones(len(Date))*11
        State = np.repeat('Berlin', len(Date))
        Gender = np.ones(len(Date))
        Age = np.repeat('A00-A04', len(Date))
        Confirmed = np.arange(1,len(Date)+1)

        df = pd.DataFrame(np.array([ID_State, State, ID_County, County, Gender, Age,  Date, Confirmed]).T,
                          columns = ['ID_State', 'State', 'ID_County', 'County', 'Gender', 'Age_RKI', 'Date', 'Confirmed'])
        Date = ['2020-01-03', '2020-01-04','2020-01-21', '2020-01-26', '2020-01-31']
        Date = np.array([np.datetime64(x) for x in Date])
        ID_County = np.ones(len(Date))*11002
        County = np.repeat('Berlin-Neu-Koeln', len(Date))
        ID_State = np.ones(len(Date))*11
        State = np.repeat('Berlin', len(Date))
        Gender = np.ones(len(Date))
        Age = np.repeat('A00-A04', len(Date))
        Confirmed = np.arange(1,len(Date)+1)

        df = df.append(pd.DataFrame(np.array([ID_State, State, ID_County, County, Gender, Age,  Date, Confirmed]).T,
                          columns = ['ID_State', 'State', 'ID_County', 'County', 'Gender', 'Age_RKI', 'Date', 'Confirmed']))

        Date = ['2020-01-02', '2020-01-03', '2020-01-04','2020-01-20', '2020-01-21', '2020-01-25', '2020-01-26', '2020-01-30', '2020-01-31']
        Date = np.array([np.datetime64(x) for x in Date])
        ID_State = np.ones(len(Date))*11000
        State = np.repeat('SK Berlin', len(Date))
        Age = np.repeat('A00-A04', len(Date))
        Confirmed = np.array([1, 4, 6, 9, 12, 16, 20, 25, 30])

        result_df = pd.DataFrame(np.array([ID_State, State, Age, Date, Confirmed]).T, columns = ['ID_County', 'County', 'Age_RKI', 'Date', 'Confirmed'])
        result_df = result_df.astype({'ID_County': 'int64', 'Confirmed': 'int64'})


        test_df = getRKIData.fuse_berlin(df)
        test_df = test_df.groupby( ['ID_County', 'County', 'Age_RKI', 'Date'])\
                             .agg({'Confirmed': sum})
        test_df = test_df.groupby(level=[1,2]).cumsum().reset_index()
        pd.testing.assert_frame_equal(test_df, result_df)

if __name__ == '__main__':
    unittest.main()