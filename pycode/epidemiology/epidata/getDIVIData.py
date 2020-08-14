import os
import sys
import pandas
from datetime import timedelta, date

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd


def adjust_data(df, date_of_data):
    # rename column 'kreis' of first date to match data of following days
    if date_of_data == date(2020, 4, 24):
        df.rename({'kreis': 'gemeindeschluessel'}, axis=1, inplace=True)

    # remove empty column
    if date_of_data == date(2020, 4, 29):
        df.drop(columns='Unnamed: 0', inplace=True)

    # add dates for data until 27.4.
    if date_of_data <= date(2020, 4, 27):
        date_str = date_of_data.strftime("%Y-%m-%d") + " 09:15:00"
        df.insert(loc=len(df.columns), column='daten_stand', value=date_str)

    # add 'bundesland' for data from 25.4.
    if date_of_data == date(2020, 4, 25):
        df.insert(loc=0, column='bundesland', value=df['gemeindeschluessel'].floordiv(1000))

    return df


def call_call_url(url_prefix, call_number, new_found=""):

    call_url = url_prefix + "/viewdocument/" + str(call_number)

    # empty data frame
    df = pandas.DataFrame()

    try:
        df = pandas.read_csv(call_url)
        if not new_found == "":
            print("New cal number found. Please copy the following to call_number_dict to increase runtime: "
                  + new_found)
    except OSError as e:
        exit_string = "ERROR: URL " + call_url + " could not be opened. " \
                      + "Hint: check your internet connection."

        sys.exit(exit_string)
    except:
        pass

    return df


def download_data_for_one_day(last_number, download_date):
    # define call numbers for dates where call number doesn't increase by 1

    call_number_dict = {date(2020, 4, 24): 3974,
                        date(2020, 5, 6): 3691,
                        date(2020, 6, 5): 3842,
                        date(2020, 6, 12): 3906,
                        date(2020, 6, 26): 3774,
                        date(2020, 6, 28): 3777,
                        date(2020, 6, 29): 3830,
                        date(2020, 6, 30): 3840,
                        date(2020, 7, 1): 3928,
                        date(2020, 7, 2): 3951,
                        date(2020, 7, 6): 3958,
                        date(2020, 7, 7): 3961,
                        date(2020, 7, 8): 3964,
                        date(2020, 7, 10): 3987,
                        date(2020, 7, 11): 4065,
                        date(2020, 7, 12): 4069,
                        date(2020, 7, 13): 4073,
                        date(2020, 7, 14): 4081,
                        date(2020, 7, 15): 4108,
                        date(2020, 7, 17): 4127,
                        date(2020, 7, 18): 4143,
                        date(2020, 7, 20): 4156,
                        date(2020, 7, 22): 4162,
                        date(2020, 7, 27): 4173,
                        date(2020, 7, 28): 4177,
                        date(2020, 7, 29): 4181,
                        date(2020, 7, 30): 4184,
                        date(2020, 8, 3): 4204,
                        date(2020, 8, 4): 4316,
                        date(2020, 8, 5): 4319,
                        date(2020, 8, 6): 4324,
                        date(2020, 8, 8): 4327,
                        date(2020, 8, 9): 4330,
                        date(2020, 8, 11): 4412,
                        date(2020, 8, 13): 4430,
                        date(2020, 8, 14): 4680
                        }

    call_date = download_date.strftime("%Y-%m-%d")

    # first links have different upload time
    time_shift_date = date(2020, 6, 5)
    if download_date < time_shift_date:
        call_time = "-09-15"
    else:
        call_time = "-12-15"

    # need extension "-2" for dates between 12.6. and 25.6.
    ext = ""
    if date(2020, 6, 12) <= download_date <= date(2020, 6, 25):
        ext = "-2"

    url_prefix = "https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" \
                 + call_date + call_time + ext

    # construction of link is different for different time periods
    # for the latest dates no call number is needed
    # number in url starts at values in call_number_dict and increases every day by 1

    sign_dict = {0:1,
                 1:-1}

    check_num_calls = 1
    if download_date in call_number_dict.keys():

        call_number = call_number_dict[download_date]

        check_num_calls = check_num_calls + 1
        df = call_call_url(url_prefix, call_number)

    else:
        call_string = ""

        # It is most likely that the difference is between 1 and 2
        for sign in range(2):
            for delta in range(1,300):

               call_number = last_number + sign_dict[sign]*delta
               if delta != 1 and delta != 2:
                  call_string = "date(" + download_date.strftime("%Y,%-m,%-d") + "): " + str(call_number) + ","

               print("Number of calls: ", check_num_calls)
               check_num_calls = check_num_calls + 1
               df = call_call_url(url_prefix, call_number,  call_string)

               if not df.empty:
                   return [call_number, df]

            # case with same call_number, which is very unlikely
            call_number = last_number
            check_num_calls = check_num_calls + 1
            df = call_call_url(url_prefix, call_number)

    return [call_number, df]


def get_divi_data(read_data=dd.defaultDict['read_data'],
                  make_plot=dd.defaultDict['make_plot'],
                    out_form=dd.defaultDict['out_form'],
                  out_folder=dd.defaultDict['out_folder'],
                  start_date=date(2020, 4, 24),
                  end_date=date.today()):

    # First csv data on 24-04-2020
    if start_date < date(2020, 4, 24):
        start_date = date(2020, 4, 24)
        print("Warning: First data available on 24-04-2020")

    directory = os.path.join(out_folder, 'Germany/')
    gd.check_dir(directory)

    filename = "FullData_DIVI"

    if read_data:
        file_in = os.path.join(directory, filename + ".json")
        # if once downloaded just read json file
        try:
            print("file_in: ", file_in)
            df = pandas.read_json(file_in)
        except ValueError:
            exit_string = "Error: The file: " + file_in + "does not exist. Call program without -r flag to get it."
            sys.exit(exit_string)
    else:
        # Get data:
        # start with empty dataframe
        df = pandas.DataFrame()

        delta = timedelta(days=1)

        last_number = 0
        while start_date <= end_date:

            [last_number, df2] = download_data_for_one_day(last_number, start_date)

            if not df2.empty:
                # data of first days needs adjustment to following data
                if start_date <= date(2020, 4, 29):
                    df2 = adjust_data(df2, start_date)
                df = df.append(df2, ignore_index=True)
                print("Success: Data of date " + start_date.strftime("%Y-%m-%d") + " has been included to dataframe")
            else:
                print("Warning: Data of date " + start_date.strftime("%Y-%m-%d") + " is not included to dataframe")

            start_date += delta

        # output data to not always download it
        if not df.empty:
            gd.write_dataframe(df, directory, filename, "json")
        else:
            exit_string = "Something went wrong, dataframe is empty."
            sys.exit(exit_string)

    # change column names
    df.rename(dd.GerEng, axis=1, inplace=True)

    df.Date = pandas.to_datetime(df.Date, format='%Y-%m-%d %H:%M')

    # insert names of  states
    df.insert(loc=0, column="State", value=df.ID_State)
    for key, value in dd.State.items():
        df.loc[df["State"] == key, ["State"]] = value

    # insert names of counties
    df.insert(loc=0, column="County", value=df.ID_County)
    for key, value in dd.County.items():
        df.loc[df["County"] == key, ["County"]] = value

    print("Available columns:", df.columns)

    # ID_County and ID_State is as defined by the "Amtlicher Gemeindeschlüssel (AGS)"
    # which is also used in the RKI data as ID_County and ID_State
    # https://de.wikipedia.org/wiki/Liste_der_Landkreise_in_Deutschland

    # The column faelle_covid_im_bundesland exits only in the data from the first day
    # The columns ICU does not exist for the 24.4.
    # and ICU_ventilated does not exist for the 24.4. and 25.4.

    # reporting_hospitals number of reporting hospitals
    # ICU is the number of covid patients in reporting hospitals
    # ICU_ventilated is the number of ventilated covid patients in reporting hospitals
    # free_ICU is the number of free ICUs in reporting hospitals
    # occupied_ICU is the number of occupied ICUs in in reporting hospitals

    # write data for counties to file

    df_counties = df[["County", "ID_County", "ICU", "ICU_ventilated", "Date"]]
    filename = "county_divi"
    gd.write_dataframe(df_counties, directory, filename, "json")

    # write data for states to file

    df_states = df.groupby(["ID_State", "State", "Date"]).agg({"ICU": sum, "ICU_ventilated": sum})
    df_states = df_states.reset_index()
    df_states.sort_index(axis=1, inplace=True)
    filename = "state_divi"
    gd.write_dataframe(df_states, directory, filename, "json")

    # write data for germany to file

    df_ger = df.groupby(["Date"]).agg({"ICU": sum, "ICU_ventilated": sum})
    df_ger = df_ger.reset_index()
    df_ger.sort_index(axis=1, inplace=True)
    filename = "germany_divi"
    gd.write_dataframe(df_ger, directory, filename, "json")


def main():
    [read_data, make_plot, out_form, out_folder] = gd.cli('Downloads data from DIVI')
    get_divi_data(read_data, make_plot, out_form, out_folder)


if __name__ == "__main__":
    main()
