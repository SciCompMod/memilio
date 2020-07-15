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
      df.drop(columns = 'Unnamed: 0', inplace=True)

   # add dates for data until 27.4.
   if date_of_data <= date(2020, 4, 27):
      date_str = date_of_data.strftime("%Y-%m-%d")+" 09:15:00"
      df.insert(loc = len(df.columns), column = 'daten_stand', value = date_str)

   # add 'bundesland' for data from 25.4.
   if date_of_data == date(2020, 4, 25):
      df = df.assign(bundesland = df['gemeindeschluessel'].floordiv(1000))

   return df

def download_data_for_one_day(download_date):

   # define time periods where construction of urls is different
   data1_start_date = date(2020, 4, 30)
   data1_end_date = date(2020, 5, 5)
   time_shift_date = date(2020, 6, 5)
   data2_start_date = date(2020, 6, 12)
   data2_end_date = date(2020, 6, 25)

   call_date = download_date.strftime("%Y-%m-%d")

   # first links have different upload time
   if download_date < time_shift_date:
      call_time = "-09-15"
   else:
      call_time = "-12-15"

   # construction of link is different for different time periods
   if download_date >= data1_start_date and download_date <= data1_end_date:
      # number in url starts at 3980 and increases every day by 1
      call_number = (download_date - data1_start_date).days  + 3980
      call_url = "https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" + call_date + call_time + "/viewdocument/" + str(call_number)
   else:
      if download_date >= data2_start_date and download_date <= data2_end_date:
         # number in url starts at 3906 and increases every day by 1
         call_number = (download_date - data2_start_date).days + 3906
         call_url = "https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" + call_date + call_time + "-2/viewdocument/" + str(call_number)

      else:
         call_url = "https://www.divi.de/divi-intensivregister-tagesreport-archiv-csv/divi-intensivregister-" + call_date + call_time + "/download"

   df = pandas.DataFrame()

   try:
      df = pandas.read_csv(call_url)
   except OSError as e:
      exit_string = "ERROR: URL " + call_url + " could not be opened."
      sys.exit(exit_string)
   except:
      print("Information: URL: ",  call_url, " not knwon")

   return df

def get_divi_data(read_data=dd.defaultDict['read_data'],
                make_plot=dd.defaultDict['make_plot'],
                out_form=dd.defaultDict['out_form'],
                out_folder=dd.defaultDict['out_folder'],
                start_date=date(2020,4,24),
                end_date=date.today()):

   # First csv data on 24-04-2020
   if start_date < date(2020,4,24):
      start_date = date(2020,4,24)
      print("Warning: First data available on 24-04-2020")

   directory = os.path.join(out_folder, 'Germany/')
   gd.check_dir(directory)

   filename = "FullData_DIVI"

   if(read_data):
       file_in = os.path.join(directory, filename+".json")
       # if once dowloaded just read json file
       try:
          df = pandas.read_json(file_in)
       except ValueError:
          exit_string = "Error: The file: " + file_in + "does not exist. Call program without -r flag to get it."
          sys.exit(exit_string)

   else:
      # Get data:
      # start with empty dataframe
      df = pandas.DataFrame()

      delta = timedelta(days=1)

      while start_date <= end_date:

         df2 = download_data_for_one_day(start_date)

         if (df2.empty != True):
            # data of first days needs adjustment to following data
            df2 = adjust_data(df2, start_date)
            df = df.append(df2, ignore_index=True)
            print("Success: Data of date " + start_date.strftime("%Y-%m-%d") + " has been included to dataframe")
         else:
            print("Warning: Data of date " + start_date.strftime("%Y-%m-%d") + " is not included to dataframe")

         start_date += delta

      # output data to not always download it
      if(df.empty != True):
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
   
   df_counties = df[["County","ID_County","ICU","ICU_ventilated","Date"]]
   filename = "county_divi"
   gd.write_dataframe(df_counties, directory, filename, "json")
   
   
   # write data for states to file
   
   df_states = df.groupby(["ID_State","State","Date"]).agg({"ICU": sum, "ICU_ventilated": sum})
   df_states = df_states.reset_index()
   df_states.sort_index(axis=1, inplace=True)
   filename = "state_divi"
   gd.write_dataframe(df_states, directory, filename, "json")
   
   
   # write data for germany to file
   
   df_ger = df.groupby(["Date"]).agg({"ICU": sum, "ICU_ventilated": sum})
   df_ger =df_ger.reset_index()
   df_ger.sort_index(axis=1, inplace=True)
   filename = "germany_divi"
   gd.write_dataframe(df_ger, directory, filename, "json")
   
def main():

   [read_data, make_plot, out_form, out_folder] = gd.cli('Downloads data from DIVI')
   get_divi_data(read_data, make_plot, out_form, out_folder)

 
if __name__ == "__main__":

    main()

