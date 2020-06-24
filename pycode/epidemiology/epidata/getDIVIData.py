import os
import sys
import pandas

import getDataIntoPandasDataFrame as gd
import defaultDict as dd

def adjust_first_data(df, start_date):
   from datetime import date

   # rename column 'kreis' of first date to match data of following days
   if start_date == date(2020, 4, 24):
      df.rename({'kreis': 'gemeindeschluessel'}, axis=1, inplace=True)

   # remove empty column
   if start_date == date(2020, 4, 29):
      df.drop(columns = 'Unnamed: 0', inplace=True)

   # add dates for data until 27.4.
   if start_date <= date(2020, 4, 27):
      date_str = start_date.strftime("%Y-%m-%d")+" 09:15:00"
      df.insert(loc = len(df.columns), column = 'daten_stand', value = date_str)

   # add 'bundesland' for data from 25.4.
   if start_date == date(2020, 4, 25):
      df = df.assign(bundesland = df['gemeindeschluessel'].floordiv(1000))

   return df

def get_divi_data(read_data=dd.defaultDict['read_data'],
                make_plot=dd.defaultDict['make_plot'],
                out_form=dd.defaultDict['out_form'],
                out_folder=dd.defaultDict['out_folder']):

   directory = os.path.join(out_folder, 'Germany/')
   gd.check_dir(directory)

   filename = "Data_Divi"

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
      # First csv data on 24-04-2020

      from datetime import timedelta, date

      start_date = date(2020, 4, 24)
      end_date = date.today()
      delta = timedelta(days=1)
      data1_end_date = date(2020, 4, 28)
      # start with empty dataframe
      df = pandas.DataFrame()

      try:
         # Read urls from file and save in a dictionary
         url_dict = {}
         url_dict_file = open("Germany/url_dictionary.txt", "r")
         for line in url_dict_file:
            date = line.split(" ")[0]
            url= line.split(" ")[1]
            url_dict.update({date: url})
         url_dict_file.close()

         # Download data with given urls
         while start_date <= end_date:

            call_url = url_dict[start_date.strftime("%Y-%m-%d")]

            df2 = pandas.DataFrame()
            df2 = gd.loadCsv('file', extension='', apiUrl = call_url)

            if (df2.empty != True):
               # data of first days needs adjustment to following data
               df2 = adjust_first_data(df2, start_date)
               df = df.append(df2, ignore_index=True)
               print("Success: Data of date " + start_date.strftime("%Y-%m-%d") + " has been included to dataframe")
            else:
               print("Warning: Data of date " + start_date.strftime("%Y-%m-%d") + " is not included to dataframe")

            start_date += delta

      except:
         # Construct urls
         # dict which cases where the next data is not +2
         date_number_dict = {
             "24-04-2020": 1580,
             "25-04-2020": 1578,
             "26-04-2020": 1579,
             "27-04-2020": 1582,
             "2020-04-28": 1585,
             "2020-04-30": 1588,
             "2020-05-06": 1601,
             "2020-05-08": 1609,
             "2020-05-09": 1610,
             "2020-05-10": 1611,
             "2020-05-11": 1614,
             "2020-05-13": 1622,
             "2020-05-14": 1623,
             "2020-05-15": 1629,
             "2020-05-16": 1628,
             "2020-05-17": 1633,
             "2020-05-18": 1634,
             "2020-05-19": 1637,
             "2020-05-20": 1640,
             "2020-05-28": 1658,
             "2020-06-04": 1674,
             "2020-06-10": 1687,
             "2020-06-13": 1694,
             "2020-06-18": 1705,
             "2020-06-19": 1710,
             "2020-06-20": 1715,
             "2020-06-21": 1716,
             "2020-06-22": 1720,
         }

         list_call_extension = ["", "-1", "-09-15", "-csv"]
         list_call_link = ["-divi-intensivregister-tagesreport-", "-divi-intensivregister-"]
         call_number=0

         while start_date <= end_date:

            # first links have different data order
            if start_date < data1_end_date:
               call_date = start_date.strftime("%d-%m-%Y")
            else:
               call_date = start_date.strftime("%Y-%m-%d")

            # get strange number of url
            # usual rule it grows by 2 from day to day
            # start point and exception is stored in dict
            try:
               call_number = date_number_dict[call_date]
            except KeyError:
               pass

            df2 = pandas.DataFrame()
            for call_link in  list_call_link:
               for call_extension in list_call_extension:
                  call_url =  'https://www.divi.de/divi-intensivregister-tagesreport-archiv/'\
                      + str(call_number) + call_link + call_date + call_extension + '/'

                  try:
                     df2 = gd.loadCsv('file',
                         extension='',
                         apiUrl = call_url)
                     # Write url to url dictionary text file
                     url_dict_file = open("Germany/url_dictionary.txt", "a")
                     url_dict_file.write(start_date.strftime("%Y-%m-%d")+" "+call_url+" \n")
                     url_dict_file.close()
                  except:
                     print("Information: URL: ",  call_url, " not knwon")

            if (df2.empty != True):
               # data of first days needs adjustment to following data
               df2 = adjust_first_data(df2, start_date)
               df = df.append(df2, ignore_index=True)
               print("Success: Data of date " + call_date + " has been included to dataframe")
            else:
               print("Warning: Data of date " + call_date + " is not included to dataframe")

            call_number = call_number + 2
            start_date += delta

      # output data to not always download it
      if(df.empty != True):
         gd.write_dataframe(df, directory, filename, "json")
      else:
         exit_string = "Something went wrong, dataframe is empty."
         sys.exit(exit_string)

   # change column names
   df.rename(dd.GerEng, axis=1, inplace=True)
   print("Available columns:", df.columns)

   df.Date = pandas.to_datetime(df.Date, format='%Y-%m-%d %H:%M')

   # ID_County is as defined by the "Amtlicher Gemeindeschlüssel (AGS)"
   # which is also used in the RKI data as ID_County
   # https://de.wikipedia.org/wiki/Liste_der_Landkreise_in_Deutschland

   # The column faelle_covid_im_bundesland exits only in the data from the first day

   # The columns falle_covid_aktuell does not exist for the 24.4.
   # and faelle_covid_aktuell_beatmet does not exist for the 24.4. and 25.4.
   # TODO: Check whether they coincide with the data from RKI

def main():

   [read_data, make_plot, out_form, out_folder] = gd.cli('Downloads data from DIVI')
   get_divi_data(read_data, make_plot, out_form, out_folder)

 
if __name__ == "__main__":

    main()

