import os
import sys
import pandas

import getDataIntoPandasDataFrame as gd
import defaultDict as dd

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

      #TODO: Write correct links to dictionary and first try to read dictionary from file.
      # Do following just if this file does not exist

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
      }

      from datetime import timedelta, date

      start_date = date(2020, 4, 24)
      end_date = date.today()
      delta = timedelta(days=1)
      data1_end_date = date(2020, 4, 28)
      # start with empty dataframe
      df = pandas.DataFrame()

      list_call_extension = ["", "-1", "-09-15", "-csv"]
      list_call_link = ["-divi-intensivregister-tagesreport-", "-divi-intensivregister-"]

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
               # output data to not always download it
               except:
                  print("Information: URL: ",  call_url, " not knwon")

         if (df2.empty != True):
            df = df.append(df2, ignore_index=True)
            print("Success: Data of date " + call_date + " has been included to dataframe")
         else:
            print("Warning: Data of date " + call_date + " is not included to dataframe")

         call_number = call_number + 2
         start_date += delta

      if(df.empty != True):
         gd.write_dataframe(df, directory, filename, "json")
      else:
         exit_string = "Something went wrong, dataframe is empty."
         sys.exit(exit_string)

   # change column names
   # TODO: Should they be translated to english and match column names from other data?
   df.rename({'anzahl_meldebereiche': 'anzahlMeldebereiche', 'faelle_covid_aktuell': 'faelleCovidAktuell',
   'faelle_covid_aktuell_beatmet': 'faelleCovidAktuellBeatmet', 'anzahl_standorte': 'anzahlStandorte',
   'betten_frei': 'bettenFrei', 'betten_belegt': 'bettenBelegt', 'daten_stand': 'datenStand'}, axis=1, inplace=True)
   print("Available columns:", df.columns)

def main():

   [read_data, make_plot, out_form, out_folder] = gd.cli('Downloads data from DIVI')
   get_divi_data(read_data, make_plot, out_form, out_folder)

 
if __name__ == "__main__":

    main()

