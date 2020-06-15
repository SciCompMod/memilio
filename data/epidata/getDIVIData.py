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
      # https://www.divi.de/images/Dokumente/Tagesdaten_Intensivregister_CSV/DIVI-Intensivregister_2020-06-03_09-15.csv
      df = gd.loadCsv('1674-divi-intensivregister-tagesreport-2020-06-04/file',
                      #apiUrl='https://www.divi.de/divi-intensivregister-tagesreport-archiv/', extension='')
                      apiUrl = 'https://www.divi.de/images/Dokumente/Tagesdaten_Intensivregister_CSV/', extension='')

      #https://www.divi.de/divi-intensivregister-tagesreport-archiv/1674-divi-intensivregister-tagesreport-2020-06-04/file

      # output data to not always download it
      gd.write_dataframe(df, directory, filename, "json")

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

