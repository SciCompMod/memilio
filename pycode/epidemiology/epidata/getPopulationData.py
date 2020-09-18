import os
import sys
import pandas
from collections import namedtuple
from epidemiology.epidata  import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd


def get_population_data(read_data=dd.defaultDict['read_data'],
                        out_form=dd.defaultDict['out_form'],
                        out_folder=dd.defaultDict['out_folder']):

   print("Warning: getpopulationdata is not working correctly. A bug workaround has been applied.")

   Data = namedtuple("Data", "filename item columns_wanted filename_out")

   d1 = Data("FullDataB", '5dc2fc92850241c3be3d704aa0945d9c_2', ["LAN_ew_RS", 'LAN_ew_GEN','LAN_ew_EWZ'], "PopulStates")
   d2 = Data("FullDataL", 'b2e6d8854d9744ca88144d30bef06a76_1', ['RS', 'GEN','EWZ'], "PopulCounties")

   #d = [d1, d2]
   d = [d1]

   directory = os.path.join(out_folder, 'Germany/')
   gd.check_dir(directory)

   for i in range(len(d)):
      get_one_data_set(read_data, out_form, directory, d[i])

def get_one_data_set(read_data, out_form, directory, d):

   if(read_data):
      # if once dowloaded just read json file
      file = os.path.join(directory, d.filename + ".json")

      try:
         df = pandas.read_json(file)

      except ValueError:
         exit_string = "Error: The file: " + file + " does not exist. Call program without -r flag to get it."
         sys.exit(exit_string)
   else:

      # Supported data formats:
      load = { 
         'csv': gd.loadCsv,
         'geojson': gd.loadGeojson
       }

      # Get data:
      df = load['csv'](d.item)

      # output data to not always download it
      gd.write_dataframe(df, directory, d.filename, "json")

   print("Available columns:", df.columns)

   # Filter data for Bundesland/Landkreis and Einwohnerzahl (EWZ)
   dfo = df[d.columns_wanted]
   dfo = dfo.rename(columns=dd.GerEng)
   gd.write_dataframe(dfo, directory, d.filename_out, out_form)


def main():

   [read_data, out_form, out_folder] = gd.cli("population")
   get_population_data(read_data, out_form, out_folder)


if __name__ == "__main__":

   main()
