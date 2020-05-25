import os
import sys
from urllib.request import urlopen
import json
import pandas
import matplotlib.pyplot as plt

from epidemiology.epidata import outputDict as od
from epidemiology.epidata import getDataIntoPandasDataFrame as gd


def loadCsv( githubUrl = 'https://github.com/datasets/covid-19/tree/master/data/', 
             CSVfile  = 'time-series-19-covid-combined.csv' ):
    """ Loads data in CSV format from John Hopkins github. (pandas DataFrame)

    This routine loads ArcGIS data sets in CSV format of the given public data 
    item ID into a pandas DataFrame and returns the DataFrame. 

    Keyword arguments:
    githubUrl -- John Hopkins github url
    CSVfile -- file name

    """
    url = githubUrl + CSVfile
    #print(url)

    df = pandas.read_csv( url )

    return df

def main():

   [read_data, make_plot, out_form, out_folder] = gd.cli('Downloads data from JH')   

   filename = os.path.join(out_folder, "FullData_JohnHopkins.json")

   if(read_data):
      # if once dowloaded just read json file
      try:
         df = pandas.read_json(filename)
      except ValueError:
         exit_string = "Error: The file: " + filename + "does not exist. Call program without -r flag to get it."
         sys.exit(exit_string)
   else:
       # Get data:
      # https://raw.githubusercontent.com/datasets/covid-19/master/data/time-series-19-covid-combined.csv
      df = loadCsv(githubUrl = 'https://raw.githubusercontent.com/datasets/covid-19/master/data/')

      # convert "Datenstand" to real date:
      #df.Datenstand = pandas.to_datetime( df.Datenstand, format='%d.%m.%Y, %H:%M Uhr').dt.tz_localize('Europe/Berlin')  

      # output data to not always download it
      df.to_json(filename)

   # Preperation for plotting/output:

   outForm = od.outForm[out_form][0]
   outFormEnd = od.outForm[out_form][1]
   outFormSpec = od.outForm[out_form][2]

   df.rename({'Country/Region': 'CountryRegion', 'Province/State': 'ProvinceState'}, axis=1, inplace=True)
   print("Available columns:", df.columns)


   ########### Coutries ##########################

   gb = df.groupby( ['CountryRegion', 'Date']).agg({"Confirmed": sum, "Recovered": sum, "Deaths": sum})

   #print(df[df.CountryRegion == "Germany"])
   #print(gb)
   #print(gb.reset_index()[gb.reset_index().CountryRegion == "Germany"])

   getattr(gb.reset_index(), outForm)(os.path.join(out_folder, "all_countries_jh" + outFormEnd), **outFormSpec)

   # Check what about external provinces. Should they be added?

   ################ Countries with given States ################################
   
   # Get all countries where States are given 
   dfD = df[~df["ProvinceState"].isnull()]

   gb = dfD.groupby( ['CountryRegion', 'ProvinceState', 'Date']).agg({"Confirmed": sum, "Recovered": sum, "Deaths": sum})

   getattr(gb.reset_index(), outForm)(os.path.join(out_folder, "all_provincestate_jh" + outFormEnd), **outFormSpec)

   #print(dfD[dfD.ProvinceState=="Saskatchewan"])
   #print(gb.reset_index()[gb.reset_index().ProvinceState=="Saskatchewan"])

   # TODO: How to handle empty values which become NaN in the beginnin but after woking on the data its just 0.0
   # One solution is to preserve them with : df['b'] = df['b'].astype(str)
   # However, what to do with the cases where after some times values occur? Do those cases exist? 

   # TODO: US more detailes 

 
if __name__ == "__main__":

    main()

