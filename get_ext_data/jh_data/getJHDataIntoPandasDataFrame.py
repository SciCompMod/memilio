import sys
from urllib.request import urlopen
import json
import pandas
import matplotlib.pyplot as plt


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

def main(get_data, read_data, make_plot):

   FullJSONData = 'FullData_JohnHopkins.json'

   if(get_data):

      # Get data:
      # https://raw.githubusercontent.com/datasets/covid-19/master/data/time-series-19-covid-combined.csv
      df = loadCsv(githubUrl = 'https://raw.githubusercontent.com/datasets/covid-19/master/data/')

      # convert "Datenstand" to real date:
      #df.Datenstand = pandas.to_datetime( df.Datenstand, format='%d.%m.%Y, %H:%M Uhr').dt.tz_localize('Europe/Berlin')  

      # output data to not always download it
      df.to_json(FullJSONData)

   elif(read_data):
      # if once dowloaded just read json file
      df = pandas.read_json(FullJSONData)

   # Preperation for plotting/output:

   df.rename({'Country/Region': 'CountryRegion', 'Province/State': 'ProvinceState'}, axis=1, inplace=True)
   print("Available columns:", df.columns)


   ########### Coutries ##########################

   gb = df.groupby( ['CountryRegion', 'Date']).agg({"Confirmed": sum, "Recovered": sum, "Deaths": sum})

   #print(df[df.CountryRegion == "Germany"])
   #print(gb)
   #print(gb.reset_index()[gb.reset_index().CountryRegion == "Germany"])

   gb.reset_index().to_json("all_countries.json", orient='records', date_format='iso')

   # Check what about external provinces. Should they be added?

   ################ Countries with given States ################################
   
   # Get all countries where States are given 
   dfD = df[~df["ProvinceState"].isnull()]

   gb = dfD.groupby( ['CountryRegion', 'ProvinceState', 'Date']).agg({"Confirmed": sum, "Recovered": sum, "Deaths": sum})

   gb.reset_index().to_json("all_provincestate.json", orient='records', date_format='iso')

   #print(dfD[dfD.ProvinceState=="Saskatchewan"])
   #print(gb.reset_index()[gb.reset_index().ProvinceState=="Saskatchewan"])

   # TODO: How to handle empty values which become NaN in the beginnin but after woking on the data its just 0.0
   # One solution is to preserve them with : df['b'] = df['b'].astype(str)
   # However, what to do with the cases where after some times values occur? Do those cases exist? 

   # TODO: US more detailes 

 
if __name__ == "__main__":

   GET_DATA = True
   READ_DATA = False
   MAKE_PLOT = True

   largv = len(sys.argv)

   if largv > 1:
      for i in range(1,largv):

          arg = sys.argv[i]

          if "READ_DATA" in arg:

             arg_split = arg.split("=")
             if len(arg_split) == 2:
                 READ_DATA = arg_split[1]
                 GET_DATA=False
             else:
                 print("Warning: your argument:", arg, "is ignored. It has to be in the form as: READ_DATA=True")

          elif "MAKE_PLOT" in arg:

             arg_split = arg.split("=")
             if len(arg_split) == 2:
                 MAKE_PLOT = arg_split[1]
             else:
                 print("Warning: your argument:", arg, "is ignored. It has to be in the form as: MAKE_PLOT=False")
          else:
             print("Warning: your argument:", arg, "is ignored.")

   main(GET_DATA, READ_DATA, MAKE_PLOT)
