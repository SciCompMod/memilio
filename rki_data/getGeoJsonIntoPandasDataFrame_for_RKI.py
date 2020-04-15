import sys
from urllib.request import urlopen
import json
import pandas
import matplotlib.pyplot as plt

GET_DATA = True
PLOT_DATA = False

def loadGeojson( itemId, apiUrl = 'https://opendata.arcgis.com/datasets/', 
                 extension = 'geojson' ):
    """ Loads ArcGIS data sets in GeoJSON format. (pandas DataFrame)

    This routine loads ArcGIS data sets in GeoJSON format of the given public 
    data item ID into a pandas DataFrame and returns the DataFrame. Trivial 
    information gets removed by JSON normalization and dropping of always 
    constant data fields.

    Keyword arguments:
    itemId -- ArcGIS public data item ID (string)
    apiUrl -- ArcGIS data sets API URL (string, default
              'https://opendata.arcgis.com/datasets/')
    extension -- Data format extension (string, default 'geojson')

    """
    url = apiUrl + itemId + '_0.' + extension

    with urlopen( url ) as res:
        data = json.loads( res.read().decode() )

    # Shape data:
    df = pandas.json_normalize( data, 'features' )

    # Make-up (removing trivial information):
    df.drop( columns = ['type', 'geometry'], inplace = True )
    df.rename( columns = lambda s: s[11:], inplace = True )

    return df

def loadCsv( itemId, apiUrl = 'https://opendata.arcgis.com/datasets/', 
                 extension = 'csv' ):
    """ Loads ArcGIS data sets in CSV format. (pandas DataFrame)

    This routine loads ArcGIS data sets in CSV format of the given public data 
    item ID into a pandas DataFrame and returns the DataFrame. 

    Keyword arguments:
    itemId -- ArcGIS public data item ID (string)
    apiUrl -- ArcGIS data sets API URL (string, default
              'https://opendata.arcgis.com/datasets/')
    extension -- Data format extension (string, default 'csv')

    """
    url = apiUrl + itemId + '_0.' + extension

    df = pandas.read_csv( url )

    return df

def main(get_data, plot_data):

   if(get_data):

      # Supported data formats:
      load = { 
         'csv': loadCsv,
         'geojson': loadGeojson
       }

      # ArcGIS public data item ID:
      itemId = 'dd4580c810204019a7b8eb3e0b329dd6'

      # Get data:
      df = load['csv'](itemId)

      # Correct Timestampes:
      # for col in [ 'Meldedatum', 'Refdatum' ]:
      #   df[col] = df[col].astype( 'datetime64' ).dt.tz_localize('Europe/Berlin')

      # convert "Datenstand" to real date:
      df.Datenstand = pandas.to_datetime( df.Datenstand, format='%d.%m.%Y, %H:%M Uhr').dt.tz_localize('Europe/Berlin')  

      # output data to not always download it
      df.to_json(r"FullData.json")

   elif(plot_data):
      df = pandas.read_json("FullData.json")
      # read json file

   # Preperation for plotting/output:

   # Correct Timestampes:
   for col in [ 'Meldedatum', 'Refdatum' ]:
      df[col] = df[col].astype( 'datetime64[ns]' ).dt.tz_localize('Europe/Berlin')

   # Use "Meldedatum" as time
   # "Refdatum" may act different to official describtion
   # on https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0,
   # be careful (big difference identified between "Refdatum" and "Meldedatum").
   dateToUse = 'Meldedatum'
   df.sort_values( dateToUse, inplace = True )

   # NeuerFall: Infected (incl. recovered) over "Meldedatum":
   gbNF = df[df.NeuerFall >= 0].groupby( dateToUse ).sum()
   gbNFcumsum = gbNF.AnzahlFall.cumsum()
   gbNFcumsum.to_json("gbNF.json")
   gbNFcumsum.plot( title = 'COVID-19 infections', grid = True, 
                               style = '-o' )
   plt.tight_layout()
   plt.show()

   # print(df.columns)

   # NeuerFall: Infected (incl. recovered) over "Meldedatum" for every state ("Bundesland"):
   gbNFst = df[df.NeuerFall >= 0].groupby( ['IdBundesland', 'Bundesland', dateToUse]).sum()

   gbNFstcumsum = gbNFst.AnzahlFall.cumsum().reset_index()
   # print(gbNFst)
   print(gbNFstcumsum)

   # output json
   gbNFstcumsum.to_json("gbNF_state.json", orient='records')

   # output nested json
   gbNFstcumsum.groupby(['IdBundesland', 'Bundesland'], as_index=False) \
               .apply(lambda x: x[[dateToUse,'AnzahlFall']].to_dict('r')) \
               .reset_index().rename(columns={0:'Dates'})\
               .to_json("gbNF_state_nested.json", orient='records')


   #gbNFstcumsum.plot( title = 'COVID-19 infections per state',x='AnzahlFall',y='dateToUse', grid = True,
   #                            style = '-o' )
   #plt.tight_layout()
   #plt.show()


   # Dead over "Meldedatum":
   gbNT = df[df.NeuerTodesfall >= 0].groupby( dateToUse ).sum()
   gbNT.AnzahlTodesfall.cumsum().plot( title = 'COVID-19 deads', grid = True,
                                    style = '-o' )
   plt.tight_layout()
   plt.show()

   # Dead by "Altersgruppe":
   gbNTAG = df[df.NeuerTodesfall >= 0].groupby( 'Altersgruppe' ).sum()
   gbNTAG.AnzahlTodesfall.plot( title = 'COVID-19 deads', grid = True, 
                             kind = 'bar' )
   plt.tight_layout()
   plt.show()


if __name__ == "__main__":


   largv = len(sys.argv)

   if largv > 1:
      for i in range(1,largv):

          arg = sys.argv[i]

          if "GET_DATA" in arg:

             arg_split = arg.split("=")
             if len(arg_split) == 2:
                 GET_DATA = arg_split[1]
             else:
                 pass

          elif "PLOT_DATA" in arg:

             arg_split = arg.split("=")
             if len(arg_split) == 2:
                 PLOT_DATA = arg_split[1]
                 GET_DATA=False
             else:
                 pass

   main(GET_DATA, PLOT_DATA)
