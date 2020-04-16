import sys
from urllib.request import urlopen
import json
import pandas
import matplotlib.pyplot as plt

GET_DATA = True
READ_DATA = False
MAKE_PLOT = True

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

def main(get_data, read_data, make_plot):

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

      # convert "Datenstand" to real date:
      df.Datenstand = pandas.to_datetime( df.Datenstand, format='%d.%m.%Y, %H:%M Uhr').dt.tz_localize('Europe/Berlin')  

      # output data to not always download it
      df.to_json(r"FullData.json")

   elif(read_data):
      # if once dowloaded just read json file
      df = pandas.read_json("FullData.json")

   # Preperation for plotting/output:

   print("Available columns:", df.columns)

   # Correct Timestampes:
   for col in [ 'Meldedatum', 'Refdatum' ]:
      df[col] = df[col].astype( 'datetime64[ns]' ).dt.tz_localize('Europe/Berlin')

   # Use "Meldedatum" as time
   # "Refdatum" may act different to official describtion
   # on https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0,
   # be careful (big difference identified between "Refdatum" and "Meldedatum").
   dateToUse = 'Meldedatum'
   df.sort_values( dateToUse, inplace = True )

   # Manipulate data to get rid of conditions: df.NeuerFall >= 0, df.NeuerTodesfall >= 0, df.NeuGenesen >=0
   # There might be a better way

   dfF = df

   dfF.loc[dfF.NeuerFall<0, ['AnzahlFall']] = 0
   dfF.loc[dfF.NeuerTodesfall<0, ['AnzahlTodesfall']] = 0
   dfF.loc[dfF.NeuGenesen<0, ['AnzahlGenesen']] = 0


   ######## Data for whole Germany all ages ##########

   # NeuerFall: Infected (incl. recovered) over "dateToUse":

   # make sum for one "dateToUse"
   gbNF = df[df.NeuerFall >= 0].groupby( dateToUse ).sum()

   # make cumulative sum of "AnzahlFall" for "dateToUse"
   gbNF_cs = gbNF.AnzahlFall.cumsum()

   # outout to json file
   gbNF_cs.to_json("infected.json")

   if(make_plot == True):
      # make plot
      gbNF_cs.plot( title = 'COVID-19 infections', grid = True, 
                               style = '-o' )
      plt.tight_layout()
      plt.show()

   # Dead over "Meldedatum":
   gbNT = df[df.NeuerTodesfall >= 0].groupby( dateToUse ).sum()
   gbNT_cs = gbNT.AnzahlTodesfall.cumsum()

   gbNT_cs.to_json("deaths.json")

   if(make_plot == True):
      gbNT_cs.plot( title = 'COVID-19 deaths', grid = True,
                                    style = '-o' )
      plt.tight_layout()
      plt.show()

      dfF.agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum}) \
         .plot( title = 'COVID-19 infections, deaths, recovered', grid = True,
                             kind = 'bar' )
      plt.tight_layout()
      plt.show()


   ############## Data for states all ages ################
   
   # NeuerFall: Infected (incl. recovered) over "dateToUse" for every state ("Bundesland"):
   gbNFst = df[df.NeuerFall >= 0].groupby( ['IdBundesland', 'Bundesland', dateToUse]).AnzahlFall.sum()

   gbNFst_cs = gbNFst.groupby(level=1).cumsum().reset_index()
  
   # print(gbNFst)
   # print(gbNFst_cs)

   # output json
   gbNFst_cs.to_json("infected_state.json", orient='records')

   # output nested json
   gbNFst_cs.groupby(['IdBundesland', 'Bundesland'], as_index=False) \
               .apply(lambda x: x[[dateToUse,'AnzahlFall']].to_dict('r')) \
               .reset_index().rename(columns={0:'Dates'})\
               .to_json("gbNF_state_nested.json", orient='records')


   # infected (incl recovered), deaths and recovered together 

   gbAllSt = dfF.groupby( ['IdBundesland', 'Bundesland', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllSt_cs = gbAllSt.groupby(level=1).cumsum().reset_index()

   # print(gbAllSt)
   # print(gbAllSt_cs)
 
   # output json
   gbAllSt_cs.to_json("all_state.json", orient='records')


   ############# Data for counties all ages ######################

   # NeuerFall: Infected (incl. recovered) over "dateToUse" for every county ("Landkreis"):
   gbNFc = df[df.NeuerFall >= 0].groupby( ['IdLandkreis', 'Landkreis', dateToUse]).sum()

   gbNFc_cs = gbNFc.groupby(level=1).AnzahlFall.cumsum().reset_index()
   
   #print(gbNFc)
   #print(gbNFc_cs)

   # output json
   gbNFc_cs.to_json("infected_county.json", orient='records')

   # infected (incl recovered), deaths and recovered together 

   gbAllC = dfF.groupby( ['IdLandkreis', 'Landkreis', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllC_cs = gbAllC.groupby(level=1).cumsum().reset_index()

   #print(gbAllC)
   #print(gbAllC_cs)

   # output json
   gbAllC_cs.to_json("all_county.json", orient='records')
   

   ######### Data whole Germany different gender ##################

   # infected (incl recovered), deaths and recovered together 

   gbAllG = dfF.groupby( ['Geschlecht', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllG_cs = gbAllG.groupby(level=0).cumsum().reset_index()

   print(gbAllG)
   print(gbAllG_cs)

   # output json
   gbAllG_cs.to_json("all_gender.json", orient='records') 

   if(make_plot == True):
      dfF.groupby( 'Geschlecht' ) \
         .agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum}) \
         . plot( title = 'COVID-19 infections, deaths, recovered', grid = True,
                             kind = 'bar' )
      plt.tight_layout()
      plt.show()

  
   ######### Data whole Germany different ages ####################

   # infected (incl recovered), deaths and recovered together 

   gbAllA = dfF.groupby( ['Altersgruppe', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllA_cs = gbAllA.groupby(level=0).cumsum().reset_index()

   print(gbAllA)
   print(gbAllA_cs)

   # output json
   gbAllA_cs.to_json("all_age.json", orient='records')

   if(make_plot == True):
      dfF.groupby( 'Altersgruppe' ) \
         .agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum}) \
         .plot( title = 'COVID-19 infections, deaths, recovered for diff ages', grid = True,
                             kind = 'bar' )
      plt.tight_layout()
      plt.show()


      # Dead by "Altersgruppe":
      gbNTAG = df[df.NeuerTodesfall >= 0].groupby( 'Altersgruppe' ).sum()
      gbNTAG.AnzahlTodesfall.plot( title = 'COVID-19 deaths', grid = True, 
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

          elif "READ_DATA" in arg:

             arg_split = arg.split("=")
             if len(arg_split) == 2:
                 READ_DATA = arg_split[1]
                 GET_DATA=False
             else:
                 pass

          elif "MAKE_PLOT" in arg:

             arg_split = arg.split("=")
             if len(arg_split) == 2:
                 MAKE_PLOT = arg_split[1]
             else:
                 pass
          else:
             print("Warning: your argument:", arg, "is ignored.")

   main(GET_DATA, READ_DATA, MAKE_PLOT)
