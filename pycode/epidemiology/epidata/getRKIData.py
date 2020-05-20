import os
import sys
import json
import pandas
import matplotlib.pyplot as plt

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import outputDict as od
from epidemiology.epidata import defaultDict as dd



def get_rki_data(read_data=dd.defaultDict['read_data'],
                 make_plot=dd.defaultDict['make_plot'],
                 out_form=dd.defaultDict['out_form'],
                 out_folder=dd.defaultDict['out_folder']):

   filename = os.path.join(out_folder, "FullData.json")

   if(read_data):
      # if once dowloaded just read json file
      #filename = os.path.join(out_folder, "FullData.json")
      
      try:
         df = pandas.read_json(filename)
      except ValueError:
         exit_string = "Error: The file: " + filename + "does not exist. Call program without -r flag to get it."
         sys.exit(exit_string)

   else:

      # Supported data formats:
      load = { 
         'csv': gd.loadCsv,
         'geojson': gd.loadGeojson
       }

      # ArcGIS public data item ID:
      itemId = 'dd4580c810204019a7b8eb3e0b329dd6_0'

      # Get data:
      df = load['csv'](itemId)

      # convert "Datenstand" to real date:
      df.Datenstand = pandas.to_datetime( df.Datenstand, format='%d.%m.%Y, %H:%M Uhr').dt.tz_localize('Europe/Berlin')  

      # output data to not always download it
      df.to_json(filename)

   
   # Preperation for plotting/output:

   outForm = od.outForm[out_form][0]
   outFormEnd = od.outForm[out_form][1]
   outFormSpec = od.outForm[out_form][2]

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
   #gbNF_cs.to_json("infected.json")
   getattr(gbNF_cs, outForm)(os.path.join(out_folder , "infected" + outFormEnd), **outFormSpec)

   if(make_plot == True):
      # make plot
      gbNF_cs.plot( title = 'COVID-19 infections', grid = True, 
                               style = '-o' )
      plt.tight_layout()
      plt.show()

   # Dead over "Meldedatum":
   gbNT = df[df.NeuerTodesfall >= 0].groupby( dateToUse ).sum()
   gbNT_cs = gbNT.AnzahlTodesfall.cumsum()

   # output
   getattr(gbNF_cs, outForm)(os.path.join(out_folder, "deaths" + outFormEnd), **outFormSpec)

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
  
   #print(gbNFst)
   #print(gbNFst_cs)

   # output
   getattr(gbNFst_cs, outForm)(os.path.join(out_folder, "infected_state" + outFormEnd), **outFormSpec)
   
   # output nested json
   gbNFst_cs.groupby(['IdBundesland', 'Bundesland'], as_index=False) \
               .apply(lambda x: x[[dateToUse,'AnzahlFall']].to_dict('r')) \
               .reset_index().rename(columns={0:'Dates'})\
               .to_json(out_folder + "gbNF_state_nested.json", orient='records')


   # infected (incl recovered), deaths and recovered together 

   gbAllSt = dfF.groupby( ['IdBundesland', 'Bundesland', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllSt_cs = gbAllSt.groupby(level=1).cumsum().reset_index()

   #print(gbAllSt)
   #print(gbAllSt.reset_index()[ (gbAllSt.reset_index().IdBundesland==16) ])
 
   # output
   getattr(gbAllSt_cs, outForm)(os.path.join(out_folder, "all_state" + outFormEnd), **outFormSpec)

   ############# Data for counties all ages ######################

   # NeuerFall: Infected (incl. recovered) over "dateToUse" for every county ("Landkreis"):
   gbNFc = df[df.NeuerFall >= 0].groupby( ['IdLandkreis', 'Landkreis', dateToUse]).sum()

   gbNFc_cs = gbNFc.groupby(level=1).AnzahlFall.cumsum().reset_index()
   
   #print(gbNFc)
   #print(gbNFc_cs)

   # output
   getattr(gbNFc_cs, outForm)( os.path.join(out_folder, "infected_county" + outFormEnd), **outFormSpec)

   # infected (incl recovered), deaths and recovered together 

   gbAllC = dfF.groupby( ['IdLandkreis', 'Landkreis', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllC_cs = gbAllC.groupby(level=1).cumsum().reset_index()

   #print(gbAllC)
   #print(gbAllC_cs)

   # output
   getattr(gbAllC_cs, outForm)(out_folder + "all_county" + outFormEnd, **outFormSpec)
   

   ######### Data whole Germany different gender ##################

   # infected (incl recovered), deaths and recovered together 

   gbAllG = dfF.groupby( ['Geschlecht', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllG_cs = gbAllG.groupby(level=0).cumsum().reset_index()

   # print(gbAllG)
   # print(gbAllG_cs)

   # output
   getattr(gbAllG_cs, outForm)(os.path.join(out_folder, "all_gender" + outFormEnd), **outFormSpec)

   if(make_plot == True):
      dfF.groupby( 'Geschlecht' ) \
         .agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum}) \
         . plot( title = 'COVID-19 infections, deaths, recovered', grid = True,
                             kind = 'bar' )
      plt.tight_layout()
      plt.show()


   ############################# Gender and State ###################################################### 

   # infected (incl recovered), deaths and recovered together 

   gbAllGState = dfF.groupby( ['IdBundesland', 'Bundesland', 'Geschlecht', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllGState_cs = gbAllGState.groupby(level=[1,2]).cumsum().reset_index()

   #print(gbAllGState)
   #print(gbAllGState_cs)

   # output
   getattr(gbAllGState_cs, outForm)(os.path.join(out_folder, "all_state_gender" + outFormEnd), **outFormSpec)

   ############# Gender and County #####################

   gbAllGCounty = dfF.groupby( ['IdLandkreis', 'Landkreis', 'Geschlecht', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllGCounty_cs = gbAllGCounty.groupby(level=[1,2]).cumsum().reset_index()

   #print(gbAllGCounty)
   #print(gbAllGCounty_cs)

   # output
   getattr(gbAllGCounty_cs, outForm)(os.path.join(out_folder, "all_county_gender" + outFormEnd), **outFormSpec)
  
   ######### Data whole Germany different ages ####################

   # infected (incl recovered), deaths and recovered together 

   gbAllA = dfF.groupby( ['Altersgruppe', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllA_cs = gbAllA.groupby(level=0).cumsum().reset_index()

   #print(gbAllA)
   #print(gbAllA_cs)

   # output
   getattr(gbAllA_cs, outForm)(os.path.join(out_folder, "all_age" + outFormEnd), **outFormSpec)

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

   ############################# Age and State ###################################################### 

   # infected (incl recovered), deaths and recovered together 

   gbAllAState = dfF.groupby( ['IdBundesland', 'Bundesland', 'Altersgruppe', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllAState_cs = gbAllAState.groupby(level=[1,2]).cumsum().reset_index()

   #print(gbAllAState.reset_index()[ (gbAllAState.reset_index().IdBundesland == 1) &  (gbAllAState.reset_index().Altersgruppe == "A00-A04") ])
   #print(gbAllAState_cs)
   #print(gbAllAState_cs[gbAllAState_cs.Altersgruppe == "A00-A04"])
   #print(gbAllAState_cs[gbAllAState_cs.IdBundesland == 1])
   #print(gbAllAState.reset_index()[ (gbAllAState.reset_index().IdBundesland == 16) & (gbAllAState_cs.Meldedatum == "2020-03-24 00:00:00+01:00") ].agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum}))
   #print(gbAllAState_cs[ (gbAllAState_cs.IdBundesland == 1) & (gbAllAState_cs.Altersgruppe == "A00-A04") ])

   # output
   getattr(gbAllAState_cs, outForm)(os.path.join(out_folder, "all_state_age" + outFormEnd), **outFormSpec)

   ############# Age and County #####################

   gbAllACounty = dfF.groupby( ['IdLandkreis', 'Landkreis', 'Altersgruppe', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllACounty_cs = gbAllACounty.groupby(level=[1,2]).cumsum().reset_index()

   #print(gbAllACounty)
   #print(gbAllACounty_cs)

   # output
   getattr(gbAllACounty_cs, outForm)(os.path.join(out_folder,  "all_county_age" + outFormEnd), **outFormSpec)


def main():

   [read_data, make_plot, out_form, out_folder] = gd.cli('Download data from RKI')

   get_rki_data(read_data, make_plot, out_form, out_folder)

if __name__ == "__main__":

   main()
