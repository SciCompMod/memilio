import sys
from urllib.request import urlopen
import json
import pandas
import matplotlib.pyplot as plt

import getDataIntoPandasDataFrame as gd
import outputDict as od

def main(get_data, read_data, make_plot, out_form):

   if(get_data):

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
   #gbNF_cs.to_json("infected.json")
   getattr(gbNF_cs, od.outForm[out_form][0])("infected" + od.outForm[out_form][1], **od.outForm[out_form][2])

   if(make_plot == True):
      # make plot
      gbNF_cs.plot( title = 'COVID-19 infections', grid = True, 
                               style = '-o' )
      plt.tight_layout()
      plt.show()

   # Dead over "Meldedatum":
   gbNT = df[df.NeuerTodesfall >= 0].groupby( dateToUse ).sum()
   gbNT_cs = gbNT.AnzahlTodesfall.cumsum()

   #gbNT_cs.to_json("deaths.json")
   getattr(gbNT_cs, od.outForm[out_form][0])("deaths" + od.outForm[out_form][1], **od.outForm[out_form][2])

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

   # output json
   #gbNFst_cs.to_json("infected_state.json", orient='records')
   getattr(gbNFst_cs, od.outForm[out_form][0])("infected_state" + od.outForm[out_form][1], **od.outForm[out_form][2])

   # output nested json
   gbNFst_cs.groupby(['IdBundesland', 'Bundesland'], as_index=False) \
               .apply(lambda x: x[[dateToUse,'AnzahlFall']].to_dict('r')) \
               .reset_index().rename(columns={0:'Dates'})\
               .to_json("gbNF_state_nested.json", orient='records')


   # infected (incl recovered), deaths and recovered together 

   gbAllSt = dfF.groupby( ['IdBundesland', 'Bundesland', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllSt_cs = gbAllSt.groupby(level=1).cumsum().reset_index()

   #print(gbAllSt)
   #print(gbAllSt.reset_index()[ (gbAllSt.reset_index().IdBundesland==16) ])
 
   # output json
   #gbAllSt_cs.to_json("all_state.json", orient='records')
   getattr(gbAllSt_cs, od.outForm[out_form][0])("all_state" + od.outForm[out_form][1], **od.outForm[out_form][2])

   ############# Data for counties all ages ######################

   # NeuerFall: Infected (incl. recovered) over "dateToUse" for every county ("Landkreis"):
   gbNFc = df[df.NeuerFall >= 0].groupby( ['IdLandkreis', 'Landkreis', dateToUse]).sum()

   gbNFc_cs = gbNFc.groupby(level=1).AnzahlFall.cumsum().reset_index()
   
   #print(gbNFc)
   #print(gbNFc_cs)

   # output json
   #gbNFc_cs.to_json("infected_county.json", orient='records')
   getattr(gbNFc_cs, od.outForm[out_form][0])("infected_county" + od.outForm[out_form][1], **od.outForm[out_form][2])

   # infected (incl recovered), deaths and recovered together 

   gbAllC = dfF.groupby( ['IdLandkreis', 'Landkreis', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllC_cs = gbAllC.groupby(level=1).cumsum().reset_index()

   #print(gbAllC)
   #print(gbAllC_cs)

   # output json
   #gbAllC_cs.to_json("all_county.json", orient='records')
   getattr(gbAllC_cs, od.outForm[out_form][0])("all_county" + od.outForm[out_form][1], **od.outForm[out_form][2])
   

   ######### Data whole Germany different gender ##################

   # infected (incl recovered), deaths and recovered together 

   gbAllG = dfF.groupby( ['Geschlecht', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllG_cs = gbAllG.groupby(level=0).cumsum().reset_index()

   # print(gbAllG)
   # print(gbAllG_cs)

   # output json
   #gbAllG_cs.to_json("all_gender.json", orient='records') 
   getattr(gbAllG_cs, od.outForm[out_form][0])("all_gender" + od.outForm[out_form][1], **od.outForm[out_form][2])

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

   # output json
   #gbAllGState_cs.to_json("all_state_gender.json", orient='records')
   getattr(gbAllGState_cs, od.outForm[out_form][0])("all_state_gender" + od.outForm[out_form][1], **od.outForm[out_form][2])

   ############# Gender and County #####################

   gbAllGCounty = dfF.groupby( ['IdLandkreis', 'Landkreis', 'Geschlecht', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllGCounty_cs = gbAllGCounty.groupby(level=[1,2]).cumsum().reset_index()

   #print(gbAllGCounty)
   #print(gbAllGCounty_cs)

   # output json
   #gbAllGCounty_cs.to_json("all_county_gender.json", orient='records')
   getattr(gbAllGCounty_cs, od.outForm[out_form][0])("all_county_gender" + od.outForm[out_form][1], **od.outForm[out_form][2])
  
   ######### Data whole Germany different ages ####################

   # infected (incl recovered), deaths and recovered together 

   gbAllA = dfF.groupby( ['Altersgruppe', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllA_cs = gbAllA.groupby(level=0).cumsum().reset_index()

   #print(gbAllA)
   #print(gbAllA_cs)

   # output json
   #gbAllA_cs.to_json("all_age.json", orient='records')
   getattr(gbAllA_cs, od.outForm[out_form][0])("all_age" + od.outForm[out_form][1], **od.outForm[out_form][2])

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

   # output json
   #gbAllAState_cs.to_json("all_state_age.json", orient='records')
   getattr(gbAllAState_cs, od.outForm[out_form][0])("all_state_age" + od.outForm[out_form][1], **od.outForm[out_form][2])

   ############# Age and County #####################

   gbAllACounty = dfF.groupby( ['IdLandkreis', 'Landkreis', 'Altersgruppe', dateToUse]).agg({"AnzahlFall": sum, "AnzahlTodesfall": sum, "AnzahlGenesen": sum})
   gbAllACounty_cs = gbAllACounty.groupby(level=[1,2]).cumsum().reset_index()

   #print(gbAllACounty)
   #print(gbAllACounty_cs)

   # output json
   #gbAllACounty_cs.to_json("all_county_age.json", orient='records')
   getattr(gbAllACounty_cs, od.outForm[out_form][0])("all_county_age" + od.outForm[out_form][1], **od.outForm[out_form][2])


if __name__ == "__main__":

   GET_DATA = True
   READ_DATA = False
   MAKE_PLOT = True
   OUT_FORM = "json"

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

          elif "OUT_FORM" in arg:

             arg_split = arg.split("=")
             if len(arg_split) == 2:
                of = arg_split[1]
                if of in ["json", "hdf5"]:
                   OUT_FORM = of
                else:
                   print("Warning: your argument:", arg, "for OUT_FORM is ignored. It has to be either hdf5 or json [default]")
             else:
                 print("Warning: your argument:", arg, "is ignored. It has to be in the form as: OUT_FORM=hdf5")

          else:
             print("Warning: your argument:", arg, "is ignored.")

   main(GET_DATA, READ_DATA, MAKE_PLOT, OUT_FORM)
