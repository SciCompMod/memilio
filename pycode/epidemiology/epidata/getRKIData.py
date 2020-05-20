import os
import sys
import json
import pandas
import numpy as np
import matplotlib.pyplot as plt

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import outputDict as od
from epidemiology.epidata import defaultDict as dd



def get_rki_data(read_data=dd.defaultDict['read_data'],
                 make_plot=dd.defaultDict['make_plot'],
                 out_form=dd.defaultDict['out_form'],
                 out_folder=dd.defaultDict['out_folder']):

   directory = os.path.join(out_folder, 'Germany/')

   if not os.path.exists(directory):
      os.makedirs(directory)

   filename = os.path.join(directory, "FullDataRKI.json")

   if(read_data):
      # if once dowloaded just read json file
      
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

      # output data to not always download it
      df.to_json(filename)


   # translate column gender from Spanish to English and standardize
   df.loc[df.Geschlecht == 'unbekannt', ['Geschlecht']] = dd.GerEng['unbekannt']
   df.loc[df.Geschlecht == 'W', ['Geschlecht']] = dd.GerEng['W']
   df.loc[df.Geschlecht == 'M', ['Geschlecht']] = dd.GerEng['M']
   df.loc[df.Altersgruppe == 'unbekannt', ['Altersgruppe']] = dd.GerEng['unbekannt']
   df.loc[df.Altersgruppe2 == 'unbekannt', ['Altersgruppe2']] = dd.GerEng['unbekannt']

   # change names of columns
   df.rename(dd.GerEng, axis=1, inplace=True)

   # Add column 'Date' with Date= Refadtum if IstErkrankungsbeginn = 1 else take Meldedatum
   df['Date'] = np.where(df['IstErkrankungsbeginn'] == 1, df['Refdatum'], df['Meldedatum'])

   # Add new column with Age with range 10 as spain data
   conditions = [
      (df[dd.GerEng['Altersgruppe2']] == '0-4') & (df[dd.GerEng['Altersgruppe2']] == '5-9'),
      (df[dd.GerEng['Altersgruppe2']] == '10-14') & (df[dd.GerEng['Altersgruppe2']] == '15-19'),
      (df[dd.GerEng['Altersgruppe2']] == '20-24') & (df[dd.GerEng['Altersgruppe2']] == '25-29'),
      (df[dd.GerEng['Altersgruppe2']] == '30-34') & (df[dd.GerEng['Altersgruppe2']] == '35-39'),
      (df[dd.GerEng['Altersgruppe2']] == '40-44') & (df[dd.GerEng['Altersgruppe2']] == '45-49'),
      (df[dd.GerEng['Altersgruppe2']] == '50-54') & (df[dd.GerEng['Altersgruppe2']] == '55-59'),
      (df[dd.GerEng['Altersgruppe2']] == '60-64') & (df[dd.GerEng['Altersgruppe2']] == '65-69'),
      (df[dd.GerEng['Altersgruppe2']] == '70-74') & (df[dd.GerEng['Altersgruppe2']] == '75-79'),
   ]

   choices = ['0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70-79']
   df['Age10'] = np.select(conditions, choices, default=dd.GerEng['unbekannt'])
   # convert "Datenstand" to real date:
   df.Datenstand = pandas.to_datetime(df.Datenstand, format='%d.%m.%Y, %H:%M Uhr').dt.tz_localize('Europe/Berlin')

   # Correct Timestampes:
   for col in [ 'Meldedatum', 'Refdatum', 'Date' ]:
      df[col] = df[col].astype( 'datetime64[ns]' ).dt.tz_localize('Europe/Berlin')

   # Be careful "Refdatum" may act different to official describtion
   # on https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0,
   # sometimes big difference identified between "Refdatum" and "Meldedatum"
   # New possibility Date is either Refdatum or Meldedatum after column 'IstErkrankungsbeginn' has been added.

   dateToUse = 'Date'
   df.sort_values( dateToUse, inplace = True )

   # Manipulate data to get rid of conditions: df.NeuerFall >= 0, df.NeuerTodesfall >= 0, df.NeuGenesen >=0
   # There might be a better way

   dfF = df

   dfF.loc[dfF.NeuerFall<0, [dd.GerEng['AnzahlFall']]] = 0
   dfF.loc[dfF.NeuerTodesfall<0, [dd.GerEng['AnzahlTodesfall']]] = 0
   dfF.loc[dfF.NeuGenesen<0, [dd.GerEng['AnzahlGenesen']]] = 0

   print("Available columns:", df.columns)

   # Preperation for plotting/output:
   outForm = od.outForm[out_form][0]
   outFormEnd = od.outForm[out_form][1]
   outFormSpec = od.outForm[out_form][2]

   ######## Data for whole Germany all ages ##########

   # NeuerFall: Infected (incl. recovered) over "dateToUse":

   # make sum for one "dateToUse"
   # old way:
   # gbNF = df[df.NeuerFall >= 0].groupby( dateToUse ).sum()
   gbNF = df[df.NeuerFall >= 0].groupby(dateToUse).agg({dd.GerEng["AnzahlFall"]: sum})

   # make cumulative sum of "AnzahlFall" for "dateToUse"
   # old way
   # gbNF_cs = gbNF.AnzahlFall.cumsum()
   gbNF_cs = gbNF.cumsum()

   # outout to json file
   getattr(gbNF_cs, outForm)(os.path.join(directory , "infected_rki" + outFormEnd), **outFormSpec)

   if(make_plot == True):
      # make plot
      gbNF_cs.plot( title = 'COVID-19 infections', grid = True, 
                               style = '-o' )
      plt.tight_layout()
      plt.show()

   # Dead over "Meldedatum":
   gbNT = df[df.NeuerTodesfall >= 0].groupby( dateToUse ).agg({dd.GerEng["AnzahlTodesfall"]: sum})
   gbNT_cs = gbNT.cumsum()

   # output
   getattr(gbNF_cs, outForm)(os.path.join(directory, "deaths_rki" + outFormEnd), **outFormSpec)

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
   #gbNFst = df[df.NeuerFall >= 0].groupby( [IdBundesland','Bundesland', dateToUse]).AnzahlFall.sum()
   gbNFst = df[df.NeuerFall >= 0].groupby( [dd.GerEng['IdBundesland'], dd.GerEng['Bundesland'], dateToUse ])\
                                 .agg({dd.GerEng["AnzahlFall"]: sum})

   gbNFst_cs = gbNFst.groupby(level=1).cumsum().reset_index()
  
   # output
   getattr(gbNFst_cs, outForm)(os.path.join(directory, "infected_state_rki" + outFormEnd), **outFormSpec)
   
   # output nested json
   # gbNFst_cs.groupby(['IdBundesland', 'Bundesland'], as_index=False) \
   #            .apply(lambda x: x[[dateToUse,'AnzahlFall']].to_dict('r')) \
   #            .reset_index().rename(columns={0:'Dates'})\
   #            .to_json(directory + "gbNF_state_nested.json", orient='records')


   # infected (incl recovered), deaths and recovered together 

   gbAllSt = dfF.groupby( [dd.GerEng['IdBundesland'], dd.GerEng['Bundesland'], dateToUse])\
                .agg({dd.GerEng["AnzahlFall"]: sum, dd.GerEng["AnzahlTodesfall"]: sum, dd.GerEng["AnzahlGenesen"]: sum})
   gbAllSt_cs = gbAllSt.groupby(level=1).cumsum().reset_index()

   # output
   getattr(gbAllSt_cs, outForm)(os.path.join(directory, "all_state_rki" + outFormEnd), **outFormSpec)

   ############# Data for counties all ages ######################

   # NeuerFall: Infected (incl. recovered) over "dateToUse" for every county ("Landkreis"):
   gbNFc = df[df.NeuerFall >= 0].groupby([dd.GerEng['IdLandkreis'], dd.GerEng['Landkreis'], dateToUse])\
                                .agg({dd.GerEng["AnzahlFall"]: sum})

   gbNFc_cs = gbNFc.groupby(level=1).cumsum().reset_index()

   # output
   getattr(gbNFc_cs, outForm)( os.path.join(directory, "infected_county_rki" + outFormEnd), **outFormSpec)

   # infected (incl recovered), deaths and recovered together 

   gbAllC = dfF.groupby( [dd.GerEng['IdLandkreis'], dd.GerEng['Landkreis'], dateToUse]).\
                agg({dd.GerEng["AnzahlFall"]: sum, dd.GerEng["AnzahlTodesfall"]: sum, dd.GerEng["AnzahlGenesen"]: sum})
   gbAllC_cs = gbAllC.groupby(level=1).cumsum().reset_index()

   # output
   getattr(gbAllC_cs, outForm)(directory + "all_county_rki" + outFormEnd, **outFormSpec)
   

   ######### Data whole Germany different gender ##################

   # infected (incl recovered), deaths and recovered together 

   gbAllG = dfF.groupby( [dd.GerEng['Geschlecht'], dateToUse])\
               .agg({dd.GerEng["AnzahlFall"]: sum, dd.GerEng["AnzahlTodesfall"]: sum, dd.GerEng["AnzahlGenesen"]: sum})
   gbAllG_cs = gbAllG.groupby(level=0).cumsum().reset_index()

   # output
   getattr(gbAllG_cs, outForm)(os.path.join(directory, "all_gender_rki" + outFormEnd), **outFormSpec)

   if(make_plot == True):
      dfF.groupby(dd.GerEng['Geschlecht'] ) \
         .agg({dd.GerEng["AnzahlFall"]: sum, dd.GerEng["AnzahlTodesfall"]: sum, dd.GerEng["AnzahlGenesen"]: sum}) \
         . plot( title = 'COVID-19 infections, deaths, recovered', grid = True,
                             kind = 'bar' )
      plt.tight_layout()
      plt.show()


   ############################# Gender and State ###################################################### 

   # infected (incl recovered), deaths and recovered together 

   gbAllGState = dfF.groupby( [dd.GerEng['IdBundesland'], dd.GerEng['Bundesland'], dd.GerEng['Geschlecht'], dateToUse])\
                    .agg({dd.GerEng["AnzahlFall"]: sum, dd.GerEng["AnzahlTodesfall"]: sum, dd.GerEng["AnzahlGenesen"]: sum})
   gbAllGState_cs = gbAllGState.groupby(level=[1,2]).cumsum().reset_index()

   # output
   getattr(gbAllGState_cs, outForm)(os.path.join(directory, "all_state_gender_rki" + outFormEnd), **outFormSpec)

   ############# Gender and County #####################

   gbAllGCounty = dfF.groupby( [dd.GerEng['IdLandkreis'], dd.GerEng['Landkreis'], dd.GerEng['Geschlecht'], dateToUse])\
                     .agg({dd.GerEng["AnzahlFall"]: sum, dd.GerEng["AnzahlTodesfall"]: sum, dd.GerEng["AnzahlGenesen"]: sum})
   gbAllGCounty_cs = gbAllGCounty.groupby(level=[1,2]).cumsum().reset_index()

   # output
   getattr(gbAllGCounty_cs, outForm)(os.path.join(directory, "all_county_gender_rki" + outFormEnd), **outFormSpec)
  
   ######### Data whole Germany different ages ####################

   # infected (incl recovered), deaths and recovered together 

   gbAllA = dfF.groupby( [dd.GerEng['Altersgruppe'], dateToUse])\
            .agg({dd.GerEng["AnzahlFall"]: sum, dd.GerEng["AnzahlTodesfall"]: sum, dd.GerEng["AnzahlGenesen"]: sum})
   gbAllA_cs = gbAllA.groupby(level=0).cumsum().reset_index()

   # output
   getattr(gbAllA_cs, outForm)(os.path.join(directory, "all_age_rki" + outFormEnd), **outFormSpec)

   if(make_plot == True):
      dfF.groupby( dd.GerEng['Altersgruppe'] ) \
         .agg({dd.GerEng["AnzahlFall"]: sum, dd.GerEng["AnzahlTodesfall"]: sum, dd.GerEng["AnzahlGenesen"]: sum}) \
         .plot( title = 'COVID-19 infections, deaths, recovered for diff ages', grid = True,
                             kind = 'bar' )
      plt.tight_layout()
      plt.show()


      # Dead by "Altersgruppe":
      gbNTAG = df[df.NeuerTodesfall >= 0].groupby( dd.GerEng['Altersgruppe'] ).sum()

      gbNTAG.AnzahlTodesfall.plot( title = 'COVID-19 deaths', grid = True, 
                             kind = 'bar' )
      plt.tight_layout()
      plt.show()

   ############################# Age and State ###################################################### 

   ##### Age_RKI #####

   # infected (incl recovered), deaths and recovered together 

   gbAllAState = dfF.groupby( [dd.GerEng['IdBundesland'], dd.GerEng['Bundesland'], dd.GerEng['Altersgruppe'], dateToUse])\
                    .agg({dd.GerEng["AnzahlFall"]: sum, dd.GerEng["AnzahlTodesfall"]: sum, dd.GerEng["AnzahlGenesen"]: sum})
   gbAllAState_cs = gbAllAState.groupby(level=[1,2]).cumsum().reset_index()

   # output
   getattr(gbAllAState_cs, outForm)(os.path.join(directory, "all_state_age_rki" + outFormEnd), **outFormSpec)

   ##### Age5 and Age10#####

   # infected (incl recovered), deaths and recovered together

   gbAllAState = dfF.groupby([dd.GerEng['IdBundesland'], dd.GerEng['Bundesland'], dd.GerEng['Altersgruppe2'], dateToUse]) \
      .agg({dd.GerEng["AnzahlFall"]: sum, dd.GerEng["AnzahlTodesfall"]: sum, dd.GerEng["AnzahlGenesen"]: sum})

   gbAllAState_cs = gbAllAState.groupby(level=[1, 2]).cumsum().reset_index()

   # output
   getattr(gbAllAState_cs, outForm)(os.path.join(directory, "all_state_age5_rki" + outFormEnd), **outFormSpec)

   ##### Age10 #####

   gbAllAState = dfF.groupby([dd.GerEng['IdBundesland'], dd.GerEng['Bundesland'], 'Age10', dateToUse]) \
      .agg({dd.GerEng["AnzahlFall"]: sum, dd.GerEng["AnzahlTodesfall"]: sum, dd.GerEng["AnzahlGenesen"]: sum})

   gbAllAState_cs = gbAllAState.groupby(level=[1, 2]).cumsum().reset_index()

   # output
   getattr(gbAllAState_cs, outForm)(os.path.join(directory, "all_state_age10_rki" + outFormEnd), **outFormSpec)



   ############# Age and County #####################

   gbAllACounty = dfF.groupby( [dd.GerEng['IdLandkreis'], dd.GerEng['Landkreis'], dd.GerEng['Altersgruppe'], dateToUse])\
                     .agg({dd.GerEng["AnzahlFall"]: sum, dd.GerEng["AnzahlTodesfall"]: sum, dd.GerEng["AnzahlGenesen"]: sum})
   gbAllACounty_cs = gbAllACounty.groupby(level=[1,2]).cumsum().reset_index()

   # output
   getattr(gbAllACounty_cs, outForm)(os.path.join(directory,  "all_county_age_rki" + outFormEnd), **outFormSpec)


def main():

   [read_data, make_plot, out_form, out_folder] = gd.cli('Download data from RKI')

   get_rki_data(read_data, make_plot, out_form, out_folder)

if __name__ == "__main__":

   main()
