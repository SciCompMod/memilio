import os
import sys
import json
import pandas
import numpy as np
import matplotlib.pyplot as plt

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd



def get_rki_data(read_data=dd.defaultDict['read_data'],
                 make_plot=dd.defaultDict['make_plot'],
                 out_form=dd.defaultDict['out_form'],
                 out_folder=dd.defaultDict['out_folder']):

   directory = os.path.join(out_folder, 'Germany/')

   gd.check_dir(directory)
   filename = os.path.join(directory, "FullDataRKI.json")

   if(read_data):
      # if once dowloaded just read json file
      
      try:
         df = pandas.read_json(filename)
      except ValueError:
         exit_string = "Error: The file: " + filename + " does not exist. Call program without -r flag to get it."
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

   # generate Test file:
   # df.head(100).to_json(os.path.join(directory, "TestDataRKI.json"))

   # store dict values in parameter to not always call dict itself
   Altersgruppe2 = dd.GerEng['Altersgruppe2']
   Altersgruppe = dd.GerEng['Altersgruppe']
   Geschlecht = dd.GerEng['Geschlecht']
   AnzahlFall = dd.GerEng['AnzahlFall']
   AnzahlGenesen = dd.GerEng['AnzahlGenesen']
   AnzahlTodesfall = dd.GerEng['AnzahlTodesfall']
   IdBundesland = dd.GerEng['IdBundesland']
   Bundesland = dd.GerEng['Bundesland']
   IdLandkreis = dd.GerEng['IdLandkreis']
   Landkreis = dd.GerEng['Landkreis']

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
      (df[Altersgruppe2] == '0-4') & (df[Altersgruppe2] == '5-9'),
      (df[Altersgruppe2] == '10-14') & (df[Altersgruppe2] == '15-19'),
      (df[Altersgruppe2] == '20-24') & (df[Altersgruppe2] == '25-29'),
      (df[Altersgruppe2] == '30-34') & (df[Altersgruppe2] == '35-39'),
      (df[Altersgruppe2] == '40-44') & (df[Altersgruppe2] == '45-49'),
      (df[Altersgruppe2] == '50-54') & (df[Altersgruppe2] == '55-59'),
      (df[Altersgruppe2] == '60-64') & (df[Altersgruppe2] == '65-69'),
      (df[Altersgruppe2] == '70-74') & (df[Altersgruppe2] == '75-79'),
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

   dfF.loc[dfF.NeuerFall<0, [AnzahlFall]] = 0
   dfF.loc[dfF.NeuerTodesfall<0, [AnzahlTodesfall]] = 0
   dfF.loc[dfF.NeuGenesen<0, [AnzahlGenesen]] = 0

   print("Available columns:", df.columns)

   ######## Data for whole Germany all ages ##########

   # NeuerFall: Infected (incl. recovered) over "dateToUse":

   # make sum for one "dateToUse"
   # old way:
   # gbNF = df[df.NeuerFall >= 0].groupby( dateToUse ).sum()
   gbNF = df[df.NeuerFall >= 0].groupby(dateToUse).agg({AnzahlFall: sum})

   # make cumulative sum of "AnzahlFall" for "dateToUse"
   # old way:
   # gbNF_cs = gbNF.AnzahlFall.cumsum()
   gbNF_cs = gbNF.cumsum()

   # outout to json file
   gd.write_dataframe(gbNF_cs, directory, "infected_rki", out_form)


   if(make_plot == True):
      # make plot
      gbNF_cs.plot( title = 'COVID-19 infections', grid = True, 
                               style = '-o' )
      plt.tight_layout()
      plt.show()

   # Dead over Date:
   gbNT = df[df.NeuerTodesfall >= 0].groupby( dateToUse ).agg({AnzahlTodesfall: sum})
   gbNT_cs = gbNT.cumsum()

   # output
   gd.write_dataframe(gbNT_cs, directory, "deaths_rki", out_form)

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
   gbNFst = df[df.NeuerFall >= 0].groupby( [IdBundesland, Bundesland, dateToUse ])\
                                 .agg({AnzahlFall: sum})

   gbNFst_cs = gbNFst.groupby(level=1).cumsum().reset_index()
  
   # output
   gd.write_dataframe(gbNFst_cs, directory, "infected_state_rki", out_form)
   
   # output nested json
   # gbNFst_cs.groupby(['IdBundesland', 'Bundesland'], as_index=False) \
   #            .apply(lambda x: x[[dateToUse,'AnzahlFall']].to_dict('r')) \
   #            .reset_index().rename(columns={0:'Dates'})\
   #            .to_json(directory + "gbNF_state_nested.json", orient='records')


   # infected (incl recovered), deaths and recovered together 

   gbAllSt = dfF.groupby( [IdBundesland, Bundesland, dateToUse])\
                .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
   gbAllSt_cs = gbAllSt.groupby(level=1).cumsum().reset_index()

   # output
   gd.write_dataframe(gbAllSt_cs, directory, "all_state_rki", out_form)

   ############# Data for counties all ages ######################

   # NeuerFall: Infected (incl. recovered) over "dateToUse" for every county ("Landkreis"):
   gbNFc = df[df.NeuerFall >= 0].groupby([IdLandkreis, Landkreis, dateToUse])\
                                .agg({AnzahlFall: sum})

   gbNFc_cs = gbNFc.groupby(level=1).cumsum().reset_index()

   # output
   gd.write_dataframe(gbNFc_cs, directory, "infected_county_rki", out_form)

   # infected (incl recovered), deaths and recovered together 

   gbAllC = dfF.groupby( [IdLandkreis, Landkreis, dateToUse]).\
                agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
   gbAllC_cs = gbAllC.groupby(level=1).cumsum().reset_index()

   # output
   gd.write_dataframe(gbAllC_cs, directory, "all_county_rki", out_form)
   

   ######### Data whole Germany different gender ##################

   # infected (incl recovered), deaths and recovered together 

   gbAllG = dfF.groupby( [Geschlecht, dateToUse])\
               .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
   gbAllG_cs = gbAllG.groupby(level=0).cumsum().reset_index()

   # output
   gd.write_dataframe(gbAllG_cs, directory, "all_gender_rki", out_form)

   if(make_plot == True):
      dfF.groupby(Geschlecht ) \
         .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum}) \
         . plot( title = 'COVID-19 infections, deaths, recovered', grid = True,
                             kind = 'bar' )
      plt.tight_layout()
      plt.show()


   ############################# Gender and State ###################################################### 

   # infected (incl recovered), deaths and recovered together 

   gbAllGState = dfF.groupby( [IdBundesland, Bundesland, Geschlecht, dateToUse])\
                    .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
   gbAllGState_cs = gbAllGState.groupby(level=[1,2]).cumsum().reset_index()

   # output
   gd.write_dataframe(gbAllGState_cs, directory, "all_state_gender_rki", out_form)

   ############# Gender and County #####################

   gbAllGCounty = dfF.groupby( [IdLandkreis, Landkreis, Geschlecht, dateToUse])\
                     .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
   gbAllGCounty_cs = gbAllGCounty.groupby(level=[1,2]).cumsum().reset_index()

   # output
   gd.write_dataframe(gbAllGCounty_cs, directory, "all_county_gender_rki", out_form)
  
   ######### Data whole Germany different ages ####################

   # infected (incl recovered), deaths and recovered together 

   gbAllA = dfF.groupby( [Altersgruppe, dateToUse])\
            .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
   gbAllA_cs = gbAllA.groupby(level=0).cumsum().reset_index()

   # output
   gd.write_dataframe(gbAllA_cs, directory, "all_age_rki", out_form)

   if(make_plot == True):
      dfF.groupby( Altersgruppe ) \
         .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum}) \
         .plot( title = 'COVID-19 infections, deaths, recovered for diff ages', grid = True,
                             kind = 'bar' )
      plt.tight_layout()
      plt.show()


      # Dead by "Altersgruppe":
      gbNTAG = df[df.NeuerTodesfall >= 0].groupby( Altersgruppe ).sum()

      gbNTAG.AnzahlTodesfall.plot( title = 'COVID-19 deaths', grid = True, 
                             kind = 'bar' )
      plt.tight_layout()
      plt.show()

   ############################# Age and State ###################################################### 

   ##### Age_RKI #####

   # infected (incl recovered), deaths and recovered together 

   gbAllAgeState = dfF.groupby( [IdBundesland, Bundesland, Altersgruppe, dateToUse])\
                    .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
   gbAllAgeState_cs = gbAllAgeState.groupby(level=[1,2]).cumsum().reset_index()

   # output
   gd.write_dataframe(gbAllAgeState_cs, directory, "all_state_age_rki", out_form)

   ##### Age5 and Age10#####

   # infected (incl recovered), deaths and recovered together

   gbAllAgeState = dfF.groupby([IdBundesland, Bundesland, dd.GerEng['Altersgruppe2'], dateToUse]) \
      .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})

   gbAllAgeState_cs = gbAllAgeState.groupby(level=[1, 2]).cumsum().reset_index()

   # output
   gd.write_dataframe(gbAllAgeState_cs, directory, "all_state_age5_rki", out_form)

   ##### Age10 #####

   gbAllAgeState = dfF.groupby([IdBundesland, Bundesland, 'Age10', dateToUse]) \
      .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})

   gbAllAgeState_cs = gbAllAgeState.groupby(level=[1, 2]).cumsum().reset_index()

   # output
   gd.write_dataframe(gbAllAgeState_cs, directory, "all_state_age10_rki", out_form)


   ############# Age and County #####################

   gbAllAgeCounty = dfF.groupby( [IdLandkreis, Landkreis, Altersgruppe, dateToUse])\
                     .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
   gbAllAgeCounty_cs = gbAllAgeCounty.groupby(level=[1,2]).cumsum().reset_index()

   # output
   gd.write_dataframe(gbAllAgeCounty_cs, directory, "all_county_age_rki", out_form)

   #### age5 ####

   gbAllAgeCounty = dfF.groupby([IdLandkreis, Landkreis, Altersgruppe2, dateToUse]) \
      .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
   gbAllAgeCounty_cs = gbAllAgeCounty.groupby(level=[1, 2]).cumsum().reset_index()

   # output
   gd.write_dataframe(gbAllAgeCounty_cs, directory, "all_county_age5_rki", out_form)

   #### age10 ####

   gbAllAgeCounty = dfF.groupby( [IdLandkreis, Landkreis, 'Age10', dateToUse])\
                     .agg({AnzahlFall: sum, AnzahlTodesfall: sum, AnzahlGenesen: sum})
   gbAllAgeCounty_cs = gbAllAgeCounty.groupby(level=[1,2]).cumsum().reset_index()

   # output
   gd.write_dataframe(gbAllAgeCounty_cs, directory, "all_county_age10_rki", out_form)


def main():

   [read_data, make_plot, out_form, out_folder] = gd.cli('Download data from RKI')

   get_rki_data(read_data, make_plot, out_form, out_folder)

if __name__ == "__main__":

   main()
