############################################################################################################
#                                                                                                          #
#        IMPORTANT NOTE: WHEN USING THIS DATA, WE HAVE TO CITE https://github.com/datadista/datasets       #
#                                                                                                          #
############################################################################################################


import sys
from urllib.request import urlopen
import json
import pandas
import matplotlib.pyplot as plt
import numpy as np

import outputDict as od

def loadCsv( githubUrl = 'https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/', 
             CSVfile  = 'nacional_covid19_rango_edad' ):
    """ Loads data in CSV format from github with Spanish Ministerio de Sanidad 
    and ISCIII (Instituto de Salud Carlos III) data. (pandas DataFrame)

    This routine loads github data sets in CSV format of the given public 
    url into a pandas DataFrame and returns the DataFrame. 

    Keyword arguments:
    githubUrl -- github url
    CSVfile -- file name

    """
    url = githubUrl + CSVfile + '.csv'
    #print(url)

    try:
        df = pandas.read_csv( url )
    except OSError as e:
        print("ERROR: URL " + url + " could not be opened.")
        df = pandas.DataFrame()
    
    
    return df


def main(get_data, read_data, make_plot, out_form):

   AgesJSONData = 'raw_spain_all_age.json'
   StatJSONData = 'raw_spain_all_state.json'

   if(get_data):
  
      # Get data:
      # https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/nacional_covid19_rango_edad.csv
      df_age = loadCsv(CSVfile  = 'nacional_covid19_rango_edad')
      
      if df_age.empty != True:
         # standardization of column titles from Spanish to English # the stupid character in front of 'fecha' is correct here.. There is a bug in the original file..
         df_age.rename({'﻿fecha': 'date', 'rango_edad': 'age', 'sexo': 'gender', 'casos_confirmados': 'confirmed', 
               'hospitalizados': 'hospitalized', 'ingresos_uci': 'icu', 'fallecidos' : 'deaths'}, axis=1, inplace=True) 

         print("Read Spanish age data from online. Available columns:", df_age.columns)

         # translate column gender from Spanish to English and standardize
         df_age.loc[df_age.gender=='ambos', ['gender']] = 'both'
         df_age.loc[df_age.gender=='mujeres', ['gender']] = 'female'
         df_age.loc[df_age.gender=='hombres', ['gender']] = 'male'
         df_age.loc[df_age.age=='80 y +', ['age']] = '80+'
         df_age.loc[df_age.age=='Total', ['age']] = 'all'

         # Correct Timestamps:
         df_age['date'] = df_age['date'].astype( 'datetime64[ns]' ).dt.tz_localize('Europe/Berlin')

         # output data to not always download it
         df_age.to_json(AgesJSONData)

      # Get data:
      # https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/ccaa_covid19_datos_isciii.csv
      df_state = loadCsv(CSVfile  = 'ccaa_covid19_datos_isciii')
      
      if df_state.empty != True:
         # standardization of column titles from Spanish to English
         df_state.rename({'Fecha': 'date', 'cod_ine': 'id_state', 'CCAA': 'state', 'Casos': 'confirmed_total', 'PCR+': 'confirmed_pcr', 
               'TestAc+': 'confirmed_anti', 'Hospitalizados': 'hospitalized', 'UCI': 'icu', 'Fallecidos' : 'deaths', 'Recuperados' : 'recovered'}, axis=1, inplace=True)

         print("Read Spanish states data from online. Available columns:", df_state.columns)

         
         # fill empty cells (nan values) with zero
         df_state.replace(np.nan, 0, inplace=True)


         # remove special characters
         df_state.loc[df_state.state=="Andalucía", ['state']] = "Andalucia"
         df_state.loc[df_state.state=="Castilla y León", ['state']] = "Castilla y Leon"
         df_state.loc[df_state.state=="Cataluña", ['state']] = "Cataluna"
         df_state.loc[df_state.state=="País Vasco", ['state']] = "Pais Vasco"
         df_state.loc[df_state.state=="Aragón", ['state']] = "Aragon"

         # Correct Timestamps:
         df_state['date'] = df_state['date'].astype( 'datetime64[ns]' ).dt.tz_localize('Europe/Berlin')

         # output data to not always download it
         df_state.to_json(StatJSONData)
         

   elif(read_data):
      # if once dowloaded just read json file

      #### ages' data
      df_age = pandas.read_json(AgesJSONData)

      print("Read from local. Available columns:", df_age.columns)

      #### states' data
      df_state = pandas.read_json(StatJSONData)

      print("Read from local. Available columns:", df_state.columns)



   # Preparation for plotting/output age specific data:

   # only consider men AND women (through information on gender away)
   df_age = df_age.loc[df_age.gender=='both']

   # write file for all age groups summed together
   df_agesum = df_age.loc[df_age.age=='all']
      # call to df_ageall.to_json("all_age.json", orient='records')
   getattr(df_agesum, od.outForm[out_form][0])("spain" + od.outForm[out_form][1], **od.outForm[out_form][2])


   # write file with information on all age groups separately
   # age_groups = ['0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+']
   df_agesep = df_age.loc[df_age.age!='all']
      # call to df_ageall.to_json("all_age.json", orient='records')
   getattr(df_agesep, od.outForm[out_form][0])("spain_all_age" + od.outForm[out_form][1], **od.outForm[out_form][2])



   # Preparation for plotting/output:

   PCR_ONLY = True # if pcr only is used 
   # dff = df_state['state'].unique()


   # if PCR only is used, the number of confirmed cases is the number of people being tested positive with PCR test
   # otherwise, the number of positive antibody tests is also taken into account
   if PCR_ONLY:
      df_state.loc[df_state.confirmed_total==0, ['confirmed_total']] = df_state.confirmed_pcr
   else:
      df_state.loc[df_state.confirmed_total==0, ['confirmed_total']] = df_state.confirmed_pcr + df_state.confirmed_anti

   states_array = df_state.state.unique()

   # output json
   # call df_state.to_json("spain_all_state.json", orient='records'), or to hdf5 alternatively
   getattr(df_state, od.outForm[out_form][0])("spain_all_state" + od.outForm[out_form][1], **od.outForm[out_form][2])




if __name__ == "__main__":

   GET_DATA = False
   READ_DATA = True
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
    
