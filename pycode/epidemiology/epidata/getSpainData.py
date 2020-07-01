############################################################################################################
#                                                                                                          #
#        IMPORTANT NOTE: WHEN USING THIS DATA, WE HAVE TO CITE https://github.com/datadista/datasets       #
#                                                                                                          #
#                                                                                                          #
#        DO NOT USE DATA FROM THE FOLLOWING REGIONS SINCE THE COLUMNS HOSPITALIZED AND ICU                 #
#        ARE NOT CORRECTLY SUMMED TO TOTAL NUMBERS ! THE SAME APPLIES TO ALL AGE DATA AT THE MOMENT !      #
#                                                                                                          #    
#               HOSPITALIZED                                   ICU                                         #
#               Castilla La Mancha (until 2020-04-11)          Castilla La Mancha (hasta 2020-04-12)       #
#               Comunidad Valenciana (hasta 2020-04-08)        Castilla y León (hasta 2020-04-17)          #
#               Madrid (hasta 2020-04-26)                      Comunidad Valenciana (hasta 2020-04-08)     #
#               Castilla y León (hasta 2020-04-06)             Galicia (hasta 2020-04-29)                  #
#               Madrid (hasta 2020-04-26)                                                                  #           
#                                                                                                          #
############################################################################################################


import os
import sys
import pandas
import numpy as np

from epidemiology.epidata import getDataIntoPandasDataFrame as gd
from epidemiology.epidata import defaultDict as dd


def get_spain_data(read_data=dd.defaultDict['read_data'],
                   make_plot=dd.defaultDict['make_plot'],
                   out_form=dd.defaultDict['out_form'],
                   out_folder=dd.defaultDict['out_folder']):

   directory = os.path.join(out_folder, 'Spain/')
   gd.check_dir(directory)


   ages_file = 'raw_spain_all_age'
   stat_file = 'raw_spain_all_state.json'

   if(read_data):

      AgesJSONData = os.path.join(directory, ages_file+'.json')
      StatJSONData = os.path.join(directory, stat_file+'.json')

      # if once dowloaded just read json file
      print("Read from local.")

      #### ages' data
      df_age = pandas.read_json(AgesJSONData)

      #### states' data
      df_state = pandas.read_json(StatJSONData)

   else:

      print("Read Spanish data from online.")

      try:
         # Get data:
         # https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/nacional_covid19_rango_edad.csv
         df_age = gd.loadCsv('nacional_covid19_rango_edad',
                           apiUrl='https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/')
      except:
         exit_string = "Files \'nacional_covid19_rango_edad\' are not available online. Check URL."
         sys.exit(exit_string)

      if df_age.empty != True:
         # output data to not always download it
         gd.write_dataframe(df_age, directory, ages_file, "json")

      try:
         # Get data:
         # https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/ccaa_covid19_datos_isciii.csv
         df_state = gd.loadCsv('ccaa_covid19_datos_isciii',
                              apiUrl='https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/')
      except:
         exit_string = "Files \'ccaa_covid19_datos_isciii\' are not available online. Check URL."
         sys.exit(exit_string)

      if df_state.empty != True:

         # output data to not always download it
         gd.write_dataframe(df_age, directory, stat_file, "json")
      else:
         exit_string = "Something went wrong, dataframe is empty!"
         sys.exit(exit_string)
         
   # Manipulate data
   if df_age.empty != True:
      # standardization of column titles from Spanish to English
      # the stupid character in front of 'fecha' is correct here. There is a bug in the original file.
      df_age.rename(dd.EsEng, axis=1, inplace=True)

      print("Available age columns:", df_age.columns)

      # translate column gender from Spanish to English and standardize
      gender = dd.EngEng['gender']
      age10 = dd.EngEng['age10']
      df_age.loc[df_age[gender] == 'ambos', [gender]] = dd.EngEng['both']
      df_age.loc[df_age[gender] == 'mujeres', [gender]] = dd.EngEng['female']
      df_age.loc[df_age[gender] == 'mujeres', [gender]] = dd.EngEng['female']
      df_age.loc[df_age[gender] == 'hombres', [gender]] = dd.EngEng['male']
      df_age.loc[df_age[age10] == '80 y +', [age10]] = dd.EngEng['80+']
      df_age.loc[df_age[age10] == '90 y +', [age10]] = dd.EngEng['90+']
      df_age.loc[df_age[age10] == 'Total', [age10]] = dd.EngEng['all']

      # Correct Timestamps:
      date = dd.EngEng['date']
      df_age[date] = df_age[date].astype('datetime64[ns]').dt.tz_localize('Europe/Berlin')

   if df_state.empty != True:
      # standardization of column titles from Spanish to English
      df_state.rename(dd.EsEng, axis=1, inplace=True)

      print("Available state columns:", df_state.columns)

      # fill empty cells (nan values) with zero
      df_state.replace(np.nan, 0, inplace=True)

      # remove special characters
      state = dd.EngEng['state']
      df_state.loc[df_state[state] == "Andalucía", [state]] = dd.EsEng["Andalucía"]
      df_state.loc[df_state[state] == "Castilla y León", [state]] = dd.EsEng["Castilla y León"]
      df_state.loc[df_state[state] == "Cataluña", [state]] = dd.EsEng["Cataluña"]
      df_state.loc[df_state[state] == "País Vasco", [state]] = dd.EsEng["País Vasco"]
      df_state.loc[df_state[state] == "Aragón", [state]] = dd.EsEng["Aragón"]

      # Correct Timestamps:
      df_state['Date'] = df_state['Date'].astype('datetime64[ns]').dt.tz_localize('Europe/Berlin')

   # only consider men AND women (through information on gender away)
   df_age = df_age.loc[ df_age[dd.EngEng["gender"]] == dd.EngEng['both'] ]

   # write file for all age groups summed together
   df_agesum = df_age.loc[df_age[dd.EngEng["age10"]] == dd.EngEng['all']]

   gd.write_dataframe(df_agesum, directory, "spain", out_form)

   # write file with information on all age groups separately
   # age_groups = ['0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+']
   df_agesep = df_age.loc[df_age[dd.EngEng["age10"]] != dd.EngEng['all']]

   gd.write_dataframe(df_agesep, directory, "spain_all_age", out_form)

   # Preparation for plotting/output:

   PCR_ONLY = False # if pcr only is used 
   # dff = df_state['state'].unique()

   # if PCR only is used, the number of confirmed cases is the number of people being tested positive with PCR test
   # otherwise, the number of positive antibody tests is also taken into account
   if PCR_ONLY:
      df_state.loc[df_state[dd.EngEng['confirmedTotal']] == 0, [dd.EngEng['confirmedTotal']]]\
         = df_state[ dd.EngEng["confirmedPcr"]]
   else:
      df_state.loc[df_state[dd.EngEng['confirmedTotal']] == 0, [dd.EngEng['confirmedTotal']]]\
         = df_state[ dd.EngEng["confirmedPcr"]] + df_state[dd.EngEng["confirmedAb"]]


   # output json
   gd.write_dataframe(df_state, directory, "spain_all_state", out_form)


def main():

    [read_data, make_plot, out_form, out_folder] = gd.cli('Download of spain data')
    get_spain_data(read_data, make_plot, out_form, out_folder)


if __name__ == "__main__":

    main()
