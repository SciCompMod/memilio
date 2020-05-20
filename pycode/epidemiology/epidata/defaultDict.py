#This dictionary ensures that in case of calling the functions
# and of calling the console scripts the default values are the same

defaultDict = {
   'read_data':False,
   'make_plot':False,
   'out_form':'json',
   'out_folder':''
}

GerEng = {
   'Geschlecht': 'Gender',
   'AnzahlFall': 'Confirmed',
   'AnzahlGenesen': 'Recovered',
   'AnzahlTodesfall': 'Deaths',
   'IdBundesland': 'ID_State',
   'Bundesland': 'State',
   'IdLandkreis': 'Id_County',
   'Landkreis': 'County',
   'Altersgruppe': 'Age_RKI',
   'Altersgruppe2': 'Age5',
   'unbekannt': 'unknown',
   'W' : 'female',
   'M' : 'male',
}

EsEng = {'fecha': 'Date',
         'Fecha': 'Date',
         'rango_edad': 'Age10',
         'sexo': 'Gender',
         'casos_confirmados': 'Confirmed',
         'Casos': 'Confirmed_total',
         'PCR+': 'Confirmed_PCR',
         'TestAc+': 'Confirmed_AB',
         'hospitalizados': 'Hospitalized',
         'Hospitalizados': 'Hospitalized',
         'fallecidos' : 'Deaths',
         'cod_ine': 'ID_State',
         'CCAA': 'State',
         'UCI': 'ICU',
         'ingresos_uci': 'ICU',
         'Fallecidos' : 'Deaths',
         'Recuperados' : 'Recovered'
         }