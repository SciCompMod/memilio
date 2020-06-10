#This dictionary ensures that in case of calling the functions
# and of calling the console scripts the default values are the same

defaultDict = {
   'read_data':False,
   'make_plot':False,
   'out_form':'json',
   'out_folder':''
}

# The following dict EngEng makes sure that for all
# languages and sources the same names are used
# Rules for keys: start with small letter,
# one word, if several words start with capital letter
# do not use underscore
EngEng = {
   'gender': 'Gender',
   'confirmed': 'Confirmed',
   'confirmedTotal': 'Confirmed_total',
   'confirmedPcr': 'Confirmed_PCR',
   'confirmedAb': 'Confirmed_AB',
   'recovered': 'Recovered',
   'deaths': 'Deaths',
   'idState': 'ID_State',
   'state': 'State',
   'idCounty': 'ID_County',
   'county': 'County',
   'ageRKI': 'Age_RKI',
   'age5': 'Age5',
   'age10': 'Age',
   'unknown': 'unknown',
   'female' : 'female',
   'male' : 'male',
   'date' : 'Date',
   'hospitalized': 'Hospitalized',
   'intensive care unit': 'ICU',
   '80+': '80+',
   '90+': '90+',
   'both': 'both',
   'all': 'all',
   'occupied_ICU': 'occupied_ICU',
   'free_ICU': 'free_ICU',
}

GerEng = {
   'Geschlecht': EngEng['gender'],
   'AnzahlFall': EngEng['confirmed'],
   'AnzahlGenesen': EngEng['recovered'],
   'AnzahlTodesfall': EngEng['deaths'],
   'IdBundesland': EngEng['idState'],
   'Bundesland': EngEng['state'],
   'IdLandkreis': EngEng['idCounty'],
   'Landkreis': EngEng['county'],
   'Altersgruppe': EngEng['ageRKI'],
   'Altersgruppe2': EngEng['age5'],
   'unbekannt': EngEng['unknown'],
   'W' : EngEng['female'],
   'M' : EngEng['male'],

   'bundesland': EngEng['idState'],
   'betten_belegt': EngEng['occupied_ICU'],
   'betten_frei': EngEng['free_ICU'],
   'daten_stand': EngEng['date'],
   'kreis': EngEng['idCounty'],
   'gemeindeschluessel': EngEng['idCounty'],
}

EsEng = {'fecha': EngEng['date'],
         'rango_edad': EngEng['age10'],
         'sexo': EngEng['gender'],
         'casos_confirmados': EngEng['confirmed'],
         'fallecidos' : EngEng['deaths'],
         'ingresos_uci': EngEng['intensive care unit'],
         'Fecha': EngEng['date'],
         'Casos': EngEng['confirmedTotal'],
         'PCR+': EngEng['confirmedPcr'],
         'TestAc+': EngEng['confirmedAb'],
         'hospitalizados': EngEng['hospitalized'],
         'Hospitalizados': EngEng['hospitalized'],
         'cod_ine': EngEng['idState'],
         'CCAA': EngEng['state'],
         'UCI': EngEng['intensive care unit'],
         'Fallecidos' : EngEng['deaths'],
         'Recuperados' : EngEng['recovered'],
         'Andalucía' : 'Andalucia',
         'Castilla y León' : 'Castilla y Leon',
         'Cataluña': 'Cataluna',
         'País Vasco' : 'Pais Vasco',
         'Aragón' : 'Aragon',
         }

