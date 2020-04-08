import sys
import requests
import json
import matplotlib.pyplot as  plt

GET_DATA = True
PLOT_DATA = True

def get_rki_data():
#   response = requests.get("https://services7.arcgis.com/mOBPykOjAyBO2ZKk/arcgis/rest/services/RKI_COVID19/FeatureServer/0/query?where=1%3D1&outFields=*&returnGeometry=false&outSR=4326&f=json")

   params = {
      "where" : "1=1",
      "outFields" : "*",
      "outSR" : "4326",
      "f" : "json",
   }

   response = requests.get("https://services7.arcgis.com/mOBPykOjAyBO2ZKk/arcgis/rest/services/RKI_COVID19/FeatureServer/0/query", params = params)
   print(response.url)

   print("response.status_code: ", response.status_code)

   if response.status_code == 200:


      # make dictionary out of it
      data = response.json()
  
      # readable output:
      text = json.dumps(data, sort_keys=True, indent=4)
      #print(text)

      
      return data
  
   else:
      pass 
 

# get data to plot first

def plot_rki_data(data):


   features = data["features"]
   print(features[0]["attributes"]["Meldedatum"])
   # attributes = features["attributes"]
   
   # go through list of dictionaries
   for feature in features:
      meldedatum = feature["attributes"]["Meldedatum"]
      print(meldedatum)
      # TODO: collect data to make statistic

   


def main(get_data, plot_data):

   if(get_data):
      data = get_rki_data();

   if(plot_data):
      plot_rki_data(data);

   print(get_data, plot_data)

   

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
             else:
                 pass

   main(GET_DATA, PLOT_DATA)
