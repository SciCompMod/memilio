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
      itemIdB = '5dc2fc92850241c3be3d704aa0945d9c_2'
      itemIdL = 'b2e6d8854d9744ca88144d30bef06a76_1'

      # Get data:
      dfB = load['csv'](itemIdB)
      dfL = load['csv'](itemIdL)

      # output data to not always download it
      dfB.to_json(r"FullDataB.json")
      dfL.to_json(r"FullDataL.json")

   elif(read_data):
      # if once dowloaded just read json file
      dfB = pandas.read_json("FullDataB.json")
      dfL = pandas.read_json("FullDataL.json")

   # Preperation for plotting/output:

   print("Available columns for states:", dfB.columns)
   print("Available columns for counties:", dfL.columns)  

   # Filter data for Bundesland/Landkreis and Einwohnerzahl (EWZ)

   dfB = dfB.rename(columns={'LAN_ew_GEN': 'Bundesland', 'LAN_ew_EWZ': 'EWZ'})
   dfBo = dfB[['Bundesland','EWZ']]
   # dfBo.to_json("PopulStates.json", orient='records')
   getattr(dfBo, od.outForm[out_form][0])("PopulStates" + od.outForm[out_form][1], **od.outForm[out_form][2])
  

   # TODO Counties are less as in RKI data. Compare it!
   dfL = dfL.rename(columns={'GEN': 'Landkreis'})
   dfLo = dfL[['Landkreis','EWZ']]
   # dfLo.to_json("PopulCounties.json", orient='records')
   getattr(dfLo, od.outForm[out_form][0])("PopulCounties" + od.outForm[out_form][1], **od.outForm[out_form][2])


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
