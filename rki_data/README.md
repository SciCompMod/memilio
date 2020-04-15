This script (without plotting instead with outputting) could be used as a cronjob to get data from rki and get every statistic we need
 
Call python3 getGeoJsonIntoPandasDataFrame_for_RKI.py

than data is downloaded and improved.
The full data set is stored in a json.
Different statistic is plotted.
One statistic is as an example also outputted as json (This could be easily done for every statistic we need and than plotted with javascript)

While running the code close one window to get the next one.


Call python3 getGeoJsonIntoPandasDataFrame_for_RKI.py PLOT_DATA=True

To use stored data to not download it all the time.

Needed python packages:
- pandas
- matplotlib
