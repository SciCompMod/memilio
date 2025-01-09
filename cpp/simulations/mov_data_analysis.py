import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from datetime import datetime, date



# txt file structure is as follows
# timestep Nr.: 0
# x,y,z
# and so on
# timestep Nr.: 1
# x,y,z
# and so on

#where x is location id, y is location type, z is amount of agents at location

# read in data
def read_data(file):
    with open(file) as f:
        lines = f.readlines()
    f.close()
    return lines

# extract data from txt file
def extract_data(lines):
    data = []
    for line in lines:
        if 'timestep' in line:
            timestep = int(line.split(' ')[-1])
            data.append([])
        else:
            data[-1].append(line.strip().split(','))
    return data

# create dataframe from data

def create_df(data):
    # we want a dataframe, data has the following structure:
    # [[x,y,z],[x,y,z],...] for each timestep
    return pd.DataFrame([[int(x[0]), x[1], int(x[2]), timestep] for timestep, timestep_data in enumerate(data) for x in timestep_data], columns=['location_id', 'location_type', 'amount', 'timestep'])
    

def print_info(df):
    #sort by amount of persons concurrently at location (amount)
    df = df.sort_values('amount', ascending=False)
    print(df.head(10))
    print(df.tail(10))


# main function
def main():
    file = '/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/cpp/simulations/loc_data_occupancy.txt'
    lines = read_data(file)
    #first line is garbage
    data = extract_data(lines[1:])
    df = create_df(data)
    print_info(df)
    return

if __name__ == '__main__':
    main()



