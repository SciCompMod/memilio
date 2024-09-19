import pandas as pd
import geopandas
import random
from memilio.epidata import geoModificationGermany as geoger
import os
import copy
from math import radians, cos, sin, asin, sqrt
import matplotlib.pyplot as plt


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    r = 6371
    return c * r


def plot_mean_distances():
    centers = os.path.join(os.path.dirname(__file__), 'centers_county.json')
    df = pd.read_json(centers).transpose()
    counties = df.index.tolist()
    state_to_county = geoger.get_stateid_to_countyids_map(
        merge_eisenach=True)
    county_to_state_map = geoger.get_countyid_to_stateid_map(
        merge_berlin=True)

    # state_names = geoger.get_state_names()
    state_names = {
        1: 'Schleswig-Holstein',
        2: 'Hamburg',
        3: 'Lower Saxony',
        4: 'Bremen',
        5: 'North Rhine-Westphalia',
        6: 'Hesse',
        7: 'Rhineland-Palatinate',
        8: 'Baden-Württemberg',
        9: 'Bavaria',
        10: 'Saarland',
        11: 'Berlin',
        12: 'Brandenburg',
        13: 'Mecklenburg-Western Pomerania',
        14: 'Saxony',
        15: 'Saxony-Anhalt',
        16: 'Thuringia',
    }

    state_names = [state_names[i] for i in range(1, 17)]

    # calculate mean distances for each county with counties in the same state
    mean_distances = {}
    for county in counties:
        coordinates_county = df.loc[county]
        state_id = county_to_state_map[county]

        # get all counties from the state
        state_counties = copy.deepcopy(state_to_county[state_id])

        # delete county from state counties
        state_counties.remove(county)

        # calculate distance to each county in the state
        distances = []
        for state_county in state_counties:
            coordinates_state_county = df.loc[state_county]

            # calculate distance between counties
            distance = haversine(coordinates_county[0], coordinates_county[1],
                                 coordinates_state_county[0], coordinates_state_county[1])
            distances.append(distance)

        # calculate mean distance
        mean_distance = sum(distances) / len(distances) if distances else 0

        # store mean distance
        mean_distances[county] = mean_distance

    # create DataFrame from mean distances
    mean_distances_df = pd.DataFrame.from_dict(
        mean_distances, orient='index', columns=['mean_distance'])

    # add column state_id
    mean_distances_df['state'] = mean_distances_df.index.map(
        county_to_state_map)

    # mean
    mean_distance_per_state = mean_distances_df.groupby('state')[
        'mean_distance'].mean()

    # weighted mean based on # counties per Bundesland
    count_per_state = mean_distances_df.groupby(
        'state')['mean_distance'].count()

    total_distance_per_state = mean_distance_per_state * count_per_state
    weighted_mean_across_states = total_distance_per_state.sum() / \
        count_per_state.sum()

    # print avg over all states
    print("Mean distance over all states: ", mean_distance_per_state.mean())
    print("Weighted mean distance over all states: ",
          weighted_mean_across_states)

    # plot
    plt.figure(figsize=(10, 8))
    mean_distance_per_state.plot(
        kind='barh', title='Mean Distance per State', xlabel='Mean Distance (km)', ylabel='State')
    plt.xlabel('Mean Distance (km)', fontsize=20)
    plt.ylabel('State', fontsize=20)
    plt.xticks(fontsize=15)
    plt.grid(axis='x')
    plt.yticks(ticks=range(len(state_names)), labels=state_names,
               fontsize=10)  # Reduce font size for better fit
    plt.tight_layout()
    plt.savefig(os.path.join(os.path.dirname(__file__),
                'mean_distance_per_state_horizontal.pdf'))


def plot_county_in_distances(radius):
    """
    Plot the number of counties within a given radius (in km) for each county
    """

    centers = os.path.join(os.path.dirname(__file__), 'centers_county.json')
    df = pd.read_json(centers).transpose()
    counties = df.index.tolist()
    state_to_county = geoger.get_stateid_to_countyids_map(
        merge_eisenach=True)
    county_to_state_map = geoger.get_countyid_to_stateid_map(
        merge_berlin=True)

    county_in_radius = {}
    for county in counties:
        coordinates_county = df.loc[county]
        state_id = county_to_state_map[county]

        # calculate distance to each county in the state
        tmp = []
        for c in counties:
            if c == county:
                continue
            coordinates_state_county = df.loc[c]

            # calculate distance between counties
            distance = haversine(coordinates_county[0], coordinates_county[1],
                                 coordinates_state_county[0], coordinates_state_county[1])
            if distance <= radius:
                tmp.append(c)

        county_in_radius[county] = tmp

    # write a txt file with the results, where each line is a county and the counties within the radius
    with open(os.path.join(os.path.dirname(__file__), f'counties_in_radius_{radius}.txt'), 'w') as f:
        for county, counties in county_in_radius.items():
            counties_str = [str(county)
                            for county in counties]  # Konvertiere int zu str
            f.write(f"{county}: {', '.join(counties_str)}\n")

    return county_in_radius


def plot_random_counties(radius):
    county_in_radius = plot_county_in_distances(radius)
    map_data = geopandas.read_file(
        os.path.join(
            os.getcwd(),
            'tools/vg2500_12-31.utm32s.shape/vg2500/VG2500_KRS.shp'))

    county_names_ids = geoger.get_county_names_and_ids()
    county_id_to_name = {item[1]: item[0] for item in county_names_ids}

    # Select 5 random counties from the county_in_radius dictionary
    random_counties = random.sample(list(county_in_radius.keys()), 10)

    # Loop over each selected county and create a plot
    for i, county_id in enumerate(random_counties):
        fig, ax = plt.subplots(figsize=(10, 10))

        # Plot all counties in a very light gray
        map_data.plot(ax=ax, color='#D3D3D3', linewidth=0.5)

        # Get the counties within the radius
        counties_within_radius = county_in_radius[county_id]

        # Find the corresponding geometries in the map data
        counties_to_highlight = map_data[map_data['ARS'].astype(
            int).isin(counties_within_radius)]

        selected_county_geom = map_data[map_data['ARS'].astype(
            int) == county_id]
        selected_county_geom.plot(ax=ax, color='blue', alpha=0.7)

        # Highlight these counties in a random color
        counties_to_highlight.plot(ax=ax, color=random.choice(
            ['red']), alpha=0.7)

        # Plot the boundaries of all counties on top
        map_data.boundary.plot(ax=ax, linewidth=1, color='black')

        # Set title and remove axis
        plt.title(
            f'Counties within a {radius} km Radius of County {county_id_to_name[county_id]}', fontsize=15)
        plt.axis('off')

        # Save the plot
        plt.savefig(os.path.join(os.path.dirname(
            __file__), f'random_counties_{radius}_{county_id}.pdf'))
        plt.close()


if __name__ == '__main__':
    # plot_mean_distances()
    # plot_county_in_distances(radius=73.5)
    plot_random_counties(radius=103)
