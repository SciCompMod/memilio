"""
@file commuter_migration_bfa.py

@brief gets data related to county migration from "Bundesagentur fuer Arbeit"
"""
import collections
import numpy as np
from epidemiology.epidata import getDataIntoPandasDataFrame as gD


def get_data(setup_dict):
    """! Gets the matrix of commuter migration patterns and all additional helper variables.
    Gets all data generated in this file

    Keyword arguments:
    @ param setup_dict dictionary with necessary values:
        'counties': Dataframe read from file 'kreise_deu.xlsx' with population data
        'num_counties': number of counties in Dataframe counties,
        'num_govregions': number of government regions in Dataframe counties,
        'path': String with datapath where migration files can be found
        'abs_tol': tolerated undetected people
        'rel_tol': relative Tolerance to undetected people
    """

    (countykey_list, countypop_list, govkey_list) = get_key_and_population_lists(setup_dict)
    (countykey2numlist, govkey2numlist) = map_keys_to_numlists(setup_dict, countykey_list, govkey_list)
    (countykey2govkey, countykey2localnumlist, gov_county_table, state_gov_table) = assign_geographical_entities(
        setup_dict,
        countykey_list, govkey_list)
    mat_commuter_migration = get_matrix_commuter_migration_patterns(setup_dict, countypop_list, govkey_list,
                                                                    countykey2numlist,
                                                                    govkey2numlist,
                                                                    countykey2govkey, countykey2localnumlist,
                                                                    gov_county_table, state_gov_table)
    return (countykey_list, countypop_list, govkey_list, countykey2numlist, govkey2numlist, gov_county_table,
            countykey2govkey, countykey2localnumlist, state_gov_table, mat_commuter_migration)


def get_key_and_population_lists(setup_dict):
    """! Get list of regional key identifiers and population sizes for counties.

    """
    # get and store all regional (county) identifiers in a list; store county populations accordingly
    # get a list of governing regions
    countykey_list = []
    countypop_list = []
    govkey_list = []
    counties = setup_dict['counties']
    for i in range(0, counties.shape[0]):
        # regional county identifieres (5 numbers)
        if len(str(counties.iloc[i][0])) == 5 and counties.iloc[i][0].isdigit():
            countykey_list.append(counties.iloc[i][0])
            countypop_list.append(counties.iloc[i][5])
            # government region keys (2 or 3 numbers)
        elif i < counties.shape[0] - 1 and len(str(counties.iloc[i][0])) < len(str(counties.iloc[i + 1][0])):
            # workaround for old gov. regions and Saxony
            if (not str(counties.iloc[i][1]).startswith('früher') and not str(counties.iloc[i][1]).startswith(
                    'Direktion')):
                # only take those keys which have less numbers than the key in the next row
                if len(str(counties.iloc[i][0])) != 4 and len(str(counties.iloc[i + 1][0])) == 5:
                    # where string length is not 4 and next key has length four
                    # these rows correspond to 'local government' regions (except for BW, RP and Saxony)
                    govkey_list.append(counties.iloc[i][0])

                elif i < counties.shape[0] - 2:
                    if len(str(counties.iloc[i][0])) == 3 and len(str(counties.iloc[i + 2][0])) == 5:
                        # workaround for BW; 'government regions' are again divided but do not appear as such in
                        # documents of the Arbeitsagentur
                        govkey_list.append(counties.iloc[i][0])

                    if len(str(counties.iloc[i][0])) == 2 and len(str(counties.iloc[i + 2][0])) == 5:
                        # workaround for RP and Saxony;
                        if str(counties.iloc[i + 1][1]).startswith('früher'):
                            # workaround for RP; 'government regions' were dissolved
                            govkey_list.append(counties.iloc[i][0])
                        elif str(counties.iloc[i + 1][1]).startswith('Direktion'):
                            # workaround for Saxony; 'Direktionsbezirke' not referred in commuter migration
                            govkey_list.append(counties.iloc[i][0])

    if len(govkey_list) != setup_dict['num_govregions']:
        print('Error. Number of government regions wrong. Having', len(govkey_list), 'instead of',
              setup_dict['num_govregions'])

    return (countykey_list, countypop_list, govkey_list)


def verify_sorted(countykey_list):
    """ verify that read countykey_list is sorted
    @param countykey_list List of county regional keys
    """
    sum_check = 0
    countykey_list_unique = np.unique(np.array(countykey_list))
    for i in range(0, len(countykey_list)):
        sum_check = int(countykey_list_unique[i]) - int(countykey_list[i])
        if sum_check > 0:
            print('Error. Input list not sorted, population per county list had to be sorted accordingly.')


def map_keys_to_numlists(setup_dict, countykey_list=None, govkey_list=None):
    """! Creates hash maps from from county regional keys and keys of its governing regions to numbered lists.

    Keyword arguments:
    @param setup_dict dictionary with necessary values.
    @param countykey_list List of county regional keys.
    @param govkey_list List of governing regions regional keys.
    """
    # create a hashmap from sorted regional identifiers (01001 - ...) to 0 - num_counties
    if countykey_list is None or govkey_list is None:
        countykey_list, _, govkey_list = get_key_and_population_lists(setup_dict)
    verify_sorted(countykey_list)

    countykey2numlist = collections.OrderedDict()
    i = 0
    for index in countykey_list:
        countykey2numlist[index] = i
        i += 1

    if i != setup_dict['num_counties']:
        print("Error. Number of counties wrong.")

    # create a hash map from sorted gov keys to local list
    govkey2numlist = collections.OrderedDict()
    i = 0
    for index in govkey_list:
        govkey2numlist[index] = i
        i += 1

    if i != setup_dict['num_govregions']:
        print("Error. Number of governing regions wrong.")

    return countykey2numlist, govkey2numlist


def assign_geographical_entities(setup_dict, countykey_list=None, govkey_list=None):
    """! Assigns counties to governing regions based on key comparison and creates list of governing regions per state.

    Only works with sorted key lists.

    Keyword arguments:
    @param setup_dict dictionary with necessary values
    @param countykey_list List of county regional keys.
    @param govkey_list List of governing regions regional keys.
    """

    if govkey_list is None or countykey_list is None:
        countykey_list, _, govkey_list = get_key_and_population_lists(setup_dict)

    verify_sorted(countykey_list)

    # Create list of government regions with lists of counties that belong to them and list of states with government
    # regions that belong to them; only works with sorted lists of keys.
    gov_county_table = []

    gov_index = 0
    col_index = 0
    col_list = []

    for i in range(0, len(countykey_list)):

        # check for belonging to currently considered government region
        if str(countykey_list[i]).startswith(str(govkey_list[gov_index])):
            col_list.append(countykey_list[i])  # add county to current government region
            col_index += 1
        # go to next government region
        if i < len(countykey_list) - 1 and (not str(countykey_list[i + 1]).startswith(str(govkey_list[gov_index]))):
            # add government region to full table
            gov_county_table.append(col_list)
            col_list = []
            gov_index += 1
            col_index = 0
    # add last government region
    gov_county_table.append(col_list)

    if len(gov_county_table) != setup_dict['num_govregions']:
        print('Error. Number of government regions wrong.')

    # create a unique hash map from county key to its government region and
    # a global key to local (in gov region) key ordering
    countykey2govkey = collections.OrderedDict()
    countykey2localnumlist = collections.OrderedDict()
    for i in range(0, len(gov_county_table)):
        for j in range(0, len(gov_county_table[i])):
            countykey2govkey[gov_county_table[i][j]] = i
            countykey2localnumlist[gov_county_table[i][j]] = j

    # create government regions list per state
    state_gov_table = []

    state_id = 1
    state_govlist_loc = []
    for i in range(0, len(govkey_list)):

        if str(int(govkey_list[i])).startswith(str(state_id)):
            state_govlist_loc.append(govkey_list[i])

        if i + 1 < len(govkey_list) and not str(int(govkey_list[i + 1])).startswith(str(state_id)):
            state_id += 1
            state_gov_table.append(state_govlist_loc)
            state_govlist_loc = []
    # add last state's list
    state_gov_table.append(state_govlist_loc)

    return countykey2govkey, countykey2localnumlist, gov_county_table, state_gov_table


def get_matrix_commuter_migration_patterns(setup_dict, countypop_list=None, govkey_list=None, countykey2numlist=None,
                                           govkey2numlist=None, countykey2govkey=None, countykey2localnumlist=None,
                                           gov_county_table=None, state_gov_table=None):
    """! Computes matrix of commuter migration patterns.

    Keyword arguments:
    @param setup_dict dictionary with necessary values
    @param countypop_list List of populations per counties.
    @param govkey_list List of governing regions regional keys.
    @param countykey2numlist Hash map from county regional keys to numbered list.
    @param govkey2numlist Hash map from governing region regional keys to numbered list.
    @param countykey2govkey Hash map from county regional keys to governing region regional keys.
    @param countykey2localnumlist Hash map from county regional keys to local numbered list (per governing region).
    @param gov_county_table Table of county regional keys per governing region.
    @param state_gov_table Table of governing region regional keys per federal state.

    """

    if countypop_list is None or govkey_list is None:
        (countykey_list, countypop_list, govkey_list) = get_key_and_population_lists(setup_dict)
    if countykey2numlist is None or govkey2numlist is None:
        (countykey2numlist, govkey2numlist) = map_keys_to_numlists(setup_dict, countykey_list, govkey_list)
    if gov_county_table is None or countykey2govkey is None or countykey2localnumlist is None or \
            state_gov_table is None:
        (countykey2govkey, countykey2localnumlist, gov_county_table, state_gov_table) = assign_geographical_entities(
            setup_dict,
            countykey_list, govkey_list)

    mat_commuter_migration = np.zeros((setup_dict['num_counties'], setup_dict['num_counties']))

    # maxium errors (of people not detected)
    max_abs_err = 0
    max_rel_err = 0

    files = []
    for n in range(1, 10):
        files.append('krpend_0' + str(n) + "_0.xlsx")
    for n in range(10, 17):
        files.append('krpend_' + str(n) + "_0.xlsx")

    n = 0

    for item in files:
        # Using the 'Einpendler' sheet to correctly distribute summed values over counties of other gov. region
        commuter_migration_file = gD.loadExcel(targetFileName=item, apiUrl=setup_dict['path'], extension='',
                                               sheet_name=3)
        # pd.read_excel(os.path.join(setup_dict['path'], item), sheet_name=3)

        counties_done = []  # counties considered as 'migration from'
        # current_row = -1  # row of matrix that belongs to county migrated from
        current_col = -1  # column of matrix that belongs to county migrated to
        checksum = 0  # sum of county migration from, to be checked against sum in document

        for i in range(0, commuter_migration_file.shape[0]):
            if (len(str(commuter_migration_file.iloc[i][0])) == 5
                    and (commuter_migration_file.iloc[i][0]).isdigit()):
                checksum = 0
                # make zero'd list of counties explicitly migrated to from county considered
                # 'implicit' migration means 'migration to' which is summed in a larger regional entity and not given in
                # detail per county
                counties_migratedfrom = []
                for j in range(0, len(gov_county_table)):
                    counties_migratedfrom.append(np.zeros(len(gov_county_table[j])))

                counties_done.append(commuter_migration_file.iloc[i][0])
                current_col = countykey2numlist[commuter_migration_file.iloc[i][0]]
                curr_county_migratedto = commuter_migration_file.iloc[i][1]
                current_key = commuter_migration_file.iloc[i][0]
                # migration to itself excluded!
                counties_migratedfrom[countykey2govkey[current_key]][countykey2localnumlist[current_key]] = 1

            if not isinstance(commuter_migration_file.iloc[i][2], float):
                # removal of nan's, regional keys are stored as strings

                # check if entry is a digit
                if (commuter_migration_file.iloc[i][2]).isdigit():
                    # explicit migration from county to county
                    if len(str(commuter_migration_file.iloc[i][2])) == 5:
                        # check if entry refers to a specific county, then set matrix value
                        current_row = countykey2numlist[commuter_migration_file.iloc[i][2]]
                        val = commuter_migration_file.iloc[i][4]
                        mat_commuter_migration[current_row, current_col] = val
                        checksum += val
                        counties_migratedfrom[countykey2govkey[commuter_migration_file.iloc[i][2]]][
                            countykey2localnumlist[commuter_migration_file.iloc[i][2]]] = 1

                    # take summed values of other REMAINING counties of government region
                    # here, some counties of the region are stated explicitly and the rest is summed
                    elif (str(commuter_migration_file.iloc[i][3]) == 'Übrige Kreise (Regierungsbezirk)' and str(
                            commuter_migration_file.iloc[i][4]).isdigit()):
                        # remove trailing zeros (dummy key w/o zeros: dummy_key_wozeros)
                        dummy_key_wozeros = str(commuter_migration_file.iloc[i][2])
                        if len(dummy_key_wozeros) > 2 and dummy_key_wozeros[2] == '0':
                            dummy_key_wozeros = dummy_key_wozeros[0:2]

                            # sum population of all counties not explicitly migrated from
                            # of the current gov region migrated from
                        dummy_pop_sum = 0
                        for k in range(0, len(gov_county_table[govkey2numlist[dummy_key_wozeros]])):
                            if counties_migratedfrom[govkey2numlist[dummy_key_wozeros]][k] < 1:
                                # get identifier (0-401) for county key
                                globindex = countykey2numlist[gov_county_table[govkey2numlist[dummy_key_wozeros]][k]]
                                # sum up
                                dummy_pop_sum += countypop_list[globindex]

                        # distribute emigration relatively to county population where migration comes from
                        for k in range(0, len(gov_county_table[govkey2numlist[dummy_key_wozeros]])):
                            if counties_migratedfrom[govkey2numlist[dummy_key_wozeros]][k] < 1:
                                # get identifier (0-401) for county key
                                globindex = countykey2numlist[gov_county_table[govkey2numlist[dummy_key_wozeros]][k]]
                                counties_migratedfrom[govkey2numlist[dummy_key_wozeros]][k] = 1

                                # set value computed relatively to county size and effective migration
                                current_row = globindex
                                val = commuter_migration_file.iloc[i][4] * countypop_list[globindex] / dummy_pop_sum
                                checksum += val
                                mat_commuter_migration[current_row, current_col] = val

                    # take summed values of ALL counties of a government region
                    # here, no single county of the region is stated explicitly, all counties are summed together
                    elif (commuter_migration_file.iloc[i][2] in govkey_list and sum(
                            counties_migratedfrom[govkey2numlist[commuter_migration_file.iloc[i][2]]]) == 0):
                        # sum population of all counties not explicitly migrated to
                        # of the current gov region migrated to
                        dummy_pop_sum = 0
                        for k in range(0, len(gov_county_table[govkey2numlist[commuter_migration_file.iloc[i][2]]])):
                            if counties_migratedfrom[govkey2numlist[commuter_migration_file.iloc[i][2]]][k] < 1:
                                # get identifier (0-401) for county key
                                globindex = countykey2numlist[gov_county_table[govkey2numlist[
                                    commuter_migration_file.iloc[i][2]]][k]]
                                # sum up
                                dummy_pop_sum += countypop_list[globindex]

                        # distribute emigration relatively to county population where migration comes from
                        for k in range(0, len(gov_county_table[govkey2numlist[commuter_migration_file.iloc[i][2]]])):
                            if counties_migratedfrom[govkey2numlist[commuter_migration_file.iloc[i][2]]][k] < 1:
                                # get identifier (0-401) for county key
                                globindex = countykey2numlist[gov_county_table[govkey2numlist[
                                    commuter_migration_file.iloc[i][2]]][k]]
                                counties_migratedfrom[govkey2numlist[commuter_migration_file.iloc[i][2]]][k] = 1

                                # set value computed relatively to county size and effective migration
                                current_row = globindex
                                val = commuter_migration_file.iloc[i][4] * countypop_list[globindex] / dummy_pop_sum
                                checksum += val
                                mat_commuter_migration[current_row, current_col] = val

                    # take summed values of other REMAINING counties of a whole Bundesland
                    # here, some counties of the Bundesland are stated explicitly and the rest is summed
                    # the first or is for the case that the right first line of the incoming people directly
                    # addresses one
                    # the latter 'or's is used if no single county nor gov region of a federal state is stated
                    # explicitly
                    # although there are existent government regions in this federal state (i.e., the state itself is
                    # not considered a governement region according to gov_list)
                    elif ((str(commuter_migration_file.iloc[i][3]) == 'Übrige Regierungsbezirke (Bundesland)' and str(
                            commuter_migration_file.iloc[i][4]).isdigit())
                          or ((commuter_migration_file.iloc[i][2]).isdigit() and str(
                                commuter_migration_file.iloc[i - 1][2]).startswith('nan'))
                          or (len(str(commuter_migration_file.iloc[i][2])) == 2 and
                              abs(float(commuter_migration_file.iloc[i][2]) - float(
                                  commuter_migration_file.iloc[i - 1][2])) == 1)
                          or (len(str(commuter_migration_file.iloc[i][2])) == 2 and
                              abs(float(commuter_migration_file.iloc[i][2]) - float(
                                  commuter_migration_file.iloc[i - 1][2])) == 2)):

                        # auxiliary key of Bundesland (key translated to int starting at zero)
                        dummy_key = int(commuter_migration_file.iloc[i][2]) - 1

                        # sum population of all counties not explicitly migrated
                        # from the current gov region migrated from
                        dummy_pop_sum = 0
                        for j in range(0, len(state_gov_table[dummy_key])):
                            # over all government regions not explicitly stated
                            gov_index = govkey2numlist[state_gov_table[dummy_key][j]]
                            for k in range(0,
                                           len(gov_county_table[gov_index])):
                                # over all counties of the considered gov region
                                if counties_migratedfrom[gov_index][k] < 1:
                                    # get identifier (0-401) for county key
                                    globindex = countykey2numlist[gov_county_table[gov_index][k]]
                                    # sum up
                                    dummy_pop_sum += countypop_list[globindex]

                        # distribute emigration relatively to county population where migration comes from
                        for j in range(0, len(
                                state_gov_table[dummy_key])):  # over all government regions not explicitly stated
                            gov_index = govkey2numlist[state_gov_table[dummy_key][j]]
                            for k in range(0,
                                           len(gov_county_table[gov_index])):
                                # over all counties of the considered gov region
                                if counties_migratedfrom[gov_index][k] < 1:
                                    # get identifier (0-401) for county key
                                    globindex = countykey2numlist[gov_county_table[gov_index][k]]
                                    counties_migratedfrom[gov_index][k] = 1

                                    # set value computed relatively to county size and effective migration
                                    current_row = globindex
                                    val = commuter_migration_file.iloc[i][4] * countypop_list[globindex] / dummy_pop_sum
                                    checksum += val
                                    mat_commuter_migration[current_row, current_col] = val

            # sum of total migration 'from'
            if str(commuter_migration_file.iloc[i][3]) == 'Einpendler aus dem Bundesgebiet':
                abs_err = abs(checksum - commuter_migration_file.iloc[i][4])
                if abs_err > max_abs_err:
                    max_abs_err = abs_err
                if abs_err / checksum > max_rel_err:
                    max_rel_err = abs_err / checksum
                if abs_err < setup_dict['abs_tol'] and abs_err / checksum < setup_dict['rel_tol']:
                    checksum = 0
                else:
                    print('Error in calculations for county ', curr_county_migratedto,
                          '\nAccumulated values:', checksum,
                          ', correct sum:', commuter_migration_file.iloc[i][4])
                    print('Absolute error:', abs_err, ', relative error:', abs_err / checksum)

        n += 1
        print('Federal state read. Progress ', n, '/ 16')
    if n != 16:
        print('Error. Files missing.')

    print('Maximum absolute error:', max_abs_err)
    print('Maximum relative error:', max_rel_err)
    return mat_commuter_migration


def main():
    """! Main program entry."""

    num_counties = 401  # number of counties
    num_govregions = 34  # number of local governing regions
    abs_tol = 100  # maximum absolute error allowed per county migration
    rel_tol = 0.01  # maximum relative error allowed per county migration
    print('Attention: You have to use your VPN access, otherwise this file is not working.')
    path = 'http://hpcagainstcorona.sc.bs.dlr.de/data/migration/'
    counties = gD.loadExcel(targetFileName='kreise_deu', apiUrl=path, extension='.xlsx', sheet_name=1)

    setup_dict = {'num_counties': num_counties,
                  'num_govregions': num_govregions,
                  'abs_tol': abs_tol,
                  'rel_tol': rel_tol,
                  'path': path,
                  'counties': counties}
    (countykey_list, countypop_list, govkey_list, countykey2numlist, govkey2numlist, gov_county_table, countykey2govkey,
     countykey2localnumlist, state_gov_table, mat_commuter_migration) = get_data(setup_dict)


if __name__ == "__main__":
    main()
