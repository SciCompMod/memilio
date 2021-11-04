# TODO: to be either removed or substantially changed before merging...
import pandas as pd
import epidemiology.epidata.defaultDict as dd
import itertools
import numpy as np

# creates a mapping from given intervals to new desired intervals
def create_intervals_mapping(from_lower_bounds, to_lower_bounds):
    """! Creates a mapping from given intervals to new desired intervals

    @param from_lower_bounds lower bounds of original intervals
    @param to_lower_bounds desired lower bounds of new intervals
    @return mapping from intervals to intervals
    """
    # compute share of all_ages intervals from population intervals
    from_to_mapping = [[] for i in range(0, len(from_lower_bounds)-1)]
    j = 0  # iterator over all age breaks
    for i in range(0, len(from_lower_bounds)-1):
        # check if lower bound is larger in lower resolved data
        if from_lower_bounds[i] >= to_lower_bounds[j]:
            # Example: min_age_pop[i]=3 and min_age_pop[i+1]=6 shall be mapped on all_ages[j]=3 and all_ages[j+1]=5
            #   Then the all ages interval from j to j+1 will obtain the share
            #       x = (all_ages[j+1] - all_ages[j]) / (min_age_pop[i+1] - min_age_pop[i])
            #   of population age group 3-6
            #   in the next step, check if 6 is larger than all_ages[j+2]=Y
            #       if no: add the remaining part 1-x to j+1
            #       if yes: compute the corresponding share and go through it iteratively
            share = 0
            # if not, the remaining share will be assigned all ages j
            while from_lower_bounds[i+1] > to_lower_bounds[j+1]:
                share += (to_lower_bounds[j+1] - to_lower_bounds[j]
                          ) / (from_lower_bounds[i+1] - from_lower_bounds[i])
                from_to_mapping[i].append([share, j])
                j += 1
            from_to_mapping[i].append([1-share, j])
            # if both upper bounds are equal, then all ages j will not get any more share from any old group
            if from_lower_bounds[i+1] == to_lower_bounds[j+1]:
                j += 1

    return from_to_mapping

mobi = 'migration'
user = 'kueh_mj/'

df1 = pd.read_csv('/home/' + user + 'memilio/data/' + mobi + '/twitter_scaled_1252.txt', sep=' ', header=None)
print("Size twitter data " + str(len(df1)) + " x " + str(len(df1.columns)))
print("Entry correct at [398,4]: " + str(df1.iloc[398,-4] == 12.52))
print('partially empty lines: ' + str(len(df1[df1.isnull().any(axis=1)])))

df2 = pd.read_csv('/home/' + user + 'memilio/data/' + mobi + '/commuter_migration_scaled_2020.txt', sep=' ', header=None)
print("Size commuter data " + str(len(df2)) + " x " + str(len(df2.columns)))
print("Entry correct at [0,1]: " + str(df2.iloc[0,1] == 710))
print('partially empty lines: ' + str(len(df2[df2.isnull().any(axis=1)])))

df3 = pd.read_json('/home/' + user + 'memilio/data/pydata/Germany/county_divi_ma7.json')
print("Size DIVI infection data " + str(len(df3)) + " x " + str(len(df3.columns)) + ", division of length by 400: " + str(len(df3)/400))
print(df3.columns)
print('partially empty lines: ' + str(len(df3[df3.isnull().any(axis=1)])))

df4 = pd.read_json('/home/' + user + 'memilio/data/pydata/Germany/all_county_age_ma7_rki.json')
print("Size RKI infection data " + str(len(df4)) + " x " + str(len(df4.columns)) + ", division of length by 400: " + str(len(df4)/400))
print(df4.columns)
print('partially empty lines (with unknown columns): ' + str(len(df4[df4.isnull().any(axis=1)])))
df4a = df4[df4.isnull().any(axis=1)]
print('partially empty lines (without unknown columns): ' + str(len(df4a[df4a.Age_RKI!='unknown'])))

df5 = pd.read_json('/home/' + user + 'memilio/data/pydata/Germany/county_current_population.json')
print("Size population data " + str(len(df5)) + " x " + str(len(df5.columns)) + ", division of length by 400: " + str(len(df5)/400))
print(df5.columns)
print('partially empty lines: ' + str(len(df5[df5.isnull().any(axis=1)])))

df6 = pd.read_json('/home/' + user + 'memilio/data/pydata/Germany/all_county_agevacc_vacc_ma7.json')

print("Size RKI vaccination data " + str(len(df6)) + " x " + str(len(df6.columns)) + ", division of length by 400: " + str(len(df6)/400/3))
print(df6.columns)
print('partially empty lines: ' + str(len(df6[df6.isnull().any(axis=1)])))

unique_age_groups_old = sorted(df6[dd.EngEng['ageRKI']].unique())
ids  = sorted(df6['ID_County'].unique())

max_age_all = 100
# get age groups separators of original vaccination table
min_age_old = []
extrapolate_agegroups = True
for age in unique_age_groups_old:
    if '-' in age:
        min_age_old.append(int(age.split('-')[0]))
    elif '+' in age:
        min_age_old.append(int(age.split('+')[0]))
    else:
        extrapolate_agegroups = False
        print("Error in provided age groups from vaccination data; "
        "can not extrapolate to infection number age groups.")
min_age_old.append(max_age_all)

population = df5

min_age_pop = []
extrapolate_agegroups = True
unique_age_groups_pop = list(population.columns)[2:]
for age in unique_age_groups_pop:
    age = age.split()[0]  # remove " years" from string
    if '-' in age:
        min_age_pop.append(int(age.split('-')[0]))
    elif '>' in age:
        min_age_pop.append(int(age.split('>')[1]))
    elif '<' in age:
        min_age_pop.append(0)
    else:
        extrapolate_agegroups = False
        print("Error in provided age groups from population data;"
        " can not extrapolate to infection number age groups.")
min_age_pop.append(max_age_all)

# new age groups, here taken from definition of RKI infection data
min_age_new = [0, 5, 15, 35, 60, 80, max_age_all]

# combine all age group breaks
min_all_ages = sorted(pd.unique(list(itertools.chain(
    min_age_old, min_age_pop, min_age_new))))

# get number of new age groups that are not vaccinated at all
j = 0
new_age_not_vacc = 0
while min_age_new[j+1] < min_age_old[0]:
    new_age_not_vacc += 1
    j += 1

# compute share of all_ages intervals from population intervals
population_to_all_ages_share = create_intervals_mapping(
    min_age_pop, min_all_ages)

# compute mappings from all ages to old and new intervals
all_ages_to_age_old_share = create_intervals_mapping(
    min_all_ages, [0] + min_age_old)
all_ages_to_age_new_share = create_intervals_mapping(
    min_all_ages, min_age_new)

# compute indices of (partially) shared intervals from old to new
age_old_to_age_new_share = create_intervals_mapping(
    [0] + min_age_old, min_age_new)
age_old_to_age_new_share[0] = []
for i in range(1, len(age_old_to_age_new_share)):
    age_old_to_age_new_share[i] = [x[1] for x in age_old_to_age_new_share[i]]

# get interval indices from all age groups that correspond to old age group
age_old_to_all_ages_indices = [[] for zz in range(0, len(min_age_old)-1)]
for i in range(0, len(unique_age_groups_old)):
    for k in range(0, len(all_ages_to_age_old_share)):
        if all_ages_to_age_old_share[k][0][1] == i + new_age_not_vacc:
            age_old_to_all_ages_indices[i].append(k)
        elif k == len(all_ages_to_age_old_share) \
            or all_ages_to_age_old_share[k][0][1] == i + new_age_not_vacc + 1:
            break

# get interval indices from all age groups that correspond to new age group
age_new_to_all_ages_indices = [[] for zz in range(0, len(min_age_new)-1)]
for i in range(0, len(min_age_new)):
    for k in range(0, len(all_ages_to_age_new_share)):
        if all_ages_to_age_new_share[k][0][1] == i:
            age_new_to_all_ages_indices[i].append(k)
        elif k == len(all_ages_to_age_new_share) \
                or all_ages_to_age_new_share[k][0][1] == i + 1:
            break

# create new data frame and add zero to all new age group columns
population_all_ages = pd.DataFrame(population[dd.EngEng['idCounty']])
for i in min_all_ages:
    population_all_ages[str(i)] = 0

# iterate over all original age groups
for i in range(0, len(population_to_all_ages_share)):
    # iterate over intervals where population shares are assigned to
    for assign_share in population_to_all_ages_share[i]:
        # assign_share[0]: share / factor, assign_share[1]: column / age group
        population_all_ages[str(min_all_ages[assign_share[1]])
                                ] += assign_share[0] * population[unique_age_groups_pop[i]]

# rename last column and save total number per county
population_all_ages.rename(
    columns={str(min_all_ages[-1]): 'Total'}, inplace=True)
# remove last entry from all ages to call remaining columns
min_all_ages = min_all_ages[:-1]
population_all_ages['Total'] = population_all_ages[[
    str(i) for i in min_all_ages]].sum(axis=1)

# TODO: a similar functionality has to be implemented as unit test
if max(
    abs(
        population[unique_age_groups_pop].sum(axis=1) -
        population_all_ages[[str(i) for i in min_all_ages]].sum(
            axis=1))) > 1e-8:
    print("ERROR")

if True:
    population_old_ages = pd.DataFrame(population[dd.EngEng['idCounty']])
    for i in range(len(age_old_to_all_ages_indices)):
        # access columns + start_age_data since county_ID (and maybe other) 
        # is in first place
        start_age_data = list(population_all_ages.columns).index('0')
        population_old_ages[unique_age_groups_old[i]] = population_all_ages.iloc[:, np.array(
            age_old_to_all_ages_indices[i])+start_age_data].sum(axis=1)

'''for i in range(len(df6)):
    temp = population_old_ages[population_old_ages['ID_County']==df6['ID_County'].values[i]]
    ratio = df6[['Vacc_partially', 'Vacc_completed', 'Vacc_refreshed']].values[i,:] / population_old_ages[df6['Age_RKI'].values[i]].values[0]
    
    
    if ratio[0] >= 1 or ratio[1] >= 1 or ratio[2] >= 1:
        print(ratio)'''
i=0
for id in ids:
    for age in unique_age_groups_old:
        temp = df6.loc[(df6.ID_County==id) & (df6.Age_RKI==age), ['Vacc_completed']].values/population_old_ages.loc[population_old_ages.ID_County==id,age].values

        temp_where = np.where(temp >= 0.97)
        if len(temp_where[0]) > 0:
            i += 1
            print(id, age)#, ':', temp[temp_where])

print(i)