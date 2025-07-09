import numpy as np

population_by_age = [
    673380, 700940, 744162, 797030, 783616, 792542, 806423, 816018, 824410, 802173,
    797364, 771006, 771740, 760061, 779030, 769457, 792341, 793349, 792901, 817999,
    843316, 834907, 857842, 893523, 952560, 968463, 994345, 1020294, 1001001, 980890,
    988288, 1015189, 1033390, 1060509, 1151887, 1139320, 1165498, 1143025, 1124922, 1092066,
    1083525, 1083458, 1099294, 1086742, 1088001, 1032845, 1014045, 999487, 976705, 947515,
    956347, 957400, 1035320, 1138063, 1173661, 1251796, 1301334, 1324870, 1351436, 1344152,
    1360738, 1339742, 1294344, 1269790, 1214279, 1169242, 1091834, 1059429, 1024594, 979854,
    949236, 904559, 890040, 859165, 842694, 797432, 711029, 651950, 550371, 468806,
    607141, 593775, 553001, 641590, 633529, 3025803
]

# Convert the population list to a numpy array
population_array = np.array(population_by_age)

# Calculate the total population
total_population = np.sum(population_array)

# definte the age groups
age_groups = [
    "0-4 years",
    "5-14 years",
    "15-34 years",
    "35-59 years",
    "60-79 years",
    "80+ years"
]
# Calculate the population by age group
population_by_age_group = {
    "0-4 years": np.sum(population_array[0:5]),
    "5-14 years": np.sum(population_array[5:15]),
    "15-34 years": np.sum(population_array[15:35]),
    "35-59 years": np.sum(population_array[35:60]),
    "60-79 years": np.sum(population_array[60:80]),
    "80+ years": np.sum(population_array[80:])
}
# Print the total population
print(f"Total population: {total_population}")
# Print the population by age group
for age_group, population in population_by_age_group.items():
    print(
        f" {age_group}: {population} ({population / total_population:.2%})")
