import json


with open('../../../data/pydata/Germany/county_current_population.json') as f :
    population_data = json.load(f)

adapted_age_groups = [(0, 4), (5, 14), (15, 34), (35, 59), (60, 79), (80,)]

adapted_population_data = []
for item in population_data:
    age_groups = []
    age_groups.append(int(item["<3 years"] + (2/3 * item["3-5 years"])))
    age_groups.append(int(1/3 * item["3-5 years"] + item["6-14 years"]))
    age_groups.append(int(item["15-17 years"] + item["18-24 years"] + item["25-29 years"] + (0.5 * item["30-39 years"])))
    age_groups.append(int(0.5 * item["30-39 years"] + item["40-49 years"] + (0.8* item["50-64 years"])))
    age_groups.append(int(0.2 * item["50-64 years"] + item["65-74 years"] + (0.43 * item[">74 years"])))
    age_groups.append(int(0.57 * item[">74 years"]))
    adapted_population_data.append({"ID_County":item["ID_County"], "Total":sum(age_groups), "0-4":age_groups[0], "5-14":age_groups[1], "15-34":age_groups[2], "35-59":age_groups[3], "60-79":age_groups[4], "80+":age_groups[5]})

json_string = json.dumps(adapted_population_data)

with open("adapted_county_current_population.json", "w") as outfile:
    outfile.write(json_string)
