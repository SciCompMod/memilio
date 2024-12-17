import requests 

header = {"Authorization": "Bearer"}

url = "http://sc-030247l:8000/"

# POST compartments

compartment_data = [
    {"name": "MildInfections",
     "description": "MildInfections",
     "tags": []
    },
    {"name": "Hospitalized",
     "description": "Hospitalized",
     "tags": []
    },
    {"name": "ICU",
     "description": "ICU",
     "tags": []
    },
    {"name": "Dead",
     "description": "Dead",
     "tags": []
    }]

# # POST all compartments
# for compartment in compartment_data:
#     post_response = requests.post(url + "compartments/", json=compartment, headers=header)

#     # post_response_json = post_response.json()
#     # print(post_response_json)

#     print(post_response.status_code)


# # GET and DELETE all compartments
# response = requests.get(url + "compartments/", headers=header)
# ids  = [compartment["id"] for compartment in response.json()]
# for id in ids: 
#     requests.delete(url + "compartments/" + id, headers = header)

group_data = [
    { "name": "age_0", "description": "0-4", "category": "age"},
    { "name": "age_1", "description": "5-14", "category": "age"},
    { "name": "age_2", "description": "15-34", "category": "age"},
    { "name": "age_3", "description": "35-59", "category": "age"},
    { "name": "age_4", "description": "60-79", "category": "age"},
    { "name": "age_5", "description": "80+", "category": "age"},
    { "name": "total", "description": "All ages", "category": "age"}
]

# for group in group_data:
#     post_response = requests.post(url + "groups/", json=group, headers=header)

intervention_data = [
    {
        "name": "School closure",
        "description": "School closure intervention",
        "tags": []
    },
    {
        "name": "Face masks & social distancing School",
        "description": "Face mask usage and social distancing measures applied at 1-25% in schools",
        "tags": []
    },
    {
        "name": "Face masks & social distancing Work",
        "description": "Face mask usage and social distancing measures applied at 1-25% in workplaces",
        "tags": []
    },
    {
        "name": "Face masks & social distancing Other",
        "description": "Face mask usage and social distancing measures applied at 1-25% in other settings",
        "tags": []
    },
    {
        "name": "Remote work",
        "description": "Implementation of remote work policies",
        "tags": []
    }
]

for intervention in intervention_data: 
    post_response = requests.post(url + "interventions/", json=intervention, headers=header)
    print(post_response.status_code)
