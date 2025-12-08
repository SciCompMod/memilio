import requests
import datetime
import json
from itertools import combinations
import numpy as np
import memilio.epidata.geoModificationGermany as gMG

# header = {'Authorization': "Bearer anythingAsPasswordIsFineCurrently"}

url = "https://zam10063.zam.kfa-juelich.de/api-dev/"


def delete_scenarios(headers):
    response = requests.get(url + "scenarios/", headers=headers)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(
            url + "scenarios/" + id, headers=headers)
        if (delete_response.status_code != 200):
            print(delete_response.status_code)
            print(delete_response.reason)

    print(requests.get(url + "scenarios/", headers=headers).json())


def delete_interventions(headers):
    response = requests.get(url + "interventions/templates/", headers=headers)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(
            url + "interventions/templates/" + id, headers=headers)
        if (delete_response.status_code != 200):
            print(delete_response.status_code)
            print(delete_response.reason)

    print(requests.get(url + "interventions/templates/", headers=headers).json())


def delete_models(headers):
    response = requests.get(url + "models/", headers=headers)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(
            url + "models/" + id, headers=headers)
        if (delete_response.status_code != 200):
            print(delete_response.status_code)
            print(delete_response.reason)

    print(requests.get(url + "models/", headers=headers).json())


def delete_parameters(headers):
    response = requests.get(url + "parameterdefinitions/", headers=headers)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(
            url + "parameterdefinitions/" + id, headers=headers)
        if (delete_response.status_code != 200):
            print(delete_response.status_code)
            print(delete_response.reason)

    print(requests.get(url + "parameterdefinitions/", headers=headers).json())


def delete_nodelists(headers):
    response = requests.get(url + "nodelists/", headers=headers)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(
            url + "nodelists/" + id, headers=headers)
        if (delete_response.status_code != 200):
            print(delete_response.status_code)
            print(delete_response.reason)
            print(delete_response.json())

    print(requests.get(url + "nodelists/", headers=headers).json())


def delete_nodes(headers):
    response = requests.get(url + "nodes/", headers=headers)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(url + "nodes/" + id, headers=headers)
        if (delete_response.status_code != 200):
            print(delete_response.status_code)
            print(delete_response.reason)

    print(requests.get(url + "nodes/", headers=headers).json())


def delete_compartments(headers):
    response = requests.get(url + "compartments/", headers=headers)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(
            url + "compartments/" + id, headers=headers)
        if (delete_response.status_code != 200):
            print(delete_response.status_code)
            print(delete_response.reason)

    print(requests.get(url + "compartments/", headers=headers).json())


def delete_groups(headers):
    response = requests.get(url + "groups/", headers=headers)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(
            url + "groups/" + id, headers=headers)
        if (delete_response.status_code != 200):
            print(delete_response.status_code)
            print(delete_response.reason)

    print(requests.get(url + "groups/", headers=headers).json())


def delete_everything_from_db(headers):
    delete_scenarios(headers)
    delete_interventions(headers)
    delete_models(headers)
    delete_parameters(headers)
    delete_nodelists(headers)
    delete_nodes(headers)
    delete_compartments(headers)
    delete_groups(headers)


def post_to_db_compartments(headers):
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

    # POST all compartments
    for compartment in compartment_data:
        post_response = requests.post(
            url + "compartments/", json=compartment, headers=headers)
        if (post_response.status_code != 200):
            print(post_response.status_code)


def post_to_db_groups(headers):
    group_data = [
        {"name": "Group1", "description": "0-4", "category": "age"},  # 0 bis 4
        {"name": "Group2", "description": "5-14", "category": "age"},  # 5-14
        {"name": "Group3", "description": "15-34", "category": "age"},  # 15-34
        {"name": "Group4", "description": "35-59", "category": "age"},  # 35-59
        {"name": "Group5", "description": "60-79", "category": "age"},  # 60 - 79
        {"name": "Group6", "description": "80+", "category": "age"},  # 80+
        {"name": "Total", "description": "Total population", "category": "age"}
    ]

    for group in group_data:
        post_response = requests.post(
            url + "groups/", json=group, headers=headers)
        if (post_response.status_code != 200):
            print(post_response.status_code)


def post_to_db_interventions(headers):
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
        post_response = requests.post(
            url + "interventions/templates/", json=intervention, headers=headers)
        if (post_response.status_code != 200):
            print(post_response.status_code)


def post_to_db_parameters(headers):
    parameter_data = [
        {"name": "TimeExposed",
         "description": "TimeExposed"},
        {"name": "TimeInfectedNoSymptoms",
         "description": "TimeInfectedNoSymptoms"},
        {"name": "TimeInfectedSymptoms",
         "description": "TimeInfectedSymptoms"},
        {"name": "TimeInfectedSevere",
         "description": "TimeInfectedSevere"},
        {"name": "TimeInfectedCritical",
         "description": "TimeInfectedCritical"},
        {"name": "TransmissionProbabilityOnContact",
         "description": "TransmissionProbabilityOnContact"},
        {"name": "RelativeTransmissionNoSymptoms",
         "description": "RelativeTransmissionNoSymptoms"},
        {"name": "RiskOfInfectionFromSymptomatic",
         "description": "RiskOfInfectionFromSymptomatic"},
        {"name": "MaxRiskOfInfectionFromSymptomatic",
         "description": "MaxRiskOfInfectionFromSymptomatic"},
        {"name": "RecoveredPerInfectedNoSymptoms",
         "description": "RecoveredPerInfectedNoSymptoms"},
        {"name": "SeverePerInfectedSymptoms",
         "description": "SeverePerInfectedSymptoms"},
        {"name": "CriticalPerSevere",
         "description": "CriticalPerSevere"},
        {"name": "DeathsPerCritical",
         "description": "DeathsPerCritical"},
        {"name": "ReducedExposedPartialImmunity",
         "description": "ReducedExposedPartialImmunity"},
        {"name": "ReducedExposedImprovedImmunity",
         "description": "ReducedExposedImprovedImmunity"},
        {"name": "ReducedInfectedSymptomsPartialImmunity",
         "description": "ReducedInfectedSymptomsPartialImmunity"},
        {"name": "ReducedInfectedSymptomsImprovedImmunity",
         "description": "ReducedInfectedSymptomsImprovedImmunity"},
        {"name": "ReducedInfectedSevereCriticalDeadPartialImmunity",
         "description": "ReducedInfectedSevereCriticalDeadPartialImmunity"},
        {"name": "ReducedInfectedSevereCriticalDeadImprovedImmunity",
         "description": "ReducedInfectedSevereCriticalDeadImprovedImmunity"},
        {"name": "ReducedTimeInfectedMild",
         "description": "ReducedTimeInfectedMild"},
        {"name": "Seasonality",
         "description": "Seasonality"}
    ]

    for parameter in parameter_data:
        post_response = requests.post(
            url + "parameterdefinitions/", json=parameter, headers=headers)
        if (post_response.status_code != 200):
            print(post_response.status_code)


def post_to_db_model(headers):
    get_compartments = requests.get(url + "compartments/", headers=headers)
    compartments = [compartment["id"]
                    for compartment in get_compartments.json()]

    get_groups = requests.get(url + "groups/", headers=headers)
    groups = [group["id"] for group in get_groups.json()]

    get_parameters = requests.get(
        url + "parameterdefinitions/", headers=headers)
    parameter_ids = [parameter["id"] for parameter in get_parameters.json()]

    model_data = {
        "name": "secirvvs",
        "description": "ODE-SECIRVVS model",
        "compartments": compartments,
        "groups": groups,
        "parameterDefinitions": parameter_ids
    }

    post_response = requests.post(
        url + "models/", json=model_data, headers=headers)
    if (post_response.status_code != 200):
        print(post_response.status_code)


def post_to_db_nodes(headers):
    nodes_init = gMG.get_county_ids()
    nodes_init.insert(0, '00000')
    # every county ID should have 5 digits
    nodes_init = [str(countyID).zfill(5) for countyID in nodes_init]
    # transform nodes_init to strings
    nodes_init = [str(countyID) for countyID in nodes_init]

    node_data = [{"nuts": f"{node}",
                  "name": f"{node}"} for node in nodes_init]

    for node in node_data:
        post_response = requests.post(
            url + "nodes/", json=node, headers=headers)
        if (post_response.status_code != 200):
            print(post_response.status_code)


def post_to_db_nodelist(headers):
    get_nodes = requests.get(url + "nodes/", headers=headers)
    nodes = [node["id"] for node in get_nodes.json()]

    nodelist_data = {
        "name": "all_counties",
        "description": "Contains all counties of Germany as well as a node for the whole of Germany.",
        "nodeIds": nodes
    }

    post_response = requests.post(
        url + "nodelists/", json=nodelist_data, headers=headers)
    if (post_response.status_code != 200):
        print(post_response.status_code)


def post_to_db_scenarios(headers, start_date_scenarios, days_simulated_scenarios=30, days_computed_casedata=10, modelparameters_entry={}, post=True):
    # Define start and end date for casedata scenario. Here, we assume that last day of the casedata scenario
    # corresponds to the first day of the simulated scenarios.
    start_date_casedata = (datetime.datetime.strptime(start_date_scenarios, "%Y-%m-%d") -
                           datetime.timedelta(days=days_computed_casedata-2)).strftime("%Y-%m-%d")
    end_date_casedata = start_date_scenarios
    # Define start and end date of simulation.
    start_date_simulation = start_date_scenarios
    end_date_simulation = (datetime.datetime.strptime(start_date_scenarios, "%Y-%m-%d") +
                           datetime.timedelta(days=days_simulated_scenarios)).strftime("%Y-%m-%d")

    # Get ids of model, nodelist and interventions.
    get_models = requests.get(url + "models/", headers=headers)
    model_id = [model["id"]
                for model in get_models.json() if model["name"] == "secirvvs"]

    get_nodelists = requests.get(url+"nodelists/", headers=headers)
    nodelist_id = [nodelist["id"] for nodelist in get_nodelists.json(
    ) if nodelist["name"] == "all_counties"]

    # Get sorted list of groups
    # TODO: Do this more carefully.
    get_groups = requests.get(url + "groups/", headers=headers)
    group_ids = []
    for i, group in enumerate(get_groups.json()):
        if group["name"] == f"Group{i+1}":
            group_ids.append(group["id"])
        elif group["name"] == f"Total":
            group_ids.append(group["id"])

    if (modelparameters_entry == {} and post):

        variantFactor = 1.4
        # Define dict that contains name of parameter and min and max values for all 6 age groups.
        parameter_dict_default = {
            "TimeExposed": [[2.67, 2.67, 2.67, 2.67, 2.67, 2.67], [4., 4., 4., 4., 4., 4.]],
            "TimeInfectedNoSymptoms": [[1.2, 1.2, 1.2, 1.2, 1.2, 1.2], [2.53, 2.53, 2.53, 2.53, 2.53, 2.53]],
            "TimeInfectedSymptoms": [[5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465], [8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085]],
            "TimeInfectedSevere": [[
                3.925, 3.925, 4.85, 6.4, 7.2, 9.], [6.075, 6.075, 7., 8.7, 9.8, 13.]],
            "TimeInfectedCritical": [[4.95, 4.95, 4.86, 14.14, 14.4, 10.], [8.95, 8.95, 8.86, 20.58, 19.8, 13.2]],
            "TransmissionProbabilityOnContact": [[
                0.02 * variantFactor, 0.05 * variantFactor, 0.05 * variantFactor,
                0.05 * variantFactor, 0.08 * variantFactor, 0.1 * variantFactor
            ], [
                0.04 * variantFactor, 0.07 * variantFactor, 0.07 * variantFactor,
                0.07 * variantFactor, 0.10 * variantFactor, 0.15 * variantFactor
            ]],
            "RelativeTransmissionNoSymptoms": [[0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]],
            "RiskOfInfectionFromSymptomatic": [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.2, 0.2, 0.2, 0.2, 0.2, 0.2]],
            "MaxRiskOfInfectionFromSymptomatic": [[0.4, 0.4, 0.4, 0.4, 0.4, 0.4], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]],
            "RecoveredPerInfectedNoSymptoms": [[
                0.2, 0.2, 0.15, 0.15, 0.15, 0.15], [
                0.3, 0.3, 0.25, 0.25, 0.25, 0.25]],
            "SeverePerInfectedSymptoms": [[
                0.006, 0.006, 0.015, 0.049, 0.15, 0.20],  [
                0.009, 0.009, 0.023, 0.074, 0.18, 0.25]],
            "CriticalPerSevere": [[0.05, 0.05, 0.05, 0.10, 0.25, 0.35], [0.10, 0.10, 0.10, 0.20, 0.35, 0.45]],
            "DeathsPerCritical": [[0.00, 0.00, 0.10, 0.10, 0.30, 0.5], [0.10, 0.10, 0.18, 0.18, 0.50, 0.7]],
            "ReducedExposedPartialImmunity": [[0.75, 0.75, 0.75, 0.75, 0.75, 0.75], [0.85, 0.85, 0.85, 0.85, 0.85, 0.85]],
            "ReducedExposedImprovedImmunity": [[0.281, 0.281, 0.281, 0.281, 0.281, 0.281], [0.381, 0.381, 0.381, 0.381, 0.381, 0.381]],
            "ReducedInfectedSymptomsPartialImmunity": [[0.6, 0.6, 0.6, 0.6, 0.6, 0.6], [0.7, 0.7, 0.7, 0.7, 0.7, 0.7]],
            "ReducedInfectedSymptomsImprovedImmunity": [[0.193, 0.193, 0.193, 0.193, 0.193, 0.193], [0.293, 0.293, 0.293, 0.293, 0.293, 0.293]],
            "ReducedInfectedSevereCriticalDeadPartialImmunity": [[0.05, 0.05, 0.05, 0.05, 0.05, 0.05], [0.15, 0.15, 0.15, 0.15, 0.15, 0.15]],
            "ReducedInfectedSevereCriticalDeadImprovedImmunity": [[0.041, 0.041, 0.041, 0.041, 0.041, 0.041], [0.141, 0.141, 0.141, 0.141, 0.141, 0.141]],
            "ReducedTimeInfectedMild": [[1., 1., 1., 1., 1., 1.], [1., 1., 1., 1., 1., 1.]],
            "Seasonality": [[0.1, 0.1, 0.1, 0.1, 0.1, 0.1], [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]]}

        # Append parameter value for age group "total". This is done by weighting the other parameter values according to
        # their share of the total population of Germany (using data from regionalstatistik).

        append_parameter_for_agegroup_total(parameter_dict_default)
        get_parameters = requests.get(
            url + "parameterdefinitions/", headers=headers)

        modelparameters_entry = []
        for parameter in get_parameters.json():
            value_entry = []
            for i, group_id in enumerate(group_ids):
                value_entry.append({
                    "groupId": f'{group_id}',
                    "valueMin": parameter_dict_default[parameter["name"]][0][i],
                    "valueMax": parameter_dict_default[parameter["name"]][1][i]
                })
            modelparameters_entry.append({
                "parameterId": f'{parameter["id"]}',
                "values": value_entry
            })

    get_interventions = requests.get(
        url + "interventions/templates/", headers=headers)

    school_closure_id = [intervention["id"] for intervention in get_interventions.json(
    ) if intervention["name"] == "School closure"]
    facemasks_school_id = [intervention["id"] for intervention in get_interventions.json(
    ) if intervention["name"] == "Face masks & social distancing School"]
    facemasks_work_id = [intervention["id"] for intervention in get_interventions.json(
    ) if intervention["name"] == "Face masks & social distancing Work"]
    facemasks_other_id = [intervention["id"] for intervention in get_interventions.json(
    ) if intervention["name"] == "Face masks & social distancing Other"]
    remote_work_id = [intervention["id"] for intervention in get_interventions.json(
    ) if intervention["name"] == "Remote work"]

    intervention_data_extended = [
        {
            "id": school_closure_id[0],
            "name": "School closure",
            "description": "School closure intervention",
            "tags": [],
            "coefficient": 1.0
        },
        {
            "id": facemasks_school_id[0],
            "name": "Face masks & social distancing School",
            "description": "Face mask usage and social distancing measures applied at 1-25% in schools",
            "tags": [],
            "coefficient": 0.25
        },
        {
            "id": facemasks_work_id[0],
            "name": "Face masks & social distancing Work",
            "description": "Face mask usage and social distancing measures applied at 1-25% in workplaces",
            "tags": [],
            "coefficient": 0.25
        },
        {
            "id": facemasks_other_id[0],
            "name": "Face masks & social distancing Other",
            "description": "Face mask usage and social distancing measures applied at 1-25% in other settings",
            "tags": [],
            "coefficient": 0.25
        },
        {
            "id": remote_work_id[0],
            "name": "Remote work",
            "description": "Implementation of remote work policies",
            "tags": [],
            "coefficient": 0.35
        }
    ]

    # Get all possible combinations of interventions
    all_combinations = []
    for r in range(1, len(intervention_data_extended) + 1):
        all_combinations.extend(combinations(intervention_data_extended, r))

    # Define Baseline scenario without any interventions and extrapolated scenario
    scenario_data = [{
        "name": "casedata",
        "description": "Extrapolierte Daten des RKI.",
        "startDate": f"{start_date_casedata}",
        "endDate": f"{end_date_casedata}",
        "modelId": model_id[0],
        "modelParameters": modelparameters_entry,
        "nodeListId": nodelist_id[0],
        "linkedInterventions": [],
        "percentiles": [50],
    },
        {
        "name": "baseline",
        "description": "Basisszenario ohne Interventionen.",
        "startDate": f"{start_date_simulation}",
        "endDate": f"{end_date_simulation}",
        "modelId": model_id[0],
        "modelParameters": modelparameters_entry,
        "nodeListId": nodelist_id[0],
        "linkedInterventions": [],
        "percentiles": [25, 50, 75],
    }]

    # Define remaining scenarios with different combinations of interventions
    for i, combination in enumerate(all_combinations):
        intervention_entry = []
        for intervention in combination:
            intervention_entry.append({
                "interventionId": intervention["id"],
                "startDate": f"{start_date_simulation}",
                "endDate": f"{end_date_simulation}",
                "coefficient": intervention["coefficient"]
            })
        scenario_data.append({
            "name": f"Scenario{i+2}",
            "description": "",
            "startDate": f"{start_date_simulation}",
            "endDate": f"{end_date_simulation}",
            "modelId": model_id[0],
            "modelParameters": modelparameters_entry,
            "nodeListId": nodelist_id[0],
            "linkedInterventions": intervention_entry,
            "percentiles": [25, 50, 75],
        })

    # Define optimal control scenario
    intervention_entry = []
    for intervention in intervention_data_extended:
        intervention_entry.append({
            "interventionId": intervention["id"],
            "startDate": f"{start_date_simulation}",
            "endDate": f"{end_date_simulation}",
            "coefficient": intervention["coefficient"]
        })
    scenario_data.append({
        "name": f"Scenario1_OptimalControl",
        "description": "",
        "startDate": f"{start_date_simulation}",
        "endDate": f"{end_date_simulation}",
        "modelId": model_id[0],
        "modelParameters": modelparameters_entry,
        "nodeListId": nodelist_id[0],
        "linkedInterventions": intervention_entry,
        "percentiles": [25, 50, 75],
    })

    if post:
        for scenario in scenario_data:
            post_response = requests.post(
                url + "scenarios/", json=scenario, headers=headers)
            if (post_response.status_code != 200):
                print(post_response.status_code)
                print(post_response.reason)

    return scenario_data


def append_parameter_for_agegroup_total(param_dict):
    # Append parameter value for age group "total". This additional parameter is computed by weighting the other
    # parameter values according to their share of the total population of Germany (using data from regionalstatistik).
    share_of_agegroup = [0.04773178, 0.09029715,
                         0.22754236, 0.34473159, 0.21830716, 0.07138996]
    for key in list(param_dict.keys()):

        param_dict[key][0].append(
            np.dot(param_dict[key][0], share_of_agegroup))
        param_dict[key][1].append(
            np.dot(param_dict[key][1], share_of_agegroup))


def post_to_db(headers, start_date_scenarios, days_simulated_scenarios, days_computed_casedata):
    # Define all necessary data and post.
    post_to_db_compartments(headers)
    post_to_db_groups(headers)
    post_to_db_interventions(headers)
    post_to_db_parameters(headers)
    post_to_db_model(headers)
    post_to_db_nodes(headers)
    post_to_db_nodelist(headers)
    post_to_db_scenarios(headers, start_date_scenarios,
                         days_simulated_scenarios, days_computed_casedata)


def main():
    """ With this script, we first delete all scenarios that are in the database and then post them again to the database.
    For the used parameters, see the functions above. """

    # Information needed for authentication.
    service_realm = "lha-a"
    client_id = "loki-front-dev"
    username = "simulation@ginger.example"
    password = "12345"

    # Request access token from idp.
    req_auth = requests.post(f'https://lokiam.de/realms/{service_realm}/protocol/openid-connect/token', data={
        'client_id': client_id, 'grant_type': 'password', 'username': username, 'password': password, 'scope': 'loki-back-audience'})

    # Raise error if response not ok.
    req_auth.raise_for_status()

    # Parse response and get access token.
    res = req_auth.json()
    token = res['access_token']

    # Set up headers with token.
    headers = {'Authorization': f'Bearer {token}', 'X-Realm': service_realm}

    print("Delete everything from db.")
    delete_everything_from_db(headers)

    print("Fill db.")
    start_date_scenarios = "2025-05-25"
    days_simulated_scenarios = 30
    days_computed_casedata = 30
    post_to_db(headers, start_date_scenarios, days_simulated_scenarios,
               days_computed_casedata)


if __name__ == "__main__":
    main()
