import requests
import datetime
import json
from itertools import combinations
import memilio.epidata.geoModificationGermany as gMG

header = {'Authorization': "Bearer anythingAsPasswordIsFineCurrently"}

url = "http://localhost:8123/"


def delete_everything_from_db():
    # GET and DELETE everything
    response = requests.get(url + "compartments/", headers=header)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(
            url + "compartments/" + id, headers=header)
        print(delete_response.status_code)

    print(requests.get(url + "compartments/", headers=header).json())

    # response = requests.get(url + "groups/", headers=header)
    # print(len(response.json()))

    response = requests.get(url + "groups/", headers=header)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(url + "groups/" + id, headers=header)
        print(delete_response.status_code)

    print(requests.get(url + "groups/", headers=header).json())

    response = requests.get(url + "nodelists/", headers=header)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(
            url + "nodelists/" + id, headers=header)
        print(delete_response.status_code)
    print(requests.get(url + "nodelists/", headers=header).json())

    response = requests.get(url + "nodes/", headers=header)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(url + "nodes/" + id, headers=header)
        print(delete_response.status_code)

    print(requests.get(url + "nodes/", headers=header).json())

    response = requests.get(url + "models/", headers=header)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(url + "models/" + id, headers=header)
        print(delete_response.status_code)

    print(requests.get(url + "models/", headers=header).json())

    response = requests.get(url + "parameterdefinitions/", headers=header)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(
            url + "parameterdefinitions/" + id, headers=header)
        print(delete_response.status_code)

    print(requests.get(url + "parameterdefinitions/", headers=header).json())

    response = requests.get(url + "interventions/templates/", headers=header)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(
            url + "interventions/templates/" + id, headers=header)
        print(delete_response.status_code)

    print(requests.get(url + "interventions/templates/", headers=header).json())

    response = requests.get(url + "scenarios/", headers=header)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(
            url + "scenarios/" + id, headers=header)
        print(delete_response.status_code)

    print(requests.get(url + "scenarios/", headers=header).json())


def post_to_db_compartments():
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
            url + "compartments/", json=compartment, headers=header)
        if (post_response.status_code != 200):
            print(post_response.status_code)


def post_to_db_groups():
    group_data = [
        {"name": "age_0", "description": "0-4", "category": "age"},
        {"name": "age_1", "description": "5-14", "category": "age"},
        {"name": "age_2", "description": "15-34", "category": "age"},
        {"name": "age_3", "description": "35-59", "category": "age"},
        {"name": "age_4", "description": "60-79", "category": "age"},
        {"name": "age_5", "description": "80+", "category": "age"}
    ]

    for group in group_data:
        post_response = requests.post(
            url + "groups/", json=group, headers=header)
        if (post_response.status_code != 200):
            print(post_response.status_code)


def post_to_db_interventions():
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
            url + "interventions/templates/", json=intervention, headers=header)
        if (post_response.status_code != 200):
            print(post_response.status_code)


def post_to_db_parameters():
    parameter_data = [
        {"name": "TimeExposed",
         "description": "TimeExposed"},
        {"name": "TimeInfectedNoSymptom",
         "description": "TimeInfectedNoSymptom"},
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
        {"name": "ReducExposedPartialImmunity",
         "description": "ReducExposedPartialImmunity"},
        {"name": "ReducExposedImprovedImmunity",
         "description": "ReducExposedImprovedImmunity"},
        {"name": "ReducInfectedSymptomsPartialImmunity",
         "description": "ReducInfectedSymptomsPartialImmunity"},
        {"name": "ReducInfectedSymptomsImprovedImmunity",
         "description": "ReducInfectedSymptomsImprovedImmunity"},
        {"name": "ReducInfectedSevereCriticalDeadPartialImmunity",
         "description": "ReducInfectedSevereCriticalDeadPartialImmunity"},
        {"name": "ReducInfectedSevereCriticalDeadImprovedImmunity",
         "description": "ReducInfectedSevereCriticalDeadImprovedImmunity"},
        {"name": "ReducTimeInfectedMild",
         "description": "ReducTimeInfectedMild"},
        {"name": "Seasonality",
         "description": "Seasonality"}
    ]

    for parameter in parameter_data:
        post_response = requests.post(
            url + "parameterdefinitions/", json=parameter, headers=header)
        if (post_response.status_code != 200):
            print(post_response.status_code)


def post_to_db_model():
    get_compartments = requests.get(url + "compartments/", headers=header)
    compartments = [compartment["id"]
                    for compartment in get_compartments.json()]

    get_groups = requests.get(url + "groups/", headers=header)
    groups = [group["id"] for group in get_groups.json()]

    get_parameters = requests.get(
        url + "parameterdefinitions/", headers=header)
    parameter_ids = [parameter["id"] for parameter in get_parameters.json()]

    model_data = {
        "name": "secirvvs",
        "description": "ODE-SECIRVVS model",
        "compartments": compartments,
        "groups": groups,
        "parameterDefinitions": parameter_ids
    }

    post_response = requests.post(
        url + "models/", json=model_data, headers=header)
    if (post_response.status_code != 200):
        print(post_response.status_code)


def post_to_db_nodes():
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
            url + "nodes/", json=node, headers=header)
        if (post_response.status_code != 200):
            print(post_response.status_code)


def post_to_db_nodelist():
    get_nodes = requests.get(url + "nodes/", headers=header)
    nodes = [node["id"] for node in get_nodes.json()]

    nodelist_data = {
        "name": "all_counties",
        "description": "Contains all counties of Germany as well as a node for the whole of Germany.",
        "nodeIds": nodes
    }

    post_response = requests.post(
        url + "nodelists/", json=nodelist_data, headers=header)
    if (post_response.status_code != 200):
        print(post_response.status_code)


def post_to_db_scenarios():
    # Define start and end date of simulation
    start_date = datetime.datetime.now().strftime("%Y-%m-%d")
    end_date = (datetime.datetime.now() +
                datetime.timedelta(days=30)).strftime("%Y-%m-%d")

    # Get ids of model, nodelist and interventions
    get_models = requests.get(url + "models/", headers=header)
    model_id = [model["id"]
                for model in get_models.json() if model["name"] == "secirvvs"]

    get_nodelists = requests.get(url+"nodelists/", headers=header)
    nodelist_id = [nodelist["id"] for nodelist in get_nodelists.json(
    ) if nodelist["name"] == "all_counties"]

    # Get sorted list of groups
    get_groups = requests.get(url + "groups/", headers=header)
    group_ids = []
    for i, group in enumerate(get_groups.json()):
        if group["name"] == f"age_{i}":
            group_ids.append(group["id"])

    variantFactor = 1.4
    # Define dict that contains name of parameter and min and max values for all 6 age groups.
    parameter_dict_default = {
        "TimeExposed": [[2.67, 2.67, 2.67, 2.67, 2.67, 2.67], [4., 4., 4., 4., 4., 4.]],
        "TimeInfectedNoSymptom": [[1.2, 1.2, 1.2, 1.2, 1.2, 1.2], [2.53, 2.53, 2.53, 2.53, 2.53, 2.53]],
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
        "ReducExposedPartialImmunity": [[0.75, 0.75, 0.75, 0.75, 0.75, 0.75], [0.85, 0.85, 0.85, 0.85, 0.85, 0.85]],
        "ReducExposedImprovedImmunity": [[0.281, 0.281, 0.281, 0.281, 0.281, 0.281], [0.381, 0.381, 0.381, 0.381, 0.381, 0.381]],
        "ReducInfectedSymptomsPartialImmunity": [[0.6, 0.6, 0.6, 0.6, 0.6, 0.6], [0.7, 0.7, 0.7, 0.7, 0.7, 0.7]],
        "ReducInfectedSymptomsImprovedImmunity": [[0.193, 0.193, 0.193, 0.193, 0.193, 0.193], [0.293, 0.293, 0.293, 0.293, 0.293, 0.293]],
        "ReducInfectedSevereCriticalDeadPartialImmunity": [[0.05, 0.05, 0.05, 0.05, 0.05, 0.05], [0.15, 0.15, 0.15, 0.15, 0.15, 0.15]],
        "ReducInfectedSevereCriticalDeadImprovedImmunity": [[0.041, 0.041, 0.041, 0.041, 0.041, 0.041], [0.141, 0.141, 0.141, 0.141, 0.141, 0.141]],
        "ReducTimeInfectedMild": [[1., 1., 1., 1., 1., 1.], [1., 1., 1., 1., 1., 1.]],
        "Seasonality": [[0.1, 0.1, 0.1, 0.1, 0.1, 0.1], [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]]}

    get_parameters = requests.get(
        url + "parameterdefinitions/", headers=header)

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
        url + "interventions/templates/", headers=header)

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
            "coefficient": 1.0.
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
        "description": "Extrapolated RKI data",
        "startDate": f"{start_date}",
        "endDate": f"{end_date}",
        "modelId": model_id[0],
        "modelParameters": modelparameters_entry,
        "nodeListId": nodelist_id[0],
        "linkedInterventions": []
    },
        {
        "name": "baseline",
        "description": "Baseline scenario without any interventions",
        "startDate": f"{start_date}",
        "endDate": f"{end_date}",
        "modelId": model_id[0],
        "modelParameters": modelparameters_entry,
        "nodeListId": nodelist_id[0],
        "linkedInterventions": []
    }]

    # Define remaining scenarios with different combinations of interventions
    for i, combination in enumerate(all_combinations):
        intervention_entry = []
        for intervention in combination:
            intervention_entry.append({
                "interventionId": intervention["id"],
                "startDate": f"{start_date}",
                "endDate": f"{end_date}",
                "coefficient": intervention["coefficient"]
            })
        scenario_data.append({
            "name": f"Scenario {i+2}",
            "description": "",
            "startDate": f"{start_date}",
            "endDate": f"{end_date}",
            "modelId": model_id[0],
            "modelParameters": modelparameters_entry,
            "nodeListId": nodelist_id[0],
            "linkedInterventions": intervention_entry
        })

    for scenario in scenario_data:
        post_response = requests.post(
            url + "scenarios/", json=scenario, headers=header)
        if (post_response.status_code != 200):
            print(post_response.status_code)


def post_to_db():
    # Define all necessary data and post.
    post_to_db_compartments()
    post_to_db_groups()
    post_to_db_interventions()
    post_to_db_parameters()
    post_to_db_model()
    post_to_db_nodes()
    post_to_db_nodelist()


def main():
    post_to_db()
    # delete_everything_from_db()


if __name__ == "__main__":
    main()
