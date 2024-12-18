import requests 
import datetime
import json
from itertools import combinations

header = {'Authorization': "Bearer anythingAsPasswordIsFineCurrently"}

url = "http://localhost:8123/"

def delete_everything_from_db():
  # GET and DELETE everything
  response = requests.get(url + "compartments/", headers=header)
  ids  = [compartment["id"] for compartment in response.json()]
  for id in ids: 
      delete_response = requests.delete(url + "compartments/" + id, headers = header)
      print(delete_response.status_code)

  print(requests.get(url + "compartments/", headers=header).json())

  # response = requests.get(url + "groups/", headers=header)
  # print(len(response.json()))

  response = requests.get(url + "groups/", headers=header)
  ids  = [compartment["id"] for compartment in response.json()]
  for id in ids: 
      delete_response=requests.delete(url + "groups/" + id, headers = header)
      print(delete_response.status_code)

  print(requests.get(url + "groups/", headers=header).json())

  response = requests.get(url + "nodelists/", headers=header)
  ids  = [compartment["id"] for compartment in response.json()]
  for id in ids: 
      delete_response=requests.delete(url + "nodelists/" + id, headers = header)
      print(delete_response.status_code)
  print(requests.get(url + "nodelists/", headers=header).json())

  response = requests.get(url + "nodes/", headers=header)
  ids  = [compartment["id"] for compartment in response.json()]
  for id in ids: 
      delete_response=requests.delete(url + "nodes/" + id, headers = header)
      print(delete_response.status_code)

  print(requests.get(url + "nodes/", headers=header).json())

  response = requests.get(url + "models/", headers=header)
  ids  = [compartment["id"] for compartment in response.json()]
  for id in ids: 
      delete_response=requests.delete(url + "models/" + id, headers = header)
      print(delete_response.status_code)

  print(requests.get(url + "models/", headers=header).json())

  response = requests.get(url + "parameterdefinitions/", headers=header)
  ids  = [compartment["id"] for compartment in response.json()]
  for id in ids: 
      delete_response=requests.delete(url + "parameterdefinitions/" + id, headers = header)
      print(delete_response.status_code)

  print(requests.get(url + "parameterdefinitions/", headers=header).json())

  response = requests.get(url + "interventions/templates/", headers=header)
  ids  = [compartment["id"] for compartment in response.json()]
  for id in ids: 
      delete_response=requests.delete(url + "interventions/templates/" + id, headers = header)
      print(delete_response.status_code)

  print(requests.get(url + "interventions/templates/", headers=header).json())

  response = requests.get(url + "scenarios/", headers=header)
  ids  = [compartment["id"] for compartment in response.json()]
  for id in ids: 
      delete_response=requests.delete(url + "scenarios/" + id, headers = header)
      print(delete_response.status_code)

  print(requests.get(url + "scenarios/", headers=header).json())


def post_to_db():
  ### Define all necessary data, post requests are at the bottom of the script. 

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
      post_response = requests.post(url + "compartments/", json=compartment, headers=header)

      # post_response_json = post_response.json()
      # print(post_response_json)

      print(post_response.status_code)



  group_data = [
      { "name": "age_0", "description": "0-4", "category": "age"},
      { "name": "age_1", "description": "5-14", "category": "age"},
      { "name": "age_2", "description": "15-34", "category": "age"},
      { "name": "age_3", "description": "35-59", "category": "age"},
      { "name": "age_4", "description": "60-79", "category": "age"},
      { "name": "age_5", "description": "80+", "category": "age"}
  ]

  for group in group_data:
      post_response = requests.post(url + "groups/", json=group, headers=header)
      print(post_response.status_code)

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
      post_response = requests.post(url + "interventions/templates/", json=intervention, headers=header)
      print(post_response.status_code)


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
      post_response = requests.post(url + "parameterdefinitions/", json=parameter, headers=header)
      print(post_response.status_code)

  get_compartments = requests.get(url + "compartments/", headers=header)
  compartments  = [compartment["id"] for compartment in get_compartments.json()]

  get_groups = requests.get(url + "groups/", headers=header)
  groups  = [group["id"] for group in get_groups.json()]

  get_parameters = requests.get(url + "parameterdefinitions/", headers=header)
  parameter_ids  = [parameter["id"] for parameter in get_parameters.json()]

  model_data = {
    "name": "secirvvs",
    "description": "ODE-SECIRVVS model",
    "compartments": compartments,
    "groups": groups,
    "parameterDefinitions": parameter_ids
  }

  post_response = requests.post(url + "models/", json=model_data, headers=header)
  print(post_response.status_code)

  nodes_init =  [
      "00000",
      "01001",
      "01002",
      "01003",
      "01004",
      "01051",
      "01053",
      "01054",
      "01055",
      "01056",
      "01057",
      "01058",
      "01059",
      "01060",
      "01061",
      "01062",
      "02000",
      "03101",
      "03102",
      "03103",
      "03151",
      "03153",
      "03154",
      "03155",
      "03157",
      "03158",
      "03159",
      "03241",
      "03251",
      "03252",
      "03254",
      "03255",
      "03256",
      "03257",
      "03351",
      "03352",
      "03353",
      "03354",
      "03355",
      "03356",
      "03357",
      "03358",
      "03359",
      "03360",
      "03361",
      "03401",
      "03402",
      "03403",
      "03404",
      "03405",
      "03451",
      "03452",
      "03453",
      "03454",
      "03455",
      "03456",
      "03457",
      "03458",
      "03459",
      "03460",
      "03461",
      "03462",
      "04011",
      "04012",
      "05111",
      "05112",
      "05113",
      "05114",
      "05116",
      "05117",
      "05119",
      "05120",
      "05122",
      "05124",
      "05154",
      "05158",
      "05162",
      "05166",
      "05170",
      "05314",
      "05315",
      "05316",
      "05334",
      "05358",
      "05362",
      "05366",
      "05370",
      "05374",
      "05378",
      "05382",
      "05512",
      "05513",
      "05515",
      "05554",
      "05558",
      "05562",
      "05566",
      "05570",
      "05711",
      "05754",
      "05758",
      "05762",
      "05766",
      "05770",
      "05774",
      "05911",
      "05913",
      "05914",
      "05915",
      "05916",
      "05954",
      "05958",
      "05962",
      "05966",
      "05970",
      "05974",
      "05978",
      "06411",
      "06412",
      "06413",
      "06414",
      "06431",
      "06432",
      "06433",
      "06434",
      "06435",
      "06436",
      "06437",
      "06438",
      "06439",
      "06440",
      "06531",
      "06532",
      "06533",
      "06534",
      "06535",
      "06611",
      "06631",
      "06632",
      "06633",
      "06634",
      "06635",
      "06636",
      "07111",
      "07131",
      "07132",
      "07133",
      "07134",
      "07135",
      "07137",
      "07138",
      "07140",
      "07141",
      "07143",
      "07211",
      "07231",
      "07232",
      "07233",
      "07235",
      "07311",
      "07312",
      "07313",
      "07314",
      "07315",
      "07316",
      "07317",
      "07318",
      "07319",
      "07320",
      "07331",
      "07332",
      "07333",
      "07334",
      "07335",
      "07336",
      "07337",
      "07338",
      "07339",
      "07340",
      "08111",
      "08115",
      "08116",
      "08117",
      "08118",
      "08119",
      "08121",
      "08125",
      "08126",
      "08127",
      "08128",
      "08135",
      "08136",
      "08211",
      "08212",
      "08215",
      "08216",
      "08221",
      "08222",
      "08225",
      "08226",
      "08231",
      "08235",
      "08236",
      "08237",
      "08311",
      "08315",
      "08316",
      "08317",
      "08325",
      "08326",
      "08327",
      "08335",
      "08336",
      "08337",
      "08415",
      "08416",
      "08417",
      "08421",
      "08425",
      "08426",
      "08435",
      "08436",
      "08437",
      "09161",
      "09162",
      "09163",
      "09171",
      "09172",
      "09173",
      "09174",
      "09175",
      "09176",
      "09177",
      "09178",
      "09179",
      "09180",
      "09181",
      "09182",
      "09183",
      "09184",
      "09185",
      "09186",
      "09187",
      "09188",
      "09189",
      "09190",
      "09261",
      "09262",
      "09263",
      "09271",
      "09272",
      "09273",
      "09274",
      "09275",
      "09276",
      "09277",
      "09278",
      "09279",
      "09361",
      "09362",
      "09363",
      "09371",
      "09372",
      "09373",
      "09374",
      "09375",
      "09376",
      "09377",
      "09461",
      "09462",
      "09463",
      "09464",
      "09471",
      "09472",
      "09473",
      "09474",
      "09475",
      "09476",
      "09477",
      "09478",
      "09479",
      "09561",
      "09562",
      "09563",
      "09564",
      "09565",
      "09571",
      "09572",
      "09573",
      "09574",
      "09575",
      "09576",
      "09577",
      "09661",
      "09662",
      "09663",
      "09671",
      "09672",
      "09673",
      "09674",
      "09675",
      "09676",
      "09677",
      "09678",
      "09679",
      "09761",
      "09762",
      "09763",
      "09764",
      "09771",
      "09772",
      "09773",
      "09774",
      "09775",
      "09776",
      "09777",
      "09778",
      "09779",
      "09780",
      "10041",
      "10042",
      "10043",
      "10044",
      "10045",
      "10046",
      "11000",
      "12051",
      "12052",
      "12053",
      "12054",
      "12060",
      "12061",
      "12062",
      "12063",
      "12064",
      "12065",
      "12066",
      "12067",
      "12068",
      "12069",
      "12070",
      "12071",
      "12072",
      "12073",
      "13003",
      "13004",
      "13071",
      "13072",
      "13073",
      "13074",
      "13075",
      "13076",
      "14511",
      "14521",
      "14522",
      "14523",
      "14524",
      "14612",
      "14625",
      "14626",
      "14627",
      "14628",
      "14713",
      "14729",
      "14730",
      "15001",
      "15002",
      "15003",
      "15081",
      "15082",
      "15083",
      "15084",
      "15085",
      "15086",
      "15087",
      "15088",
      "15089",
      "15090",
      "15091",
      "16051",
      "16052",
      "16053",
      "16054",
      "16055",
      "16061",
      "16062",
      "16063",
      "16064",
      "16065",
      "16066",
      "16067",
      "16068",
      "16069",
      "16070",
      "16071",
      "16072",
      "16073",
      "16074",
      "16075",
      "16076",
      "16077"
    ]

  node_data = [{"nuts": f"{node}",
                "name": f"{node}"} for node in nodes_init]

  for node in node_data: 
      post_response = requests.post(url + "nodes/", json=node, headers=header)
      print(post_response.status_code)

  get_nodes = requests.get(url + "nodes/", headers=header)
  nodes  = [node["id"] for node in get_nodes.json()]

  nodelist_data = {
    "name": "all_counties",
    "description": "Contains all counties of Germany as well as a node for the whole of Germany.",
    "nodeIds": nodes
  }

  post_response = requests.post(url + "nodelists/", json=nodelist_data, headers=header)
  print(post_response.status_code)


  ### Scenario 
  # Define start and end date of simulation
  start_date = datetime.datetime.now().strftime("%Y-%m-%d")
  end_date = (datetime.datetime.now() + datetime.timedelta(days=30)).strftime("%Y-%m-%d")

  # Get ids of model, nodelist and interventions
  get_models = requests.get(url + "models/", headers=header)
  model_id =[model["id"] for model in get_models.json() if model["name"]=="secirvvs"]

  get_nodelists = requests.get(url+"nodelists/", headers = header)
  nodelist_id = [nodelist["id"] for nodelist in get_nodelists.json() if nodelist["name"]=="all_counties"]


  # Get sorted list of groups
  get_groups = requests.get(url + "groups/", headers=header)
  group_ids = []
  for i, group in enumerate(get_groups.json()):
    if group["name"] == f"age_{i}":
      group_ids.append(group["id"])

  variantFactor = 1.4
  # Define dict that contains name of parameter and min and max values for all 6 age groups. 
  parameter_dict_default = {
      "TimeExposed": [[2.67,2.67,2.67,2.67,2.67,2.67], [4.,4.,4.,4.,4.,4.]],
      "TimeInfectedNoSymptom":[[1.2,1.2,1.2,1.2,1.2,1.2],[2.53,2.53,2.53,2.53,2.53,2.53]],
      "TimeInfectedSymptoms": [[5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465],[8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085]], 
      "TimeInfectedSevere": [[
                  3.925, 3.925, 4.85, 6.4, 7.2, 9.], [6.075, 6.075, 7., 8.7, 9.8, 13.]],
      "TimeInfectedCritical": [[4.95, 4.95, 4.86, 14.14, 14.4, 10.], [8.95, 8.95, 8.86, 20.58, 19.8, 13.2]],
      "TransmissionProbabilityOnContact": [[
                  0.02 * variantFactor, 0.05 * variantFactor, 0.05 * variantFactor,
                  0.05 * variantFactor, 0.08 * variantFactor, 0.1 * variantFactor
              ], [
                  0.04 * variantFactor, 0.07 * variantFactor, 0.07 * variantFactor,
                  0.07 * variantFactor, 0.10 * variantFactor, 0.15 * variantFactor
              ] ],
      "RelativeTransmissionNoSymptoms": [[0.5,0.5,0.5,0.5,0.5,0.5],[0.5,0.5,0.5,0.5,0.5,0.5]],
      "RiskOfInfectionFromSymptomatic": [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.2, 0.2, 0.2, 0.2, 0.2, 0.2]],
      "MaxRiskOfInfectionFromSymptomatic": [[0.4,0.4,0.4,0.4,0.4,0.4],[0.5,0.5,0.5,0.5,0.5,0.5]],
      "RecoveredPerInfectedNoSymptoms": [[
                  0.2, 0.2, 0.15, 0.15, 0.15, 0.15], [
                  0.3, 0.3, 0.25, 0.25, 0.25, 0.25]],
      "SeverePerInfectedSymptoms": [[
                  0.006, 0.006, 0.015, 0.049, 0.15, 0.20],  [
                  0.009, 0.009, 0.023, 0.074, 0.18, 0.25]],
      "CriticalPerSevere" : [[0.05, 0.05, 0.05, 0.10, 0.25, 0.35], [0.10, 0.10, 0.10, 0.20, 0.35, 0.45]],
      "DeathsPerCritical": [[0.00, 0.00, 0.10, 0.10, 0.30, 0.5], [0.10, 0.10, 0.18, 0.18, 0.50, 0.7]],
      "ReducExposedPartialImmunity": [[0.75,0.75,0.75,0.75,0.75,0.75], [0.85,0.85,0.85,0.85,0.85,0.85]],
      "ReducExposedImprovedImmunity": [[0.281,0.281,0.281,0.281,0.281,0.281], [0.381,0.381,0.381,0.381,0.381,0.381]],
      "ReducInfectedSymptomsPartialImmunity": [[0.6,0.6,0.6,0.6,0.6,0.6], [0.7,0.7,0.7,0.7,0.7,0.7]],
      "ReducInfectedSymptomsImprovedImmunity": [[0.193,0.193,0.193,0.193,0.193,0.193], [0.293,0.293,0.293,0.293,0.293,0.293]],
      "ReducInfectedSevereCriticalDeadPartialImmunity": [[0.05,0.05,0.05,0.05,0.05,0.05], [0.15,0.15,0.15,0.15,0.15,0.15]],
      "ReducInfectedSevereCriticalDeadImprovedImmunity": [[0.041,0.041,0.041,0.041,0.041,0.041],[0.141,0.141,0.141,0.141,0.141,0.141]],
      "ReducTimeInfectedMild": [[1.,1.,1.,1.,1.,1.], [1.,1.,1.,1.,1.,1.]],
      "Seasonality": [[0.1,0.1,0.1, 0.1, 0.1, 0.1],[0.3, 0.3, 0.3, 0.3, 0.3, 0.3]]}

  # get_parameters = requests.get(url + "parameterdefinitions/", headers=header)

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



  get_interventions = requests.get(url + "interventions/templates/", headers=header)

  school_closure_id = [intervention["id"] for intervention in get_interventions.json() if intervention["name"]=="School closure"]
  facemasks_school_id = [intervention["id"] for intervention in get_interventions.json() if intervention["name"]=="Face masks & social distancing School"]
  facemasks_work_id = [intervention["id"] for intervention in get_interventions.json() if intervention["name"]=="Face masks & social distancing Work"]
  facemasks_other_id = [intervention["id"] for intervention in get_interventions.json() if intervention["name"]=="Face masks & social distancing Other"]
  remote_work_id = [intervention["id"] for intervention in get_interventions.json() if intervention["name"]=="Remote work"]

  intervention_data_extended = [
      {
          "id": school_closure_id[0],
          "name": "School closure",
          "description": "School closure intervention",
          "tags": [],
          "coefficient" : 0.
      },
      {
          "id": facemasks_school_id[0],
          "name": "Face masks & social distancing School",
          "description": "Face mask usage and social distancing measures applied at 1-25% in schools",
          "tags": [],
          "coefficient" : 0.25
      },
      {
          "id": facemasks_work_id[0],
          "name": "Face masks & social distancing Work",
          "description": "Face mask usage and social distancing measures applied at 1-25% in workplaces",
          "tags": [],
          "coefficient" : 0.25
      },
      {
          "id": facemasks_other_id[0],
          "name": "Face masks & social distancing Other",
          "description": "Face mask usage and social distancing measures applied at 1-25% in other settings",
          "tags": [],
          "coefficient" : 0.25
      },
      {
          "id": remote_work_id[0],
          "name": "Remote work",
          "description": "Implementation of remote work policies",
          "tags": [],
          "coefficient" : 0.35
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

  # with open("scenarios.json", "w") as outfile: 
  #     json.dump(scenario_data, outfile)
    
  for scenario in scenario_data: 
      post_response = requests.post(url + "scenarios/", json=scenario, headers=header)
      if (post_response.status_code != 200):
        print(post_response.status_code)
  print()

def main():
  post_to_db()
  # delete_everything_from_db()

if __name__ == "__main__":
    main()

