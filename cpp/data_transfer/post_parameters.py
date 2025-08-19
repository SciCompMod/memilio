import requests
import datetime
import os
import csv
import json
from itertools import combinations
import numpy as np
import memilio.epidata.geoModificationGermany as gMG

# header = {'Authorization': "Bearer anythingAsPasswordIsFineCurrently"}


def delete_scenarios(url, headers):

    response = requests.get(url + "scenarios/", headers=headers)
    ids = [compartment["id"] for compartment in response.json()]
    for id in ids:
        delete_response = requests.delete(
            url + "scenarios/" + id, headers=headers)
        if (delete_response.status_code != 200):
            print(delete_response.status_code)
            print(delete_response.reason)
        else:
            print("Deleting scenarios has been successful.")

    print(requests.get(url + "scenarios/", headers=headers).json())


def read_mcmc_data(data_dir):
    """!Reads MCMC data from CSV files in a specified directory and organizes them into a dictionary.

    The dictionary has the groups as keys. The values are the dates which are obtained from the file names and
    the actual data from the CSV files.
    """
    data = {}
    for file in os.listdir(data_dir):
        # only consider files that start with pipeline_mcmc_summary_0
        if file.startswith("pipeline_mcmc_summary_0"):
            group, date = file.replace(
                "pipeline_mcmc_summary_0_", "").split("_")
            # get date from file name and read data from file. Then add data to dictionary.
            date = date.split(".")[0]
            with open(os.path.join(data_dir, file), "r") as f:
                reader = csv.DictReader(f)
                csv_data = [row for row in reader]
                if group not in data:
                    data[group] = []
                data[group].append({"date": date, "data": csv_data})
    return data


def get_mcmc_data_at_t(t, data_dir):
    """!Returns the MCMC data for a specific date t from the specified directory.
    """
    data = read_mcmc_data(data_dir)
    data_at_t = {}
    for group in data:
        for entry in data[group]:
            if entry["date"] == t:
                data_at_t[group] = entry["data"]
    if not data_at_t:
        print("No data found for date " + t)
    return data_at_t


def get_mcmc_model_params(url, t, data_dir, headers):
    """!Reads MCMC data for a specific date t and creates an entry for the model parameters.
    """
    modelparameters_entry = []
    data_t = get_mcmc_data_at_t(t, data_dir)
    parameter_list = requests.get(
        url + "parameterdefinitions/", headers=headers).json()

    get_groups = requests.get(url + "groups/", headers=headers)
    group_ids = []
    for i, group in enumerate(get_groups.json()):
        if group["name"] == f"age_{i}":
            group_ids.append(group["id"])

    # list all parameters that are in data_t
    data_t_parameters = []
    for group in data_t:
        for entry in data_t[group]:
            # ignore all entires starting with "sim" or "ChangePoint"
            if entry["Name"].startswith("sim") or entry["Name"].startswith("ChangePoint") or entry["Name"].startswith("ContactPattern"):
                continue
            if entry["Name"] not in data_t_parameters:
                data_t_parameters.append(entry["Name"])
    # append Seasonality as this is not estimated using mcmc
    data_t_parameters.append('Seasonality')

    for parameter in data_t_parameters:
        # seach the parameter in the parameter_list
        parameter_id = None
        for entry in parameter_list:
            if entry["name"] == parameter:
                parameter_id = entry["id"]
        # if parameter is not found in parameter_list, print error message
        if parameter_id is None:
            print("Parameter " + parameter + " not found in parameter list.")
            continue

        param_val = []

        # For 'Seasonality use default parameters, as this is not estimated by mcmc
        seasonality_default = [[0.1, 0.1, 0.1, 0.1, 0.1, 0.1], [
            0.3, 0.3, 0.3, 0.3, 0.3, 0.3]]
        if parameter == 'Seasonality':
            for group in range(len(data_t)):
                param_val.append(
                    [seasonality_default[0][group], seasonality_default[1][group]])

        # use data from mcmc
        else:
            for group in data_t:
                for entry in data_t[group]:
                    if entry["Name"] == parameter:
                        param_val.append(
                            [float(entry["hdi_25%"]), float(entry["hdi_75%"])])
            if len(param_val) != len(data_t):
                print("Parameter " + parameter +
                      " has not the correct number of groups.")
                continue

        # Append parameter value for age group "total". This is done by weighting the other parameter values according to
        # their share of the total population of Germany (using data from regionalstatistik).
        share_of_agegroup = [0.04773178, 0.09029715,
                             0.22754236, 0.34473159, 0.21830716, 0.07138996]
        min_vector = [values[0] for values in param_val]
        max_vector = [values[1] for values in param_val]
        param_val.append([np.dot(share_of_agegroup, min_vector),
                         np.dot(share_of_agegroup, max_vector)])

        # Get sorted list of groups
        get_groups = requests.get(url + "groups/", headers=headers)
        group_ids = []
        for i, group in enumerate(get_groups.json()):
            if group["name"] == f"age_{i}":
                group_ids.append(group["id"])

        value_entry = []
        for group_index, group_id in enumerate(group_ids):
            value_entry.append({
                "groupId": f'{group_id}',
                "valueMin": param_val[group_index][0],
                "valueMax": param_val[group_index][1]
            })
        modelparameters_entry.append({
            "parameterId": parameter_id,
            "values": value_entry
        })
    return modelparameters_entry


def post_mcmc_parameters(url, t, data_dir, service_realm="", client_id="", username="", password=""):

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

    modelparameters_entry = []
    data_t = get_mcmc_data_at_t(t, data_dir)
    parameter_list = requests.get(
        url + "parameterdefinitions/", headers=headers).json()

    # Get mcmc parameters and put it into the right format for the API.
    parameters_mcmc = get_mcmc_model_params(
        url, t, data_dir, headers)

    # # list all scenarios
    scenarios = requests.get(url + "scenarios/", headers=headers).json()

    new_scenarios = []

    for scenario in scenarios:
        # load full set of scenario data
        scenario_data = requests.get(
            url + "scenarios/" + scenario['id'], headers=headers)

        if (scenario_data.status_code != 200):
            print(scenario_data.status_code)
            print(scenario_data.reason)

        # Upload the model parameters to the scenario
        scenario_data.json()['modelParameters'] = parameters_mcmc

        new_scenarios.append(scenario_data.json())

    #     # put
    #     # response = requests.put(
    #     #     url + f"scenarios/{scenario['id']}",
    #     #     headers=header,
    #     #     json=scenario_data
    #     # )

    #     # patch
    #     # updated_data = {
    #     #     "modelParameters": parameters_mcmc
    #     # }
    #     # response = requests.patch(
    #     #     url + f"scenarios/{scenario['id']}",
    #     #     headers=header,
    #     #     json=updated_data
    #     # )

    # delete scenarios with old parameters
    delete_scenarios(url, headers)

    # Create new headers for posting scenarios.
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

    # post scenarios with updated parameters
    for scenario in new_scenarios:
        response = requests.post(
            url + "scenarios/",
            headers=headers,
            json=scenario
        )

        if response.status_code == 200:
            print(f"Scenario {scenario['id']} updated successfully.")
        else:
            print(
                f"Failed to update scenario {scenario['id']}. Status code: {response.status_code}, Response: {response.text}")


def main():
    url = "https://zam10063.zam.kfa-juelich.de/api-new/"
    cwd = os.getcwd()
    mcmc_dir = os.path.join(cwd, "mcmc data")
    # date_today = datetime.datetime.now().strftime("%Y-%m-%d")
    date_today = "2025-02-11"

    # Information needed for authentication.
    service_realm = ""
    client_id = ""
    username = ""
    password = ""

    # # Request access token from idp.
    # req_auth = requests.post(f'https://lokiam.de/realms/{service_realm}/protocol/openid-connect/token', data={
    #     'client_id': client_id, 'grant_type': 'password', 'username': username, 'password': password, 'scope': 'loki-back-audience'})

    # # Raise error if response not ok.
    # req_auth.raise_for_status()

    # # Parse response and get access token.
    # res = req_auth.json()
    # token = res['access_token']

    # # Set up headers with token.
    # headers = {'Authorization': f'Bearer {token}', 'X-Realm': service_realm}

    # Update scenarios with mcmc parameters.
    post_mcmc_parameters(url, date_today, mcmc_dir, service_realm=service_realm,
                         client_id=client_id, username=username, password=password)


if __name__ == "__main__":
    main()
