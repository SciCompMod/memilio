import requests
import os
from zipfile import ZipFile
from pathlib import PurePath

import logging
import time

log = logging.getLogger(__file__)

# header = {'Authorization': "Bearer anythingAsPasswordIsFineCurrently"}

def read_scenario_description(path_to_scenario_description):
    """! Zips scenario for put to backend.

    @param[in] path_to_scenario_description Path where zips will be stored.
    """
    print(path_to_scenario_description)
    log.info(path_to_scenario_description)

    with open(os.path.join(path_to_scenario_description, "scenario_description.txt"), "r") as text_file:
        description_string = text_file.read()
    
    return description_string

def write_zip(path_to_saved_zips, zipped_name, percentiles=['p50'], case_data=False):
    """! Zips scenario for put to backend.

    @param[in] path_to_saved_zips Path where zips will be stored.
    @param[in] zipped_name name of zip. ".zip" is added to name.
    @param[in] percentiles. List[str]
    @param[in] case_data Flag to distinguish between casedata and scenarios.
    """
    print(path_to_saved_zips)
    log.info(path_to_saved_zips)

    zipfile = os.path.join(path_to_saved_zips, zipped_name + ".zip")
    with ZipFile(zipfile, "w") as zip_object:
        for percentile in percentiles:
            if case_data:
                # The results from the casedata scenario are stored in a directory above path_to_saved_zips which is
                # why we are adapting the path here. (We go up twice because path ends with "/".)
                path_parts = PurePath(path_to_saved_zips).parts[0:-1]
                real_path_to_data = os.path.join(*path_parts)
                print(real_path_to_data)
                log.info(real_path_to_data)
            else:
                # Scenario results are stored  in subfolder of path_to_saved_zips specified by scenario id and name.
                real_path_to_data = os.path.join(
                    path_to_saved_zips, zipped_name, percentile)

            # Write results in zip.
            zip_object.write(os.path.join(
                real_path_to_data, "Results.h5"), f"{percentile}/Results.h5")
            zip_object.write(os.path.join(
                real_path_to_data, "Results_sum.h5"), f"{percentile}/Results_sum.h5")

    return zipfile

def put_scenario(scenario_id, zip_file, url, delay, service_realm, client_id, username, password, scenario_description, max_tries=3, delay_check=5):
    # https://stackoverflow.com/questions/18208109/python-requests-post-a-zip-file-with-multipart-form-data
    headers = define_headers(service_realm, client_id, username, password)

    try_iteration = 0
    upload_scenario_successful = False
    upload_description_successful = False
    if scenario_description is None:
        upload_description_successful = True
    while (try_iteration < max_tries) and not (upload_scenario_successful and upload_description_successful):
        print(f'Waiting for {1+delay*try_iteration} seconds before putting scenario')
        time.sleep(1+try_iteration*delay)
        try_iteration += 1

        if upload_scenario_successful is False:
            with open(zip_file, 'rb') as fileobj:
                put_response = requests.put(url + "scenarios/" + scenario_id, headers=headers,
                                            files={"file": (zip_file, fileobj)})
            print(
                f'Put HTTP response code for scenario {scenario_id} was {put_response.status_code}, reason was {put_response.reason}.')

            if put_response.status_code != 200:
                print(put_response.text)
        
        if upload_description_successful is False:
            put_response = requests.put(url + "scenarios/" + scenario_id + "/description", headers=headers,
                                        data={"description": scenario_description})

            print(
                f'Put HTTP response code for description of scenario {scenario_id} was {put_response.status_code}, reason was {put_response.reason}.')

            if put_response.status_code != 200:
                print(put_response.text)

        time.sleep(delay_check)
        get_scenario_response = requests.get(
            url + "scenarios/" + scenario_id, headers=headers)
        if get_scenario_response.status_code != 200:
            print(get_scenario_response)
            continue
        else:
            get_scenario_response = get_scenario_response.json()
        
        # Check upload of scenario
        if (upload_scenario_successful is False) and (get_scenario_response["timestampSimulated"] is not None):
            print(
                f'Upload of scenario {get_scenario_response["id"]} was successful, timestampSimulated is not None.')
            upload_scenario_successful = True
        else:
            print(
                f'Upload of scenario {get_scenario_response["id"]} was not successful, timestampSimulated is None.')
            if try_iteration < max_tries:
                print(f'Retrying Upload of scenario {get_scenario_response["id"]}.')

        # Check upload of description
        if (upload_description_successful is False) and (get_scenario_response["description"] == scenario_description):
            print(
                f'Upload of description of scenario {get_scenario_response["id"]} was successful, timestampSimulated is not None.')
            upload_description_successful = True
        else:
            print(
                f'Upload of description of scenario {get_scenario_response["id"]} was not successful.')
            if try_iteration < max_tries:
                print(f'Retrying Upload of scenario {get_scenario_response["id"]}.')



def put_scenarios(path_to_scenario_results, url, delay, service_realm, client_id, username, password):
    """ Puts scenarios into database.

    @param[in] path_to_scenario_results Directory from where we can access simulation results and where the zips for 
        uploading will be stored. We assume that the casedata results lie in the parent folder of this directory, i.e. 
        the files Results.h5 and Result_sum.h5 lie in the parent directory. For each scenario, we expect the results 
        to be stored in a subfolder of the directory 'path_to_scenario_results' specified by the scenario id and name. 
        In this folder we have subfolders for the percentiles and for each percentile, we expect the files Results.h5 
        and Result_sum.h5. 
    @param[in] url URL of API.
    """
    print(f'Uploading scenarios to {url} from {path_to_scenario_results}')
    headers = define_headers(service_realm, client_id, username, password)
    scenarios = requests.get(url + "scenarios/", headers=headers)
    print(f'scenarios response code was {scenarios.status_code}')
    print(f'scenarios response was {scenarios.content}')
    scenarios = scenarios.json()
    for scenario in scenarios:
        scenario_description = None
        if scenario["name"] == "casedata":
            # if not os.path.exists(os.path.join(path_to_scenario_results+"Results.h5")):
            #     log.error(
            #         f'Path {os.path.join(path_to_scenario_results+"Results.h5")} does not exist')
            #     continue

            percentiles = ['p50']
            zip_file = write_zip(path_to_saved_zips=path_to_scenario_results,
                                 zipped_name=f"{scenario['name']}_{scenario['id']}", percentiles=percentiles,
                                 case_data=True)
            
        # Scenario of optimal control
        elif "OptimalControl" in scenario["name"]:
            scenario_path = path_to_scenario_results + \
                f"{scenario['name']}_{scenario['id']}/"
            if not os.path.exists(scenario_path):
                log.error(f'Path {scenario_path} does not exist')
                continue

            percentiles = ['p25', 'p50', 'p75']
            zip_file = write_zip(path_to_saved_zips=path_to_scenario_results,
                                 zipped_name=f"{scenario['name']}_{scenario['id']}", percentiles=percentiles,
                                 case_data=False)

            scenario_description = read_scenario_description(path_to_scenario_results)

        else:
            scenario_path = path_to_scenario_results + \
                f"{scenario['name']}_{scenario['id']}/"
            if not os.path.exists(scenario_path):
                log.error(f'Path {scenario_path} does not exist')
                continue

            percentiles = ['p25', 'p50', 'p75']
            zip_file = write_zip(path_to_saved_zips=path_to_scenario_results,
                                 zipped_name=f"{scenario['name']}_{scenario['id']}", percentiles=percentiles,
                                 case_data=False)

        print(f'Uploading {scenario["name"]} with id {scenario["id"]}')
        put_scenario(scenario['id'], zip_file=zip_file, url=url, delay=delay,
                     service_realm=service_realm, client_id=client_id, username=username, password=password, scenario_description=scenario_description)


def define_headers(service_realm, client_id, username, password):
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

    return headers


def main():
    url = 'https://zam10063.zam.kfa-juelich.de/api-dev/'

    # Delay for time between upload of scenarios.
    delay = 420

    # Information needed for authentication.
    service_realm = ""
    client_id = ""
    username = ""
    password = ""

    # put_scenarios(os.path.join('/home/jadebeck/', "results_osecirvvs/"), url)
    put_scenarios(os.path.join(
        os.getcwd(), "results_retros_test/"), url, delay, service_realm, client_id, username, password)


if __name__ == '__main__':
    main()
