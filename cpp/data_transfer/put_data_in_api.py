import requests
import os
from zipfile import ZipFile
from pathlib import PurePath

import logging
import time

log = logging.getLogger(__file__)

header = {'Authorization': "Bearer anythingAsPasswordIsFineCurrently"}


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


def put_scenario(scenario_id, zip_file, url, delay):
    # https://stackoverflow.com/questions/18208109/python-requests-post-a-zip-file-with-multipart-form-data
    fileobj = open(zip_file, 'rb')
    put_response = requests.put(url + "scenarios/" + scenario_id, headers=header,
                                files={"file": (zip_file, fileobj)})
    print(
        f'Put HTTP response code for scenario {scenario_id} was {put_response.status_code}, reason was {put_response.reason}.')
    
    if put_response.status_code!=200:
        print(put_response.text)

    print(f'Waiting for {delay} seconds before uploading next scenario')
    time.sleep(delay)

    get_scenario_response = requests.get(
        url + "scenarios/" + scenario_id, headers=header).json()
    if get_scenario_response["timestampSimulated"] != None:
        print(
            f'Upload of scenario {get_scenario_response["id"]} was successful, timestampSimulated is not None.')
    else:
        print(
            f'Upload of scenario {get_scenario_response["id"]} was not successful, timestampSimulated is None.')


def put_scenarios(path_to_scenario_results, url, delay=420):
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
    scenarios = requests.get(url + "scenarios/", headers=header)
    print(f'scenarios response code was {scenarios.status_code}')
    print(f'scenarios response was {scenarios.content}')
    scenarios = scenarios.json()
    for scenario in scenarios:
        if scenario["name"] == "casedata":
            if not os.path.exists(path_to_scenario_results):
                log.error(
                    f'Path {path_to_scenario_results} does not exist')
                continue

            percentiles = ['p50']
            zip_file = write_zip(path_to_saved_zips=path_to_scenario_results,
                                 zipped_name=f"{scenario['name']}_{scenario['id']}", percentiles=percentiles,
                                 case_data=True)

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
        put_scenario(scenario['id'], zip_file=zip_file, url=url, delay=delay)


def main():
    url = 'https://zam10063.zam.kfa-juelich.de/api-new/'
    # put_scenarios(os.path.join('/home/jadebeck/', "results_osecirvvs/"), url)
    put_scenarios(os.path.join(os.getcwd(), "results_retros_adapted_deaths/"), url)


if __name__ == '__main__':
    main()
