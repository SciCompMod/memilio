import requests
import os
from zipfile import ZipFile

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
    log.info(path_to_saved_zips)

    zipfile = os.path.join(path_to_saved_zips, zipped_name + ".zip")
    with ZipFile(zipfile, "w") as zip_object:
        for percentile in percentiles:
            if case_data:
                # The results from the casedata scenario are stored in a directory above path_to_saved_zips which is
                # why we are adapting the path here. (We go up twice because path ends with "/".)
                real_path_to_data = os.path.dirname(
                    os.path.dirname(path_to_saved_zips))

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


def put_scenario(scenario_id, zip_file, url):
    # https://stackoverflow.com/questions/18208109/python-requests-post-a-zip-file-with-multipart-form-data
    fileobj = open(zip_file, 'rb')
    put_response = requests.put(url + "scenarios/" + scenario_id, headers=header,
                                files={"file": (zip_file, fileobj)})
    print(
        f'Put HTTP response code for scenario {scenario_id} was {put_response.status_code}')


def put_scenarios(path_to_scenario_results, url, delay=10):
    """ Puts scenarios into database.

    @param[in] path_to_scenario_results Directory from where we can access simulation results and where the zips for 
        uploading will be stored. The casedata results lie in the parent folder of this directory; the scenario results 
        lie in a subfolder of this folder specified by the scenario id and name.
    @param[in] url URL of API.
    @param[in] delay Delay between uploads in seconds. 
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
                f"/{scenario['name']}_{scenario['id']}/"
            if not os.path.exists(scenario_path):
                log.error(f'Path {scenario_path} does not exist')
                continue

            percentiles = ['p25', 'p50', 'p75']
            zip_file = write_zip(path_to_saved_zips=path_to_scenario_results,
                                 zipped_name=f"{scenario['name']}_{scenario['id']}", percentiles=percentiles,
                                 case_data=False)

        print(f'uploading {scenario["id"]}')
        put_scenario(scenario['id'], zip_file=zip_file, url=url)
        print(f'uploaded {scenario["id"]}')
        print(f'waiting for  {delay} seconds before uploading next scenario')
        time.sleep(delay)


def main():
    url = 'https://zam10063.zam.kfa-juelich.de/api-new/'
    put_scenarios(os.path.join('/home/jadebeck/', "results_osecirvvs/"), url)


if __name__ == '__main__':
    main()
