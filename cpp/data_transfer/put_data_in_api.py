import requests
import os
from zipfile import ZipFile

import logging
import time

log = logging.getLogger(__file__)

header = {'Authorization': "Bearer anythingAsPasswordIsFineCurrently"}


def write_zip(path_to_data, zipped_name, percentiles=['p50'], case_data=False):
    """! Zips scenario for put to backend.

    @param[in] path_to_data Path to data.
    @param[in] zipped_name name of zip. ".zip" is added to name.
    @param[in] percentiles. List[str]
    """
    log.info(path_to_data)
    zipfile = os.path.join(path_to_data, zipped_name + ".zip")
    with ZipFile(zipfile, "w") as zip_object:
        for percentile in percentiles:
            if not case_data:
                real_path_to_data = os.path.join(path_to_data, zipped_name, percentile)
            else:
                real_path_to_data = os.path.join(path_to_data, zipped_name)
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
    print(f'Put HTTP response code for scenario {scenario_id} was {put_response.status_code}')


def put_scenarios(path_to_scenario_results, url, delay=10):
    print(f'Uploading scenarios to {url} from {path_to_scenario_results}')
    scenarios = requests.get(url + "scenarios/", headers=header)
    print(f'scenarios response code was {scenarios.status_code}')
    print(f'scenarios response was {scenarios.content}')
    scenarios = scenarios.json()
    for scenario in scenarios:
        real_path_to_scenario_results = path_to_scenario_results  # + f"{scenario['name']}_{scenario['id']}/"
        scenario_path = real_path_to_scenario_results + f"/{scenario['name']}_{scenario['id']}/"
        if not os.path.exists(scenario_path):
            log.error(f'Path {scenario_path} does not exist')
            continue
        if scenario["name"] == "casedata":
            percentiles = ['p50']
            zip_file = write_zip(path_to_data=real_path_to_scenario_results,
                                 zipped_name=f"{scenario['name']}_{scenario['id']}", percentiles=percentiles,
                                 case_data=True)
        else:
            percentiles = ['p25', 'p50', 'p75']
            zip_file = write_zip(path_to_data=real_path_to_scenario_results,
                                 zipped_name=f"{scenario['name']}_{scenario['id']}", percentiles=percentiles,
                                 case_data=False)
        # parse id
        # path_to_data = os.path.join(path_to_scenario_results)

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
