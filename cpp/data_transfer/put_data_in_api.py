import requests
import os
from zipfile import ZipFile

header = {'Authorization': "Bearer anythingAsPasswordIsFineCurrently"}


def write_zip(path_to_data, zipped_name, percentiles=['p50'], case_data=False):
    """! Zips scenario for put to backend.

    @param[in] path_to_data Path to data.
    @param[in] zipped_name name of zip. ".zip" is added to name.
    @param[in] percentiles. List[str]
    """
    zipfile = os.path.join(path_to_data, zipped_name + ".zip", "w")
    with ZipFile(zipfile) as zip_object:
        for percentile in percentiles:
            if not case_data:
                real_path_to_data = path_to_data + percentile
            else:
                real_path_to_data = path_to_data
            zip_object.write(os.path.join(
                real_path_to_data, "Results.h5"), f"{percentile}/Results.h5")
            zip_object.write(os.path.join(
                real_path_to_data, "Results_sum.h5"), f"{percentile}/Results_sum.h5")
    return zipfile


def put_scenario(scenario_id, zip_file, url):
    # https://stackoverflow.com/questions/18208109/python-requests-post-a-zip-file-with-multipart-form-data
    fileobj = open(zip_file, 'rb')
    put_response = requests.put(url + "scenarios/" + scenario_id, headers=header,
                                files={"archive", zip_file, fileobj})
    print(put_response.status_code)

def put_scenarios(path_to_scenario_results, url):
    scenarios = requests.get(url + "scenarios/", headers=header).json()
    for scenario in scenarios:
        # parse id
        path_to_data = os.path.join(path_to_scenario_results, f"scenario['name']_scenario['id']")
        zip_file = write_zip(path_to_data=path_to_data)
        put_scenario(scenario['id'], zip_file=zip_file, url=url)
        
