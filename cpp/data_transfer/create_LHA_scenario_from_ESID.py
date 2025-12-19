import os
import subprocess
import json
import requests
import datetime
import tempfile
import time
import pandas as pd
import concurrent.futures


class LHASimulation:
    def __init__(self, data_dir, results_dir, build_dir, run_data_url, headers):
        self.data_dir = data_dir
        self.results_dir = results_dir
        self.build_dir = build_dir
        self.run_data_url = run_data_url
        self.headers = headers
        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir)

    def preprocess_lha_data(self, lha_data_filename):
        file = self.data_dir + f"/Germany/pydata/{lha_data_filename}.csv"
        df = pd.read_csv(file, sep=",")
        # TODO: Adjust this to new dataset structure, probably first 5 columns that need to be deleted.
        # Drop first five columns as they contain informaion about the data set and are not needed in for the initialization of the ODE model
        df.drop(df.columns[[0, 1, 2, 3, 4]], axis=1)
        df.to_csv(file, sep=";")

    def _process_scenario(self, scenario, num_runs):
        scenario_data_run = None
        print("Processing scenario: ", scenario['name'])
        try:
            max_retries = 100
            retry_delay = 2

            for attempt in range(max_retries):
                try:
                    scenario_data_run = requests.get(
                        self.run_data_url + "scenarios/" + scenario['id'], headers=self.headers).json()
                    break
                except requests.exceptions.RequestException as e:
                    if attempt < max_retries - 1:
                        print(
                            f"Retry {attempt+1}/{max_retries} for scenario {scenario['name']}. Error: {e}")
                        time.sleep(retry_delay * (2**attempt))
                    else:
                        raise

            temp_dir = tempfile.TemporaryDirectory()

            scenario_id = scenario["id"]

            # Get scenario data and list of parameters from scenario information and stroe this in temp_dir.
            # We can access these files during the simulation, afterwards they will be deleted.
            # Then we get the scenario data and save it as json.
            response_scenario = requests.get(
                self.run_data_url + f"scenarios/{scenario_id}/", headers=self.headers)
            with open(os.path.join(temp_dir.name, "scenario_data_run_lha.json"), 'w') as f:
                json.dump(response_scenario.json(), f, ensure_ascii=False)

            # Get the parameter list. This is the same for all scenarios.
            response_params = requests.get(
                self.run_data_url + f"parameterdefinitions/", headers=self.headers)
            with open(os.path.join(temp_dir.name, "parameter_list_lha.json"), 'w') as f:
                json.dump(response_params.json(), f, ensure_ascii=False)

            # Get start_date of scenario
            start_date = datetime.datetime.strptime(
                scenario['startDate'], '%Y-%m-%d').date()

            num_days_sim = (datetime.datetime.strptime(
                scenario['endDate'],
                "%Y-%m-%d") - datetime.datetime.strptime(
                scenario['startDate'],
                "%Y-%m-%d")).days

            # Get LHA ids for which there is a scenario in database. For now, this is set artificially to the ID of Cologne.
            # TODO: This needs to be adapted so that we obtain the node id of the LHA from the creatorID in the scenario information.
            lha_ids = [5314]

            # Run scenario for each LHA ID separately as we assume that LHAs don't share data for now.
            for lha_id in lha_ids:
                print(f"Computing scenario for LHA {lha_id}")

                # Check if data is available for considered LHA.

                # If yes, preprocess dataset and give LHA id to run function.
                # TODO: Check if LHA data is available and extract to data_dir.
                # TODO: Get dataset id and write into description
                filename = f"lha_{lha_id}_observations_{start_date.year:04d}-{start_date.month:02d}-{start_date.day:02d}"
                # Preprocess LHA data, i.e. remove columns that describe dataset and do not contain data from LHA
                # TODO: This function needs to be adapted to final dataset structure.
                # preprocess_lha_data(data_dir, filename)

                # Run scenario

                result_dir_lha = f"{self.results_dir}/lha_{lha_id}"
                year = start_date.year
                month = start_date.month
                day = start_date.day

                run = subprocess.run(
                    [
                        f"{self.build_dir}/bin/simulate_lha_scenario {self.data_dir} {result_dir_lha} {temp_dir.name} {str(lha_id)} {str(year)} {str(month)} {str(day)} {str(num_days_sim)} {str(num_runs)}"
                    ],
                    shell=True,
                    check=True,
                    capture_output=True,
                    text=True,
                )
                print("STDOUT:")
                print(run.stdout)
                print("STDERR:")
                print(run.stderr)

                temp_dir.cleanup()

        except Exception as e:
            print(
                f"Error processing LHA scenario {scenario.get('id', 'N/A')}: {e}")
            return f"Failed to process LHA scenario {scenario.get('id', 'N/A')}: {e}"

    def run(self, num_runs=10, max_workers=10):

        # list all scenarios
        scenarios = requests.get(
            self.run_data_url + "scenarios/", headers=self.headers).json()

        # Only use scenarios that are optimal control
        scenarios[:] = (scenario for i, scenario in enumerate(
            scenarios) if "LHA" in scenario["name"])

        print(
            f"Processing {len(scenarios)} scenarios with {max_workers} workers.")

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(self._process_scenario, scenario, num_runs)
                       for scenario in scenarios]

            for future in concurrent.futures.as_completed(futures):
                try:
                    result = future.result()
                    print(result)
                except Exception as exc:
                    print(f'Exception: {exc}')

        print("Finished LHA simulation.")


def main():

    cwd = os.getcwd()
    run_data_url = "https://zam10063.zam.kfa-juelich.de/api-dev/"

    service_realm = ""
    client_id = ""
    username = ""
    password = ""

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

    sim = LHASimulation(data_dir=os.path.join(cwd, "data"), results_dir=os.path.join(
        cwd, "results_lhas"), build_dir=os.path.join(cwd, "build"), run_data_url=run_data_url, headers=headers)
    sim.run(num_runs=3, max_workers=1)


if __name__ == "__main__":
    main()
