import requests

header = {'Authorization': "Bearer anythingAsPasswordIsFineCurrently"}

url = "http://localhost:8123/"

# For test: Get an arbitrary scenario id and upload arbitrary data
scenarios = requests.get(
            url + "scenarios/",  headers=header)

scenario_ids = [scenario["id"] for scenario in scenarios.json()]


scenario_id = scenario_ids[0]

put_response = requests.put(url + "scenarios/" + scenario_id, headers=header)

print(put_response.status_code)