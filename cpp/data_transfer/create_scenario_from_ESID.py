# TODO
#
# for each scenario we get one scenario.json:
# - read in scenario.json
# - read in model.json (sevirvvs.json)
# - set parameters from scenario.json
#      - from model.json: for each group, check if group is of category AgeGroup and accumulate to find out number of AgeGroups
#      - for each parameterID check corresponding model parameter
#      - then check groupID for corresponding model group
#      - set value(s)
# - set NPIs from scenario.json
#      - for each interventionID check intervention name and description
#      - apply intervention to the model with values
# - create graph model
# - initialize model as usual with data from RKI or other data
# - run simulation
# - post results to ESID backend with new API
#
import json
import time
import datetime
import numpy as np
from enum import Enum
import os
import requests
from itertools import combinations


import memilio.simulation as mio
import memilio.simulation.osecirvvs as osecirvvs


class Location(Enum):
    Home = 0
    School = 1
    Work = 2
    Other = 3


class Intervention(Enum):
    Home = 0
    SchoolClosure = 1
    HomeOffice = 2
    GatheringBanFacilitiesClosure = 3
    PhysicalDistanceAndMasks = 4
    SeniorAwareness = 5


class InterventionLevel(Enum):
    Main = 0
    PhysicalDistanceAndMasks = 1
    SeniorAwareness = 2
    Holidays = 3


def get_invervention_list(url, header):
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
            "coefficient": 0.
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

    return all_combinations


class Simulation:

    def __init__(self, data_dir, results_dir, run_data_url):
        self.num_groups = 6
        self.data_dir = data_dir
        self.results_dir = results_dir
        self.szenario_data = ""
        self.intervention_list = []
        self.parameter_list = []
        self.run_data_url = run_data_url
        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir)

    def get_parameter_values(self, parameters, parameter_name):
        parameter = next(
            (entry for entry in parameters if entry['name'] == parameter_name), None)
        if parameter:
            min_value = parameter['values'][0]['valueMin']
            max_value = parameter['values'][0]['valueMax']
            return min_value, max_value
        print(f"Parameter {parameter_name} not found in parameters")
        return None, None

    def set_covid_parameters(self, model):
        # read parameters given from the scenario
        # TODO: Fix when the scenario data is fixed
        parameters = self.scenario_data['modelParameters']

        # Create a mapping from parameterId to name
        id_to_name = {entry['id']: entry['name']
                      for entry in self.parameter_list}

        # Add the corresponding name to each entry in parameters
        # get min, max values for each parameter
        parameter_values = {}
        for parameter in parameters:
            parameter['name'] = id_to_name.get(
                parameter['parameterId'], "Unknown")
            param_name = parameter['name']
            min_value, max_value = self.get_parameter_values(
                parameters, param_name)

            # Store values in the dictionary with dynamically generated keys
            parameter_values[f"{param_name}Min"] = min_value
            parameter_values[f"{param_name}Max"] = max_value

        def array_assign_uniform_distribution(param, min, max, num_groups=6):
            if isinstance(
                    min, (int, float)) and isinstance(
                        max, (int, float)):
                min = min * np.ones(num_groups)
                max = max * np.ones(num_groups)
            elif not (isinstance(min, (list)) and isinstance(
                    max, (list))):
                raise TypeError(
                    "Invalid type for parameter 'min' or 'max. \
                            Expected a scalar or a list. Must be the same for both.")
            for i in range(num_groups):
                param[mio.AgeGroup(i)] = mio.UncertainValue(
                    0.5 * (max[i] + min[i]))
                param[mio.AgeGroup(i)].set_distribution(
                    mio.ParameterDistributionUniform(min[i], max[i]))

        # if the parameter is not found in the scenario data, we use the default values
        # times
        timeExposedMin = parameter_values.get("TimeExposedMin", 2.67)
        timeExposedMax = parameter_values.get("TimeExposedMax", 4.)
        timeInfectedNoSymptomsMin = parameter_values.get(
            "TimeInfectedNoSymptomsMin", 1.2)
        timeInfectedNoSymptomsMax = parameter_values.get(
            "TimeInfectedNoSymptomsMax", 2.53)
        timeInfectedSymptomsMin = parameter_values.get(
            "TimeInfectedSymptomsMin", [
                5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465])
        timeInfectedSymptomsMax = parameter_values.get(
            "TimeInfectedSymptomsMax", [
                8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085])
        timeInfectedSevereMin = parameter_values.get(
            "TimeInfectedSevereMin", [
                3.925, 3.925, 4.85, 6.4, 7.2, 9.])
        timeInfectedSevereMax = parameter_values.get(
            "TimeInfectedSevereMax", [6.075, 6.075, 7., 8.7, 9.8, 13.])
        timeInfectedCriticalMin = parameter_values.get(
            "TimeInfectedCriticalMin", [4.95, 4.95, 4.86, 14.14, 14.4, 10.])
        timeInfectedCriticalMax = parameter_values.get(
            "TimeInfectedCriticalMax", [8.95, 8.95, 8.86, 20.58, 19.8, 13.2])

        array_assign_uniform_distribution(
            model.parameters.TimeExposed, timeExposedMin, timeExposedMax)

        array_assign_uniform_distribution(
            model.parameters.TimeInfectedNoSymptoms, timeInfectedNoSymptomsMin,
            timeInfectedNoSymptomsMax)

        array_assign_uniform_distribution(
            model.parameters.TimeInfectedSymptoms, timeInfectedSymptomsMin,
            timeInfectedSymptomsMax)

        array_assign_uniform_distribution(
            model.parameters.TimeInfectedSevere, timeInfectedSevereMin,
            timeInfectedSevereMax)

        array_assign_uniform_distribution(
            model.parameters.TimeInfectedCritical, timeInfectedCriticalMin,
            timeInfectedCriticalMax)

        # probabilities
        variantFactor = 1.4
        transmissionProbabilityOnContactMin = parameter_values.get(
            "TransmissionProbabilityOnContactMin",
            [
                0.02 * variantFactor, 0.05 * variantFactor, 0.05 * variantFactor,
                0.05 * variantFactor, 0.08 * variantFactor, 0.1 * variantFactor
            ]
        )
        transmissionProbabilityOnContactMax = parameter_values.get(
            "TransmissionProbabilityOnContactMax",
            [
                0.04 * variantFactor, 0.07 * variantFactor, 0.07 * variantFactor,
                0.07 * variantFactor, 0.10 * variantFactor, 0.15 * variantFactor
            ]
        )
        relativeTransmissionNoSymptomsMin = parameter_values.get(
            "RelativeTransmissionNoSymptomsMin", 0.5
        )
        relativeTransmissionNoSymptomsMax = parameter_values.get(
            "RelativeTransmissionNoSymptomsMax", 0.5
        )

        # The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
        # depends on incidence and test and trace capacity
        riskOfInfectionFromSymptomaticMin = parameter_values.get(
            "RiskOfInfectionFromSymptomaticMin", 0.0
        )
        riskOfInfectionFromSymptomaticMax = parameter_values.get(
            "RiskOfInfectionFromSymptomaticMax", 0.2
        )
        maxRiskOfInfectionFromSymptomaticMin = parameter_values.get(
            "MaxRiskOfInfectionFromSymptomaticMin", 0.4
        )
        maxRiskOfInfectionFromSymptomaticMax = parameter_values.get(
            "MaxRiskOfInfectionFromSymptomaticMax", 0.5
        )
        recoveredPerInfectedNoSymptomsMin = parameter_values.get(
            "RecoveredPerInfectedNoSymptomsMin", [
                0.2, 0.2, 0.15, 0.15, 0.15, 0.15]
        )
        recoveredPerInfectedNoSymptomsMax = parameter_values.get(
            "RecoveredPerInfectedNoSymptomsMax", [
                0.3, 0.3, 0.25, 0.25, 0.25, 0.25]
        )
        severePerInfectedSymptomsMin = parameter_values.get(
            "SeverePerInfectedSymptomsMin", [
                0.006, 0.006, 0.015, 0.049, 0.15, 0.20]
        )
        severePerInfectedSymptomsMax = parameter_values.get(
            "SeverePerInfectedSymptomsMax", [
                0.009, 0.009, 0.023, 0.074, 0.18, 0.25]
        )
        criticalPerSevereMin = parameter_values.get(
            "CriticalPerSevereMin", [0.05, 0.05, 0.05, 0.10, 0.25, 0.35]
        )
        criticalPerSevereMax = parameter_values.get(
            "CriticalPerSevereMax", [0.10, 0.10, 0.10, 0.20, 0.35, 0.45]
        )
        deathsPerCriticalMin = parameter_values.get(
            "DeathsPerCriticalMin", [0.00, 0.00, 0.10, 0.10, 0.30, 0.5]
        )
        deathsPerCriticalMax = parameter_values.get(
            "DeathsPerCriticalMax", [0.10, 0.10, 0.18, 0.18, 0.50, 0.7]
        )
        reducExposedPartialImmunityMin = parameter_values.get(
            "ReducedExposedPartialImmunityMin", 0.75
        )
        reducExposedPartialImmunityMax = parameter_values.get(
            "ReducedExposedPartialImmunityMax", 0.85
        )
        reducExposedImprovedImmunityMin = parameter_values.get(
            "ReducedExposedImprovedImmunityMin", 0.281
        )
        reducExposedImprovedImmunityMax = parameter_values.get(
            "ReducedExposedImprovedImmunityMax", 0.381
        )
        reducInfectedSymptomsPartialImmunityMin = parameter_values.get(
            "ReducedInfectedSymptomsPartialImmunityMin", 0.6
        )
        reducInfectedSymptomsPartialImmunityMax = parameter_values.get(
            "ReducedInfectedSymptomsPartialImmunityMax", 0.7
        )
        reducInfectedSymptomsImprovedImmunityMin = parameter_values.get(
            "ReducedInfectedSymptomsImprovedImmunityMin", 0.193
        )
        reducInfectedSymptomsImprovedImmunityMax = parameter_values.get(
            "ReducedInfectedSymptomsImprovedImmunityMax", 0.293
        )
        reducInfectedSevereCriticalDeadPartialImmunityMin = parameter_values.get(
            "ReducedInfectedSevereCriticalDeadPartialImmunityMin", 0.05
        )
        reducInfectedSevereCriticalDeadPartialImmunityMax = parameter_values.get(
            "ReducedInfectedSevereCriticalDeadPartialImmunityMax", 0.15
        )
        reducInfectedSevereCriticalDeadImprovedImmunityMin = parameter_values.get(
            "ReducedInfectedSevereCriticalDeadImprovedImmunityMin", 0.041
        )
        reducInfectedSevereCriticalDeadImprovedImmunityMax = parameter_values.get(
            "ReducedInfectedSevereCriticalDeadImprovedImmunityMax", 0.141
        )
        reducTimeInfectedMild = parameter_values.get(
            "ReducedTimeInfectedMild", 1.0)

        array_assign_uniform_distribution(
            model.parameters.TransmissionProbabilityOnContact,
            transmissionProbabilityOnContactMin,
            transmissionProbabilityOnContactMax)

        array_assign_uniform_distribution(
            model.parameters.RelativeTransmissionNoSymptoms,
            relativeTransmissionNoSymptomsMin,
            relativeTransmissionNoSymptomsMax)

        array_assign_uniform_distribution(
            model.parameters.RiskOfInfectionFromSymptomatic,
            riskOfInfectionFromSymptomaticMin,
            riskOfInfectionFromSymptomaticMax)

        array_assign_uniform_distribution(
            model.parameters.MaxRiskOfInfectionFromSymptomatic,
            maxRiskOfInfectionFromSymptomaticMin,
            maxRiskOfInfectionFromSymptomaticMax)

        array_assign_uniform_distribution(
            model.parameters.RecoveredPerInfectedNoSymptoms,
            recoveredPerInfectedNoSymptomsMin,
            recoveredPerInfectedNoSymptomsMax)

        array_assign_uniform_distribution(
            model.parameters.SeverePerInfectedSymptoms,
            severePerInfectedSymptomsMin,
            severePerInfectedSymptomsMax)

        array_assign_uniform_distribution(
            model.parameters.CriticalPerSevere,
            criticalPerSevereMin,
            criticalPerSevereMax)

        array_assign_uniform_distribution(
            model.parameters.DeathsPerCritical,
            deathsPerCriticalMin,
            deathsPerCriticalMax)

        array_assign_uniform_distribution(
            model.parameters.ReducExposedPartialImmunity,
            reducExposedPartialImmunityMin,
            reducExposedPartialImmunityMax)

        array_assign_uniform_distribution(
            model.parameters.ReducExposedImprovedImmunity,
            reducExposedImprovedImmunityMin,
            reducExposedImprovedImmunityMax)

        array_assign_uniform_distribution(
            model.parameters.ReducInfectedSymptomsPartialImmunity,
            reducInfectedSymptomsPartialImmunityMin,
            reducInfectedSymptomsPartialImmunityMax)

        array_assign_uniform_distribution(
            model.parameters.ReducInfectedSymptomsImprovedImmunity,
            reducInfectedSymptomsImprovedImmunityMin,
            reducInfectedSymptomsImprovedImmunityMax)

        array_assign_uniform_distribution(
            model.parameters.ReducInfectedSevereCriticalDeadPartialImmunity,
            reducInfectedSevereCriticalDeadPartialImmunityMin,
            reducInfectedSevereCriticalDeadPartialImmunityMax)

        array_assign_uniform_distribution(
            model.parameters.ReducInfectedSevereCriticalDeadImprovedImmunity,
            reducInfectedSevereCriticalDeadImprovedImmunityMin,
            reducInfectedSevereCriticalDeadImprovedImmunityMax)

        array_assign_uniform_distribution(
            model.parameters.ReducTimeInfectedMild,
            reducTimeInfectedMild,
            reducTimeInfectedMild)

        # start day is set to the n-th day of the year
        start_date = datetime.datetime.strptime(
            self.scenario_data['startDate'], '%Y-%m-%d').date()
        start_day = mio.Date(
            start_date.year, start_date.month, start_date.day).day_in_year

        model.parameters.StartDay = parameter_values.get(
            "StartDay", start_day)

        model.parameters.Seasonality = mio.UncertainValue(
            0.5 * (parameter_values.get(
                "SeasonalityMax", 0.3
            ) + parameter_values.get(
                "SeasonalityMin", 0.1
            )))
        model.parameters.Seasonality.set_distribution(
            mio.ParameterDistributionUniform(parameter_values.get(
                "SeasonalityMin", 0.1
            ), parameter_values.get(
                "SeasonalityMax", 0.3
            )))

        model.parameters.StartDayNewVariant = parameter_values.get(
            "StartDayNewVariant", mio.Date(
                start_date.year, start_date.month, start_date.day).day_in_year
        )

    def set_contact_matrices(self, model):
        contact_matrices = mio.ContactMatrixGroup(
            len(list(Location)), self.num_groups)
        locations = ["home", "school_pf_eig", "work", "other"]

        for i, location in enumerate(locations):
            baseline_file = os.path.join(
                self.data_dir, "contacts", "baseline_" + location + ".txt")
            contact_matrices[i] = mio.ContactMatrix(
                mio.read_mobility_plain(baseline_file),
                np.zeros((self.num_groups, self.num_groups))
            )
        model.parameters.ContactPatterns.cont_freq_mat = contact_matrices

    def set_npis(self, params):
        # set up for damping
        contacts = params.ContactPatterns
        dampings = contacts.dampings

        group_weights_all = np.ones(self.num_groups)
        group_weights_seniors = np.vectorize(
            lambda i: 1.0 if i == 5 else (0.5 if i == 4 else 0.0))(
            np.arange(self.num_groups))

        loc_home = Location.Home.value
        loc_school = Location.School.value
        loc_work = Location.Work.value
        loc_other = Location.Other.value

        lvl_main = InterventionLevel.Main.value
        lvl_pd_and_masks = InterventionLevel.PhysicalDistanceAndMasks.value
        lvl_seniors = InterventionLevel.SeniorAwareness.value
        lvl_holidays = InterventionLevel.Holidays.value

        typ_home = Intervention.Home.value
        typ_school = Intervention.SchoolClosure.value
        typ_homeoffice = Intervention.HomeOffice.value
        typ_gathering = Intervention.GatheringBanFacilitiesClosure.value
        typ_distance = Intervention.PhysicalDistanceAndMasks.value
        typ_senior = Intervention.SeniorAwareness.value

        def damping_helper(
                t, min, max, damping_level, type, location,
                group_weights=group_weights_all):
            v = mio.UncertainValue(0.5 * (max + min))
            v.set_distribution(mio.ParameterDistributionUniform(min, max))
            return mio.DampingSampling(
                value=v,
                level=damping_level,
                type=type,
                time=t,
                matrix_indices=location,
                group_weights=group_weights)

        def school_closure(t, min, max):
            return damping_helper(
                t, min, max, lvl_main, typ_school, [loc_school])

        def home_office(t, min, max):
            return damping_helper(
                t, min, max, lvl_main, typ_homeoffice, [loc_work])

        def physical_distancing_school(t, min, max):
            return damping_helper(
                t, min, max, lvl_pd_and_masks, typ_distance, [loc_school])

        def physical_distancing_work(t, min, max):
            return damping_helper(
                t, min, max, lvl_pd_and_masks, typ_distance, [loc_work])

        def physical_distancing_other(t, min, max):
            return damping_helper(
                t, min, max, lvl_pd_and_masks, typ_distance, [loc_other])

        # read interventions given from the scenario
        interventions_scenario = self.scenario_data['linkedInterventions']

        for intervention in interventions_scenario:
            # search intervention in self.intervention_list
            intervention_data = next(
                (entry for entry in self.intervention_list if entry['id'] == intervention['interventionId']), None)

            if not intervention_data:
                print(
                    f"Intervention {intervention['interventionId']} not found in intervention list")
                continue

            start_date_day = 0

            if intervention_data['name'] == 'School closure':
                dampings.append(school_closure(
                    start_date_day, intervention['coefficient'], intervention['coefficient']))
            elif intervention_data['name'] == 'Remote work':
                dampings.append(home_office(
                    start_date_day, intervention['coefficient'], intervention['coefficient']))
            elif intervention_data['name'] == "Face masks & social distancing School":
                dampings.append(physical_distancing_school(
                    start_date_day, intervention['coefficient'], intervention['coefficient']))
            elif intervention_data['name'] == "Face masks & social distancing Work":
                dampings.append(physical_distancing_work(
                    start_date_day, intervention['coefficient'], intervention['coefficient']))
            elif intervention_data['name'] == "Face masks & social distancing Other":
                dampings.append(physical_distancing_other(
                    start_date_day, intervention['coefficient'], intervention['coefficient']))
            else:
                print(
                    f"Intervention {intervention_data['name']} not implemented yet!")

        params.ContactPatterns.dampings = dampings

    def get_graph(self, extrapolate):
        model = osecirvvs.Model(self.num_groups)
        self.set_covid_parameters(model)
        self.set_contact_matrices(model)
        self.set_npis(model.parameters)

        graph = osecirvvs.ModelGraph()

        scaling_factor_infected = np.ones(self.num_groups)
        scaling_factor_icu = 1.0
        tnt_capacity_factor = 1.43 / 100000.

        path_population_data = os.path.join(
            self.data_dir, "pydata", "Germany",
            "county_current_population.json")

        # get start date in mio.Date format
        start_date = datetime.datetime.strptime(
            self.scenario_data['startDate'], '%Y-%m-%d').date()
        start_date = mio.Date(
            start_date.year, start_date.month, start_date.day)

        # get nums between start and end date
        end_date = datetime.datetime.strptime(
            self.scenario_data['endDate'], '%Y-%m-%d').date()
        end_date = mio.Date(
            end_date.year, end_date.month, end_date.day)

        osecirvvs.set_nodes(
            model.parameters,
            start_date,
            end_date, self.data_dir,
            path_population_data, True, graph, scaling_factor_infected,
            scaling_factor_icu, tnt_capacity_factor, end_date - start_date - 1, extrapolate)

        osecirvvs.set_edges(
            self.data_dir, graph, len(Location))

        return graph

    def run(self, num_days_sim, num_runs=10):
        mio.set_log_level(mio.LogLevel.Warning)

        header = {'Authorization': "Bearer anythingAsPasswordIsFineCurrently"}

        # list all scenarios
        scenarios = requests.get(
            self.run_data_url + "scenarios/", headers=header).json()

        self.parameter_list = requests.get(
            self.run_data_url + "parameterdefinitions/", headers=header).json()

        # read intervention list
        self.intervention_list = requests.get(
            self.run_data_url + "interventions/templates/", headers=header).json()

        # TODO: Assuming all parameters are equal for each scenarios.
        # Therefore, we only need to build the graph once?

        for scenario in scenarios:
            extrapolate = False
            if scenario['name'] == 'casedata':
                extrapolate = True

            # load full set of scenario data
            self.scenario_data = requests.get(
                self.run_data_url + "scenarios/" + scenario['id'], headers=header).json()

            # for testing overwrite startDate and endDate
            # self.scenario_data['startDate'] = "2024-12-22"
            # self.scenario_data['endDate'] = "2024-01-01"

            graph = self.get_graph(extrapolate)

            # in the casedata scenario, we are just interested in the extrapolated data. Therefore, we skip the simulation.
            if extrapolate:
                continue

            study = osecirvvs.ParameterStudy(
                graph, 0., num_days_sim, 0.5, num_runs)
            ensemble = study.run(False)

            ensemble_results = []
            ensemble_params = []
            for run in range(num_runs):
                graph_run = ensemble[run]
                ensemble_results.append(
                    osecirvvs.interpolate_simulation_result(graph_run))
                ensemble_params.append(
                    [graph_run.get_node(node_indx).property.model
                     for node_indx in range(graph.num_nodes)])

            node_ids = [graph.get_node(i).id for i in range(graph.num_nodes)]

            save_percentiles = True
            save_single_runs = False

            res_dir_scenario = os.path.join(
                self.results_dir,
                f'{scenario["name"]}_{scenario["id"]}'
            )

            # create directory if it does not exist
            if not os.path.exists(res_dir_scenario):
                os.makedirs(res_dir_scenario)

            osecirvvs.save_results(
                ensemble_results, ensemble_params, node_ids, res_dir_scenario,
                save_single_runs, save_percentiles, num_days_sim, True)


if __name__ == "__main__":
    cwd = os.getcwd()
    run_data_url = "http://localhost:8123/"
    # with open(scenario_data_path) as f:
    #     scenario_data = json.load(f)
    mcmc_dir = os.path.join(cwd, "mcmc data")
    date_today = '2025-01-13'
    sim = Simulation(
        data_dir=os.path.join(cwd, "data"),
        results_dir=os.path.join(cwd, "results_osecirvvs"), run_data_url=run_data_url)
    sim.run(num_days_sim=30, num_runs=3)
