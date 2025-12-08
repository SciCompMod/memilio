import json
import time
import datetime
import numpy as np
from enum import Enum
import os
import requests
import pandas as pd
import tempfile
import concurrent.futures
import subprocess


import memilio.simulation as mio
import memilio.simulation.osecirvvs as osecirvvs
from memilio.epidata import geoModificationGermany as geoger

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


class Simulation:

    def __init__(self, data_dir, results_dir, parameter_list, intervention_list, df_interventions):
        self.num_groups = 6
        self.data_dir = data_dir
        self.results_dir = results_dir
        self.intervention_list = []
        self.parameter_list = []

        node_ids = geoger.get_county_ids(True, True)
        path_vacc_data = os.path.join(
            self.data_dir, "Germany", "pydata", "vacc_county_ageinf_ma7.json")
        path_case_data = os.path.join(
            self.data_dir, "Germany", "pydata", "cases_all_county_age_ma7.json")
        path_population_data = os.path.join(
            self.data_dir, "Germany", "pydata", "county_current_population.json")
        self.vacc_data = osecirvvs.read_vaccination_data(path_vacc_data)
        self.case_data = osecirvvs.read_confirmed_cases_data(path_case_data)
        self.population_data = osecirvvs.read_population_data(
            path_population_data, node_ids)
        
        self.parameter_list = parameter_list
        self.intervention_list = intervention_list
        self.df_interventions = df_interventions

    def _process_scenario(self, scenario_data_run, num_runs):
        print("Processing scenario: ", scenario_data_run['name'])
        try:
            num_days_sim = (datetime.datetime.strptime(
                scenario_data_run['endDate'], "%Y-%m-%d") - datetime.datetime.strptime(scenario_data_run['startDate'], "%Y-%m-%d")).days

            # for testing
            # scenario_data_run['startDate'] = "2024-12-22"
            # scenario_data_run['endDate'] = "2024-01-01"

            extrapolate = False
            graph = self.get_graph(extrapolate, scenario_data_run)

            if extrapolate:
                return f"Processed casedata scenario {scenario_data_run['id']}"

            print(
                f"Starting simulation for scenario: {scenario_data_run['id']} ({scenario_data_run['name']})")
            study = osecirvvs.ParameterStudy(
                graph, 0., num_days_sim, 0.5, num_runs)
            ensemble = study.run(False)

            ensemble_results = []
            ensemble_params = []
            # ensemble_edges = []
            for run in range(num_runs):
                graph_run = ensemble[run]
                ensemble_results.append(
                    osecirvvs.interpolate_simulation_result(graph_run))
                ensemble_params.append(
                    [graph_run.get_node(node_indx).property.model
                     for node_indx in range(graph.num_nodes)])
                # ensemble_edges.append(
                #     [graph_run.get_edge(edge_indx).property.mobility_results
                #      for edge_indx in range(graph.num_edges)])

            node_ids = [graph.get_node(i).id for i in range(graph.num_nodes)]

            save_percentiles = True
            save_single_runs = False

            # save edge results
            # self.save_results_edges(graph, ensemble_edges, num_days_sim)

            osecirvvs.save_results(
                ensemble_results, ensemble_params, node_ids, self.results_dir,
                save_single_runs, save_percentiles, num_days_sim, True)
            return f"Successfully processed scenario {scenario_data_run['id']}"

        except Exception as e:
            print(
                f"Error processing scenario {scenario_data_run.get('id', 'N/A')}: {e}")
            return f"Failed to process scenario {scenario_data_run.get('id', 'N/A')}: {e}"

    def get_parameter_values(self, parameters, parameter_name):
        parameter = next(
            (entry for entry in parameters if entry['name'] == parameter_name), None)
        if parameter:
            min_value = parameter['values'][0]['valueMin']
            max_value = parameter['values'][0]['valueMax']
            return min_value, max_value
        print(f"Parameter {parameter_name} not found in parameters")
        return None, None

    def set_covid_parameters(self, model, scenario_data):
        parameters = scenario_data['modelParameters']

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
            scenario_data['startDate'], '%Y-%m-%d').date()
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

    def set_npis(self, params, scenario_data):
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
        interventions_scenario = scenario_data['linkedInterventions']
        coefficients = {}

        for intervention in interventions_scenario:
            # search intervention in self.intervention_list
            intervention_data = next(
                (entry for entry in self.intervention_list if entry['id'] == intervention['interventionId']), None)

            if not intervention_data:
                print(
                    f"Intervention {intervention['interventionId']} not found in intervention list")
                continue
            
            coefficients[intervention_data['name']] = intervention['coefficient']
        
        for row in self.df_interventions.rows:

            if 'SchoolClosure' in self.df_interventions.columns:
                coefficient = row.get('SchoolClosure') * coefficients['School closure']
                dampings.append(school_closure(
                    row.get('Time'), coefficient, coefficient))
            if 'HomeOffice' in self.df_interventions.columns:
                coefficient = row.get('HomeOffice') * coefficients['Remote work']
                dampings.append(home_office(
                    row.get('Time'), coefficient, coefficient))
            if 'PhysicalDistancingSchool' in self.df_interventions.columns:
                coefficient = row.get('PhysicalDistancingSchool') * coefficients['Face masks & social distancing School']
                dampings.append(physical_distancing_school(
                    row.get('Time'), coefficient, coefficient))
            if 'PhysicalDistancingWork' in self.df_interventions.columns:
                coefficient = row.get('PhysicalDistancingWork') * coefficients['Face masks & social distancing Work']
                dampings.append(physical_distancing_work(
                    row.get('Time'), coefficient, coefficient))
            if 'PhysicalDistancingOther' in self.df_interventions.columns:
                coefficient = row.get('PhysicalDistancingOther') * coefficients['Face masks & social distancing Other']
                dampings.append(physical_distancing_other(
                    row.get('Time'), coefficient, coefficient))
                
        params.ContactPatterns.dampings = dampings

    def get_edge_indices(self, graph):
        edge_indices = []
        pop_object = graph.get_node(0).property.populations
        infection_states_mild = [
            osecirvvs.InfectionState.ExposedNaive,
            osecirvvs.InfectionState.ExposedPartialImmunity,
            osecirvvs.InfectionState.ExposedImprovedImmunity,
            osecirvvs.InfectionState.InfectedNoSymptomsNaive,
            osecirvvs.InfectionState.InfectedNoSymptomsPartialImmunity,
            osecirvvs.InfectionState.InfectedNoSymptomsImprovedImmunity,
            osecirvvs.InfectionState.InfectedNoSymptomsNaiveConfirmed,
            osecirvvs.InfectionState.InfectedNoSymptomsPartialImmunityConfirmed,
            osecirvvs.InfectionState.InfectedNoSymptomsImprovedImmunityConfirmed,
            osecirvvs.InfectionState.InfectedSymptomsNaive,
            osecirvvs.InfectionState.InfectedSymptomsPartialImmunity,
            osecirvvs.InfectionState.InfectedSymptomsImprovedImmunity,
            osecirvvs.InfectionState.InfectedSymptomsNaiveConfirmed,
            osecirvvs.InfectionState.InfectedSymptomsPartialImmunityConfirmed,
            osecirvvs.InfectionState.InfectedSymptomsImprovedImmunityConfirmed,
        ]

        indicies_mild = []

        for age_group in range(self.num_groups):
            for infection_state in infection_states_mild:
                indicies_mild.append(pop_object.get_flat_index(osecirvvs.MultiIndex_PopulationsArray(
                    mio.AgeGroup(age_group), infection_state)))

        edge_indices.append(indicies_mild)

        return edge_indices

    def get_graph(self, extrapolate, scenario_data):
        model = osecirvvs.Model(self.num_groups)
        self.set_covid_parameters(model, scenario_data)
        self.set_contact_matrices(model)
        self.set_npis(model.parameters, scenario_data)

        graph = osecirvvs.ModelGraph()

        scaling_factor_infected = np.ones(self.num_groups)
        scaling_factor_icu = 1.0
        tnt_capacity_factor = 1.43 / 100000.

        path_population_data = os.path.join(
            self.data_dir, "Germany", "pydata",
            "county_current_population.json")

        # get start date in mio.Date format
        start_date = datetime.datetime.strptime(
            scenario_data['startDate'], '%Y-%m-%d').date()
        start_date = mio.Date(
            start_date.year, start_date.month, start_date.day)

        # get nums between start and end date
        end_date = datetime.datetime.strptime(
            scenario_data['endDate'], '%Y-%m-%d').date()
        end_date = mio.Date(
            end_date.year, end_date.month, end_date.day)

        osecirvvs.set_nodes_cached(
            model.parameters, start_date, end_date, self.data_dir, path_population_data,
            True, graph, scaling_factor_infected, scaling_factor_icu, tnt_capacity_factor,
            end_date - start_date - 1, extrapolate, True, self.case_data, self.population_data,
            self.vacc_data)

        # get indices for the edges
        # edge_indices = self.get_edge_indices(graph)
        edge_indices = []

        osecirvvs.set_edges(
            self.data_dir, graph, len(Location), edge_indices)

        return graph

    def save_results_edges(self, graph, ensemble_edges, num_days_sim):
        edge_indx_pair = {}
        for edge_indx in range(graph.num_edges):
            edge = graph.get_edge(edge_indx)
            edge_indx_pair[edge_indx] = (
                edge.start_node_idx, edge.end_node_idx)

        # save edges with one file per day
        # convert time series into numpy arrays
        for run_indx in range(len(ensemble_edges)):
            for edge_indx in range(len(ensemble_edges[0])):
                ensemble_edges[run_indx][edge_indx] = ensemble_edges[run_indx][edge_indx].as_ndarray(
                )

        percentiles = [25, 50, 75]
        # calc percentiles for each edge
        ensemble_edges_percentiles = []
        for edge_indx in range(len(ensemble_edges[0])):
            edge_percentiles = []
            for percentile in percentiles:
                edge_percentile = np.percentile(
                    [ensemble_edges[run_indx][edge_indx]
                     for run_indx in range(len(ensemble_edges))],
                    percentile, axis=0)
                edge_percentiles.append(edge_percentile)
            ensemble_edges_percentiles.append(edge_percentiles)

        # create one file per day
        for day in range(1, num_days_sim + 1):
            day_data = []
            for edge_indx in range(len(ensemble_edges_percentiles)):
                start_node_indx, end_node_indx = edge_indx_pair[edge_indx]
                start_node_id = graph.get_node(start_node_indx).id
                end_node_id = graph.get_node(end_node_indx).id
                for percentile_indx, percentile in enumerate(percentiles):
                    # since we use dt=0.5, we need to multiply the day by 2 to get the correct index
                    day_indx = (day * 2) - 1
                    day_data.append({
                        "day": day,
                        "start_node": start_node_id,
                        "end_node": end_node_id,
                        "percentile": percentile,
                        "mild_infected": int(ensemble_edges_percentiles[edge_indx][percentile_indx][1][day_indx]),
                        "total": int(ensemble_edges_percentiles[edge_indx][percentile_indx][2][day_indx])
                    })
            # save to csv
            df = pd.DataFrame(day_data)
            df.to_csv(
                f"{self.results_dir}/edges_day_{day}.csv", index=False)

    def run(self, scenario_data_run, num_runs=10):
        mio.set_log_level(mio.LogLevel.Warning)

        # self.header = {'Authorization': "Bearer anythingAsPasswordIsFineCurrently"}

        self._process_scenario(scenario_data_run, num_runs)

        print("Finished simulating results.")

class OptimalControlSimulation:

    def __init__(self, data_dir, results_dir, build_dir, run_data_url, headers):
        self.data_dir = data_dir
        self.results_dir = results_dir
        self.build_dir = build_dir
        self.run_data_url = run_data_url
        self.headers = headers
        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir)
    
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
                    # Run with local data
                    # with open(os.path.join(self.data_dir, "scenario_data_run.json"), "r") as file:
                    #     scenario_data_run = json.load(file)
                    break
                except requests.exceptions.RequestException as e:
                    if attempt < max_retries - 1:
                        print(
                            f"Retry {attempt+1}/{max_retries} for scenario {scenario['name']}. Error: {e}")
                        time.sleep(retry_delay * (2**attempt))
                    else:
                        raise
            
            temp_dir = tempfile.TemporaryDirectory()

            with open(os.path.join(temp_dir.name, "scenario_data_run.json"), "w") as outfile:
                json.dump(scenario_data_run, outfile)
            with open(os.path.join(temp_dir.name, "parameter_list.json"), "w") as outfile:
                json.dump(self.parameter_list, outfile)
            with open(os.path.join(temp_dir.name, "intervention_list.json"), "w") as outfile:
                json.dump(self.intervention_list, outfile)

            optimal_control_result_filename = "OptimalControl_output.txt"
            # run optimal control in c++
            run = subprocess.run(
                [
                    '{build_directory}/bin/simulate_optimal_control_scenario -DataDirectory \\"{temp_directory}\\" -OutputFileName \\"{filename}\\" -ConstraintInfectedCases {constraint:.0f}'.format(
                        build_directory = self.build_dir,
                        temp_directory = temp_dir.name, 
                        filename = optimal_control_result_filename, 
                        constraint = 1e5)
                ],
                shell=True,
                check=True,
                capture_output=False,
                text=True,
            )

            df_interventions = pd.read_csv(os.path.join(temp_dir.name, optimal_control_result_filename))
            
            res_dir_scenario = os.path.join(
                self.results_dir,
                f'{scenario_data_run["name"]}_{scenario_data_run["id"]}'
            )

            if not os.path.exists(res_dir_scenario):
                os.makedirs(res_dir_scenario)
                
            self.save_optimal_controls_description(res_dir_scenario, df_interventions)

            sim = Simulation(
                data_dir=self.data_dir,
                results_dir=res_dir_scenario,
                parameter_list=self.parameter_list, 
                intervention_list=self.intervention_list, 
                df_interventions=df_interventions)
            sim.run(scenario_data_run, num_runs=num_runs)

            temp_dir.cleanup()

        except Exception as e:
            print(
                f"Error processing optimal conntrol scenario {scenario.get('id', 'N/A')}: {e}")
            return f"Failed to process optimal conntrol scenario {scenario.get('id', 'N/A')}: {e}"

    def save_optimal_controls_description(self, res_dir_scenario, df_interventions):
        
        description = "Optimal control values found at:"

        for column_name, values in df_interventions.items():
            column_str = ' ' + column_name + ': ['
            column_str += ", ".join(str(x) for x in values)
            column_str += "];"
            description += column_str

        with open(os.path.join(res_dir_scenario, "scenario_description.txt"), "w") as text_file:
            text_file.write(description)


    def run(self, num_runs=10, max_workers=10):
        mio.set_log_level(mio.LogLevel.Warning)

        # header = {'Authorization': "Bearer anythingAsPasswordIsFineCurrently"}

        # list all scenarios
        scenarios = requests.get(
            self.run_data_url + "scenarios/", headers=self.headers).json()
        
        # Only use scenarios that are optimal control 
        scenarios[:] = (scenario for i, scenario in enumerate(scenarios) if "OptimalControl" in scenario["name"])

        self.parameter_list = requests.get(
            self.run_data_url + "parameterdefinitions/", headers=self.headers).json()

        # read intervention list
        self.intervention_list = requests.get(
            self.run_data_url + "interventions/templates/", headers=self.headers).json()
        
        # Run with local data
        # with open(os.path.join(self.data_dir, "scenarios.json"), "rb") as file:
        #     scenarios = json.load(file)
        #     scenarios[:] = [scenarios[0]]
        # with open(os.path.join(self.data_dir, "parameter_list.json"), "rb") as file:
        #     self.parameter_list = json.load(file)
        # with open(os.path.join(self.data_dir, "intervention_list.json"), "rb") as file:
        #     self.intervention_list = json.load(file)

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

        print("Finished optimization.")


if __name__ == "__main__":
    cwd = os.getcwd()
    run_data_url = "https://zam10063.zam.kfa-juelich.de/api-dev/"

    # Information needed for authentication.
    # service_realm = ""
    # client_id = ""
    # username = ""
    # password = ""

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
    headers = {"username": "MyVeryOwnHost"}

    run_data_url = "https://zam10063.zam.kfa-juelich.de/api-dev/"
    optsim = OptimalControlSimulation(
        data_dir=os.path.join(cwd, "data"),
        results_dir=os.path.join(cwd, "results_osecirvvs_optimization"),
        build_dir=os.path.join(cwd, "cpp/build"),
        run_data_url=run_data_url,
        headers=headers)  
    optsim.run(3, 1)