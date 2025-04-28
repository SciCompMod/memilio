#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Henrik Zunker
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
import numpy as np
import datetime
import os
import memilio.simulation as mio
import memilio.simulation.osecir as osecir
import memilio.plot.createGIF as mp

from enum import Enum
from memilio.simulation.osecir import (Model, Simulation,
                                       interpolate_simulation_result)


class Location(Enum):
    """ """
    Home = 0
    School = 1
    Work = 2
    Other = 3


class Intervention(Enum):
    """ """
    Home = 0
    SchoolClosure = 1
    HomeOffice = 2
    GatheringBanFacilitiesClosure = 3
    PhysicalDistanceAndMasks = 4
    SeniorAwareness = 5


class InterventionLevel(Enum):
    """ """
    Main = 0
    PhysicalDistanceAndMasks = 1
    SeniorAwareness = 2
    Holidays = 3


class Simulation:
    """ """

    def __init__(self, data_dir, start_date, results_dir):
        self.num_groups = 6
        self.data_dir = data_dir
        self.start_date = start_date
        self.results_dir = results_dir
        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir)

    def set_covid_parameters(self, model):
        """

        :param model: 

        """
        def array_assign_uniform_distribution(param, min, max, num_groups=6):
            """

            :param param: 
            :param min: 
            :param max: 
            :param num_groups:  (Default value = 6)

            """
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

        timeExposedMin = 2.67
        timeExposedMax = 4.
        timeInfectedNoSymptomsMin = 1.2
        timeInfectedNoSymptomsMax = 2.53
        timeInfectedSymptomsMin = [
            5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465]
        timeInfectedSymptomsMax = [8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085]
        timeInfectedSevereMin = [3.925, 3.925, 4.85, 6.4, 7.2, 9.]
        timeInfectedSevereMax = [6.075, 6.075, 7., 8.7, 9.8, 13.]
        timeInfectedCriticalMin = [4.95, 4.95, 4.86, 14.14, 14.4, 10.]
        timeInfectedCriticalMax = [8.95, 8.95, 8.86, 20.58, 19.8, 13.2]

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
        transmissionProbabilityOnContactMin = [
            0.02, 0.05, 0.05, 0.05, 0.08, 0.15]
        transmissionProbabilityOnContactMax = [
            0.04, 0.07, 0.07, 0.07, 0.10, 0.20]
        relativeTransmissionNoSymptomsMin = 1
        relativeTransmissionNoSymptomsMax = 1

        # The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
        # depends on incidence and test and trace capacity
        riskOfInfectionFromSymptomaticMin = 0.1
        riskOfInfectionFromSymptomaticMax = 0.3
        maxRiskOfInfectionFromSymptomaticMin = 0.3
        maxRiskOfInfectionFromSymptomaticMax = 0.5
        recoveredPerInfectedNoSymptomsMin = [0.2, 0.2, 0.15, 0.15, 0.15, 0.15]
        recoveredPerInfectedNoSymptomsMax = [0.3, 0.3, 0.25, 0.25, 0.25, 0.25]
        severePerInfectedSymptomsMin = [0.006, 0.006, 0.015, 0.049, 0.15, 0.20]
        severePerInfectedSymptomsMax = [0.009, 0.009, 0.023, 0.074, 0.18, 0.25]
        criticalPerSevereMin = [0.05, 0.05, 0.05, 0.10, 0.25, 0.35]
        criticalPerSevereMax = [0.10, 0.10, 0.10, 0.20, 0.35, 0.45]
        deathsPerCriticalMin = [0.00, 0.00, 0.10, 0.10, 0.30, 0.5]
        deathsPerCriticalMax = [0.10, 0.10, 0.18, 0.18, 0.50, 0.7]

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

        # start day is set to the n-th day of the year
        model.parameters.StartDay = self.start_date.timetuple().tm_yday

        model.parameters.Seasonality.value = 0.2

        seasonality_min = 0.1
        seasonality_max = 0.3
        model.parameters.Seasonality = mio.UncertainValue(
            0.5 * (seasonality_max + seasonality_min))
        model.parameters.Seasonality.set_distribution(
            mio.ParameterDistributionUniform(seasonality_min, seasonality_max))

    def set_contact_matrices(self, model):
        """

        :param model: 

        """
        contact_matrices = mio.ContactMatrixGroup(
            len(list(Location)), self.num_groups)
        locations = ["home", "school_pf_eig", "work", "other"]

        for i, location in enumerate(locations):
            baseline_file = os.path.join(
                self.data_dir, "Germany", "contacts", "baseline_" + location + ".txt")
            minimum_file = os.path.join(
                self.data_dir, "Germany", "contacts", "minimum_" + location + ".txt")
            contact_matrices[i] = mio.ContactMatrix(
                mio.read_mobility_plain(baseline_file),
                mio.read_mobility_plain(minimum_file)
            )
        model.parameters.ContactPatterns.cont_freq_mat = contact_matrices

    def set_npis(self, params, end_date):
        """

        :param params: 
        :param end_date: 

        """
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
            """

            :param t: 
            :param min: 
            :param max: 
            :param damping_level: 
            :param type: 
            :param location: 
            :param group_weights:  (Default value = group_weights_all)

            """
            v = mio.UncertainValue(0.5 * (max + min))
            v.set_distribution(mio.ParameterDistributionUniform(min, max))
            return mio.DampingSampling(
                value=v,
                level=damping_level,
                type=type,
                time=t,
                matrix_indices=location,
                group_weights=group_weights)

        def contacts_at_home(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_main, typ_home, [loc_home])

        def school_closure(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_main, typ_school, [loc_school])

        def home_office(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_main, typ_homeoffice, [loc_work])

        def social_events(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_main, typ_gathering, [loc_other])

        def social_events_work(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_main, typ_gathering, [loc_work])

        def physical_distancing_home_school(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_pd_and_masks, typ_distance, [loc_home, loc_school])

        def physical_distancing_work_other(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_pd_and_masks, typ_distance, [loc_work, loc_other])

        def senior_awareness(t, min, max):
            """

            :param t: 
            :param min: 
            :param max: 

            """
            return damping_helper(
                t, min, max, lvl_seniors, typ_senior, [loc_home, loc_other],
                group_weights_seniors)

        # SPRING 2020 LOCKDOWN SCENARIO
        start_spring_date = datetime.date(
            2020, 3, 18)

        if start_spring_date < end_date:
            start_spring = (start_spring_date - self.start_date).days
            dampings.append(contacts_at_home(start_spring, 0.6, 0.8))
            dampings.append(school_closure(start_spring, 1.0, 1.0))
            dampings.append(home_office(start_spring, 0.2, 0.3))
            dampings.append(social_events(start_spring, 0.6, 0.8))
            dampings.append(social_events_work(start_spring, 0.1, 0.2))
            dampings.append(
                physical_distancing_home_school(
                    start_spring, 0.4, 0.6))
            dampings.append(
                physical_distancing_work_other(
                    start_spring, 0.4, 0.6))
            dampings.append(senior_awareness(start_spring, 0.0, 0.0))

        # SUMMER 2020 SCENARIO
        start_summer_date = datetime.date(year=2020, month=5, day=15)
        if start_summer_date < end_date:
            start_summer = (start_summer_date - self.start_date).days
            school_reopen_time = (datetime.date(
                2020, 6, 15) - self.start_date).days
            dampings.append(contacts_at_home(start_summer, 0.0, 0.2))
            # schools partially reopened
            dampings.append(school_closure(start_summer, 0.5, 0.5))
            # school fully reopened
            dampings.append(school_closure(school_reopen_time, 0.0, 0.0))
            dampings.append(home_office(start_summer, 0.2, 0.3))
            dampings.append(social_events(start_summer, 0.0, 0.2))
            dampings.append(social_events_work(start_summer, 0.0, 0.05))
            dampings.append(
                physical_distancing_home_school(
                    start_summer, 0.0, 0.2))
            dampings.append(
                physical_distancing_work_other(
                    start_summer, 0.0, 0.2))
            dampings.append(senior_awareness(start_summer, 0.0, 0.0))

            # autumn enforced attention
            start_autumn_date = datetime.date(year=2020, month=10, day=1)
            if start_autumn_date < end_date:
                start_autumn = (start_autumn_date - self.start_date).days
                dampings.append(contacts_at_home(start_autumn, 0.2, 0.4))
                dampings.append(
                    physical_distancing_home_school(
                        start_autumn, 0.2, 0.4))
                dampings.append(
                    physical_distancing_work_other(
                        start_autumn, 0.2, 0.4))

            # autumn lockdown light
            start_autumn_lockdown_date = datetime.date(2020, 11, 1)
            if (start_autumn_lockdown_date < end_date):
                start_autumn_lockdown = (
                    start_autumn_lockdown_date - self.start_date).days
                dampings.append(contacts_at_home(
                    start_autumn_lockdown, 0.4, 0.6))
                dampings.append(school_closure(
                    start_autumn_lockdown, 0.0, 0.0))
                dampings.append(home_office(start_autumn_lockdown, 0.2, 0.3))
                dampings.append(social_events(start_autumn_lockdown, 0.6, 0.8))
                dampings.append(social_events_work(
                    start_autumn_lockdown, 0.0, 0.1))
                dampings.append(
                    physical_distancing_home_school(
                        start_autumn_lockdown, 0.2, 0.4))
                dampings.append(
                    physical_distancing_work_other(
                        start_autumn_lockdown, 0.4, 0.6))
                dampings.append(senior_awareness(
                    start_autumn_lockdown, 0.0, 0.0))

            # winter lockdown
            start_winter_lockdown_date = datetime.date(2020, 12, 16)
            if (start_winter_lockdown_date < end_date):
                min = 0.6
                max = 0.8  # for strictest scenario: 0.8 - 1.0
                start_winter_lockdown = (
                    start_winter_lockdown_date - self.start_date).days
                dampings.append(contacts_at_home(
                    start_winter_lockdown, min, max))
                dampings.append(school_closure(
                    start_winter_lockdown, 1.0, 1.0))
                dampings.append(home_office(start_winter_lockdown, 0.2, 0.3))
                dampings.append(social_events(start_winter_lockdown, min, max))
                dampings.append(social_events_work(
                    start_winter_lockdown, 0.1, 0.2))
                dampings.append(
                    physical_distancing_home_school(
                        start_winter_lockdown, 0.2, 0.4))
                dampings.append(
                    physical_distancing_work_other(
                        start_winter_lockdown, min, max))
                dampings.append(senior_awareness(
                    start_winter_lockdown, 0.0, 0.0))

                # relaxing of restrictions over christmas days
                xmas_date = datetime.date(2020, 12, 24)
                xmas = (xmas_date - self.start_date).days
                dampings.append(contacts_at_home(xmas, 0.0, 0.0))
                dampings.append(home_office(xmas, 0.4, 0.5))
                dampings.append(social_events(xmas, 0.4, 0.6))
                dampings.append(
                    physical_distancing_home_school(
                        xmas, 0.0, 0.0))
                dampings.append(physical_distancing_work_other(xmas, 0.4, 0.6))

                # after christmas
                after_xmas_date = datetime.date(2020, 12, 27)
                after_xmas = (after_xmas_date - self.start_date).days
                dampings.append(contacts_at_home(after_xmas, min, max))
                dampings.append(home_office(after_xmas, 0.2, 0.3))
                dampings.append(social_events(after_xmas, 0.6, 0.8))
                dampings.append(
                    physical_distancing_home_school(after_xmas, 0.2, 0.4))
                dampings.append(
                    physical_distancing_work_other(after_xmas, min, max))

            # local dynamic NPIs
            dynamic_npis = params.DynamicNPIsInfectedSymptoms
            local_npis = []
            # increased from [0.4, 0.6] in Nov
            local_npis.append(contacts_at_home(0, 0.6, 0.8))
            local_npis.append(school_closure(0, 0.25, 0.25))  # see paper
            local_npis.append(home_office(0, 0.2, 0.3))
            local_npis.append(social_events(0, 0.6, 0.8))
            local_npis.append(social_events_work(0, 0.1, 0.2))
            local_npis.append(physical_distancing_home_school(0, 0.6, 0.8))
            local_npis.append(physical_distancing_work_other(0, 0.6, 0.8))
            local_npis.append(senior_awareness(0, 0.0, 0.0))

            dynamic_npis.interval = 3.0
            dynamic_npis.duration = 14.0
            dynamic_npis.base_value = 100000
            dynamic_npis.set_threshold(200.0, local_npis)

            # school holidays(holiday periods are set per node, see set_nodes)
            contacts.school_holiday_damping = damping_helper(
                0, 1.0, 1.0, lvl_holidays, typ_school, [loc_school])
            contacts.dampings = dampings

    def get_graph(self, end_date):
        """

        :param end_date: 

        """
        model = Model(self.num_groups)
        self.set_covid_parameters(model)
        self.set_contact_matrices(model)
        self.set_npis(model.parameters, end_date)

        graph = osecir.ModelGraph()

        scaling_factor_infected = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
        scaling_factor_icu = 1.0
        tnt_capacity_factor = 7.5 / 100000.

        data_dir_Germany = os.path.join(self.data_dir, "Germany")
        mobility_data_file = os.path.join(
            data_dir_Germany, "mobility", "commuter_mobility_2022.txt")
        pydata_dir = os.path.join(data_dir_Germany, "pydata")

        path_population_data = os.path.join(pydata_dir,
                                            "county_current_population.json")

        mio.osecir.set_nodes(
            model.parameters,
            mio.Date(self.start_date.year,
                     self.start_date.month, self.start_date.day),
            mio.Date(end_date.year,
                     end_date.month, end_date.day), pydata_dir,
            path_population_data, True, graph, scaling_factor_infected,
            scaling_factor_icu, tnt_capacity_factor, 0, False)

        mio.osecir.set_edges(
            mobility_data_file, graph, len(Location))

        return graph

    def run(self, num_days_sim, num_runs=10, save_graph=True, create_gif=True):
        """

        :param num_days_sim: 
        :param num_runs:  (Default value = 10)
        :param save_graph:  (Default value = True)
        :param create_gif:  (Default value = True)

        """
        mio.set_log_level(mio.LogLevel.Warning)
        end_date = self.start_date + datetime.timedelta(days=num_days_sim)

        graph = self.get_graph(end_date)

        if save_graph:
            path_graph = os.path.join(self.results_dir, "graph")
            if not os.path.exists(path_graph):
                os.makedirs(path_graph)
            osecir.write_graph(graph, path_graph)

        study = osecir.ParameterStudy(
            graph, 0., num_days_sim, 0.5, num_runs)
        ensemble = study.run()

        ensemble_results = []
        ensemble_params = []
        for run in range(num_runs):
            graph_run = ensemble[run]
            ensemble_results.append(interpolate_simulation_result(graph_run))
            ensemble_params.append(
                [graph_run.get_node(node_indx).property.model
                 for node_indx in range(graph.num_nodes)])

        node_ids = [graph.get_node(i).id for i in range(graph.num_nodes)]

        save_percentiles = True
        save_single_runs = False

        osecir.save_results(
            ensemble_results, ensemble_params, node_ids, self.results_dir,
            save_single_runs, save_percentiles)
        if create_gif:
            # any compartments in the model (see InfectionStates)
            compartments = [c for c in range(1, 8)]
            mp.create_gif_map_plot(
                self.results_dir + "/p75", self.results_dir, compartments)
        return 0


if __name__ == "__main__":
    file_path = os.path.dirname(os.path.abspath(__file__))
    sim = Simulation(
        data_dir=os.path.join(file_path, "../../../data"),
        start_date=datetime.date(year=2020, month=12, day=12),
        results_dir=os.path.join(file_path, "../../../results_osecir"))
    num_days_sim = 50
    sim.run(num_days_sim, num_runs=2)
