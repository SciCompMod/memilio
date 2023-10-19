import numpy as np
import datetime
import copy
import os.path as path
import memilio.simulation as mio
import memilio.simulation.secir as secir

from enum import Enum, auto
from memilio.simulation.secir import AgeGroup, Index_InfectionState
from memilio.simulation.secir import InfectionState as State
from memilio.simulation.secir import (Model, Simulation,
                                      interpolate_simulation_result, simulate)


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

    def __init__(self, data_dir, start_date):
        self.num_groups = 6
        self.data_dir = data_dir
        self.start_date = start_date

    def set_covid_parameters(self, model):
        def array_assign_uniform_distribution(param, min, max, num_groups=6):
            if isinstance(
                    min, (int, float)) and isinstance(
                        max, (int, float)):
                min = min * np.ones(num_groups)
                max = max * np.ones(num_groups)
            elif not (isinstance(min, (list, tuple)) and isinstance(
                    max, (list, tuple))):
                raise TypeError(
                    "Invalid type for parameter 'min' or 'max. \
                            Expected a scalar or a vector. Must be the same for both.")
            for i in range(num_groups):
                param[secir.AgeGroup(i)] = mio.UncertainValue(
                    0.5 * (max[i] + min[i]))
                param[secir.AgeGroup(i)].set_distribution(
                    mio.ParameterDistributionUniform(min[i], max[i]))

        t_incubation = 5.2
        t_serial_interval_min = 0.5 * 2.67 + 0.5 * 5.2
        t_serial_interval_max = 0.5 * 4.00 + 0.5 * 5.2
        timeInfectedSymptomsMin = [
            5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465]
        timeInfectedSymptomsMax = [8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085]
        timeInfectedSevereMin = [3.925, 3.925, 4.85, 6.4, 7.2, 9.]
        timeInfectedSevereMax = [6.075, 6.075, 7., 8.7, 9.8, 13.]
        timeInfectedCriticalMin = [4.95, 4.95, 4.86, 14.14, 14.4, 10.]
        timeInfectedCriticalMax = [8.95, 8.95, 8.86, 20.58, 19.8, 13.2]

        array_assign_uniform_distribution(
            model.parameters.IncubationTime, t_incubation, t_incubation)

        array_assign_uniform_distribution(
            model.parameters.SerialInterval, t_serial_interval_min,
            t_serial_interval_max)

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

        model.parameters.StartDay = (
            self.start_date -
            datetime.date(year=self.start_date.year, month=1, day=1)).days + 1
        model.parameters.Seasonality.value = 0.2

        seasonality_min = 0.1
        seasonality_max = 0.3
        model.parameters.Seasonality = mio.UncertainValue(
            0.5 * (seasonality_max + seasonality_min))
        model.parameters.Seasonality.set_distribution(
            mio.ParameterDistributionUniform(seasonality_min, seasonality_max))

    def set_contact_matrices(self, model):
        contact_matrices = mio.ContactMatrixGroup(
            len(list(Location)), self.num_groups)
        contact_matrices[0] = mio.ContactMatrix(
            mio.secir.read_mobility_plain(
                path.join(self.data_dir, "contacts", "baseline_home.txt")),
            mio.secir.read_mobility_plain(
                path.join(self.data_dir, "contacts", "minimum_home.txt")))
        contact_matrices[1] = mio.ContactMatrix(
            mio.secir.read_mobility_plain(
                path.join(
                    self.data_dir, "contacts", "baseline_school_pf_eig.txt")),
            mio.secir.read_mobility_plain(
                path.join(
                    self.data_dir, "contacts", "minimum_school_pf_eig.txt")))
        contact_matrices[2] = mio.ContactMatrix(
            mio.secir.read_mobility_plain(
                path.join(self.data_dir, "contacts", "baseline_work.txt")),
            mio.secir.read_mobility_plain(
                path.join(self.data_dir, "contacts", "minimum_work.txt")))
        contact_matrices[3] = mio.ContactMatrix(
            mio.secir.read_mobility_plain(
                path.join(self.data_dir, "contacts", "baseline_other.txt")),
            mio.secir.read_mobility_plain(
                path.join(self.data_dir, "contacts", "minimum_other.txt")))

        model.parameters.ContactPatterns.cont_freq_mat = contact_matrices

    def set_npis(self, params, end_date):
        contacts = params.ContactPatterns
        dampings = contacts.dampings

        group_weights_all = np.ones(self.num_groups)
        group_weights_seniors = np.vectorize(
            lambda i: 1.0 if i == 5 else (0.5 if i == 4 else 0.0))(
            np.arange(self.num_groups))

        loc_h = Location.Home.value
        loc_s = Location.School.value
        loc_w = Location.Work.value
        loc_o = Location.Other.value

        lvl_m = InterventionLevel.Main.value
        lvl_p = InterventionLevel.PhysicalDistanceAndMasks.value
        lvl_s = InterventionLevel.SeniorAwareness.value
        lvl_h = InterventionLevel.Holidays.value

        typ_h = Intervention.Home.value
        typ_s = Intervention.SchoolClosure.value
        typ_ho = Intervention.HomeOffice.value
        typ_g = Intervention.GatheringBanFacilitiesClosure.value
        typ_p = Intervention.PhysicalDistanceAndMasks.value
        typ_se = Intervention.SeniorAwareness.value

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

        def contacts_at_home(t, min, max):
            return damping_helper(
                t, min, max, lvl_m, typ_h, [loc_h])

        def school_closure(t, min, max):
            return damping_helper(
                t, min, max, lvl_m, typ_s, [loc_s])

        def home_office(t, min, max):
            return damping_helper(
                t, min, max, lvl_m, typ_ho, [loc_w])

        def social_events(t, min, max):
            return damping_helper(
                t, min, max, lvl_m, typ_g, [loc_o])

        def social_events_work(t, min, max):
            return damping_helper(
                t, min, max, lvl_m, typ_g, [loc_w])

        def physical_distancing_home_school(t, min, max):
            return damping_helper(
                t, min, max, lvl_p, typ_p, [loc_h, loc_s])

        def physical_distancing_work_other(t, min, max):
            return damping_helper(
                t, min, max, lvl_p, typ_p, [loc_w, loc_o])

        def senior_awareness(t, min, max):
            return damping_helper(
                t, min, max, lvl_s, typ_se, [loc_h, loc_o],
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
                0, 1.0, 1.0, lvl_h, typ_s, [loc_s])
            contacts.dampings = dampings

    def get_graph(self, end_date):
        model = Model(self.num_groups)
        self.set_covid_parameters(model)
        self.set_contact_matrices(model)
        self.set_npis(model.parameters, end_date)

        graph = secir.ModelGraph()

        scaling_factor_infected = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
        scaling_factor_icu = 1.0
        tnt_capacity_factor = 7.5 / 100000.

        path_population_data = path.join(
            self.data_dir, "pydata", "Germany",
            "county_current_population.json")

        s_day = mio.Date(self.start_date.year,
                         self.start_date.month, self.start_date.day)
        e_day = mio.Date(end_date.year,
                         end_date.month, end_date.day)
        mio.secir.set_nodes(
            model.parameters, s_day, e_day, self.data_dir,
            path_population_data, True, graph, scaling_factor_infected,
            scaling_factor_icu, tnt_capacity_factor, 0, False)

        mio.secir.set_edges(
            self.data_dir, graph, len(Location))

        return graph

    def run(self, num_days_sim, num_runs=10):
        mio.set_log_level(mio.LogLevel.Warning)
        end_date = self.start_date + datetime.timedelta(days=num_days_sim)

        graph = self.get_graph(end_date)

        def find_indices_of_true_values(boolean_array):
            true_indices = np.where(boolean_array)[0]
            return true_indices

        def print_param(parameters):
            for i in range(6):
                A0 = secir.AgeGroup(i)
                print("age group " + str(i))
                print(parameters.RelativeTransmissionNoSymptoms[A0].value)
                print(parameters.TransmissionProbabilityOnContact[A0].value)
                print(parameters.RecoveredPerInfectedNoSymptoms[A0].value)
                print(parameters.RiskOfInfectionFromSymptomatic[A0].value)
                print(parameters.SeverePerInfectedSymptoms[A0].value)
                print(parameters.CriticalPerSevere[A0].value)
                print(parameters.DeathsPerCritical[A0].value)
                print(parameters.MaxRiskOfInfectionFromSymptomatic[A0].value)
                print(parameters.IncubationTime[A0].value)
                print(parameters.TimeInfectedSymptoms[A0].value)
                print(parameters.SerialInterval[A0].value)
                print(parameters.TimeInfectedSevere[A0].value)
                print(parameters.TimeInfectedCritical[A0].value)

        # find_indices_of_true_values([np.any(np.isnan(graph.get_node(i).property.result.get_value(0))) for i in range(400)])

        def handle_result(graph, run_idx):
            handle_result.interpolated = interpolate_simulation_result(graph)
            print("run " + str(run_idx))

        mio.secir.write_graph(graph, "graph_python")

        study = secir.ParameterStudy(
            graph, 0., num_days_sim, 0.5, num_runs)
        study.run(handle_result)
        last_result = handle_result.interpolated
        return [ts.as_ndarray() for ts in last_result]


if __name__ == "__main__":
    file_path = path.dirname(path.abspath(__file__))
    sim = Simulation(
        data_dir=path.join(file_path, "../../../data"),
        start_date=datetime.date(year=2020, month=12, day=12))
    num_days_sim = 30
    sim.run(num_days_sim, num_runs=2)
