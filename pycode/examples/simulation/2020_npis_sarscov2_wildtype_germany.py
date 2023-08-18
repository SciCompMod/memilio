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
    Work = 1
    School = 2
    Other = 3


class Intervention(Enum):
    Home = auto()
    SchoolClosure = auto()
    HomeOffice = auto()
    GatheringBanFacilitiesClosure = auto()
    PhysicalDistanceAndMasks = auto()
    SeniorAwareness = auto()


class Simulation:

    def __init__(self, data_dir):

        self.start_date = mio.Date(2020, 5, 15)
        start_day = (datetime.date(year=2020, month=5,
                                   day=15) - datetime.date(year=2020, month=1, day=1)).days

        self.end_date = mio.Date(2020, 9, 1)

        num_groups = 6

        model = Model(num_groups)

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

        contact_matrices = mio.ContactMatrixGroup(
            num_groups, len(list(Location)))
        contact_matrices[0] = mio.ContactMatrix(
            mio.secir.read_mobility_plain(
                path.join(data_dir, "contacts", "baseline_home.txt")),
            mio.secir.read_mobility_plain(
                path.join(data_dir, "contacts", "minimum_home.txt")))
        contact_matrices[1] = mio.ContactMatrix(
            mio.secir.read_mobility_plain(
                path.join(
                    data_dir, "contacts", "baseline_school_pf_eig.txt")),
            mio.secir.read_mobility_plain(
                path.join(
                    data_dir, "contacts", "minimum_school_pf_eig.txt")))
        contact_matrices[2] = mio.ContactMatrix(
            mio.secir.read_mobility_plain(
                path.join(data_dir, "contacts", "baseline_work.txt")),
            mio.secir.read_mobility_plain(
                path.join(data_dir, "contacts", "minimum_work.txt")))
        contact_matrices[3] = mio.ContactMatrix(
            mio.secir.read_mobility_plain(
                path.join(data_dir, "contacts", "baseline_other.txt")),
            mio.secir.read_mobility_plain(
                path.join(data_dir, "contacts", "minimum_other.txt")))

        params = model.parameters
        params.ContactPatterns.cont_freq_mat = contact_matrices
        params.StartDay = start_day
        params.Seasonality.value = 0.2

        graph = secir.ModelGraph()

        scaling_factor_infected = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
        scaling_factor_icu = 1.0
        tnt_capacity_factor = 7.5 / 100000.

        path_population_data = path.join(
            data_dir, "pydata", "Germany", "county_current_population.json")

        mio.secir.set_nodes(
            params, self.start_date, self.end_date, data_dir,
            path_population_data, True, graph, scaling_factor_infected,
            scaling_factor_icu, tnt_capacity_factor, 0, False)

        mio.secir.set_edges(
            data_dir, graph, len(Location))

        self.graph = graph

    def run2(self, num_runs=10):
        mio.set_log_level(mio.LogLevel.Off)

        # calculate tmax
        tmax = (
            datetime.date(
                self.end_date.year, self.end_date.month, self.end_date.day) -
            datetime.date(
                self.start_date.year, self.start_date.month, self.start_date.
                day)).days

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
            group = secir.AgeGroup(0)
            print("run {} with infection rate {:.2G}".format(handle_result.c, graph.get_node(
                0).property.model.parameters.TransmissionProbabilityOnContact[group].value))
            print("compartments at t = {}:".format(
                graph.get_node(0).property.result.get_time(0)))
            print(graph.get_node(0).property.result.get_value(0))
            print("compartments at t = {}:".format(
                graph.get_node(0).property.result.get_last_time()))
            print(graph.get_node(0).property.result.get_last_value())
            if (any(np.any(np.isnan(graph.get_node(i).property.result.get_value(1))) for i in range(400))):
                indx_nan = find_indices_of_true_values(
                    [np.any(
                        np.isnan(
                            graph.get_node(i).property.result.get_value(1)))
                     for i in range(400)])
                print_param(graph.get_node(
                    indx_nan[0]).property.model.parameters)
            handle_result.c += 1
        handle_result.c = 0

        study = secir.ParameterStudy(
            self.graph, 0., tmax, 0.5, num_runs)
        study.run(handle_result)

        # self.last_result = handle_result.interpolated
        return [ts.as_ndarray() for ts in self.last_result]

        # self._param_study.run(dampings, dampings_school,
        #                       local_npis, handle_result)
        # self.last_result = handle_result.interpolated
        # return ts

    def run(self, params=[0.1, 0.5, 0.25, 0.1, 0.05, 0.1, 0.1, 0.0]):
        # todo
        dampings = DampingHelper()

        #######
        # spring lockdown
        #######
        spring_lockdown_date = datetime.date(year=2020, month=3, day=18)
        spring_lockdown_day = (spring_lockdown_date - self._start_date).days

        # home
        dampings.add(0.6, Intervention.Home.value, 1,
                     spring_lockdown_day, [Location.Home.value])
        # school
        dampings.add(1.0, Intervention.SchoolClosure.value, 1,
                     spring_lockdown_day, [Location.School.value])
        # homeoffice, work closures
        dampings.add(0.2, Intervention.HomeOffice.value, 1,
                     spring_lockdown_day, [Location.Work.value])
        # gatherings, events
        dampings.add(
            0.6, Intervention.GatheringBanFacilitiesClosure.value, 1,
            spring_lockdown_day, [Location.Other.value])
        dampings.add(
            0.1, Intervention.GatheringBanFacilitiesClosure.value, 1,
            spring_lockdown_day, [Location.Work.value])
        # distancing, masks
        dampings.add(
            0.4, Intervention.PhysicalDistanceAndMasks.value, 2,
            spring_lockdown_day, [Location.Home.value, Location.School.value])
        dampings.add(
            0.4, Intervention.PhysicalDistanceAndMasks.value, 2,
            spring_lockdown_day, [Location.Work.value, Location.Other.value])
        # senior protection
        dampings.add(0.0, Intervention.SeniorAwareness.value, 3, spring_lockdown_day, [
                     Location.Home.value, Location.Other.value], [0.0, 0.0, 0.0, 0.0, 0.5, 1.0])

        #######
        # SUMMER
        #######
        school_reopen_day = (datetime.date(
            2020, 6, 15) - self._start_date).days
        dampings.add(0.0, Intervention.SchoolClosure.value, 1,
                     school_reopen_day, [Location.School.value])

        start_summer_day = (datetime.date(2020, 5, 15) - self._start_date).days
        # home
        dampings.add(params[0], Intervention.Home.value, 1,
                     start_summer_day, [Location.Home.value])
        # school
        dampings.add(
            params[1],
            Intervention.SchoolClosure.value, 1, start_summer_day,
            [Location.School.value])
        # homeoffice, work closures
        dampings.add(params[2], Intervention.HomeOffice.value,
                     1, start_summer_day, [Location.Work.value])
        # gatherings, events
        dampings.add(
            params[3],
            Intervention.GatheringBanFacilitiesClosure.value, 1,
            start_summer_day, [Location.Other.value])
        dampings.add(
            params[4],
            Intervention.GatheringBanFacilitiesClosure.value, 1,
            start_summer_day, [Location.Work.value])
        # distancing, masks
        dampings.add(
            params[5],
            Intervention.PhysicalDistanceAndMasks.value, 2, start_summer_day,
            [Location.Home.value, Location.School.value])
        dampings.add(
            params[6],
            Intervention.PhysicalDistanceAndMasks.value, 2, start_summer_day,
            [Location.Work.value, Location.Other.value])
        # senior protection
        dampings.add(params[7], Intervention.SeniorAwareness.value, 3, start_summer_day, [
                     Location.Home.value, Location.Other.value], [0.0, 0.0, 0.0, 0.0, 0.5, 1.0])

        #######
        # school holidays
        #######
        dampings_school = DampingHelper()
        holiday_day = 0  # dynamic
        dampings_school.add(
            1.0, Intervention.SchoolClosure.value, 5, holiday_day,
            [Location.School.value])

        #######
        # dynamic lockdown
        #######
        local_npi_day = 0  # dynamic
        local_npi_1 = DampingHelper()

        # home
        local_npi_1.add(0.2, Intervention.Home.value, 1,
                        local_npi_day, [Location.Home.value])
        # school
        # local_npi_1.add(0.0, Intervetion.SchoolClosure.value, 1, local_npi_day, [Location.School.value])
        # homeoffice, work closures
        local_npi_1.add(0.2, Intervention.HomeOffice.value, 1,
                        local_npi_day, [Location.Work.value])
        # gatherings, events
        local_npi_1.add(
            0.2, Intervention.GatheringBanFacilitiesClosure.value, 1,
            local_npi_day, [Location.Other.value])
        local_npi_1.add(
            0.0, Intervention.GatheringBanFacilitiesClosure.value, 1,
            local_npi_day, [Location.Work.value])
        # distancing, masks
        local_npi_1.add(
            0.2, Intervention.PhysicalDistanceAndMasks.value, 2, local_npi_day,
            [Location.Home.value, Location.School.value])
        local_npi_1.add(
            0.2, Intervention.PhysicalDistanceAndMasks.value, 2, local_npi_day,
            [Location.Work.value, Location.Other.value])
        # senior protection
        local_npi_1.add(0.0, Intervention.SeniorAwareness.value, 3, local_npi_day, [
                        Location.Home.value, Location.Other.value], [0.0, 0.0, 0.0, 0.0, 0.5, 1.0])

        local_npi_2 = DampingHelper()
        # home
        local_npi_2.add(0.6, Intervention.Home.value, 1,
                        local_npi_day, [Location.Home.value])
        # school
        local_npi_2.add(1.0, Intervention.SchoolClosure.value,
                        1, local_npi_day, [Location.School.value])
        # homeoffice, work closures
        local_npi_2.add(0.2, Intervention.HomeOffice.value, 1,
                        local_npi_day, [Location.Work.value])
        # gatherings, events
        local_npi_2.add(
            0.6, Intervention.GatheringBanFacilitiesClosure.value, 1,
            local_npi_day, [Location.Other.value])
        local_npi_2.add(
            0.1, Intervention.GatheringBanFacilitiesClosure.value, 1,
            local_npi_day, [Location.Work.value])
        # distancing, masks
        local_npi_2.add(
            0.6, Intervention.PhysicalDistanceAndMasks.value, 2, local_npi_day,
            [Location.Home.value, Location.School.value])
        local_npi_2.add(
            0.6, Intervention.PhysicalDistanceAndMasks.value, 2, local_npi_day,
            [Location.Work.value, Location.Other.value])
        # senior protection
        local_npi_2.add(0.0, Intervention.SeniorAwareness.value, 3, local_npi_day, [
                        Location.Home.value, Location.Other.value], [0.0, 0.0, 0.0, 0.0, 0.5, 1.0])

        local_npis = [local_npi_1, local_npi_2]

        def handle_result(graph):
            handle_result.interpolated = interpolate_simulation_result(graph)
        self._param_study.run(dampings, dampings_school,
                              local_npis, handle_result)
        self.last_result = handle_result.interpolated
        return [ts.as_ndarray() for ts in self.last_result]


if __name__ == "__main__":
    # TODO: get abs path
    output_path = path.dirname(path.abspath(__file__))
    sim = Simulation(
        data_dir=path.join(output_path, "../../../data"))
    sim.run2()
