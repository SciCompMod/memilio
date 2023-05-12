import numpy as np
import datetime
import copy
import os.path as path
import memilio.simulation as mio
from enum import Enum, auto

from memilio.simulation import ContactMatrix, Damping, UncertainContactMatrix
import memilio.simulation.secir as secir
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

    def array_assign_uniform_distribution()

    def __init__(self, data_dir):
        t0 = 0.0
        tmax = 45.0
        dt = 0.5

        self._start_date = datetime.date(year=2020, month=5, day=15)
        start_day = (
            self._start_date - datetime.date(year=2020, month=1, day=1)).days

        num_groups = 6

        model = Model(num_groups)

        t_incubation = 5.2
        t_serial_interval_min = 0.5 * 2.67 + 0.5 * 5.2
        t_serial_interval_max = 0.5 * 4.00 + 0.5 * 5.2
        timeInfectedSymptomsMin = {
            5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465}
        timeInfectedSymptomsMax = {8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085}
        timeInfectedSevereMin = {3.925, 3.925, 4.85, 6.4, 7.2, 9.}
        timeInfectedSevereMax = {6.075, 6.075, 7., 8.7, 9.8, 13.}
        timeInfectedCriticalMin = {4.95, 4.95, 4.86, 14.14, 14.4, 10.}
        timeInfectedCriticalMax = {8.95, 8.95, 8.86, 20.58, 19.8, 13.2}

        model.parameters.IncubationTime = mio.UncertainValue(
            0.5 * (t_incubation + t_incubation))
        model.parameters.IncubationTime.set_distribution(
            mio.ParameterDistributionUniform(t_incubation, t_incubation))

        model.parameters.SerialInterval = mio.UncertainValue(
            0.5 * (t_serial_interval_min + t_serial_interval_max))
        model.parameters.SerialInterval.set_distribution(
            mio.ParameterDistributionUniform(
                t_serial_interval_min, t_serial_interval_max))

        model.parameters.TimeInfectedSymptoms = mio.UncertainValue(
            0.5 * (timeInfectedSymptomsMin + timeInfectedSymptomsMax))
        model.parameters.TimeInfectedSymptoms.set_distribution(
            mio.ParameterDistributionUniform(
                timeInfectedSymptomsMin, timeInfectedSymptomsMax))

        model.parameters.TimeInfectedSevere = mio.UncertainValue(
            0.5 * (timeInfectedSevereMin + timeInfectedSevereMax))
        model.parameters.TimeInfectedSevere.set_distribution(
            mio.ParameterDistributionUniform(
                timeInfectedSevereMin, timeInfectedSevereMax))

        model.parameters.TimeInfectedCritical = mio.UncertainValue(
            0.5 * (timeInfectedCriticalMin + timeInfectedCriticalMax))
        model.parameters.TimeInfectedCritical.set_distribution(
            mio.ParameterDistributionUniform(
                timeInfectedCriticalMin, timeInfectedCriticalMax))

        # probabilities
        transmissionProbabilityOnContactMin = {
            0.02, 0.05, 0.05, 0.05, 0.08, 0.15}
        transmissionProbabilityOnContactMax = {
            0.04, 0.07, 0.07, 0.07, 0.10, 0.20}
        relativeTransmissionNoSymptomsMin = 1
        relativeTransmissionNoSymptomsMax = 1
        # The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
        #  depends on incidence and test and trace capacity
        riskOfInfectionFromSymptomaticMin = 0.1
        riskOfInfectionFromSymptomaticMax = 0.3
        maxRiskOfInfectionFromSymptomaticMin = 0.3
        maxRiskOfInfectionFromSymptomaticMax = 0.5
        recoveredPerInfectedNoSymptomsMin = {0.2, 0.2, 0.15, 0.15, 0.15, 0.15}
        recoveredPerInfectedNoSymptomsMax = {0.3, 0.3, 0.25, 0.25, 0.25, 0.25}
        severePerInfectedSymptomsMin = {0.006, 0.006, 0.015, 0.049, 0.15, 0.20}
        severePerInfectedSymptomsMax = {0.009, 0.009, 0.023, 0.074, 0.18, 0.25}
        criticalPerSevereMin = {0.05, 0.05, 0.05, 0.10, 0.25, 0.35}
        criticalPerSevereMax = {0.10, 0.10, 0.10, 0.20, 0.35, 0.45}
        deathsPerCriticalMin = {0.00, 0.00, 0.10, 0.10, 0.30, 0.5}
        deathsPerCriticalMax = {0.10, 0.10, 0.18, 0.18, 0.50, 0.7}

        model.parameters.TransmissionProbabilityOnContact = mio.UncertainValue(
            0.5 * (transmissionProbabilityOnContactMin + transmissionProbabilityOnContactMax))
        model.parameters.TransmissionProbabilityOnContact.set_distribution(
            mio.ParameterDistributionUniform(
                transmissionProbabilityOnContactMin,
                transmissionProbabilityOnContactMax))

        model.parameters.RelativeTransmissionNoSymptoms = mio.UncertainValue(
            0.5 * (relativeTransmissionNoSymptomsMin + relativeTransmissionNoSymptomsMax))
        model.parameters.RelativeTransmissionNoSymptoms.set_distribution(
            mio.ParameterDistributionUniform(
                relativeTransmissionNoSymptomsMin,
                relativeTransmissionNoSymptomsMax))

        model.parameters.RiskOfInfectionFromSymptomatic = mio.UncertainValue(
            0.5 * (riskOfInfectionFromSymptomaticMin + riskOfInfectionFromSymptomaticMax))
        model.parameters.RiskOfInfectionFromSymptomatic.set_distribution(
            mio.ParameterDistributionUniform(
                riskOfInfectionFromSymptomaticMin,
                riskOfInfectionFromSymptomaticMax))

        model.parameters.MaxRiskOfInfectionFromSymptomatic = mio.UncertainValue(
            0.5 * (maxRiskOfInfectionFromSymptomaticMin + maxRiskOfInfectionFromSymptomaticMax))
        model.parameters.MaxRiskOfInfectionFromSymptomatic.set_distribution(
            mio.ParameterDistributionUniform(
                maxRiskOfInfectionFromSymptomaticMin,
                maxRiskOfInfectionFromSymptomaticMax))

        model.parameters.RecoveredPerInfectedNoSymptoms = mio.UncertainValue(
            0.5 * (recoveredPerInfectedNoSymptomsMin + recoveredPerInfectedNoSymptomsMax))
        model.parameters.RecoveredPerInfectedNoSymptoms.set_distribution(
            mio.ParameterDistributionUniform(
                recoveredPerInfectedNoSymptomsMin,
                recoveredPerInfectedNoSymptomsMax))

        model.parameters.SeverePerInfectedSymptoms = mio.UncertainValue(
            0.5 * (severePerInfectedSymptomsMin + severePerInfectedSymptomsMax))
        model.parameters.SeverePerInfectedSymptoms.set_distribution(
            mio.ParameterDistributionUniform(
                severePerInfectedSymptomsMin,
                severePerInfectedSymptomsMax))

        model.parameters.CriticalPerSevere = mio.UncertainValue(
            0.5 * (criticalPerSevereMin + criticalPerSevereMax))
        model.parameters.CriticalPerSevere.set_distribution(
            mio.ParameterDistributionUniform(
                criticalPerSevereMin,
                criticalPerSevereMax))

        model.parameters.DeathsPerCritical = mio.UncertainValue(
            0.5 * (deathsPerCriticalMin + deathsPerCriticalMax))
        model.parameters.DeathsPerCritical.set_distribution(
            mio.ParameterDistributionUniform(
                deathsPerCriticalMin,
                deathsPerCriticalMax))

        contact_matrices = mio.ContactMatrixGroup(
            num_groups, len(list(Location)))
        contact_matrices[0] = mio.ContactMatrix(
            np.loadtxt(
                path.join(data_dir, "contacts", "baseline_home.txt")),
            np.loadtxt(
                path.join(data_dir, "contacts", "minimum_home.txt")))
        contact_matrices[1] = mio.ContactMatrix(
            np.loadtxt(
                path.join(
                    data_dir, "contacts", "baseline_school_pf_eig.txt")),
            np.loadtxt(
                path.join(
                    data_dir, "contacts", "minimum_school_pf_eig.txt")))
        contact_matrices[2] = mio.ContactMatrix(
            np.loadtxt(
                path.join(data_dir, "contacts", "baseline_work.txt")),
            np.loadtxt(
                path.join(data_dir, "contacts", "minimum_work.txt")))
        contact_matrices[3] = mio.ContactMatrix(
            np.loadtxt(
                path.join(data_dir, "contacts", "baseline_other.txt")),
            np.loadtxt(
                path.join(data_dir, "contacts", "minimum_other.txt")))

        params = model.parameters
        params.get_contact_patterns().cont_freq_mat = contact_matrices
        params.set_start_day(start_day)
        params.set_seasonality(0.2)

        # read data
        county_ids = mio.get_county_ids(
            path.join(data_dir, "pydata", "Germany"))

        commuter = np.loadtxt(
            path.join(data_dir, "migration", "commuter_migration_scaled.txt"))
        twitter = np.loadtxt(
            path.join(data_dir, "migration", "twitter_scaled_1252.txt"))

        num_nodes = commuter.shape[0]
        params_nodes = []
        for i in range(num_nodes):
            params_nodes.append(copy.copy(params))

        graph = secir.ModelGraph()

        scaling_factor_infected = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
        scaling_factor_icu = 1.0
        tnt_capacity_factor = 7.5 / 100000.
        migrating_compartments = [State.Susceptible, State.Exposed,
                                  State.InfectedNoSymptoms, State.InfectedSymptoms, State.Recovered]

        # TODO: Call with templates
        read_function_nodes = mio.read_input_data_county()
        read_function_edges = np.loadtxt()
        node_id_function = mio.get_node_ids()

        # TODO: How to call io functions
    #     set_node_function = mio.set_nodes<mio::osecir::TestAndTraceCapacity, mio::osecir::ContactPatterns, mio::osecir::Model,
    #                    mio::MigrationParameters, mio::osecir::Parameters, decltype(read_function_nodes),
    #                    decltype(node_id_function)>;
    #      =
    #     mio::set_edges<ContactLocation, mio::osecir::Model, mio::MigrationParameters, mio::MigrationCoefficientGroup,
    #                    mio::osecir::InfectionState, decltype(read_function_edges)>;
    # BOOST_OUTCOME_TRY(
    #     set_node_function(params, start_date, end_date, data_dir,
    #                       mio::path_join((data_dir / "pydata" / "Germany").string(), "county_current_population.json"),
    #                       true, params_graph, read_function_nodes, node_id_function, scaling_factor_infected,
    #                       scaling_factor_icu, tnt_capacity_factor, 0, false));
    # BOOST_OUTCOME_TRY(set_edge_function(data_dir, params_graph, migrating_compartments, contact_locations.size(),
    #                                     read_function_edges));

        params_nodes = read_population_data_county(
            path.join(data_dir, "pydata", "Germany"),
            (self._start_date.year, self._start_date.month, self._start_date.day),
            county_ids,
            [2.5], 1.0,
            params_nodes)

        # set counties
        for i in range(num_nodes):
            params_nodes[i].set_icu_capacity(1e200)
            tnt_capacity = params_nodes[i].populations.get_total(
            ) * 7.5 / 100000.
            params_nodes[i].set_test_and_trace_capacity(tnt_capacity)
            params_nodes[i].apply_constraints()
            graph.add_node(params_nodes[i])

        # set migration between counties
        for i in range(num_nodes):
            for j in range(num_nodes):
                num_coeffs = num_groups * 8
                mig_params = mio.MigrationParameters(
                    num_coeffs, len(list(Location)))

                populations = graph.get_node(i).populations

                # commuters
                min_commuter_age_idx = 2
                # Nur Altersgruppen 3-5 (idx 2-4) pendeln zur Arbeit
                max_commuter_age_idx = 4
                working_population = 0.0
                for age_idx in range(min_commuter_age_idx,
                                     max_commuter_age_idx):
                    working_population += populations.get_group_total(
                        0, age_idx) * (0.33 if age_idx == 4 else 1.0)
                commuter_coeff = commuter[i, j] / working_population
                for age_idx in range(min_commuter_age_idx,
                                     max_commuter_age_idx):
                    for compartment in [
                            State.Susceptible, State.Exposed,
                            State.InfectedNoSymptoms, State.InfectedSymptoms,
                            State.Recovered]:
                        idx = populations.get_flat_index(
                            [age_idx, int(compartment)])
                        mig_params.coefficients[Location.Work.value].baseline[idx] = commuter_coeff * (
                            0.33 if age_idx == 4 else 1.0)

                # twitter
                total_population = populations.get_total()
                twitter_coeff = twitter[i, j] / total_population
                for age_idx in range(min_commuter_age_idx,
                                     max_commuter_age_idx):
                    for compartment in [
                            State.Susceptible, State.Exposed,
                            State.InfectedNoSymptoms, State.InfectedSymptoms,
                            State.Recovered]:
                        idx = populations.get_flat_index(
                            [age_idx, compartment])
                        mig_params.coefficients[Location.Other.value].baseline[idx] = twitter_coeff * (
                            0.33 if age_idx == 4 else 1.0)

                # conservative balancing of ICU
                mig_params.icu_balance.abs_fill_threshold = 0.3
                mig_params.icu_balance.fill_difference_threshold = 0.1
                mig_params.icu_balance.coefficient = 0.0005

                if commuter_coeff > 4e-5 or twitter_coeff > 1e-5:
                    graph.add_edge(i, j, mig_params)

        self._param_study = ParameterStudy(graph, t0, tmax, dt, num_runs=1)

        self._params_nodes = params_nodes
        self._t0 = t0
        self._tmax = tmax
        self._data_dir = data_dir
        self.last_result = None

    def cmp(self):
        county_ids = get_county_ids(
            path.join(self._data_dir, "pydata", "Germany"))
        ts = [TimeSeries(6*8) for id in county_ids]
        # very expensive!
        for t in range(int(self._t0), int(self._tmax) + 1):
            date = self._start_date + datetime.timedelta(days=(t - self._t0))
            params_nodes = read_population_data_county(
                path.join(self._data_dir, "pydata", "Germany"),
                (date.year, date.month, date.day),
                county_ids,
                [2.5], 1.0,
                self._params_nodes)
            for i in range(len(params_nodes)):
                ts[i].add_time_point(
                    t, params_nodes[i].populations.get_compartments())
        return ts

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
    sim = Simulation(
        data_dir=path.join(
            path.dirname(path.abspath(__file__)),
            "..", "..", "data"))
    print(sim.run())
