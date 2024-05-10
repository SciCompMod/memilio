/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#include "matchers.h"
#include "temp_file_register.h"
#include "test_data_dir.h"
#include "memilio/data/analyze_result.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/epidemiology/damping_sampling.h"
#include "memilio/epidemiology/simulation_day.h"
#include "memilio/io/io.h"
#include "memilio/io/result_io.h"
#include "memilio/mobility/graph.h"
#include "memilio/utils/stl_util.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/uncertain_value.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/parameter_space.h"
#include "ode_secirvvs/parameters.h"
#include "ode_secirvvs/parameters_io.h"
#include "ode_secirvvs/analyze_result.h"

#include "gtest/gtest.h"
#include "gmock/gmock-matchers.h"
#include <algorithm>
#include <iterator>
#include <limits>

TEST(TestOdeSECIRVVS, simulateDefault)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.1;

    mio::osecirvvs::Model model(1);
    model.populations.set_total(10);
    model.populations.set_difference_from_total({(mio::AgeGroup)0, mio::osecirvvs::InfectionState::SusceptibleNaive},
                                                10);
    model.parameters.get<mio::osecirvvs::DailyPartialVaccination>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirvvs::DailyPartialVaccination>().array().setConstant(0);
    model.parameters.get<mio::osecirvvs::DailyFullVaccination>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirvvs::DailyFullVaccination>().array().setConstant(0);
    model.parameters.get<mio::osecirvvs::DailyBoosterVaccination>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirvvs::DailyBoosterVaccination>().array().setConstant(0);
    mio::TimeSeries<double> result = simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

TEST(TestOdeSECIRVVS, overflow_vaccinations)
{
    const double t0   = 0;
    const double tmax = 1;
    const double dt   = 1;

    // init simple model
    mio::osecirvvs::Model model(1);
    model.populations[{mio::AgeGroup(0), mio::osecirvvs::InfectionState::SusceptibleNaive}]            = 10.;
    model.populations[{mio::AgeGroup(0), mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}]  = 10.;
    model.populations[{mio::AgeGroup(0), mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}] = 10.;
    model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>()[mio::AgeGroup(0)]         = 0.0;

    // set vaccination rates higher than total population for each layer
    const size_t daily_vaccinations = 100;
    model.parameters.get<mio::osecirvvs::DailyPartialVaccination>().resize(mio::SimulationDay(tmax + 1));
    model.parameters.get<mio::osecirvvs::DailyFullVaccination>().resize(mio::SimulationDay(tmax + 1));
    model.parameters.get<mio::osecirvvs::DailyBoosterVaccination>().resize(mio::SimulationDay(tmax + 1));
    for (size_t i = 0; i <= tmax; ++i) {
        auto num_vaccinations = static_cast<double>((i + 1) * daily_vaccinations);
        model.parameters.get<mio::osecirvvs::DailyPartialVaccination>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
        model.parameters.get<mio::osecirvvs::DailyFullVaccination>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
        model.parameters.get<mio::osecirvvs::DailyBoosterVaccination>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
    }

    // simulate one step with explicit Euler
    auto integrator = std::make_shared<mio::EulerIntegratorCore>();
    auto result     = simulate_flows(t0, tmax, dt, model, integrator);

    // get the flow indices for each type of vaccination and also the indices of the susceptible compartments
    auto flow_indx_partial_vaccination =
        model.get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptibleNaive,
                                  mio::osecirvvs::InfectionState::TemporaryImmunPartialImmunity>({mio::AgeGroup(0)});
    auto flow_indx_full_vaccination =
        model.get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptiblePartialImmunity,
                                  mio::osecirvvs::InfectionState::TemporaryImmunImprovedImmunity>({mio::AgeGroup(0)});
    auto flow_indx_booster_vaccination =
        model.get_flat_flow_index<mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity,
                                  mio::osecirvvs::InfectionState::TemporaryImmunImprovedImmunity>({mio::AgeGroup(0)});
    auto indx_S_naive =
        model.populations.get_flat_index({mio::AgeGroup(0), mio::osecirvvs::InfectionState::SusceptibleNaive});
    auto indx_S_partial = model.populations.get_flat_index(
        {mio::AgeGroup(0), mio::osecirvvs::InfectionState::SusceptiblePartialImmunity});
    auto indx_S_improved = model.populations.get_flat_index(
        {mio::AgeGroup(0), mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity});

    // check that the number of vaccinated people is never higher than the total number of susceptible people
    EXPECT_NEAR(result[1].get_last_value()[flow_indx_partial_vaccination], result[0].get_value(0)[indx_S_naive], 1e-10);
    EXPECT_NEAR(result[1].get_last_value()[flow_indx_full_vaccination], result[0].get_value(0)[indx_S_partial], 1e-10);
    EXPECT_NEAR(result[1].get_last_value()[flow_indx_booster_vaccination], result[0].get_value(0)[indx_S_improved],
                1e-10);
}

TEST(TestOdeSECIRVVS, smooth_vaccination_rate)
{
    const ScalarType tmax = 2.;

    // init simple model
    mio::osecirvvs::Model model(1);
    auto& daily_vaccinations = model.parameters.get<mio::osecirvvs::DailyPartialVaccination>();
    daily_vaccinations.resize(mio::SimulationDay(tmax + 1));

    daily_vaccinations[{mio::AgeGroup(0), mio::SimulationDay(0)}] = 0;
    daily_vaccinations[{mio::AgeGroup(0), mio::SimulationDay(1)}] = 10;
    daily_vaccinations[{mio::AgeGroup(0), mio::SimulationDay(2)}] = 110;

    const auto eps1 = 0.15;

    // test when t is out of the range
    Eigen::VectorXd result = model.vaccinations_at(5, daily_vaccinations, eps1);
    EXPECT_EQ(result.size(), 1);
    EXPECT_NEAR(result[0], 0, 1e-12);

    // test when t i below the lower bound
    result = model.vaccinations_at(0.5, daily_vaccinations, eps1);
    EXPECT_EQ(result.size(), 1);
    EXPECT_NEAR(result[0], 10, 1e-12);

    result = model.vaccinations_at(1.5, daily_vaccinations, eps1);
    EXPECT_EQ(result.size(), 1);
    EXPECT_NEAR(result[0], 100, 1e-12);

    // test when t is withing the range of the smoothing
    result = model.vaccinations_at(0.85, daily_vaccinations, eps1);
    EXPECT_NEAR(result[0], 10.0, 1e-12);

    result = model.vaccinations_at(0.90, daily_vaccinations, eps1);
    EXPECT_NEAR(result[0], 32.5, 1e-12);

    result = model.vaccinations_at(0.95, daily_vaccinations, eps1);
    EXPECT_NEAR(result[0], 77.5, 1e-12);

    result = model.vaccinations_at(1.0, daily_vaccinations, eps1);
    EXPECT_NEAR(result[0], 100.0, 1e-12);

    // Test also with a different epsilon
    const auto eps2 = 0.4;
    result          = model.vaccinations_at(0.6, daily_vaccinations, eps2);
    EXPECT_NEAR(result[0], 10.0, 1e-12);

    result = model.vaccinations_at(0.8, daily_vaccinations, eps2);
    EXPECT_NEAR(result[0], 55, 1e-12);

    result = model.vaccinations_at(1., daily_vaccinations, eps2);
    EXPECT_NEAR(result[0], 100.0, 1e-12);
}

void assign_uniform_distribution(mio::UncertainValue& p, double min, double max, bool set_invalid_initial_value)
{
    auto invalid_initial = max == 0 ? 1.0 : max * 1.001;
    auto valid_initial   = (max + min) * 0.5;
    auto initial         = set_invalid_initial_value ? invalid_initial : valid_initial;
    p                    = mio::UncertainValue(initial);
    p.set_distribution(mio::ParameterDistributionUniform(min, max));
}

template <size_t N>
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue, mio::AgeGroup>& array,
                                       const double (&min)[N], const double (&max)[N], bool set_invalid_initial_value)
{
    for (auto i = mio::AgeGroup(0); i < array.size<mio::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min[size_t(i)], max[size_t(i)], set_invalid_initial_value);
    }
}

void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue, mio::AgeGroup>& array, double min,
                                       double max, bool set_invalid_initial_value)
{
    for (auto i = mio::AgeGroup(0); i < array.size<mio::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min, max, set_invalid_initial_value);
    }
}

void set_synthetic_population_data(mio::osecirvvs::Model::Populations& populations, bool set_invalid_initial_value)
{
    for (mio::AgeGroup i = 0; i < mio::get<mio::AgeGroup>(populations.size()); i++) {
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::ExposedNaive}], 10, 20,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::ExposedImprovedImmunity}], 10, 21,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::ExposedPartialImmunity}], 10, 22,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}], 1, 10,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}],
                                    2, 10, set_invalid_initial_value);
        assign_uniform_distribution(
            populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}], 3, 10,
            set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}], 4, 10,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}],
                                    5, 10, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}],
                                    6, 10, set_invalid_initial_value);

        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}],
                                    5, 11, set_invalid_initial_value);
        assign_uniform_distribution(
            populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}], 5, 12,
            set_invalid_initial_value);
        assign_uniform_distribution(
            populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}], 5, 13,
            set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}], 5,
                                    14, set_invalid_initial_value);
        assign_uniform_distribution(
            populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}], 5, 15,
            set_invalid_initial_value);
        assign_uniform_distribution(
            populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}], 5, 16,
            set_invalid_initial_value);

        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedSevereNaive}], 1, 2,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedSevereImprovedImmunity}], 1,
                                    3, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedSeverePartialImmunity}], 1,
                                    4, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalNaive}], 1, 5,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalPartialImmunity}],
                                    1, 6, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalImprovedImmunity}],
                                    1, 7, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}], 200,
                                    300, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}], 200,
                                    400, set_invalid_initial_value);
        populations.set_difference_from_group_total<mio::AgeGroup>(
            {i, mio::osecirvvs::InfectionState::SusceptibleNaive}, 1000);
    }
}

void set_demographic_parameters(mio::osecirvvs::Model::ParameterSet& parameters, bool set_invalid_initial_value)
{
    assign_uniform_distribution(parameters.get<mio::osecirvvs::ICUCapacity>(), 20, 50, set_invalid_initial_value);
    assign_uniform_distribution(parameters.get<mio::osecirvvs::TestAndTraceCapacity>(), 100, 200,
                                set_invalid_initial_value);
    parameters.get<mio::osecirvvs::DailyPartialVaccination>().resize(mio::SimulationDay(size_t(1000)));
    parameters.get<mio::osecirvvs::DailyPartialVaccination>().array().setConstant(5);
    parameters.get<mio::osecirvvs::DailyFullVaccination>().resize(mio::SimulationDay(size_t(1000)));
    parameters.get<mio::osecirvvs::DailyFullVaccination>().array().setConstant(3);
    parameters.get<mio::osecirvvs::DailyBoosterVaccination>().resize(mio::SimulationDay(size_t(1000)));
    parameters.get<mio::osecirvvs::DailyBoosterVaccination>().array().setConstant(3);
}

void set_contact_parameters(mio::osecirvvs::Model::ParameterSet& parameters, bool set_invalid_initial_value)
{
    auto& contacts       = parameters.get<mio::osecirvvs::ContactPatterns>();
    auto& contact_matrix = contacts.get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(0.5);
    contact_matrix[0].get_baseline().diagonal().setConstant(5.0);
    contact_matrix[0].add_damping(0.3, mio::SimulationTime(5.0));

    auto& npis      = parameters.get<mio::osecirvvs::DynamicNPIsInfectedSymptoms>();
    auto npi_groups = Eigen::VectorXd::Ones(contact_matrix[0].get_num_groups());
    auto npi_value  = mio::UncertainValue(0.5);
    assign_uniform_distribution(npi_value, 0.25, 0.75, set_invalid_initial_value);
    npis.set_threshold(10.0, {mio::DampingSampling(npi_value, mio::DampingLevel(0), mio::DampingType(0),
                                                   mio::SimulationTime(0), {0}, npi_groups)});
    npis.set_base_value(100'000);
    npis.set_interval(mio::SimulationTime(3.0));
    npis.set_duration(mio::SimulationTime(14.0));
    parameters.get_end_dynamic_npis() = 10.0; //required for dynamic NPIs to have effect in this model
}

void set_covid_parameters(mio::osecirvvs::Model::ParameterSet& params, bool set_invalid_initial_value)
{
    //times
    const double timeExposedMin            = 2.67;
    const double timeExposedMax            = 4.;
    const double timeInfectedNoSymptomsMin = 1.2;
    const double timeInfectedNoSymptomsMax = 2.53;
    const double timeInfectedSymptomsMin[] = {5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465};
    const double timeInfectedSymptomsMax[] = {8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085};
    const double timeInfectedSevereMin[]   = {3.925, 3.925, 4.85, 6.4, 7.2, 9.};
    const double timeInfectedSevereMax[]   = {6.075, 6.075, 7., 8.7, 9.8, 13.};
    const double timeInfectedCriticalMin[] = {4.95, 4.95, 4.86, 14.14, 14.4, 10.};
    const double timeInfectedCriticalMax[] = {8.95, 8.95, 8.86, 20.58, 19.8, 13.2};

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeExposed>(), timeExposedMin, timeExposedMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedNoSymptoms>(), timeInfectedNoSymptomsMin,
                                      timeInfectedNoSymptomsMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSymptoms>(), timeInfectedSymptomsMin,
                                      timeInfectedSymptomsMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSevere>(), timeInfectedSevereMin,
                                      timeInfectedSevereMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedCritical>(), timeInfectedCriticalMin,
                                      timeInfectedCriticalMax, set_invalid_initial_value);

    //probabilities
    double fac_variant                                 = 1.4;
    const double transmissionProbabilityOnContactMin[] = {0.02 * fac_variant, 0.05 * fac_variant, 0.05 * fac_variant,
                                                          0.05 * fac_variant, 0.08 * fac_variant, 0.1 * fac_variant};

    const double transmissionProbabilityOnContactMax[] = {0.04 * fac_variant, 0.07 * fac_variant, 0.07 * fac_variant,
                                                          0.07 * fac_variant, 0.10 * fac_variant, 0.15 * fac_variant};
    const double relativeTransmissionNoSymptomsMin     = 0.5;
    const double relativeTransmissionNoSymptomsMax     = 0.5;
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    const double riskOfInfectionFromSymptomaticMin    = 0.0;
    const double riskOfInfectionFromSymptomaticMax    = 0.2;
    const double maxRiskOfInfectionFromSymptomaticMin = 0.4;
    const double maxRiskOfInfectionFromSymptomaticMax = 0.5;
    const double recoveredPerInfectedNoSymptomsMin[]  = {0.2, 0.2, 0.15, 0.15, 0.15, 0.15};
    const double recoveredPerInfectedNoSymptomsMax[]  = {0.3, 0.3, 0.25, 0.25, 0.25, 0.25};
    const double severePerInfectedSymptomsMin[]       = {0.006, 0.006, 0.015, 0.049, 0.15, 0.20};
    const double severePerInfectedSymptomsMax[]       = {0.009, 0.009, 0.023, 0.074, 0.18, 0.25};
    const double criticalPerSevereMin[]               = {0.05, 0.05, 0.05, 0.10, 0.25, 0.35};
    const double criticalPerSevereMax[]               = {0.10, 0.10, 0.10, 0.20, 0.35, 0.45};
    const double deathsPerCriticalMin[]               = {0.00, 0.00, 0.10, 0.10, 0.30, 0.5};
    const double deathsPerCriticalMax[]               = {0.10, 0.10, 0.18, 0.18, 0.50, 0.7};

    const double reducExposedPartialImmunityMin                     = 0.75;
    const double reducExposedPartialImmunityMax                     = 0.85;
    const double reducExposedImprovedImmunityMin                    = 0.281;
    const double reducExposedImprovedImmunityMax                    = 0.381;
    const double reducInfectedSymptomsPartialImmunityMin            = 0.6;
    const double reducInfectedSymptomsPartialImmunityMax            = 0.7;
    const double reducInfectedSymptomsImprovedImmunityMin           = 0.193;
    const double reducInfectedSymptomsImprovedImmunityMax           = 0.293;
    const double reducInfectedSevereCriticalDeadPartialImmunityMin  = 0.05;
    const double reducInfectedSevereCriticalDeadPartialImmunityMax  = 0.15;
    const double reducInfectedSevereCriticalDeadImprovedImmunityMin = 0.041;
    const double reducInfectedSevereCriticalDeadImprovedImmunityMax = 0.141;
    const double reducTimeInfectedMildMin                           = 0.8;
    const double reducTimeInfectedMildMax                           = 1.0;

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TransmissionProbabilityOnContact>(),
                                      transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RelativeTransmissionNoSymptoms>(),
                                      relativeTransmissionNoSymptomsMin, relativeTransmissionNoSymptomsMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic>(),
                                      riskOfInfectionFromSymptomaticMin, riskOfInfectionFromSymptomaticMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic>(),
                                      maxRiskOfInfectionFromSymptomaticMin, maxRiskOfInfectionFromSymptomaticMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms>(),
                                      recoveredPerInfectedNoSymptomsMin, recoveredPerInfectedNoSymptomsMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::SeverePerInfectedSymptoms>(),
                                      severePerInfectedSymptomsMin, severePerInfectedSymptomsMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::CriticalPerSevere>(), criticalPerSevereMin,
                                      criticalPerSevereMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::DeathsPerCritical>(), deathsPerCriticalMin,
                                      deathsPerCriticalMax, set_invalid_initial_value);

    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedPartialImmunity>(),
                                      reducExposedPartialImmunityMin, reducExposedPartialImmunityMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedImprovedImmunity>(),
                                      reducExposedImprovedImmunityMin, reducExposedImprovedImmunityMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity>(),
                                      reducInfectedSymptomsPartialImmunityMin, reducInfectedSymptomsPartialImmunityMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity>(),
                                      reducInfectedSymptomsImprovedImmunityMin,
                                      reducInfectedSymptomsImprovedImmunityMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity>(),
                                      reducInfectedSevereCriticalDeadPartialImmunityMin,
                                      reducInfectedSevereCriticalDeadPartialImmunityMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity>(),
                                      reducInfectedSevereCriticalDeadImprovedImmunityMin,
                                      reducInfectedSevereCriticalDeadImprovedImmunityMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducTimeInfectedMild>(), reducTimeInfectedMildMin,
                                      reducTimeInfectedMildMax, set_invalid_initial_value);

    //sasonality
    const double seasonality_min = 0.1;
    const double seasonality_max = 0.3;

    assign_uniform_distribution(params.get<mio::osecirvvs::Seasonality>(), seasonality_min, seasonality_max,
                                set_invalid_initial_value);
}

mio::osecirvvs::Model make_model(int num_age_groups, bool set_invalid_initial_value = false)
{
    assert(num_age_groups <= 6 && "Provide more values in functions above to test more age groups.");
    mio::osecirvvs::Model model(num_age_groups);
    set_covid_parameters(model.parameters, set_invalid_initial_value);
    set_synthetic_population_data(model.populations, set_invalid_initial_value);
    set_demographic_parameters(model.parameters, set_invalid_initial_value);
    set_contact_parameters(model.parameters, set_invalid_initial_value);
    model.parameters.apply_constraints();
    return model;
}

TEST(TestOdeSECIRVVS, draw_sample_graph)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters> graph;

    auto num_age_groups = 6;
    //create model with invalid initials so the test fails if no sampling is done
    graph.add_node(0, make_model(num_age_groups, /*set_invalid_initial_value*/ true));
    graph.add_node(0, make_model(num_age_groups, /*set_invalid_initial_value*/ true));
    graph.add_edge(0, 1, Eigen::VectorXd::Constant(num_age_groups, num_age_groups));

    auto sampled_graph = mio::osecirvvs::draw_sample(graph);

    ASSERT_EQ(sampled_graph.nodes().size(), graph.nodes().size());
    ASSERT_EQ(sampled_graph.edges().size(), graph.edges().size());

    // spot check for sampling
    auto& parameters0          = sampled_graph.nodes()[0].property.parameters;
    auto& populations0         = sampled_graph.nodes()[0].property.populations;
    auto& timeInfectedCritical = parameters0.get<mio::osecirvvs::TimeInfectedCritical>()[mio::AgeGroup(1)];
    ASSERT_GE(double(timeInfectedCritical), 4.95);
    ASSERT_LE(double(timeInfectedCritical), 8.95);
    auto& param_exp_factor = parameters0.get<mio::osecirvvs::ReducExposedPartialImmunity>()[mio::AgeGroup(0)];
    ASSERT_GE(double(param_exp_factor), 0.75);
    ASSERT_LE(double(param_exp_factor), 0.85);
    auto& compartment_inf =
        populations0[{mio::AgeGroup(2), mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}];
    ASSERT_GE(double(compartment_inf), 5.0);
    ASSERT_LE(double(compartment_inf), 10.0);
    auto& npi_value =
        parameters0.get<mio::osecirvvs::DynamicNPIsInfectedSymptoms>().get_thresholds()[0].second[0].get_value();
    ASSERT_GE(double(npi_value), 0.25);
    ASSERT_LE(double(npi_value), 0.75);

    // special cases
    ASSERT_NEAR(populations0.get_total(), 1000 * num_age_groups, 1e-2);
    ASSERT_TRUE((parameters0.get<mio::osecirvvs::InfectiousnessNewVariant>().array(),
                 parameters0.get<mio::osecirvvs::TransmissionProbabilityOnContact>().array() * 1.0) //using high variant
                    .all());

    // spot check for parameters that should be equal or different between nodes
    auto& parameters1  = sampled_graph.nodes()[1].property.parameters;
    auto& populations1 = sampled_graph.nodes()[1].property.populations;
    ASSERT_EQ(parameters1.get<mio::osecirvvs::DynamicNPIsInfectedSymptoms>().get_thresholds()[0].second[0].get_value(),
              parameters0.get<mio::osecirvvs::DynamicNPIsInfectedSymptoms>().get_thresholds()[0].second[0].get_value());
    ASSERT_TRUE((parameters1.get<mio::osecirvvs::TimeInfectedSymptoms>().array() ==
                 parameters0.get<mio::osecirvvs::TimeInfectedSymptoms>().array())
                    .all());
    //these could fail in very(!) rare cases if they are randomly sampled to the same value
    ASSERT_NE(parameters1.get<mio::osecirvvs::ICUCapacity>(), parameters0.get<mio::osecirvvs::ICUCapacity>())
        << "Failure might be spurious, check RNG seeds.";
    ASSERT_FALSE((populations1.array() == populations0.array()).all()) << "Failure might be spurious, check RNG seeds.";
}

TEST(TestOdeSECIRVVS, draw_sample_model)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    auto num_age_groups = 6;
    auto model          = make_model(num_age_groups, /*set_invalid_initial_value*/ true);
    mio::osecirvvs::draw_sample(model);

    // spot check for sampling
    auto& parameters           = model.parameters;
    auto& populations          = model.populations;
    auto& timeInfectedCritical = parameters.get<mio::osecirvvs::TimeInfectedCritical>()[mio::AgeGroup(1)];
    ASSERT_GE(double(timeInfectedCritical), 4.95);
    ASSERT_LE(double(timeInfectedCritical), 8.95);
    auto& param_exp_factor = parameters.get<mio::osecirvvs::ReducExposedPartialImmunity>()[mio::AgeGroup(0)];
    ASSERT_GE(double(param_exp_factor), 0.75);
    ASSERT_LE(double(param_exp_factor), 0.85);
    auto& compartment_inf =
        populations[{mio::AgeGroup(2), mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}];
    ASSERT_GE(double(compartment_inf), 5.0);
    ASSERT_LE(double(compartment_inf), 10.0);

    // special cases
    ASSERT_NEAR(populations.get_total(), 1000 * num_age_groups, 1e-2);
    ASSERT_TRUE((parameters.get<mio::osecirvvs::InfectiousnessNewVariant>().array(),
                 parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>().array() * 1.0)
                    .all());
}

TEST(TestOdeSECIRVVS, checkPopulationConservation)
{
    auto num_age_groups = 2;
    auto model          = make_model(num_age_groups);
    auto num_days       = 30;

    auto result = mio::osecirvvs::simulate(0, num_days, 0.1, model);

    double num_persons = 0.0;
    for (auto i = 0; i < result.get_last_value().size(); i++) {
        EXPECT_GE(result.get_last_value()[i], -1e-3);
        num_persons += result.get_last_value()[i];
    }
    EXPECT_NEAR(num_persons, model.populations.get_total(), 1e-10);
}

#if defined(MEMILIO_HAS_HDF5) && defined(MEMILIO_HAS_JSONCPP)

TEST(TestOdeSECIRVVS, read_confirmed_cases)
{
    auto num_age_groups = 6; //reading data requires RKI data age groups
    auto model          = std::vector<mio::osecirvvs::Model>({make_model(num_age_groups)});
    std::vector<int> region{1002};
    auto path = mio::path_join(TEST_DATA_DIR, "pydata/Germany/cases_all_county_age_ma7.json");
    std::vector<std::vector<int>> t_Exposed(1);
    std::vector<std::vector<int>> t_InfectedNoSymptoms(1);
    std::vector<std::vector<int>> t_InfectedSymptoms(1);
    std::vector<std::vector<int>> t_InfectedSevere(1);
    std::vector<std::vector<int>> t_InfectedCritical(1);
    std::vector<std::vector<int>> t_imm_interval1(1);

    std::vector<std::vector<double>> mu_C_R(1);
    std::vector<std::vector<double>> mu_I_H(1);
    std::vector<std::vector<double>> mu_H_U(1);

    std::vector<std::vector<double>> num_InfectedSymptoms(1);
    std::vector<std::vector<double>> num_death(1);
    std::vector<std::vector<double>> num_timm(1);
    std::vector<std::vector<double>> num_Exposed(1);
    std::vector<std::vector<double>> num_InfectedNoSymptoms(1);
    std::vector<std::vector<double>> num_InfectedSevere(1);
    std::vector<std::vector<double>> num_icu(1);

    num_InfectedSymptoms[0]   = std::vector<double>(num_age_groups, 0.0);
    num_death[0]              = std::vector<double>(num_age_groups, 0.0);
    num_timm[0]               = std::vector<double>(num_age_groups, 0.0);
    num_Exposed[0]            = std::vector<double>(num_age_groups, 0.0);
    num_InfectedNoSymptoms[0] = std::vector<double>(num_age_groups, 0.0);
    num_InfectedSevere[0]     = std::vector<double>(num_age_groups, 0.0);
    num_icu[0]                = std::vector<double>(num_age_groups, 0.0);
    for (size_t group = 0; group < static_cast<size_t>(num_age_groups); group++) {

        t_Exposed[0].push_back(static_cast<int>(
            std::round(model[0].parameters.template get<mio::osecirvvs::TimeExposed>()[(mio::AgeGroup)group])));
        t_InfectedNoSymptoms[0].push_back(static_cast<int>(std::round(
            model[0].parameters.template get<mio::osecirvvs::TimeInfectedNoSymptoms>()[(mio::AgeGroup)group])));
        t_InfectedSymptoms[0].push_back(static_cast<int>(std::round(
            model[0].parameters.template get<mio::osecirvvs::TimeInfectedSymptoms>()[(mio::AgeGroup)group])));
        t_InfectedSevere[0].push_back(static_cast<int>(
            std::round(model[0].parameters.template get<mio::osecirvvs::TimeInfectedSevere>()[(mio::AgeGroup)group])));
        t_InfectedCritical[0].push_back(static_cast<int>(std::round(
            model[0].parameters.template get<mio::osecirvvs::TimeInfectedCritical>()[(mio::AgeGroup)group])));
        t_imm_interval1[0].push_back(static_cast<int>(std::round(
            model[0].parameters.template get<mio::osecirvvs::TimeTemporaryImmunityPI>()[(mio::AgeGroup)group])));

        mu_C_R[0].push_back(
            model[0].parameters.template get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)group]);
        mu_I_H[0].push_back(
            model[0].parameters.template get<mio::osecirvvs::SeverePerInfectedSymptoms>()[(mio::AgeGroup)group]);
        mu_H_U[0].push_back(
            model[0].parameters.template get<mio::osecirvvs::CriticalPerSevere>()[(mio::AgeGroup)group]);
    }

    auto read = mio::osecirvvs::details::read_confirmed_cases_data(
        path, region, {2020, 12, 01}, num_Exposed, num_InfectedNoSymptoms, num_InfectedSymptoms, num_InfectedSevere,
        num_icu, num_death, num_timm, t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
        t_InfectedCritical, t_imm_interval1, mu_C_R, mu_I_H, mu_H_U, std::vector<double>(size_t(num_age_groups), 1.0));

    ASSERT_THAT(read, IsSuccess());
}

TEST(TestOdeSECIRVVS, read_data)
{
    auto num_age_groups = 6; //reading data requires RKI data age groups
    auto model1         = std::vector<mio::osecirvvs::Model>({make_model(num_age_groups)});
    auto model2         = std::vector<mio::osecirvvs::Model>({make_model(num_age_groups)});
    auto model3         = std::vector<mio::osecirvvs::Model>({make_model(num_age_groups)});

    auto read_result1 = mio::osecirvvs::read_input_data_county(
        model1, {2020, 12, 01}, {1002}, std::vector<double>(size_t(num_age_groups), 1.0), 1.0, TEST_DATA_DIR, 10);

    auto read_result2 = mio::osecirvvs::read_input_data(
        model2, {2020, 12, 01}, {1002}, std::vector<double>(size_t(num_age_groups), 1.0), 1.0, TEST_DATA_DIR, 10);

    auto read_result_district = mio::osecirvvs::read_input_data(model3, {2020, 12, 01}, {1002},
                                                                std::vector<double>(size_t(num_age_groups), 1.0), 1.0,
                                                                mio::path_join(TEST_DATA_DIR, "pydata/District"), 10);

    ASSERT_THAT(read_result1, IsSuccess());
    ASSERT_THAT(read_result2, IsSuccess());
    ASSERT_THAT(read_result_district, IsSuccess());

    // values were generated by the tested function; can only check stability of the implementation, not correctness
    auto expected_values =
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecirvvs::InfectionState::Count)) << 410.72, 3630.45,
         0.213004, 3.19834, 0.840323, 0.0946686, 1.42149, 0.373477, 0, 0, 0, 0.558774, 4.4528, 0.955134, 0, 0, 0,
         0.0103231, 0.00287609, 0.00365536, 0.0297471, 3.5, 4, 2893.07, 0, 0, 0, 2652.78, 714.319, 763.188, 8032.87,
         0.49701, 7.4628, 1.96075, 0.319507, 4.79751, 1.26048, 0, 0, 0, 0.93129, 8.36257, 1.79379, 0, 0, 0, 0.0132265,
         0.00302746, 0.00384774, 0.0297471, 3.5, 4, 5916.15, 0, 0, 0, 3658.05, 799.805, 5554.16, 43483, 5.32163,
         43.3153, 9.75737, 2.7505, 22.3877, 5.04313, 0, 0, 0, 8.84909, 39.8328, 7.32558, 0, 0, 0, 0.56467, 0.071343,
         0.0777408, 0.0753594, 3.5, 4, 21861.8, 0, 0, 0, 2589.83, 845.38, 6617.66, 48070, 5.33096, 40.6794, 9.01336,
         2.83472, 21.6311, 4.79282, 0, 0, 0, 8.55241, 37.2832, 6.74428, 0, 0, 0, 1.78101, 0.218787, 0.234499, 0.487853,
         3.5, 4, 24154, 0, 0, 0, 3367.68, 774.244, 1528.03, 21414.7, 0.528786, 8.62793, 2.62254, 0.329244, 5.37211,
         1.6329, 0, 0, 0, 0.939561, 8.52246, 2.1149, 0, 0, 0, 0.856992, 0.223414, 0.328498, 2.61775, 3.5, 4, 16069.9, 0,
         0, 0, 4038.53, 836.776, 144.225, 1979.35, 0.22575, 9.11334, 5.90344, 0.0707574, 2.85642, 1.85033, 0, 0, 0,
         0.0889777, 2.09765, 1.10935, 0, 0, 0, 0.0681386, 0.046887, 0.146922, 4.75954, 3.5, 4, 8288.72, 0, 0, 0,
         4447.91, 823.157)
            .finished();

    ASSERT_THAT(print_wrap(model1[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(model2[0].populations.array().cast<double>()), 1e-5, 1e-5));

    ASSERT_THAT(print_wrap(model1[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(model3[0].populations.array().cast<double>()), 1e-5, 1e-5));

    ASSERT_THAT(print_wrap(model1[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));
    ASSERT_THAT(print_wrap(model2[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));
    ASSERT_THAT(print_wrap(model3[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));

    // some more tests which are actually not necessary but can help if the previous tests fails or needs to get changed
    for (mio::AgeGroup i = 0; i < model1[0].parameters.get_num_groups(); i++) {
        EXPECT_GE(double(model1[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}]), 0);
        EXPECT_GE(double(model2[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}]), 0);
        EXPECT_GE(double(model3[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}]), 0);
        EXPECT_GE(double(model1[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}]), 0);
        EXPECT_GE(double(model2[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}]), 0);
        EXPECT_GE(double(model3[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}]), 0);
        EXPECT_GE(double(model1[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(double(model2[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(double(model3[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(double(model1[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(double(model2[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(double(model3[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(
            double(model1[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}]), 0);
        EXPECT_GE(
            double(model2[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}]), 0);
        EXPECT_GE(
            double(model3[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}]), 0);
        EXPECT_GE(double(model1[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}]),
                  0);
        EXPECT_GE(double(model2[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}]),
                  0);
        EXPECT_GE(double(model3[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}]),
                  0);
        EXPECT_GE(double(model2[0].populations[{i, mio::osecirvvs::InfectionState::TemporaryImmunPartialImmunity}]), 0);
        EXPECT_GE(double(model3[0].populations[{i, mio::osecirvvs::InfectionState::TemporaryImmunImprovedImmunity}]),
                  0);

        // currently dead and confirmed after commuting compartments are initialized as zero
        EXPECT_EQ(double(model1[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}]),
                  0);
        EXPECT_EQ(double(model2[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}]),
                  0);
        EXPECT_EQ(double(model3[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}]),
                  0);
        EXPECT_EQ(double(model1[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}]),
                  0);
        EXPECT_EQ(double(model2[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}]),
                  0);
        EXPECT_EQ(double(model3[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}]),
                  0);
        EXPECT_EQ(
            double(
                model1[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model2[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model3[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model1[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model2[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model3[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(model1[0]
                       .populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(model2[0]
                       .populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(model3[0]
                       .populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model1[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model2[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model3[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(double(model1[0].populations[{i, mio::osecirvvs::InfectionState::DeadNaive}]), 0);
        EXPECT_EQ(double(model2[0].populations[{i, mio::osecirvvs::InfectionState::DeadNaive}]), 0);
        EXPECT_EQ(double(model3[0].populations[{i, mio::osecirvvs::InfectionState::DeadNaive}]), 0);
        EXPECT_EQ(double(model1[0].populations[{i, mio::osecirvvs::InfectionState::DeadPartialImmunity}]), 0);
        EXPECT_EQ(double(model2[0].populations[{i, mio::osecirvvs::InfectionState::DeadPartialImmunity}]), 0);
        EXPECT_EQ(double(model3[0].populations[{i, mio::osecirvvs::InfectionState::DeadPartialImmunity}]), 0);
        EXPECT_EQ(double(model1[0].populations[{i, mio::osecirvvs::InfectionState::DeadImprovedImmunity}]), 0);
        EXPECT_EQ(double(model2[0].populations[{i, mio::osecirvvs::InfectionState::DeadImprovedImmunity}]), 0);
        EXPECT_EQ(double(model3[0].populations[{i, mio::osecirvvs::InfectionState::DeadImprovedImmunity}]), 0);
    }
}

TEST(TestOdeSECIRVVS, export_time_series_init)
{
    TempFileRegister temp_file_register;
    auto tmp_results_dir = temp_file_register.get_unique_path();
    ASSERT_THAT(mio::create_directory(tmp_results_dir), IsSuccess());

    auto num_age_groups = 6; // Data to be read requires RKI confirmed cases data age groups
    auto model          = make_model(num_age_groups);

    // Test exporting time series
    ASSERT_THAT(mio::osecirvvs::export_input_data_county_timeseries(
                    std::vector<mio::osecirvvs::Model>{model}, tmp_results_dir, {0}, {2020, 12, 01},
                    std::vector<double>(size_t(num_age_groups), 1.0), 1.0, 2,
                    mio::path_join(TEST_DATA_DIR, "county_divi_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "cases_all_county_age_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "county_current_population.json"), true,
                    mio::path_join(TEST_DATA_DIR, "vacc_county_ageinf_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "immunity_population.txt")),
                IsSuccess());

    auto data_extrapolated = mio::read_result(mio::path_join(tmp_results_dir, "Results_rki.h5"));
    ASSERT_THAT(data_extrapolated, IsSuccess());

    // Values were generated by the tested function export_input_data_county_timeseries;
    // can only check stability of the implementation, not correctness
    auto expected_results =
        mio::read_result(mio::path_join(TEST_DATA_DIR, "export_time_series_initialization_osecirts.h5")).value();

    ASSERT_THAT(print_wrap(data_extrapolated.value()[0].get_groups().matrix()),
                MatrixNear(print_wrap(expected_results[0].get_groups().matrix()), 1e-5, 1e-5));
}

// Model initialization should return same start values as export time series on that day
TEST(TestOdeSECIRVVS, model_initialization)
{
    auto num_age_groups = 6; // Data to be read requires RKI confirmed cases data age groups
    auto model          = make_model(num_age_groups);
    // Vector assignment necessary as read_input_data_county changes model
    auto model_vector = std::vector<mio::osecirvvs::Model>{model};

    ASSERT_THAT(mio::osecirvvs::read_input_data_county(model_vector, {2020, 12, 01}, {0},
                                                       std::vector<double>(size_t(num_age_groups), 1.0), 1.0,
                                                       TEST_DATA_DIR, 2, false),
                IsSuccess());

    // Values from data/export_time_series_init_osecirvvs.h5, for reading in comparison
    // operator for return of mio::read_result and model population needed.
    auto expected_values =
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecirvvs::InfectionState::Count)) << 138748, 2.11327e+06,
         0.213004, 3.19834, 0.840323, 0.0946686, 1.42149, 0.373477, 0, 0, 0, 0.558774, 4.4528, 0.955134, 0, 0, 0,
         0.0103231, 0.00287609, 0.00365536, 0.0297471, 3.5, 4, 1.21334e+06, 0, 0, 0, 2652.78, 714.319, 309974,
         4.7235e+06, 0.49701, 7.4628, 1.96075, 0.319507, 4.79751, 1.26048, 0, 0, 0, 0.93129, 8.36257, 1.79379, 0, 0, 0,
         0.0132265, 0.00302746, 0.00384774, 0.0297471, 3.5, 4, 2.71151e+06, 0, 0, 0, 3658.05, 799.805, 1.44128e+06,
         1.19121e+07, 5.32163, 43.3153, 9.75737, 2.7505, 22.3877, 5.04313, 0, 0, 0, 8.84909, 39.8328, 7.32558, 0, 0, 0,
         0.56467, 0.071343, 0.0777408, 0.0753594, 3.5, 4, 5.86047e+06, 0, 0, 0, 2589.83, 845.38, 2.40269e+06,
         1.86176e+07, 5.33096, 40.6794, 9.01336, 2.83472, 21.6311, 4.79282, 0, 0, 0, 8.55241, 37.2832, 6.74428, 0, 0, 0,
         1.78101, 0.218787, 0.234499, 0.487853, 3.5, 4, 9.00942e+06, 0, 0, 0, 3367.68, 774.244, 578004, 9.57445e+06,
         0.528786, 8.62793, 2.62254, 0.329244, 5.37211, 1.6329, 0, 0, 0, 0.939561, 8.52246, 2.1149, 0, 0, 0, 0.856992,
         0.223414, 0.328498, 2.61775, 3.5, 4, 6.35731e+06, 0, 0, 0, 4038.53, 836.776, 61810.3, 2.53029e+06, 0.22575,
         9.11334, 5.90344, 0.0707574, 2.85642, 1.85033, 0, 0, 0, 0.0889777, 2.09765, 1.10935, 0, 0, 0, 0.0681386,
         0.046887, 0.146922, 4.75954, 3.5, 4, 3.58492e+06, 0, 0, 0, 4447.91, 823.157)
            .finished();

    ASSERT_THAT(print_wrap(model_vector[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));
}

TEST(TestOdeSECIRVVS, run_simulation)
{
    auto num_age_groups = 3;
    auto model          = make_model(num_age_groups);
    auto num_days       = 30;

    auto result = mio::osecirvvs::simulate(0, num_days, 0.1, model);
    result = mio::interpolate_simulation_result(result); // Reduce influence of time steps chosen by the integrator.

    // Load result of a previous run; only tests stability, not correctness.
    auto expected_result =
        mio::read_result(mio::path_join(TEST_DATA_DIR, "results_osecirts.h5")).value()[0].get_groups();

    ASSERT_THAT(print_wrap(result.matrix()), MatrixNear(print_wrap(expected_result.matrix()), 1e-5, 1e-5));
}

#endif

TEST(TestOdeSECIRVVS, parameter_percentiles)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    //build small graph
    auto model = make_model(5);
    auto graph = mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters>();
    graph.add_node(0, model);

    //sample a few times
    auto sampled_graphs = std::vector<mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters>>();
    std::generate_n(std::back_inserter(sampled_graphs), 10, [&graph]() {
        return mio::osecirvvs::draw_sample(graph);
    });

    //extract nodes from graph
    auto sampled_nodes = std::vector<std::vector<mio::osecirvvs::Model>>();
    std::transform(sampled_graphs.begin(), sampled_graphs.end(), std::back_inserter(sampled_nodes), [](auto&& g) {
        auto models = std::vector<mio::osecirvvs::Model>();
        std::transform(g.nodes().begin(), g.nodes().end(), std::back_inserter(models), [](auto&& n) {
            return n.property;
        });
        return models;
    });

    //compute percentiles
    auto percentile_params = mio::osecirvvs::ensemble_params_percentile(sampled_nodes, 0.6)[0].parameters;

    //spot check parameters
    auto p       = double(percentile_params.get<mio::osecirvvs::ReducTimeInfectedMild>()[mio::AgeGroup(2)]);
    auto samples = std::vector<double>();
    std::transform(sampled_nodes.begin(), sampled_nodes.end(), std::back_inserter(samples),
                   [](const std::vector<mio::osecirvvs::Model>& nodes) {
                       return nodes[0].parameters.get<mio::osecirvvs::ReducTimeInfectedMild>()[mio::AgeGroup(2)];
                   });

    std::nth_element(samples.begin(), samples.begin() + 6, samples.end());
    ASSERT_THAT(p, samples[6]);

    p       = double(percentile_params.get<mio::osecirvvs::TimeExposed>()[mio::AgeGroup(2)]);
    samples = std::vector<double>();
    std::transform(sampled_nodes.begin(), sampled_nodes.end(), std::back_inserter(samples),
                   [](const std::vector<mio::osecirvvs::Model>& nodes) {
                       return nodes[0].parameters.get<mio::osecirvvs::TimeExposed>()[mio::AgeGroup(2)];
                   });
    std::nth_element(samples.begin(), samples.begin() + 6, samples.end());
    ASSERT_THAT(p, samples[6]);
}

TEST(TestOdeSECIRVVS, get_infections_relative)
{
    auto model = make_model(2);
    auto sim   = mio::osecirvvs::Simulation<>(model);
    auto y     = sim.get_result()[0];

    auto relative_infections = get_infections_relative(sim, 0.0, y);

    // see model population init to obtain sum 105=2*(7+7.5+8+9.5+10+10.5)
    ASSERT_DOUBLE_EQ(relative_infections, 105 / model.populations.get_total());
}

TEST(TestOdeSECIRVVS, get_migration_factors)
{
    auto num_age_groups = 2;
    auto model          = make_model(num_age_groups);
    auto sim            = mio::osecirvvs::Simulation<>(model);
    auto y              = sim.get_result()[0];

    auto migration_factors = get_migration_factors(sim, 0.0, y);

    auto expected_values = (Eigen::VectorXd(Eigen::Index(mio::osecirvvs::InfectionState::Count) * num_age_groups) << 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
                               .finished();
    ASSERT_THAT(print_wrap(migration_factors), MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));
}

TEST(TestOdeSECIRVVS, test_commuters)
{
    auto model                                      = make_model(2);
    auto migration_factor                           = 0.1;
    auto non_detection_factor                       = 0.3;
    model.parameters.get_start_commuter_detection() = 0.0;
    model.parameters.get_end_commuter_detection()   = 20.0;
    model.parameters.get_commuter_nondetection()    = non_detection_factor;
    auto sim                                        = mio::osecirvvs::Simulation<>(model);
    auto before_testing                             = sim.get_result().get_last_value().eval();
    auto migrated                                   = (sim.get_result().get_last_value() * migration_factor).eval();
    auto migrated_tested                            = migrated.eval();

    test_commuters(sim, migrated_tested, 0.0);

    ASSERT_NEAR(migrated_tested[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsNaive)],
                migrated[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsNaive)] * non_detection_factor,
                1e-5);
    ASSERT_NEAR(
        sim.get_result().get_last_value()[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed)],
        before_testing[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed)] +
            migrated[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsNaive)] * (1 - non_detection_factor),
        1e-5);
    ASSERT_NEAR(migrated_tested[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity)],
                migrated[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity)] *
                    non_detection_factor,
                1e-5);
    ASSERT_NEAR(
        sim.get_result()
            .get_last_value()[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed)],
        before_testing[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed)] +
            migrated[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity)] *
                (1 - non_detection_factor),
        1e-5);
    ASSERT_NEAR(migrated_tested[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity)],
                migrated[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity)] *
                    non_detection_factor,
                1e-5);
    ASSERT_NEAR(
        sim.get_result()
            .get_last_value()[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed)],
        before_testing[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed)] +
            migrated[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity)] *
                (1 - non_detection_factor),
        1e-5);

    ASSERT_NEAR(migrated_tested[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive)],
                migrated[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive)] * non_detection_factor,
                1e-5);
    ASSERT_NEAR(sim.get_result()
                    .get_last_value()[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed)],
                before_testing[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed)] +
                    migrated[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive)] *
                        (1 - non_detection_factor),
                1e-5);
    ASSERT_NEAR(migrated_tested[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity)],
                migrated[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity)] *
                    non_detection_factor,
                1e-5);
    ASSERT_NEAR(
        sim.get_result()
            .get_last_value()[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed)],
        before_testing[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed)] +
            migrated[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity)] *
                (1 - non_detection_factor),
        1e-5);
    ASSERT_NEAR(migrated_tested[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity)],
                migrated[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity)] *
                    non_detection_factor,
                1e-5);
    ASSERT_NEAR(
        sim.get_result().get_last_value()[Eigen::Index(
            mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed)],
        before_testing[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed)] +
            migrated[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity)] *
                (1 - non_detection_factor),
        1e-5);
}

TEST(TestOdeSECIRVVS, check_constraints_parameters)
{
    auto model = mio::osecirvvs::Model(1);
    ASSERT_EQ(model.parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::osecirvvs::Seasonality>(-0.2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::Seasonality>(0.2);
    model.parameters.set<mio::osecirvvs::ICUCapacity>(-2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ICUCapacity>(2);
    model.parameters.set<mio::osecirvvs::TimeExposed>(-2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeExposed>(2);
    model.parameters.set<mio::osecirvvs::TimeInfectedNoSymptoms>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeInfectedNoSymptoms>(5);
    model.parameters.set<mio::osecirvvs::TimeInfectedSymptoms>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeInfectedSymptoms>(2);
    model.parameters.set<mio::osecirvvs::TimeInfectedSevere>(-1);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeInfectedSevere>(2);
    model.parameters.set<mio::osecirvvs::TimeInfectedCritical>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeInfectedCritical>(10.);
    model.parameters.set<mio::osecirvvs::TimeTemporaryImmunityPI>(0.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeTemporaryImmunityPI>(90.);
    model.parameters.set<mio::osecirvvs::TimeTemporaryImmunityII>(-20.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeTemporaryImmunityII>(90.);
    model.parameters.set<mio::osecirvvs::TimeWaningPartialImmunity>(0.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeWaningPartialImmunity>(100.);
    model.parameters.set<mio::osecirvvs::TimeWaningImprovedImmunity>(0.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeWaningImprovedImmunity>(200);
    model.parameters.set<mio::osecirvvs::TransmissionProbabilityOnContact>(2.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TransmissionProbabilityOnContact>(0.5);
    model.parameters.set<mio::osecirvvs::RelativeTransmissionNoSymptoms>(-1.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::RelativeTransmissionNoSymptoms>(0.5);
    model.parameters.set<mio::osecirvvs::RecoveredPerInfectedNoSymptoms>(3.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::RecoveredPerInfectedNoSymptoms>(0.5);
    model.parameters.set<mio::osecirvvs::RiskOfInfectionFromSymptomatic>(-0.8);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::RiskOfInfectionFromSymptomatic>(0.5);
    model.parameters.set<mio::osecirvvs::SeverePerInfectedSymptoms>(-0.1);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::SeverePerInfectedSymptoms>(0.5);
    model.parameters.set<mio::osecirvvs::CriticalPerSevere>(-1.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::CriticalPerSevere>(0.5);
    model.parameters.set<mio::osecirvvs::DeathsPerCritical>(1.1);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::DeathsPerCritical>(0.5);
    model.parameters.set<mio::osecirvvs::VaccinationGap>(0.2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::VaccinationGap>(2);
    model.parameters.set<mio::osecirvvs::DaysUntilEffectivePartialVaccination>(-2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::DaysUntilEffectivePartialVaccination>(7);
    model.parameters.set<mio::osecirvvs::DaysUntilEffectiveImprovedVaccination>(-0.2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::DaysUntilEffectiveImprovedVaccination>(7);
    model.parameters.set<mio::osecirvvs::DaysUntilEffectiveBoosterImmunity>(-2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::DaysUntilEffectiveBoosterImmunity>(7);
    model.parameters.set<mio::osecirvvs::ReducExposedPartialImmunity>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducExposedPartialImmunity>(0.5);
    model.parameters.set<mio::osecirvvs::ReducExposedImprovedImmunity>(-0.2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducExposedImprovedImmunity>(0.5);
    model.parameters.set<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity>(0.5);
    model.parameters.set<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity>(0.);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity>(0.5);
    model.parameters.set<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity>(-4);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity>(0.5);
    model.parameters.set<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity>(-4);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity>(0.5);
    model.parameters.set<mio::osecirvvs::ReducTimeInfectedMild>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducTimeInfectedMild>(1);
    model.parameters.set<mio::osecirvvs::InfectiousnessNewVariant>(-4);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSECIRVVS, apply_constraints_parameters)
{
    const double tol_times = 1e-1;
    auto model             = mio::osecirvvs::Model(1);
    auto indx_agegroup     = mio::AgeGroup(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    model.parameters.set<mio::osecirvvs::Seasonality>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::Seasonality>(), 0);

    model.parameters.set<mio::osecirvvs::ICUCapacity>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::Seasonality>(), 0);

    model.parameters.set<mio::osecirvvs::TimeExposed>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeExposed>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TimeInfectedNoSymptoms>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedNoSymptoms>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TimeInfectedSymptoms>(1e-10);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedSymptoms>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TimeInfectedSevere>(-1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedSevere>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TimeInfectedCritical>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedCritical>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TimeTemporaryImmunityPI>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeTemporaryImmunityPI>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TimeTemporaryImmunityII>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeTemporaryImmunityII>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TimeWaningPartialImmunity>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeWaningPartialImmunity>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TimeWaningImprovedImmunity>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeWaningImprovedImmunity>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TransmissionProbabilityOnContact>(2.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR(model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>()[indx_agegroup], 0.0, 1e-14);

    model.parameters.set<mio::osecirvvs::RelativeTransmissionNoSymptoms>(-1.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::RelativeTransmissionNoSymptoms>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::RecoveredPerInfectedNoSymptoms>(3.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::RiskOfInfectionFromSymptomatic>(-0.8);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::SeverePerInfectedSymptoms>(-0.1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::SeverePerInfectedSymptoms>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::CriticalPerSevere>(-1.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::CriticalPerSevere>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::DeathsPerCritical>(1.1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::DeathsPerCritical>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::VaccinationGap>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::VaccinationGap>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::DaysUntilEffectivePartialVaccination>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::DaysUntilEffectivePartialVaccination>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::DaysUntilEffectiveImprovedVaccination>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::DaysUntilEffectiveImprovedVaccination>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::DaysUntilEffectiveBoosterImmunity>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::DaysUntilEffectiveBoosterImmunity>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::ReducExposedPartialImmunity>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::ReducExposedPartialImmunity>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirvvs::ReducExposedImprovedImmunity>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::ReducExposedImprovedImmunity>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity>(0.);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity>(-4);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity>(-4);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity>()[indx_agegroup],
              1);

    model.parameters.set<mio::osecirvvs::ReducTimeInfectedMild>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::ReducTimeInfectedMild>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirvvs::InfectiousnessNewVariant>(-4);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::InfectiousnessNewVariant>()[indx_agegroup], 1);

    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSECIRVVS, apply_variant_function)
{
    auto model = mio::osecirvvs::Model(1);
    model.parameters.set<mio::osecirvvs::TransmissionProbabilityOnContact>(0.2);

    model.parameters.set<mio::osecirvvs::StartDay>(0);
    model.parameters.set<mio::osecirvvs::StartDayNewVariant>(10);
    model.parameters.set<mio::osecirvvs::InfectiousnessNewVariant>(2.0);
    auto sim = mio::osecirvvs::Simulation<>(model);

    // test that the transmission probability is not changed due to calling the advance function
    sim.advance(0.01);
    EXPECT_NEAR(sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>()[mio::AgeGroup(0)],
                0.2, 1e-10);

    // test if the transmission probability is set to the correct value after applying the variant.
    // The referece values are calculated using equation (36) in doi.org/10.1371/journal.pcbi.1010054
    auto base_infectiousness = sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>();

    // however, the parameter should stay unchanged if the new variant is not present in the population.
    sim.apply_variant(0, base_infectiousness);
    EXPECT_NEAR(sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>()[mio::AgeGroup(0)],
                0.2, 1e-10);

    sim.apply_variant(9, base_infectiousness);
    EXPECT_NEAR(sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>()[mio::AgeGroup(0)],
                0.2, 1e-10);

    sim.apply_variant(10, base_infectiousness);
    EXPECT_NEAR(sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>()[mio::AgeGroup(0)],
                0.99 * base_infectiousness[mio::AgeGroup(0)] +
                    0.01 * base_infectiousness[mio::AgeGroup(0)] *
                        sim.get_model().parameters.get<mio::osecirvvs::InfectiousnessNewVariant>()[mio::AgeGroup(0)],
                1e-10);

    sim.apply_variant(45, base_infectiousness);
    EXPECT_NEAR(sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>()[mio::AgeGroup(0)],
                0.68 * base_infectiousness[mio::AgeGroup(0)] +
                    0.32 * base_infectiousness[mio::AgeGroup(0)] *
                        sim.get_model().parameters.get<mio::osecirvvs::InfectiousnessNewVariant>()[mio::AgeGroup(0)],
                1e-10);

    sim.apply_variant(1000, base_infectiousness);
    EXPECT_NEAR(sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>()[mio::AgeGroup(0)],
                0.4, 1e-10);
}