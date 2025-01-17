/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#include "ode_secirts/infection_state.h"
#include "ode_secirts/model.h"
#include "ode_secirts/parameter_space.h"
#include "ode_secirts/parameters.h"
#include "ode_secirts/parameters_io.h"
#include "ode_secirts/analyze_result.h"

#include "gtest/gtest.h"
#include "gmock/gmock-matchers.h"
#include <algorithm>
#include <iterator>
#include <limits>

const mio::osecirts::Model<double>& osecirts_testing_model()
{
    static mio::osecirts::Model<double> model(1);
    model.populations.array().setConstant(1);
    auto nb_groups = model.parameters.get_num_groups();

    for (mio::AgeGroup i = 0; i < nb_groups; i++) {

        // parameters
        //times
        model.parameters.get<mio::osecirts::TimeExposed<double>>()[i]                = 1.0;
        model.parameters.get<mio::osecirts::TimeInfectedNoSymptoms<double>>()[i]     = 1.0;
        model.parameters.get<mio::osecirts::TimeInfectedSymptoms<double>>()[i]       = 1.0;
        model.parameters.get<mio::osecirts::TimeInfectedSevere<double>>()[i]         = 1.0;
        model.parameters.get<mio::osecirts::TimeInfectedCritical<double>>()[i]       = 1.0;
        model.parameters.get<mio::osecirts::TimeTemporaryImmunityPI<double>>()[i]    = 1.0;
        model.parameters.get<mio::osecirts::TimeTemporaryImmunityII<double>>()[i]    = 1.0;
        model.parameters.get<mio::osecirts::TimeWaningPartialImmunity<double>>()[i]  = 1.0;
        model.parameters.get<mio::osecirts::TimeWaningImprovedImmunity<double>>()[i] = 1.0;

        //probabilities
        model.parameters.get<mio::osecirts::TransmissionProbabilityOnContact<double>>()[i]  = 0.5;
        model.parameters.get<mio::osecirts::RelativeTransmissionNoSymptoms<double>>()[i]    = 0.5;
        model.parameters.get<mio::osecirts::RiskOfInfectionFromSymptomatic<double>>()[i]    = 0.5;
        model.parameters.get<mio::osecirts::MaxRiskOfInfectionFromSymptomatic<double>>()[i] = 0.5;
        model.parameters.get<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>()[i]    = 0.5;
        model.parameters.get<mio::osecirts::SeverePerInfectedSymptoms<double>>()[i]         = 0.5;
        model.parameters.get<mio::osecirts::CriticalPerSevere<double>>()[i]                 = 0.5;
        model.parameters.get<mio::osecirts::DeathsPerCritical<double>>()[i]                 = 0.5;

        model.parameters.get<mio::osecirts::ReducExposedPartialImmunity<double>>()[i]                     = 0.5;
        model.parameters.get<mio::osecirts::ReducExposedImprovedImmunity<double>>()[i]                    = 0.5;
        model.parameters.get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<double>>()[i]            = 0.5;
        model.parameters.get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<double>>()[i]           = 0.5;
        model.parameters.get<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[i]  = 0.5;
        model.parameters.get<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[i] = 0.5;
        model.parameters.get<mio::osecirts::ReducTimeInfectedMild<double>>()[i]                           = 1.0;
    }

    model.parameters.get<mio::osecirts::ICUCapacity<double>>()          = 10;
    model.parameters.get<mio::osecirts::TestAndTraceCapacity<double>>() = 1.0;
    const size_t daily_vaccinations                                     = 1;
    const size_t num_days                                               = 5;
    model.parameters.get<mio::osecirts::DailyPartialVaccinations<double>>().resize(mio::SimulationDay(num_days));
    model.parameters.get<mio::osecirts::DailyFullVaccinations<double>>().resize(mio::SimulationDay(num_days));
    model.parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>().resize(mio::SimulationDay(num_days));
    for (size_t i = 0; i < num_days; ++i) {
        for (mio::AgeGroup j = 0; j < nb_groups; ++j) {
            auto num_vaccinations = static_cast<double>(i * daily_vaccinations);
            model.parameters.get<mio::osecirts::DailyPartialVaccinations<double>>()[{j, mio::SimulationDay(i)}] =
                num_vaccinations;
            model.parameters.get<mio::osecirts::DailyFullVaccinations<double>>()[{j, mio::SimulationDay(i)}] =
                num_vaccinations;
            model.parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>()[{j, mio::SimulationDay(i)}] =
                num_vaccinations;
        }
    }

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecirts::ContactPatterns<double>>();
    const double cont_freq                  = 1;
    const double fact                       = 1.0 / (double)(size_t)nb_groups;
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    return model;
}

TEST(TestOdeSECIRTS, get_flows)
{
    auto& model = osecirts_testing_model();

    Eigen::VectorXd y     = Eigen::VectorXd::Constant(model.get_initial_values().size(), 1.0);
    Eigen::VectorXd flows = Eigen::VectorXd::Zero(model.get_initial_flows().size());

    model.get_flows(y, y, 0.0, flows);

    std::vector<double> expected_values = {0.09375,  0.90625,  1,   0.5, 0.5, 0.5, 0.5,      0.5,      0.5, 0.5, 0.5,
                                           0.5,      0.5,      0,   0.5, 0.5, 1,   0.046875, 0.953125, 1,   0.5, 0.5,
                                           0.5,      0.5,      0.5, 0.5, 0.5, 0.5, 0.5,      0.5,      0,   0.5, 0.5,
                                           0.046875, 0.953125, 1,   0.5, 0.5, 0.5, 0.5,      0.5,      0.5, 0.5, 0.5,
                                           0.5,      0.5,      0,   0.5, 0.5, 0,   1,        1,        1};

    EXPECT_EQ(flows.size(), 53);

    // Compare expected with actual values
    for (size_t i = 0; i < expected_values.size(); i++) {
        EXPECT_NEAR(flows(i), expected_values[i], 1e-10);
    }
}

TEST(TestOdeSECIRTS, Simulation)
{
    auto integrator = std::make_shared<mio::EulerIntegratorCore<double>>();
    auto sim        = mio::osecirts::Simulation<double>(osecirts_testing_model(), 0, 1);
    sim.set_integrator(integrator);
    sim.advance(1);

    EXPECT_EQ(sim.get_result().get_num_time_points(), 2); // stores initial value and single step

    Eigen::VectorXd expected_result(29);
    expected_result << 1.0, 2.0, 0.09375, 0.046875, 0.046875, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5,
        0.5, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.0, 1.5, 1.5, 1.5, 2.90625, 7.90625;
    EXPECT_EQ(sim.get_result().get_last_value(), expected_result);
}

using FlowSim = mio::osecirts::Simulation<double, mio::FlowSimulation<double, mio::osecirts::Model<double>>>;

TEST(TestOdeSECIRTS, FlowSimulation)
{
    auto integrator = std::make_shared<mio::EulerIntegratorCore<double>>();
    FlowSim sim(osecirts_testing_model(), 0, 1);

    sim.set_integrator(integrator);
    sim.advance(1);

    EXPECT_EQ(sim.get_result().get_num_time_points(), 2); // stores initial value and single step

    Eigen::VectorXd expected_result(29);
    expected_result << 1.0, 2.0, 0.09375, 0.046875, 0.046875, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5,
        0.5, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.0, 1.5, 1.5, 1.5, 2.90625, 7.90625;
    EXPECT_EQ(sim.get_result().get_last_value(), expected_result);

    Eigen::VectorXd expected_flows(53);
    expected_flows << 0.09375, 0.90625, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 1.0,
        0.046875, 0.953125, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.046875, 0.953125,
        1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 1.0, 1.0, 1.0;
    EXPECT_EQ(sim.get_flows().get_last_value(), expected_flows);
}

TEST(TestOdeSECIRTS, simulateDefault)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.1;

    mio::osecirts::Model<double> model(1);
    model.populations.set_total(10);
    model.populations.set_difference_from_total({(mio::AgeGroup)0, mio::osecirts::InfectionState::SusceptibleNaive},
                                                10);
    model.parameters.get<mio::osecirts::DailyPartialVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirts::DailyPartialVaccinations<double>>().array().setConstant(0);
    model.parameters.get<mio::osecirts::DailyFullVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirts::DailyFullVaccinations<double>>().array().setConstant(0);
    model.parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>().array().setConstant(0);
    auto result = mio::simulate<double, mio::osecirts::Model<double>>(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

// Test model initialization with total population of 0 and ensure get_flows returns no NaN values
TEST(TestOdeSECIRTS, population_zero_no_nan)
{
    // initialize simple model with total population 0
    mio::osecirts::Model<double> model(1);
    model.populations.set_total(0.0);

    model.parameters.get<mio::osecirts::DailyPartialVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirts::DailyPartialVaccinations<double>>().array().setConstant(0);
    model.parameters.get<mio::osecirts::DailyFullVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirts::DailyFullVaccinations<double>>().array().setConstant(0);
    model.parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>().array().setConstant(0);

    // call the get_flows function
    auto dydt_default = Eigen::VectorXd(model.get_initial_flows().size());
    dydt_default.setZero();
    auto y0 = model.get_initial_values();
    model.get_flows(y0, y0, 0, dydt_default);

    // check that there are now NaN values in dydt_default
    for (int i = 0; i < dydt_default.size(); i++) {
        EXPECT_FALSE(std::isnan(dydt_default[i]));
    }
}

TEST(TestOdeSECIRTS, overflow_vaccinations)
{
    const double t0   = 0;
    const double tmax = 1;
    const double dt   = 1;

    // init simple model
    mio::osecirts::Model<double> model(1);
    model.populations[{mio::AgeGroup(0), mio::osecirts::InfectionState::SusceptibleNaive}]            = 10.;
    model.populations[{mio::AgeGroup(0), mio::osecirts::InfectionState::SusceptiblePartialImmunity}]  = 10.;
    model.populations[{mio::AgeGroup(0), mio::osecirts::InfectionState::SusceptibleImprovedImmunity}] = 10.;
    model.parameters.get<mio::osecirts::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)] = 0.0;

    // set vaccination rates higher than total population for each layer
    const size_t daily_vaccinations = 100;
    model.parameters.get<mio::osecirts::DailyPartialVaccinations<double>>().resize(
        mio::SimulationDay(static_cast<size_t>(tmax + 1)));
    model.parameters.get<mio::osecirts::DailyFullVaccinations<double>>().resize(
        mio::SimulationDay(static_cast<size_t>(tmax + 1)));
    model.parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>().resize(
        mio::SimulationDay(static_cast<size_t>(tmax + 1)));
    for (size_t i = 0; i <= tmax; ++i) {
        auto num_vaccinations = static_cast<double>((i + 1) * daily_vaccinations);
        model.parameters
            .get<mio::osecirts::DailyPartialVaccinations<double>>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
        model.parameters
            .get<mio::osecirts::DailyFullVaccinations<double>>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
        model.parameters
            .get<mio::osecirts::DailyBoosterVaccinations<double>>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
    }

    // simulate one step with explicit Euler
    auto integrator = std::make_shared<mio::EulerIntegratorCore<double>>();
    auto result     = mio::simulate_flows<double, mio::osecirts::Model<double>>(t0, tmax, dt, model, integrator);

    // get the flow indices for each type of vaccination and also the indices of the susceptible compartments
    auto flow_indx_partial_vaccination =
        model.get_flat_flow_index<mio::osecirts::InfectionState::SusceptibleNaive,
                                  mio::osecirts::InfectionState::TemporaryImmunePartialImmunity>({mio::AgeGroup(0)});
    auto flow_indx_full_vaccination =
        model.get_flat_flow_index<mio::osecirts::InfectionState::SusceptiblePartialImmunity,
                                  mio::osecirts::InfectionState::TemporaryImmuneImprovedImmunity>({mio::AgeGroup(0)});
    auto flow_indx_booster_vaccination =
        model.get_flat_flow_index<mio::osecirts::InfectionState::SusceptibleImprovedImmunity,
                                  mio::osecirts::InfectionState::TemporaryImmuneImprovedImmunity>({mio::AgeGroup(0)});
    auto indx_S_naive =
        model.populations.get_flat_index({mio::AgeGroup(0), mio::osecirts::InfectionState::SusceptibleNaive});
    auto indx_S_partial =
        model.populations.get_flat_index({mio::AgeGroup(0), mio::osecirts::InfectionState::SusceptiblePartialImmunity});
    auto indx_S_improved = model.populations.get_flat_index(
        {mio::AgeGroup(0), mio::osecirts::InfectionState::SusceptibleImprovedImmunity});

    // check that the number of vaccinated people is never higher than the total number of susceptible people
    EXPECT_NEAR(result[1].get_last_value()[flow_indx_partial_vaccination], result[0].get_value(0)[indx_S_naive], 1e-10);
    EXPECT_NEAR(result[1].get_last_value()[flow_indx_full_vaccination], result[0].get_value(0)[indx_S_partial], 1e-10);
    EXPECT_NEAR(result[1].get_last_value()[flow_indx_booster_vaccination], result[0].get_value(0)[indx_S_improved],
                1e-10);
}

TEST(TestOdeSECIRTS, smooth_vaccination_rate)
{
    const ScalarType tmax = 2.;

    // init simple model
    mio::osecirts::Model model(1);
    auto& daily_vaccinations = model.parameters.get<mio::osecirts::DailyPartialVaccinations<double>>();
    daily_vaccinations.resize(mio::SimulationDay(static_cast<size_t>(tmax + 1)));

    daily_vaccinations[{mio::AgeGroup(0), mio::SimulationDay(0)}] = 0;
    daily_vaccinations[{mio::AgeGroup(0), mio::SimulationDay(1)}] = 10;
    daily_vaccinations[{mio::AgeGroup(0), mio::SimulationDay(2)}] = 110;

    const auto eps1 = 0.15;

    // test when t is out of bounds
    Eigen::VectorXd result = model.vaccinations_at(5.5, daily_vaccinations, eps1);
    EXPECT_EQ(result.size(), 1);
    EXPECT_NEAR(result[0], 0, 1e-12);

    // test when ub + 1 is out of bounds
    result = model.vaccinations_at(1.85, daily_vaccinations, eps1);
    EXPECT_EQ(result.size(), 1);
    EXPECT_NEAR(result[0], 100, 1e-12);

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

namespace
{
void assign_uniform_distribution(mio::UncertainValue<double>& p, double min, double max, bool set_invalid_initial_value)
{
    auto invalid_initial = max == 0 ? 1.0 : max * 1.001;
    auto valid_initial   = (max + min) * 0.5;
    auto initial         = set_invalid_initial_value ? invalid_initial : valid_initial;
    p                    = mio::UncertainValue(initial);
    p.set_distribution(mio::ParameterDistributionUniform(min, max));
}

template <size_t N>
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue<double>, mio::AgeGroup>& array,
                                       const double (&min)[N], const double (&max)[N], bool set_invalid_initial_value)
{
    for (auto i = mio::AgeGroup(0); i < array.size<mio::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min[size_t(i)], max[size_t(i)], set_invalid_initial_value);
    }
}

void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue<double>, mio::AgeGroup>& array,
                                       double min, double max, bool set_invalid_initial_value)
{
    for (auto i = mio::AgeGroup(0); i < array.size<mio::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min, max, set_invalid_initial_value);
    }
}
} // namespace

void set_synthetic_population_data(mio::osecirts::Model<double>::Populations& populations,
                                   bool set_invalid_initial_value)
{
    for (mio::AgeGroup i = 0; i < mio::get<mio::AgeGroup>(populations.size()); i++) {
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::ExposedNaive}], 10, 20,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::ExposedImprovedImmunity}], 10, 21,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::ExposedPartialImmunity}], 10, 22,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsNaive}], 1, 10,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunity}],
                                    2, 10, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunity}],
                                    3, 10, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedSymptomsNaive}], 4, 10,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity}], 5,
                                    10, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunity}],
                                    6, 10, set_invalid_initial_value);

        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsNaiveConfirmed}],
                                    5, 11, set_invalid_initial_value);
        assign_uniform_distribution(
            populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}], 5, 12,
            set_invalid_initial_value);
        assign_uniform_distribution(
            populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}], 5, 13,
            set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedSymptomsNaiveConfirmed}], 5,
                                    14, set_invalid_initial_value);
        assign_uniform_distribution(
            populations[{i, mio::osecirts::InfectionState::InfectedSymptomsPartialImmunityConfirmed}], 5, 15,
            set_invalid_initial_value);
        assign_uniform_distribution(
            populations[{i, mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}], 5, 16,
            set_invalid_initial_value);

        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedSevereNaive}], 1, 2,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedSevereImprovedImmunity}], 1,
                                    3, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedSeverePartialImmunity}], 1,
                                    4, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedCriticalNaive}], 1, 5,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedCriticalPartialImmunity}], 1,
                                    6, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::InfectedCriticalImprovedImmunity}],
                                    1, 7, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::SusceptibleImprovedImmunity}], 200,
                                    300, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirts::InfectionState::SusceptiblePartialImmunity}], 200,
                                    400, set_invalid_initial_value);
        populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecirts::InfectionState::SusceptibleNaive},
                                                                   1000);
    }
}

void set_demographic_parameters(mio::osecirts::Model<double>::ParameterSet& parameters, bool set_invalid_initial_value)
{
    assign_uniform_distribution(parameters.get<mio::osecirts::ICUCapacity<double>>(), 20, 50,
                                set_invalid_initial_value);
    assign_uniform_distribution(parameters.get<mio::osecirts::TestAndTraceCapacity<double>>(), 100, 200,
                                set_invalid_initial_value);
    parameters.get<mio::osecirts::DailyPartialVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    parameters.get<mio::osecirts::DailyPartialVaccinations<double>>().array().setConstant(5);
    parameters.get<mio::osecirts::DailyFullVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    parameters.get<mio::osecirts::DailyFullVaccinations<double>>().array().setConstant(3);
    parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>().array().setConstant(3);
}

void set_contact_parameters(mio::osecirts::Model<double>::ParameterSet& parameters, bool set_invalid_initial_value)
{
    auto& contacts       = parameters.get<mio::osecirts::ContactPatterns<double>>();
    auto& contact_matrix = contacts.get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(0.5);
    contact_matrix[0].get_baseline().diagonal().setConstant(5.0);
    contact_matrix[0].add_damping(0.3, mio::SimulationTime(5.0));

    auto& npis      = parameters.get<mio::osecirts::DynamicNPIsInfectedSymptoms<double>>();
    auto npi_groups = Eigen::VectorXd::Ones(contact_matrix[0].get_num_groups());
    auto npi_value  = mio::UncertainValue(0.5);
    assign_uniform_distribution(npi_value, 0.25, 0.75, set_invalid_initial_value);
    npis.set_threshold(10.0, {mio::DampingSampling<double>(npi_value, mio::DampingLevel(0), mio::DampingType(0),
                                                           mio::SimulationTime(0), {0}, npi_groups)});
    npis.set_base_value(100'000);
    npis.set_interval(mio::SimulationTime(3.0));
    npis.set_duration(mio::SimulationTime(14.0));
    parameters.get_end_dynamic_npis() = 10.0; //required for dynamic NPIs to have effect in this model
}

void set_covid_parameters(mio::osecirts::Model<double>::ParameterSet& params, bool set_invalid_initial_value)
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

    array_assign_uniform_distribution(params.get<mio::osecirts::TimeExposed<double>>(), timeExposedMin, timeExposedMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::TimeInfectedNoSymptoms<double>>(),
                                      timeInfectedNoSymptomsMin, timeInfectedNoSymptomsMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::TimeInfectedSymptoms<double>>(),
                                      timeInfectedSymptomsMin, timeInfectedSymptomsMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::TimeInfectedSevere<double>>(), timeInfectedSevereMin,
                                      timeInfectedSevereMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::TimeInfectedCritical<double>>(),
                                      timeInfectedCriticalMin, timeInfectedCriticalMax, set_invalid_initial_value);

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

    array_assign_uniform_distribution(params.get<mio::osecirts::TransmissionProbabilityOnContact<double>>(),
                                      transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::RelativeTransmissionNoSymptoms<double>>(),
                                      relativeTransmissionNoSymptomsMin, relativeTransmissionNoSymptomsMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::RiskOfInfectionFromSymptomatic<double>>(),
                                      riskOfInfectionFromSymptomaticMin, riskOfInfectionFromSymptomaticMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::MaxRiskOfInfectionFromSymptomatic<double>>(),
                                      maxRiskOfInfectionFromSymptomaticMin, maxRiskOfInfectionFromSymptomaticMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>(),
                                      recoveredPerInfectedNoSymptomsMin, recoveredPerInfectedNoSymptomsMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::SeverePerInfectedSymptoms<double>>(),
                                      severePerInfectedSymptomsMin, severePerInfectedSymptomsMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::CriticalPerSevere<double>>(), criticalPerSevereMin,
                                      criticalPerSevereMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::DeathsPerCritical<double>>(), deathsPerCriticalMin,
                                      deathsPerCriticalMax, set_invalid_initial_value);

    array_assign_uniform_distribution(params.get<mio::osecirts::ReducExposedPartialImmunity<double>>(),
                                      reducExposedPartialImmunityMin, reducExposedPartialImmunityMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::ReducExposedImprovedImmunity<double>>(),
                                      reducExposedImprovedImmunityMin, reducExposedImprovedImmunityMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<double>>(),
                                      reducInfectedSymptomsPartialImmunityMin, reducInfectedSymptomsPartialImmunityMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<double>>(),
                                      reducInfectedSymptomsImprovedImmunityMin,
                                      reducInfectedSymptomsImprovedImmunityMax, set_invalid_initial_value);
    array_assign_uniform_distribution(
        params.get<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(),
        reducInfectedSevereCriticalDeadPartialImmunityMin, reducInfectedSevereCriticalDeadPartialImmunityMax,
        set_invalid_initial_value);
    array_assign_uniform_distribution(
        params.get<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(),
        reducInfectedSevereCriticalDeadImprovedImmunityMin, reducInfectedSevereCriticalDeadImprovedImmunityMax,
        set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirts::ReducTimeInfectedMild<double>>(),
                                      reducTimeInfectedMildMin, reducTimeInfectedMildMax, set_invalid_initial_value);

    //sasonality
    const double seasonality_min = 0.1;
    const double seasonality_max = 0.3;

    assign_uniform_distribution(params.get<mio::osecirts::Seasonality<double>>(), seasonality_min, seasonality_max,
                                set_invalid_initial_value);
}

namespace
{
mio::osecirts::Model<double> make_model(int num_age_groups, bool set_invalid_initial_value = false)
{
    assert(num_age_groups <= 6 && "Provide more values in functions above to test more age groups.");
    mio::osecirts::Model model(num_age_groups);
    set_covid_parameters(model.parameters, set_invalid_initial_value);
    set_synthetic_population_data(model.populations, set_invalid_initial_value);
    set_demographic_parameters(model.parameters, set_invalid_initial_value);
    set_contact_parameters(model.parameters, set_invalid_initial_value);
    model.parameters.apply_constraints();
    return model;
}
} // namespace

TEST(TestOdeSECIRTS, draw_sample_graph)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    mio::Graph<mio::osecirts::Model<double>, mio::MobilityParameters<double>> graph;

    auto num_age_groups = 6;
    //create model with invalid initials so the test fails if no sampling is done
    graph.add_node(0, make_model(num_age_groups, /*set_invalid_initial_value*/ true));
    graph.add_node(0, make_model(num_age_groups, /*set_invalid_initial_value*/ true));
    graph.add_edge(0, 1, Eigen::VectorXd::Constant(num_age_groups, num_age_groups));

    auto sampled_graph = mio::osecirts::draw_sample(graph);

    ASSERT_EQ(sampled_graph.nodes().size(), graph.nodes().size());
    ASSERT_EQ(sampled_graph.edges().size(), graph.edges().size());

    // spot check for sampling
    auto& parameters0          = sampled_graph.nodes()[0].property.parameters;
    auto& populations0         = sampled_graph.nodes()[0].property.populations;
    auto& timeInfectedCritical = parameters0.get<mio::osecirts::TimeInfectedCritical<double>>()[mio::AgeGroup(1)];
    ASSERT_GE(double(timeInfectedCritical), 4.95);
    ASSERT_LE(double(timeInfectedCritical), 8.95);
    auto& param_exp_factor = parameters0.get<mio::osecirts::ReducExposedPartialImmunity<double>>()[mio::AgeGroup(0)];
    ASSERT_GE(double(param_exp_factor), 0.75);
    ASSERT_LE(double(param_exp_factor), 0.85);
    auto& compartment_inf =
        populations0[{mio::AgeGroup(2), mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity}];
    ASSERT_GE(double(compartment_inf), 5.0);
    ASSERT_LE(double(compartment_inf), 10.0);
    auto& npi_value =
        parameters0.get<mio::osecirts::DynamicNPIsInfectedSymptoms<double>>().get_thresholds()[0].second[0].get_value();
    ASSERT_GE(double(npi_value), 0.25);
    ASSERT_LE(double(npi_value), 0.75);

    // special cases
    ASSERT_NEAR(populations0.get_total(), 1000 * num_age_groups, 1e-2);
    ASSERT_TRUE(
        (parameters0.get<mio::osecirts::InfectiousnessNewVariant<double>>().array(),
         parameters0.get<mio::osecirts::TransmissionProbabilityOnContact<double>>().array() * 1.0) //using high variant
            .all());

    // spot check for parameters that should be equal or different between nodes
    auto& parameters1  = sampled_graph.nodes()[1].property.parameters;
    auto& populations1 = sampled_graph.nodes()[1].property.populations;
    ASSERT_EQ(
        parameters1.get<mio::osecirts::DynamicNPIsInfectedSymptoms<double>>().get_thresholds()[0].second[0].get_value(),
        parameters0.get<mio::osecirts::DynamicNPIsInfectedSymptoms<double>>()
            .get_thresholds()[0]
            .second[0]
            .get_value());
    ASSERT_TRUE((parameters1.get<mio::osecirts::TimeInfectedSymptoms<double>>().array() ==
                 parameters0.get<mio::osecirts::TimeInfectedSymptoms<double>>().array())
                    .all());
    //these could fail in very(!) rare cases if they are randomly sampled to the same value
    ASSERT_NE(parameters1.get<mio::osecirts::ICUCapacity<double>>(),
              parameters0.get<mio::osecirts::ICUCapacity<double>>())
        << "Failure might be spurious, check RNG seeds.";
    ASSERT_FALSE((populations1.array() == populations0.array()).all()) << "Failure might be spurious, check RNG seeds.";
}

TEST(TestOdeSECIRTS, draw_sample_model)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    auto num_age_groups = 6;
    auto model          = make_model(num_age_groups, /*set_invalid_initial_value*/ true);
    mio::osecirts::draw_sample(model);

    // spot check for sampling
    auto& parameters           = model.parameters;
    auto& populations          = model.populations;
    auto& timeInfectedCritical = parameters.get<mio::osecirts::TimeInfectedCritical<double>>()[mio::AgeGroup(1)];
    ASSERT_GE(double(timeInfectedCritical), 4.95);
    ASSERT_LE(double(timeInfectedCritical), 8.95);
    auto& param_exp_factor = parameters.get<mio::osecirts::ReducExposedPartialImmunity<double>>()[mio::AgeGroup(0)];
    ASSERT_GE(double(param_exp_factor), 0.75);
    ASSERT_LE(double(param_exp_factor), 0.85);
    auto& compartment_inf =
        populations[{mio::AgeGroup(2), mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity}];
    ASSERT_GE(double(compartment_inf), 5.0);
    ASSERT_LE(double(compartment_inf), 10.0);

    // special cases
    ASSERT_NEAR(populations.get_total(), 1000 * num_age_groups, 1e-2);
    ASSERT_TRUE((parameters.get<mio::osecirts::InfectiousnessNewVariant<double>>().array(),
                 parameters.get<mio::osecirts::TransmissionProbabilityOnContact<double>>().array() * 1.0)
                    .all());
}

TEST(TestOdeSECIRTS, checkPopulationConservation)
{
    auto num_age_groups = 2;
    auto model          = make_model(num_age_groups);
    auto num_days       = 30;

    auto result = mio::osecirts::simulate<double>(0, num_days, 0.1, model);

    double num_persons = 0.0;
    for (auto i = 0; i < result.get_last_value().size(); i++) {
        EXPECT_GE(result.get_last_value()[i], -1e-3);
        num_persons += result.get_last_value()[i];
    }
    EXPECT_NEAR(num_persons, model.populations.get_total(), 1e-10);
}

#if defined(MEMILIO_HAS_HDF5) && defined(MEMILIO_HAS_JSONCPP)

TEST(TestOdeSECIRTS, read_confirmed_cases)
{
    auto num_age_groups = 6; //reading data requires RKI data age groups
    auto model          = std::vector<mio::osecirts::Model<double>>({make_model(num_age_groups)});
    const std::vector<int> region{1002};
    auto path = mio::path_join(TEST_DATA_DIR, "pydata/Germany/cases_all_county_age_ma7.json");

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

    ASSERT_THAT(mio::osecirts::details::read_confirmed_cases_data(
                    path, region, {2020, 12, 01}, num_Exposed, num_InfectedNoSymptoms, num_InfectedSymptoms,
                    num_InfectedSevere, num_icu, num_death, num_timm, model,
                    std::vector<double>(size_t(num_age_groups), 1.0), 0),
                IsSuccess());

    // read again with invalid date
    ASSERT_THAT(mio::osecirts::details::read_confirmed_cases_data(
                    path, region, {3020, 12, 01}, num_Exposed, num_InfectedNoSymptoms, num_InfectedSymptoms,
                    num_InfectedSevere, num_icu, num_death, num_timm, model,
                    std::vector<double>(size_t(num_age_groups), 1.0), 0),
                IsFailure(mio::StatusCode::OutOfRange));

    // call the compute function with empty case data
    const std::vector<mio::ConfirmedCasesDataEntry> empty_case_data;
    ASSERT_THAT(mio::osecirts::details::compute_confirmed_cases_data(
                    empty_case_data, num_Exposed, num_InfectedNoSymptoms, num_InfectedSymptoms, num_InfectedSevere,
                    num_icu, num_death, num_timm, region, {2020, 12, 01}, model,
                    std::vector<double>(size_t(num_age_groups), 1.0), 0),
                IsFailure(mio::StatusCode::InvalidValue));
}

TEST(TestOdeSECIRTS, set_divi_data_invalid_dates)
{
    mio::set_log_level(mio::LogLevel::off);
    auto model = mio::osecirts::Model<double>(1);
    model.populations.array().setConstant(1);
    auto model_vector = std::vector<mio::osecirts::Model<double>>{model};

    // Test with date before DIVI dataset was available.
    EXPECT_THAT(mio::osecirts::details::set_divi_data(model_vector, "", {1001}, {2019, 12, 01}, 1.0), IsSuccess());
    // Assure that populations is the same as before.
    EXPECT_THAT(print_wrap(model_vector[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(model.populations.array().cast<double>()), 1e-10, 1e-10));

    // Test with data after DIVI dataset was no longer updated.
    EXPECT_THAT(mio::osecirts::details::set_divi_data(model_vector, "", {1001}, {2025, 12, 01}, 1.0), IsSuccess());
    EXPECT_THAT(print_wrap(model_vector[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(model.populations.array().cast<double>()), 1e-10, 1e-10));

    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSECIRTS, set_confirmed_cases_data_with_ICU)
{
    auto num_age_groups = 6;
    auto model          = mio::osecirts::Model<double>(num_age_groups);
    model.populations.array().setConstant(1);

    const std::vector<std::vector<double>> immunity_population = {{0.04, 0.04, 0.075, 0.08, 0.035, 0.01},
                                                                  {0.61, 0.61, 0.62, 0.62, 0.58, 0.41},
                                                                  {0.35, 0.35, 0.305, 0.3, 0.385, 0.58}};

    // read case data
    auto case_data =
        mio::read_confirmed_cases_data(mio::path_join(TEST_DATA_DIR, "cases_all_county_age_ma7.json")).value();

    // Change dates of the case data so that no ICU data is available at that time.
    const auto t0 = mio::Date(2025, 1, 1);
    auto day_add  = 0;
    for (auto& entry : case_data) {
        entry.date = offset_date_by_days(t0, day_add);
        day_add++;
    }

    // ICU occupancy before function is called
    auto ICU_before = std::vector<double>(size_t(num_age_groups) * 3, 0.0);
    for (auto age_group = mio::AgeGroup(0); age_group < (mio::AgeGroup)num_age_groups; age_group++) {
        ICU_before[(size_t)age_group] =
            model.populations[{age_group, mio::osecirts::InfectionState::InfectedCriticalNaive}].value();
        ICU_before[(size_t)age_group + 1] =
            model.populations[{age_group, mio::osecirts::InfectionState::InfectedCriticalPartialImmunity}].value();
        ICU_before[(size_t)age_group + 2] =
            model.populations[{age_group, mio::osecirts::InfectionState::InfectedCriticalImprovedImmunity}].value();
    }

    // get day in mid of the data
    auto mid_day = case_data[(size_t)case_data.size() / 2].date;

    auto model_vector       = std::vector<mio::osecirts::Model<double>>{model};
    auto scaling_factor_inf = std::vector<double>(size_t(model.parameters.get_num_groups()), 1.0);
    EXPECT_THAT(mio::osecirts::details::set_confirmed_cases_data(model_vector, case_data, {1002}, mid_day,
                                                                 scaling_factor_inf, immunity_population),
                IsSuccess());

    // Get new setted ICU compartment
    auto ICU_after = std::vector<double>(size_t(num_age_groups) * 3, 0.0);
    for (auto age_group = mio::AgeGroup(0); age_group < (mio::AgeGroup)num_age_groups; age_group++) {
        ICU_after[(size_t)age_group] =
            model_vector[0].populations[{age_group, mio::osecirts::InfectionState::InfectedCriticalNaive}].value();
        ICU_after[(size_t)age_group + 1] =
            model_vector[0]
                .populations[{age_group, mio::osecirts::InfectionState::InfectedCriticalPartialImmunity}]
                .value();
        ICU_after[(size_t)age_group + 2] =
            model_vector[0]
                .populations[{age_group, mio::osecirts::InfectionState::InfectedCriticalImprovedImmunity}]
                .value();
    }
    // Test if ICU was changed
    EXPECT_NE(ICU_before, ICU_after);
}

TEST(TestOdeSECIRTS, read_data)
{
    auto num_age_groups = 6; //reading data requires RKI data age groups
    auto model1         = std::vector<mio::osecirts::Model<double>>({make_model(num_age_groups)});
    auto model2         = std::vector<mio::osecirts::Model<double>>({make_model(num_age_groups)});
    auto model3         = std::vector<mio::osecirts::Model<double>>({make_model(num_age_groups)});

    const std::vector<std::vector<double>> immunity_population = {{0.04, 0.04, 0.075, 0.08, 0.035, 0.01},
                                                                  {0.61, 0.61, 0.62, 0.62, 0.58, 0.41},
                                                                  {0.35, 0.35, 0.305, 0.3, 0.385, 0.58}};

    auto read_result1 = mio::osecirts::read_input_data_county(model1, {2020, 12, 01}, {1002},
                                                              std::vector<double>(size_t(num_age_groups), 1.0), 1.0,
                                                              TEST_DATA_DIR, 10, immunity_population);

    auto read_result2 =
        mio::osecirts::read_input_data(model2, {2020, 12, 01}, {1002}, std::vector<double>(size_t(num_age_groups), 1.0),
                                       1.0, TEST_DATA_DIR, 10, immunity_population);

    auto read_result_district =
        mio::osecirts::read_input_data(model3, {2020, 12, 01}, {1002}, std::vector<double>(size_t(num_age_groups), 1.0),
                                       1.0, mio::path_join(TEST_DATA_DIR, "pydata/District"), 10, immunity_population);

    ASSERT_THAT(read_result1, IsSuccess());
    ASSERT_THAT(read_result2, IsSuccess());
    ASSERT_THAT(read_result_district, IsSuccess());

    // values were generated by the tested function; can only check stability of the implementation, not correctness
    auto expected_values =
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecirts::InfectionState::Count)) << 411.933, 6282.58,
         0.213004, 3.19834, 0.840323, 0.0946686, 1.42149, 0.373477, 0, 0, 0, 0.558774, 4.4528, 0.955134, 0, 0, 0,
         0.0103231, 0.00287609, 0.00365536, 0.0297471, 3.5, 4, 3606.1, 0, 0, 0, 0.649663, 0.0771142, 766.409, 11687.7,
         0.49701, 7.4628, 1.96075, 0.319507, 4.79751, 1.26048, 0, 0, 0, 0.93129, 8.36257, 1.79379, 0, 0, 0, 0.0132265,
         0.00302746, 0.00384774, 0.0297471, 3.5, 4, 6712.35, 0, 0, 0, 3.24832, 0.385571, 5568.96, 46059.3, 5.32163,
         43.3153, 9.75737, 2.7505, 22.3877, 5.04313, 0, 0, 0, 8.84909, 39.8328, 7.32558, 0, 0, 0, 0.56467, 0.071343,
         0.0777408, 0.0753594, 3.5, 4, 22691, 0, 0, 0, 13.4975, 1.37363, 6631.47, 51423.9, 5.33096, 40.6794, 9.01336,
         2.83472, 21.6311, 4.79282, 0, 0, 0, 8.55241, 37.2832, 6.74428, 0, 0, 0, 1.78101, 0.218787, 0.234499, 0.487853,
         3.5, 4, 24913, 0, 0, 0, 13.8504, 1.38643, 1532.28, 25449.8, 0.528786, 8.62793, 2.62254, 0.329244, 5.37211,
         1.6329, 0, 0, 0, 0.939561, 8.52246, 2.1149, 0, 0, 0, 0.856992, 0.223414, 0.328498, 2.61775, 3.5, 4, 16901.9, 0,
         0, 0, 3.38605, 0.46498, 151.979, 6426.73, 0.22575, 9.11334, 5.90344, 0.0707574, 2.85642, 1.85033, 0, 0, 0,
         0.0889777, 2.09765, 1.10935, 0, 0, 0, 0.0681386, 0.046887, 0.146922, 4.75954, 3.5, 4, 9103.97, 0, 0, 0,
         0.530478, 0.155246)
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
        EXPECT_GE(double(model1[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsNaive}]), 0);
        EXPECT_GE(double(model2[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsNaive}]), 0);
        EXPECT_GE(double(model3[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsNaive}]), 0);
        EXPECT_GE(double(model1[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsNaive}]), 0);
        EXPECT_GE(double(model2[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsNaive}]), 0);
        EXPECT_GE(double(model3[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsNaive}]), 0);
        EXPECT_GE(double(model1[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(double(model2[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(double(model3[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(double(model1[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(double(model2[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(double(model3[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(double(model1[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunity}]),
                  0);
        EXPECT_GE(double(model2[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunity}]),
                  0);
        EXPECT_GE(double(model3[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunity}]),
                  0);
        EXPECT_GE(double(model1[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunity}]),
                  0);
        EXPECT_GE(double(model2[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunity}]),
                  0);
        EXPECT_GE(double(model3[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunity}]),
                  0);
        EXPECT_GE(double(model2[0].populations[{i, mio::osecirts::InfectionState::TemporaryImmunePartialImmunity}]), 0);
        EXPECT_GE(double(model3[0].populations[{i, mio::osecirts::InfectionState::TemporaryImmuneImprovedImmunity}]),
                  0);

        // currently dead and confirmed after commuting compartments are initialized as zero
        EXPECT_EQ(double(model1[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsNaiveConfirmed}]),
                  0);
        EXPECT_EQ(double(model2[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsNaiveConfirmed}]),
                  0);
        EXPECT_EQ(double(model3[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsNaiveConfirmed}]),
                  0);
        EXPECT_EQ(double(model1[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsNaiveConfirmed}]), 0);
        EXPECT_EQ(double(model2[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsNaiveConfirmed}]), 0);
        EXPECT_EQ(double(model3[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsNaiveConfirmed}]), 0);
        EXPECT_EQ(
            double(
                model1[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model2[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model3[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(model1[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(model2[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(model3[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model1[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model2[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model3[0].populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model1[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model2[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model3[0].populations[{i, mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(double(model1[0].populations[{i, mio::osecirts::InfectionState::DeadNaive}]), 0);
        EXPECT_EQ(double(model2[0].populations[{i, mio::osecirts::InfectionState::DeadNaive}]), 0);
        EXPECT_EQ(double(model3[0].populations[{i, mio::osecirts::InfectionState::DeadNaive}]), 0);
        EXPECT_EQ(double(model1[0].populations[{i, mio::osecirts::InfectionState::DeadPartialImmunity}]), 0);
        EXPECT_EQ(double(model2[0].populations[{i, mio::osecirts::InfectionState::DeadPartialImmunity}]), 0);
        EXPECT_EQ(double(model3[0].populations[{i, mio::osecirts::InfectionState::DeadPartialImmunity}]), 0);
        EXPECT_EQ(double(model1[0].populations[{i, mio::osecirts::InfectionState::DeadImprovedImmunity}]), 0);
        EXPECT_EQ(double(model2[0].populations[{i, mio::osecirts::InfectionState::DeadImprovedImmunity}]), 0);
        EXPECT_EQ(double(model3[0].populations[{i, mio::osecirts::InfectionState::DeadImprovedImmunity}]), 0);
    }
}

TEST(TestOdeSECIRTS, export_time_series_init)
{
    TempFileRegister temp_file_register;
    auto tmp_results_dir = temp_file_register.get_unique_path();
    ASSERT_THAT(mio::create_directory(tmp_results_dir), IsSuccess());

    auto num_age_groups = 6; // Data to be read requires RKI confirmed cases data age groups
    auto model          = make_model(num_age_groups);

    const std::vector<std::vector<double>> immunity_population = {{0.04, 0.04, 0.075, 0.08, 0.035, 0.01},
                                                                  {0.61, 0.61, 0.62, 0.62, 0.58, 0.41},
                                                                  {0.35, 0.35, 0.305, 0.3, 0.385, 0.58}};

    // Test exporting time series
    ASSERT_THAT(mio::osecirts::export_input_data_county_timeseries(
                    std::vector<mio::osecirts::Model<double>>{model}, tmp_results_dir, {0}, {2020, 12, 01},
                    std::vector<double>(size_t(num_age_groups), 1.0), 1.0, 2,
                    mio::path_join(TEST_DATA_DIR, "county_divi_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "cases_all_county_age_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "county_current_population.json"), immunity_population,
                    mio::path_join(TEST_DATA_DIR, "vacc_county_ageinf_ma7.json")),
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
TEST(TestOdeSECIRTS, model_initialization)
{
    auto num_age_groups = 6; // Data to be read requires RKI confirmed cases data age groups
    auto model          = make_model(num_age_groups);
    // Vector assignment necessary as read_input_data_county changes model
    auto model_vector = std::vector<mio::osecirts::Model<double>>{model};

    const std::vector<std::vector<double>> immunity_population = {{0.04, 0.04, 0.075, 0.08, 0.035, 0.01},
                                                                  {0.61, 0.61, 0.62, 0.62, 0.58, 0.41},
                                                                  {0.35, 0.35, 0.305, 0.3, 0.385, 0.58}};

    ASSERT_THAT(mio::osecirts::read_input_data_county(model_vector, {2020, 12, 01}, {0},
                                                      std::vector<double>(size_t(num_age_groups), 1.0), 1.0,
                                                      TEST_DATA_DIR, 2, immunity_population, false),
                IsSuccess());

    // Values from data/export_time_series_init_osecirts.h5, for reading in comparison
    // operator for return of mio::read_result and model population needed.
    auto expected_values =
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecirts::InfectionState::Count)) << 138749, 2.11592e+06,
         0.213004, 3.19834, 0.840323, 0.0946686, 1.42149, 0.373477, 0, 0, 0, 0.558774, 4.4528, 0.955134, 0, 0, 0,
         0.0103231, 0.00287609, 0.00365536, 0.0297471, 3.5, 4, 1.21406e+06, 0, 0, 0, 0.649663, 0.0771142, 309977,
         4.72715e+06, 0.49701, 7.4628, 1.96075, 0.319507, 4.79751, 1.26048, 0, 0, 0, 0.93129, 8.36257, 1.79379, 0, 0, 0,
         0.0132265, 0.00302746, 0.00384774, 0.0297471, 3.5, 4, 2.71231e+06, 0, 0, 0, 3.24832, 0.385571, 1.44129e+06,
         1.19147e+07, 5.32163, 43.3153, 9.75737, 2.7505, 22.3877, 5.04313, 0, 0, 0, 8.84909, 39.8328, 7.32558, 0, 0, 0,
         0.56467, 0.071343, 0.0777408, 0.0753594, 3.5, 4, 5.8613e+06, 0, 0, 0, 13.4975, 1.37363, 2.4027e+06, 1.8621e+07,
         5.33096, 40.6794, 9.01336, 2.83472, 21.6311, 4.79282, 0, 0, 0, 8.55241, 37.2832, 6.74428, 0, 0, 0, 1.78101,
         0.218787, 0.234499, 0.487853, 3.5, 4, 9.01018e+06, 0, 0, 0, 13.8504, 1.38643, 578008, 9.57848e+06, 0.528786,
         8.62793, 2.62254, 0.329244, 5.37211, 1.6329, 0, 0, 0, 0.939561, 8.52246, 2.1149, 0, 0, 0, 0.856992, 0.223414,
         0.328498, 2.61775, 3.5, 4, 6.35814e+06, 0, 0, 0, 3.38605, 0.46498, 61818.1, 2.53474e+06, 0.22575, 9.11334,
         5.90344, 0.0707574, 2.85642, 1.85033, 0, 0, 0, 0.0889777, 2.09765, 1.10935, 0, 0, 0, 0.0681386, 0.046887,
         0.146922, 4.75954, 3.5, 4, 3.58574e+06, 0, 0, 0, 0.530478, 0.155246)
            .finished();

    ASSERT_THAT(print_wrap(model_vector[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));
}

#endif

TEST(TestOdeSECIRTS, parameter_percentiles)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    //build small graph
    auto model = make_model(5);
    auto graph = mio::Graph<mio::osecirts::Model<double>, mio::MobilityParameters<double>>();
    graph.add_node(0, model);

    //sample a few times
    auto sampled_graphs = std::vector<mio::Graph<mio::osecirts::Model<double>, mio::MobilityParameters<double>>>();
    std::generate_n(std::back_inserter(sampled_graphs), 10, [&graph]() {
        return mio::osecirts::draw_sample(graph);
    });

    //extract nodes from graph
    auto sampled_nodes = std::vector<std::vector<mio::osecirts::Model<double>>>();
    std::transform(sampled_graphs.begin(), sampled_graphs.end(), std::back_inserter(sampled_nodes), [](auto&& g) {
        auto models = std::vector<mio::osecirts::Model<double>>();
        std::transform(g.nodes().begin(), g.nodes().end(), std::back_inserter(models), [](auto&& n) {
            return n.property;
        });
        return models;
    });

    //compute percentiles
    auto percentile_params = mio::osecirts::ensemble_params_percentile(sampled_nodes, 0.6)[0].parameters;

    //spot check parameters
    auto p       = double(percentile_params.get<mio::osecirts::ReducTimeInfectedMild<double>>()[mio::AgeGroup(2)]);
    auto samples = std::vector<double>();
    std::transform(sampled_nodes.begin(), sampled_nodes.end(), std::back_inserter(samples),
                   [](const std::vector<mio::osecirts::Model<double>>& nodes) {
                       return nodes[0].parameters.get<mio::osecirts::ReducTimeInfectedMild<double>>()[mio::AgeGroup(2)];
                   });

    std::nth_element(samples.begin(), samples.begin() + 6, samples.end());
    ASSERT_THAT(p, samples[6]);

    p       = double(percentile_params.get<mio::osecirts::TimeExposed<double>>()[mio::AgeGroup(2)]);
    samples = std::vector<double>();
    std::transform(sampled_nodes.begin(), sampled_nodes.end(), std::back_inserter(samples),
                   [](const std::vector<mio::osecirts::Model<double>>& nodes) {
                       return nodes[0].parameters.get<mio::osecirts::TimeExposed<double>>()[mio::AgeGroup(2)];
                   });
    std::nth_element(samples.begin(), samples.begin() + 6, samples.end());
    ASSERT_THAT(p, samples[6]);
}

TEST(TestOdeSECIRTS, get_infections_relative)
{
    auto model = make_model(2);
    auto sim   = mio::osecirts::Simulation<double>(model);
    auto y     = sim.get_result()[0];

    auto relative_infections = mio::osecirts::get_infections_relative<double>(sim, 0.0, y);

    // see model population init to obtain sum 105=2*(7+7.5+8+9.5+10+10.5)
    ASSERT_DOUBLE_EQ(relative_infections, 105 / model.populations.get_total());
}

TEST(TestOdeSECIRTS, get_migration_factors)
{
    auto num_age_groups = 2;
    auto model          = make_model(num_age_groups);
    auto sim            = mio::osecirts::Simulation<double>(model);
    auto y              = sim.get_result()[0];

    auto migration_factors = mio::osecirts::get_migration_factors<double>(sim, 0.0, y);

    auto expected_values = (Eigen::VectorXd(Eigen::Index(mio::osecirts::InfectionState::Count) * num_age_groups) << 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
                               .finished();
    ASSERT_THAT(print_wrap(migration_factors), MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));
}

TEST(TestOdeSECIRTS, test_commuters)
{
    auto model                                      = make_model(2);
    auto migration_factor                           = 0.1;
    auto non_detection_factor                       = 0.3;
    model.parameters.get_start_commuter_detection() = 0.0;
    model.parameters.get_end_commuter_detection()   = 20.0;
    model.parameters.get_commuter_nondetection()    = non_detection_factor;
    auto sim                                        = mio::osecirts::Simulation<>(model);
    auto before_testing                             = sim.get_result().get_last_value().eval();
    auto migrated                                   = (sim.get_result().get_last_value() * migration_factor).eval();
    auto migrated_tested                            = migrated.eval();

    mio::osecirts::test_commuters<double>(sim, migrated_tested, 0.0);

    ASSERT_NEAR(migrated_tested[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsNaive)],
                migrated[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsNaive)] * non_detection_factor,
                1e-5);
    ASSERT_NEAR(
        sim.get_result().get_last_value()[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsNaiveConfirmed)],
        before_testing[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsNaiveConfirmed)] +
            migrated[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsNaive)] * (1 - non_detection_factor),
        1e-5);
    ASSERT_NEAR(migrated_tested[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity)],
                migrated[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity)] *
                    non_detection_factor,
                1e-5);
    ASSERT_NEAR(
        sim.get_result()
            .get_last_value()[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsPartialImmunityConfirmed)],
        before_testing[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsPartialImmunityConfirmed)] +
            migrated[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity)] *
                (1 - non_detection_factor),
        1e-5);
    ASSERT_NEAR(migrated_tested[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunity)],
                migrated[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunity)] *
                    non_detection_factor,
                1e-5);
    ASSERT_NEAR(
        sim.get_result()
            .get_last_value()[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunityConfirmed)],
        before_testing[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunityConfirmed)] +
            migrated[Eigen::Index(mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunity)] *
                (1 - non_detection_factor),
        1e-5);

    ASSERT_NEAR(migrated_tested[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsNaive)],
                migrated[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsNaive)] * non_detection_factor,
                1e-5);
    ASSERT_NEAR(sim.get_result()
                    .get_last_value()[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsNaiveConfirmed)],
                before_testing[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsNaiveConfirmed)] +
                    migrated[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsNaive)] *
                        (1 - non_detection_factor),
                1e-5);
    ASSERT_NEAR(migrated_tested[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunity)],
                migrated[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunity)] *
                    non_detection_factor,
                1e-5);
    ASSERT_NEAR(
        sim.get_result()
            .get_last_value()[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed)],
        before_testing[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed)] +
            migrated[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunity)] *
                (1 - non_detection_factor),
        1e-5);
    ASSERT_NEAR(migrated_tested[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunity)],
                migrated[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunity)] *
                    non_detection_factor,
                1e-5);
    ASSERT_NEAR(
        sim.get_result()
            .get_last_value()[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed)],
        before_testing[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed)] +
            migrated[Eigen::Index(mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunity)] *
                (1 - non_detection_factor),
        1e-5);
}

TEST(TestOdeSECIRTS, check_constraints_parameters)
{
    auto model = mio::osecirts::Model<double>(1);
    ASSERT_EQ(model.parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::osecirts::Seasonality<double>>(-0.2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::Seasonality<double>>(0.2);
    model.parameters.set<mio::osecirts::ICUCapacity<double>>(-2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::ICUCapacity<double>>(2);
    model.parameters.set<mio::osecirts::TestAndTraceCapacity<double>>(-1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::TestAndTraceCapacity<double>>(1);
    model.parameters.set<mio::osecirts::TestAndTraceCapacityMaxRiskNoSymptoms<double>>(-1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::TestAndTraceCapacityMaxRiskNoSymptoms<double>>(1);
    model.parameters.set<mio::osecirts::TestAndTraceCapacityMaxRiskSymptoms<double>>(-1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::TestAndTraceCapacityMaxRiskSymptoms<double>>(1);
    model.parameters.set<mio::osecirts::TimeExposed<double>>(-2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::TimeExposed<double>>(2);
    model.parameters.set<mio::osecirts::TimeInfectedNoSymptoms<double>>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::TimeInfectedNoSymptoms<double>>(5);
    model.parameters.set<mio::osecirts::TimeInfectedSymptoms<double>>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::TimeInfectedSymptoms<double>>(2);
    model.parameters.set<mio::osecirts::TimeInfectedSevere<double>>(-1);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::TimeInfectedSevere<double>>(2);
    model.parameters.set<mio::osecirts::TimeInfectedCritical<double>>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::TimeInfectedCritical<double>>(10.);
    model.parameters.set<mio::osecirts::TimeTemporaryImmunityPI<double>>(0.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::TimeTemporaryImmunityPI<double>>(90.);
    model.parameters.set<mio::osecirts::TimeTemporaryImmunityII<double>>(-20.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::TimeTemporaryImmunityII<double>>(90.);
    model.parameters.set<mio::osecirts::TimeWaningPartialImmunity<double>>(0.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::TimeWaningPartialImmunity<double>>(100.);
    model.parameters.set<mio::osecirts::TimeWaningImprovedImmunity<double>>(0.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::TimeWaningImprovedImmunity<double>>(200);
    model.parameters.set<mio::osecirts::TransmissionProbabilityOnContact<double>>(2.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::TransmissionProbabilityOnContact<double>>(0.5);
    model.parameters.set<mio::osecirts::RelativeTransmissionNoSymptoms<double>>(-1.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::RelativeTransmissionNoSymptoms<double>>(0.5);
    model.parameters.set<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>(3.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>(0.5);
    model.parameters.set<mio::osecirts::RiskOfInfectionFromSymptomatic<double>>(-0.8);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::RiskOfInfectionFromSymptomatic<double>>(0.5);
    model.parameters.set<mio::osecirts::SeverePerInfectedSymptoms<double>>(-0.1);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::SeverePerInfectedSymptoms<double>>(0.5);
    model.parameters.set<mio::osecirts::CriticalPerSevere<double>>(-1.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::CriticalPerSevere<double>>(0.5);
    model.parameters.set<mio::osecirts::DeathsPerCritical<double>>(1.1);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::DeathsPerCritical<double>>(0.5);
    model.parameters.set<mio::osecirts::DaysUntilEffectivePartialVaccination<double>>(-2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::DaysUntilEffectivePartialVaccination<double>>(7);
    model.parameters.set<mio::osecirts::DaysUntilEffectiveImprovedVaccination<double>>(-0.2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::DaysUntilEffectiveImprovedVaccination<double>>(7);
    model.parameters.set<mio::osecirts::DaysUntilEffectiveBoosterImmunity<double>>(-2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::DaysUntilEffectiveBoosterImmunity<double>>(7);
    model.parameters.set<mio::osecirts::ReducExposedPartialImmunity<double>>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::ReducExposedPartialImmunity<double>>(0.5);
    model.parameters.set<mio::osecirts::ReducExposedImprovedImmunity<double>>(-0.2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::ReducExposedImprovedImmunity<double>>(0.5);
    model.parameters.set<mio::osecirts::ReducInfectedSymptomsPartialImmunity<double>>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::ReducInfectedSymptomsPartialImmunity<double>>(0.5);
    model.parameters.set<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<double>>(0.);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<double>>(0.5);
    model.parameters.set<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(-4);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(0.5);
    model.parameters.set<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(-4);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(0.5);
    model.parameters.set<mio::osecirts::ReducTimeInfectedMild<double>>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::ReducTimeInfectedMild<double>>(1);
    model.parameters.set<mio::osecirts::InfectiousnessNewVariant<double>>(-4);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirts::InfectiousnessNewVariant<double>>(2.0);
    ASSERT_EQ(model.parameters.check_constraints(), 0);
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSECIRTS, apply_constraints_parameters)
{
    const double tol_times = 1e-1;
    auto model             = mio::osecirts::Model(1);
    auto indx_agegroup     = mio::AgeGroup(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    model.parameters.set<mio::osecirts::Seasonality<double>>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::Seasonality<double>>(), 0);

    model.parameters.set<mio::osecirts::ICUCapacity<double>>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::Seasonality<double>>(), 0);

    model.parameters.set<mio::osecirts::TestAndTraceCapacity<double>>(-1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::TestAndTraceCapacity<double>>(), 0);

    model.parameters.set<mio::osecirts::TestAndTraceCapacityMaxRiskNoSymptoms<double>>(-1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::TestAndTraceCapacityMaxRiskNoSymptoms<double>>(), 0);

    model.parameters.set<mio::osecirts::TestAndTraceCapacityMaxRiskSymptoms<double>>(-1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::TestAndTraceCapacityMaxRiskSymptoms<double>>(), 0);

    model.parameters.set<mio::osecirts::TimeExposed<double>>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::TimeExposed<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirts::TimeInfectedNoSymptoms<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::TimeInfectedNoSymptoms<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirts::TimeInfectedSymptoms<double>>(1e-10);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::TimeInfectedSymptoms<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirts::TimeInfectedSevere<double>>(-1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::TimeInfectedSevere<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirts::TimeInfectedCritical<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::TimeInfectedCritical<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirts::TimeTemporaryImmunityPI<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::TimeTemporaryImmunityPI<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirts::TimeTemporaryImmunityII<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::TimeTemporaryImmunityII<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirts::TimeWaningPartialImmunity<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::TimeWaningPartialImmunity<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirts::TimeWaningImprovedImmunity<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::TimeWaningImprovedImmunity<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirts::TransmissionProbabilityOnContact<double>>(2.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR(model.parameters.get<mio::osecirts::TransmissionProbabilityOnContact<double>>()[indx_agegroup], 0.0,
                1e-14);

    model.parameters.set<mio::osecirts::RelativeTransmissionNoSymptoms<double>>(-1.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::RelativeTransmissionNoSymptoms<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>(3.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirts::RiskOfInfectionFromSymptomatic<double>>(-0.8);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::RiskOfInfectionFromSymptomatic<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirts::SeverePerInfectedSymptoms<double>>(-0.1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::SeverePerInfectedSymptoms<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirts::CriticalPerSevere<double>>(-1.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::CriticalPerSevere<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirts::DeathsPerCritical<double>>(1.1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::DeathsPerCritical<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirts::DaysUntilEffectivePartialVaccination<double>>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::DaysUntilEffectivePartialVaccination<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirts::DaysUntilEffectiveImprovedVaccination<double>>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::DaysUntilEffectiveImprovedVaccination<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirts::DaysUntilEffectiveBoosterImmunity<double>>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::DaysUntilEffectiveBoosterImmunity<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirts::ReducExposedPartialImmunity<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::ReducExposedPartialImmunity<double>>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirts::ReducExposedImprovedImmunity<double>>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::ReducExposedImprovedImmunity<double>>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirts::ReducInfectedSymptomsPartialImmunity<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<double>>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<double>>(0.);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<double>>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(-4);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(
        model.parameters.get<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[indx_agegroup],
        1);

    model.parameters.set<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(-4);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(
        model.parameters.get<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[indx_agegroup],
        1);

    model.parameters.set<mio::osecirts::ReducTimeInfectedMild<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::ReducTimeInfectedMild<double>>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirts::InfectiousnessNewVariant<double>>(-4);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirts::InfectiousnessNewVariant<double>>()[indx_agegroup], 1);

    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSECIRTS, apply_variant_function)
{
    auto model = mio::osecirts::Model<double>(1);
    model.parameters.set<mio::osecirts::TransmissionProbabilityOnContact<double>>(0.2);

    model.parameters.set<mio::osecirts::StartDay>(0);
    model.parameters.set<mio::osecirts::StartDayNewVariant>(10);
    model.parameters.set<mio::osecirts::InfectiousnessNewVariant<double>>(2.0);

    // set vaccinations
    const size_t daily_vaccinations = 1;
    const size_t num_days           = 5;
    model.parameters.get<mio::osecirts::DailyPartialVaccinations<double>>().resize(mio::SimulationDay(num_days));
    model.parameters.get<mio::osecirts::DailyFullVaccinations<double>>().resize(mio::SimulationDay(num_days));
    model.parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>().resize(mio::SimulationDay(num_days));
    for (size_t i = 0; i < num_days; ++i) {
        auto num_vaccinations = static_cast<double>(i * daily_vaccinations);
        model.parameters
            .get<mio::osecirts::DailyPartialVaccinations<double>>()[{mio::AgeGroup(0), mio::SimulationDay(i)}] =
            num_vaccinations;
        model.parameters
            .get<mio::osecirts::DailyFullVaccinations<double>>()[{mio::AgeGroup(0), mio::SimulationDay(i)}] =
            num_vaccinations;
        model.parameters
            .get<mio::osecirts::DailyBoosterVaccinations<double>>()[{mio::AgeGroup(0), mio::SimulationDay(i)}] =
            num_vaccinations;
    }

    auto sim = mio::osecirts::Simulation<double>(model);

    // test that the transmission probability is not changed due to calling the advance function
    sim.advance(0.01);
    EXPECT_NEAR(
        sim.get_model().parameters.get<mio::osecirts::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)],
        0.2, 1e-10);

    // test if the transmission probability is set to the correct value after applying the variant.
    // The referece values are calculated using equation (36) in doi.org/10.1371/journal.pcbi.1010054
    auto base_infectiousness =
        sim.get_model().parameters.get<mio::osecirts::TransmissionProbabilityOnContact<double>>();

    // however, the parameter should stay unchanged if the new variant is not present in the population.
    sim.apply_variant(0, base_infectiousness);
    EXPECT_NEAR(
        sim.get_model().parameters.get<mio::osecirts::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)],
        0.2, 1e-10);

    sim.apply_variant(9, base_infectiousness);
    EXPECT_NEAR(
        sim.get_model().parameters.get<mio::osecirts::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)],
        0.2, 1e-10);

    sim.apply_variant(10, base_infectiousness);
    EXPECT_NEAR(
        sim.get_model().parameters.get<mio::osecirts::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)],
        0.99 * base_infectiousness[mio::AgeGroup(0)] +
            0.01 * base_infectiousness[mio::AgeGroup(0)] *
                sim.get_model().parameters.get<mio::osecirts::InfectiousnessNewVariant<double>>()[mio::AgeGroup(0)],
        1e-10);

    sim.apply_variant(45, base_infectiousness);
    EXPECT_NEAR(
        sim.get_model().parameters.get<mio::osecirts::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)],
        0.68 * base_infectiousness[mio::AgeGroup(0)] +
            0.32 * base_infectiousness[mio::AgeGroup(0)] *
                sim.get_model().parameters.get<mio::osecirts::InfectiousnessNewVariant<double>>()[mio::AgeGroup(0)],
        1e-10);

    sim.apply_variant(1000, base_infectiousness);
    EXPECT_NEAR(
        sim.get_model().parameters.get<mio::osecirts::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)],
        0.4, 1e-10);
}
