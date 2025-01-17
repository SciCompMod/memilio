/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "gmock/gmock.h"
#include <algorithm>
#include <iterator>
#include <limits>

TEST(TestOdeSECIRVVS, simulateDefault)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.1;

    mio::osecirvvs::Model<double> model(1);
    model.populations.set_total(10);
    model.populations.set_difference_from_total({(mio::AgeGroup)0, mio::osecirvvs::InfectionState::SusceptibleNaive},
                                                10);
    model.parameters.get<mio::osecirvvs::DailyPartialVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirvvs::DailyPartialVaccinations<double>>().array().setConstant(0);
    model.parameters.get<mio::osecirvvs::DailyFullVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirvvs::DailyFullVaccinations<double>>().array().setConstant(0);
    mio::TimeSeries<double> result = simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

TEST(TestOdeSECIRVVS, reduceToSecirAndCompareWithPreviousRun)
{
    // double t0   = 0;
    // double tmax = 50;

    mio::osecirvvs::Model<double> model(1);

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10;

    model.populations.set_total(nb_total_t0);
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::ExposedNaive}]                     = nb_exp_t0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::ExposedImprovedImmunity}]          = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::ExposedPartialImmunity}]           = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}]          = nb_car_t0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}] = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}] = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}] =
        0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}] = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] =
        0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}]           = nb_inf_t0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}]  = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}] = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}] = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}]         = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}] =
        0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedSevereNaive}]             = nb_hosp_t0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedSevereImprovedImmunity}]  = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedSeverePartialImmunity}]   = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedCriticalNaive}]           = nb_icu_t0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedCriticalPartialImmunity}] = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::InfectedCriticalImprovedImmunity}] = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}]      = nb_rec_t0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}]       = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::DeadNaive}]                        = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::DeadPartialImmunity}]              = 0;
    model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::DeadImprovedImmunity}]             = 0;
    model.populations.set_difference_from_total({(mio::AgeGroup)0, mio::osecirvvs::InfectionState::SusceptibleNaive},
                                                nb_total_t0);

    model.parameters.get<mio::osecirvvs::ICUCapacity<double>>()          = 10000;
    model.parameters.get<mio::osecirvvs::TestAndTraceCapacity<double>>() = 10000;
    model.parameters.get<mio::osecirvvs::DailyPartialVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirvvs::DailyPartialVaccinations<double>>().array().setConstant(0);
    model.parameters.get<mio::osecirvvs::DailyFullVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirvvs::DailyFullVaccinations<double>>().array().setConstant(0);

    auto& contacts       = model.parameters.get<mio::osecirvvs::ContactPatterns<double>>();
    auto& contact_matrix = contacts.get_cont_freq_mat();
    contact_matrix[0]    = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    //times
    model.parameters.get<mio::osecirvvs::TimeExposed<double>>()[mio::AgeGroup(0)]            = 3.2;
    model.parameters.get<mio::osecirvvs::TimeInfectedNoSymptoms<double>>()[mio::AgeGroup(0)] = 2.;
    model.parameters.get<mio::osecirvvs::TimeInfectedSymptoms<double>>()[mio::AgeGroup(0)]   = 5;
    model.parameters.get<mio::osecirvvs::TimeInfectedSevere<double>>()[mio::AgeGroup(0)]     = 10;
    model.parameters.get<mio::osecirvvs::TimeInfectedCritical<double>>()[mio::AgeGroup(0)]   = 8;

    //probabilities
    model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)] = 0.05;
    model.parameters.get<mio::osecirvvs::RelativeTransmissionNoSymptoms<double>>()[mio::AgeGroup(0)]   = 1;
    model.parameters.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>()[mio::AgeGroup(0)]   = 0.25;
    model.parameters.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>()[mio::AgeGroup(0)]   = 0.09;
    model.parameters.get<mio::osecirvvs::SeverePerInfectedSymptoms<double>>()[mio::AgeGroup(0)]        = 0.2;
    model.parameters.get<mio::osecirvvs::CriticalPerSevere<double>>()[mio::AgeGroup(0)]                = 0.25;
    model.parameters.get<mio::osecirvvs::DeathsPerCritical<double>>()[mio::AgeGroup(0)]                = 0.3;

    // TODO: Reduction not possible like this, division by zero!
    model.parameters.get<mio::osecirvvs::ReducExposedPartialImmunity<double>>()[mio::AgeGroup(0)]           = 1.0;
    model.parameters.get<mio::osecirvvs::ReducExposedImprovedImmunity<double>>()[mio::AgeGroup(0)]          = 1.0;
    model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>()[mio::AgeGroup(0)]  = 1.0;
    model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>()[mio::AgeGroup(0)] = 0;
    model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[mio::AgeGroup(0)] =
        0;
    model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[mio::AgeGroup(0)] =
        0;
    model.parameters.get<mio::osecirvvs::ReducTimeInfectedMild<double>>()[mio::AgeGroup(0)] = 1;

    model.parameters.get<mio::osecirvvs::Seasonality<double>>() = 0.2;

    mio::set_log_level(mio::LogLevel::err);
    model.apply_constraints();
    mio::set_log_level(mio::LogLevel::warn);
    // TODO: gets stuck by division by zero!!
    // auto integrator = std::make_shared<mio::RKIntegratorCore>();
    // integrator->set_dt_min(0.3);
    // integrator->set_dt_max(1.0);
    // integrator->set_rel_tolerance(1e-4);
    // integrator->set_abs_tolerance(1e-1);
    // mio::TimeSeries<double> secihurd = simulate(t0, tmax, 0.1, model, integrator);

    // auto compare = load_test_data_csv<double>("secihurd-compare.csv");

    // ASSERT_EQ(compare.size(), static_cast<size_t>(secihurd.get_num_time_points()));
    // for (size_t i = 0; i < compare.size(); i++) {
    //     ASSERT_EQ(compare[i].size(), static_cast<size_t>(secihurd.get_num_elements()) + 1) << "at row " << i;
    //     EXPECT_NEAR(secihurd.get_time(i), compare[i][0], 1e-10) << "at row " << i;
    //     for (size_t j = 1; j < compare[i].size(); j++) {
    //         // TODO: extract naive compartments
    //         EXPECT_NEAR(secihurd.get_value(i)[j - 1], compare[i][j], 1e-10) << " at row " << i;
    //     }
    // }
}

void assign_uniform_distribution(mio::UncertainValue<double>& p, double min, double max, bool set_invalid_initial_value)
{
    auto invalid_initial = max == 0 ? 1.0 : max * 1.001;
    auto valid_initial   = (max + min) * 0.5;
    auto initial         = set_invalid_initial_value ? invalid_initial : valid_initial;
    p                    = mio::UncertainValue<double>(initial);
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

void set_synthetic_population_data(mio::osecirvvs::Model<double>::Populations& populations,
                                   bool set_invalid_initial_value)
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

void set_demographic_parameters(mio::osecirvvs::Model<double>::ParameterSet& parameters, bool set_invalid_initial_value)
{
    assign_uniform_distribution(parameters.get<mio::osecirvvs::ICUCapacity<double>>(), 20, 50,
                                set_invalid_initial_value);
    parameters.get<mio::osecirvvs::DailyPartialVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    parameters.get<mio::osecirvvs::DailyPartialVaccinations<double>>().array().setConstant(5);
    parameters.get<mio::osecirvvs::DailyFullVaccinations<double>>().resize(mio::SimulationDay(size_t(1000)));
    parameters.get<mio::osecirvvs::DailyFullVaccinations<double>>().array().setConstant(3);
}

void set_contact_parameters(mio::osecirvvs::Model<double>::ParameterSet& parameters, bool set_invalid_initial_value)
{
    auto& contacts       = parameters.get<mio::osecirvvs::ContactPatterns<double>>();
    auto& contact_matrix = contacts.get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(0.5);
    contact_matrix[0].get_baseline().diagonal().setConstant(5.0);
    contact_matrix[0].add_damping(0.3, mio::SimulationTime(5.0));

    auto& npis      = parameters.get<mio::osecirvvs::DynamicNPIsInfectedSymptoms<double>>();
    auto npi_groups = Eigen::VectorXd::Ones(contact_matrix[0].get_num_groups());
    auto npi_value  = mio::UncertainValue<double>(0.5);
    assign_uniform_distribution(npi_value, 0.25, 0.75, set_invalid_initial_value);
    npis.set_threshold(10.0, {mio::DampingSampling<double>(npi_value, mio::DampingLevel(0), mio::DampingType(0),
                                                           mio::SimulationTime(0), {0}, npi_groups)});
    npis.set_base_value(100'000);
    npis.set_interval(mio::SimulationTime(3.0));
    npis.set_duration(mio::SimulationTime(14.0));
    parameters.get_end_dynamic_npis() = 10.0; //required for dynamic NPIs to have effect in this model
    parameters.template get<mio::osecirvvs::DynamicNPIsImplementationDelay<double>>() = 7;
}

void set_covid_parameters(mio::osecirvvs::Model<double>::ParameterSet& params, bool set_invalid_initial_value)
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

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeExposed<double>>(), timeExposedMin, timeExposedMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedNoSymptoms<double>>(),
                                      timeInfectedNoSymptomsMin, timeInfectedNoSymptomsMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSymptoms<double>>(),
                                      timeInfectedSymptomsMin, timeInfectedSymptomsMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSevere<double>>(), timeInfectedSevereMin,
                                      timeInfectedSevereMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedCritical<double>>(),
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

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>(),
                                      transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RelativeTransmissionNoSymptoms<double>>(),
                                      relativeTransmissionNoSymptomsMin, relativeTransmissionNoSymptomsMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>(),
                                      riskOfInfectionFromSymptomaticMin, riskOfInfectionFromSymptomaticMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<double>>(),
                                      maxRiskOfInfectionFromSymptomaticMin, maxRiskOfInfectionFromSymptomaticMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>(),
                                      recoveredPerInfectedNoSymptomsMin, recoveredPerInfectedNoSymptomsMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::SeverePerInfectedSymptoms<double>>(),
                                      severePerInfectedSymptomsMin, severePerInfectedSymptomsMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::CriticalPerSevere<double>>(), criticalPerSevereMin,
                                      criticalPerSevereMax, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::DeathsPerCritical<double>>(), deathsPerCriticalMin,
                                      deathsPerCriticalMax, set_invalid_initial_value);

    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedPartialImmunity<double>>(),
                                      reducExposedPartialImmunityMin, reducExposedPartialImmunityMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedImprovedImmunity<double>>(),
                                      reducExposedImprovedImmunityMin, reducExposedImprovedImmunityMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>(),
                                      reducInfectedSymptomsPartialImmunityMin, reducInfectedSymptomsPartialImmunityMax,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>(),
                                      reducInfectedSymptomsImprovedImmunityMin,
                                      reducInfectedSymptomsImprovedImmunityMax, set_invalid_initial_value);
    array_assign_uniform_distribution(
        params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(),
        reducInfectedSevereCriticalDeadPartialImmunityMin, reducInfectedSevereCriticalDeadPartialImmunityMax,
        set_invalid_initial_value);
    array_assign_uniform_distribution(
        params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(),
        reducInfectedSevereCriticalDeadImprovedImmunityMin, reducInfectedSevereCriticalDeadImprovedImmunityMax,
        set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducTimeInfectedMild<double>>(),
                                      reducTimeInfectedMildMin, reducTimeInfectedMildMax, set_invalid_initial_value);

    //sasonality
    const double seasonality_min = 0.1;
    const double seasonality_max = 0.3;

    assign_uniform_distribution(params.get<mio::osecirvvs::Seasonality<double>>(), seasonality_min, seasonality_max,
                                set_invalid_initial_value);
}

mio::osecirvvs::Model<double> make_model(int num_age_groups, bool set_invalid_initial_value = false)
{
    assert(num_age_groups <= 6 && "Provide more values in functions above to test more age groups.");
    mio::osecirvvs::Model<double> model(num_age_groups);
    set_covid_parameters(model.parameters, set_invalid_initial_value);
    set_synthetic_population_data(model.populations, set_invalid_initial_value);
    set_demographic_parameters(model.parameters, set_invalid_initial_value);
    set_contact_parameters(model.parameters, set_invalid_initial_value);
    model.parameters.apply_constraints();
    return model;
}

TEST(TestOdeSECIRVVS, draw_sample)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>> graph;

    auto num_age_groups = 6;
    //create model with invalid initials so the test fails if no sampling is done
    graph.add_node(0, make_model(num_age_groups, /*set_invalid_initial_value*/ true));
    graph.add_node(0, make_model(num_age_groups, /*set_invalid_initial_value*/ true));
    graph.add_edge(0, 1, Eigen::VectorXd::Constant(num_age_groups, num_age_groups));

    auto sampled_graph = mio::osecirvvs::draw_sample(graph, true);

    ASSERT_EQ(sampled_graph.nodes().size(), graph.nodes().size());
    ASSERT_EQ(sampled_graph.edges().size(), graph.edges().size());

    // spot check for sampling
    auto& parameters0          = sampled_graph.nodes()[0].property.parameters;
    auto& populations0         = sampled_graph.nodes()[0].property.populations;
    auto& timeInfectedCritical = parameters0.get<mio::osecirvvs::TimeInfectedCritical<double>>()[mio::AgeGroup(1)];
    ASSERT_GE(double(timeInfectedCritical), 4.95);
    ASSERT_LE(double(timeInfectedCritical), 8.95);
    auto& param_exp_factor = parameters0.get<mio::osecirvvs::ReducExposedPartialImmunity<double>>()[mio::AgeGroup(0)];
    ASSERT_GE(double(param_exp_factor), 0.75);
    ASSERT_LE(double(param_exp_factor), 0.85);
    auto& compartment_inf =
        populations0[{mio::AgeGroup(2), mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}];
    ASSERT_GE(double(compartment_inf), 5.0);
    ASSERT_LE(double(compartment_inf), 10.0);
    auto& npi_value = parameters0.get<mio::osecirvvs::DynamicNPIsInfectedSymptoms<double>>()
                          .get_thresholds()[0]
                          .second[0]
                          .get_value();
    ASSERT_GE(double(npi_value), 0.25);
    ASSERT_LE(double(npi_value), 0.75);

    // special cases
    ASSERT_NEAR(populations0.get_total(), 1000 * num_age_groups, 1e-2);
    ASSERT_TRUE(
        (parameters0.get<mio::osecirvvs::InfectiousnessNewVariant<double>>().array(),
         parameters0.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>().array() * 1.0) //using high variant
            .all());

    // spot check for parameters that should be equal or different between nodes
    auto& parameters1  = sampled_graph.nodes()[1].property.parameters;
    auto& populations1 = sampled_graph.nodes()[1].property.populations;
    ASSERT_EQ(parameters1.get<mio::osecirvvs::DynamicNPIsInfectedSymptoms<double>>()
                  .get_thresholds()[0]
                  .second[0]
                  .get_value(),
              parameters0.get<mio::osecirvvs::DynamicNPIsInfectedSymptoms<double>>()
                  .get_thresholds()[0]
                  .second[0]
                  .get_value());
    ASSERT_TRUE((parameters1.get<mio::osecirvvs::TimeInfectedSymptoms<double>>().array() ==
                 parameters0.get<mio::osecirvvs::TimeInfectedSymptoms<double>>().array())
                    .all());
    //these could fail in very(!) rare cases if they are randomly sampled to the same value
    ASSERT_NE(parameters1.get<mio::osecirvvs::ICUCapacity<double>>(),
              parameters0.get<mio::osecirvvs::ICUCapacity<double>>())
        << "Failure might be spurious, check RNG seeds.";
    ASSERT_FALSE((populations1.array() == populations0.array()).all()) << "Failure might be spurious, check RNG seeds.";
}

TEST(TestOdeSECIRVVS, checkPopulationConservation)
{
    auto num_age_groups = 2;
    auto model          = make_model(num_age_groups);
    auto num_days       = 30;

    auto result = mio::osecirvvs::simulate<double>(0, num_days, 0.1, model);

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
    auto model          = std::vector<mio::osecirvvs::Model<double>>({make_model(num_age_groups)});
    std::vector<int> region{1002};
    auto path = mio::path_join(TEST_DATA_DIR, "pydata/Germany/cases_all_county_age_ma7.json");
    std::vector<std::vector<int>> t_Exposed(1);
    std::vector<std::vector<int>> t_InfectedNoSymptoms(1);
    std::vector<std::vector<int>> t_InfectedSymptoms(1);
    std::vector<std::vector<int>> t_InfectedSevere(1);
    std::vector<std::vector<int>> t_InfectedCritical(1);

    std::vector<std::vector<double>> mu_C_R(1);
    std::vector<std::vector<double>> mu_I_H(1);
    std::vector<std::vector<double>> mu_H_U(1);

    std::vector<std::vector<double>> num_InfectedSymptoms(1);
    std::vector<std::vector<double>> num_death(1);
    std::vector<std::vector<double>> num_rec(1);
    std::vector<std::vector<double>> num_Exposed(1);
    std::vector<std::vector<double>> num_InfectedNoSymptoms(1);
    std::vector<std::vector<double>> num_InfectedSevere(1);
    std::vector<std::vector<double>> num_icu(1);

    num_InfectedSymptoms[0]   = std::vector<double>(num_age_groups, 0.0);
    num_death[0]              = std::vector<double>(num_age_groups, 0.0);
    num_rec[0]                = std::vector<double>(num_age_groups, 0.0);
    num_Exposed[0]            = std::vector<double>(num_age_groups, 0.0);
    num_InfectedNoSymptoms[0] = std::vector<double>(num_age_groups, 0.0);
    num_InfectedSevere[0]     = std::vector<double>(num_age_groups, 0.0);
    num_icu[0]                = std::vector<double>(num_age_groups, 0.0);
    for (size_t group = 0; group < static_cast<size_t>(num_age_groups); group++) {

        t_Exposed[0].push_back(static_cast<int>(
            std::round(model[0].parameters.template get<mio::osecirvvs::TimeExposed<double>>()[(mio::AgeGroup)group])));
        t_InfectedNoSymptoms[0].push_back(static_cast<int>(std::round(
            model[0].parameters.template get<mio::osecirvvs::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)group])));
        t_InfectedSymptoms[0].push_back(static_cast<int>(std::round(
            model[0].parameters.template get<mio::osecirvvs::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)group])));
        t_InfectedSevere[0].push_back(static_cast<int>(std::round(
            model[0].parameters.template get<mio::osecirvvs::TimeInfectedSevere<double>>()[(mio::AgeGroup)group])));
        t_InfectedCritical[0].push_back(static_cast<int>(std::round(
            model[0].parameters.template get<mio::osecirvvs::TimeInfectedCritical<double>>()[(mio::AgeGroup)group])));

        mu_C_R[0].push_back(model[0].parameters.template get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>()[(
            mio::AgeGroup)group]);
        mu_I_H[0].push_back(
            model[0]
                .parameters.template get<mio::osecirvvs::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)group]);
        mu_H_U[0].push_back(
            model[0].parameters.template get<mio::osecirvvs::CriticalPerSevere<double>>()[(mio::AgeGroup)group]);
    }

    auto read = mio::osecirvvs::details::read_confirmed_cases_data(
        path, region, {2020, 12, 01}, num_Exposed, num_InfectedNoSymptoms, num_InfectedSymptoms, num_InfectedSevere,
        num_icu, num_death, num_rec, t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere,
        t_InfectedCritical, mu_C_R, mu_I_H, mu_H_U, std::vector<double>(size_t(num_age_groups), 1.0));

    ASSERT_THAT(read, IsSuccess());
}

TEST(TestOdeSECIRVVS, set_divi_data_invalid_dates)
{
    mio::set_log_level(mio::LogLevel::off);
    auto model = mio::osecirvvs::Model<double>(1);
    model.populations.array().setConstant(1);
    auto model_vector = std::vector<mio::osecirvvs::Model<double>>{model};

    // Test with date before DIVI dataset was available.
    EXPECT_THAT(mio::osecirvvs::details::set_divi_data(model_vector, "", {1001}, {2019, 12, 01}, 1.0), IsSuccess());
    // Assure that populations is the same as before.
    EXPECT_THAT(print_wrap(model_vector[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(model.populations.array().cast<double>()), 1e-10, 1e-10));

    // Test with data after DIVI dataset was no longer updated.
    EXPECT_THAT(mio::osecirvvs::details::set_divi_data(model_vector, "", {1001}, {2025, 12, 01}, 1.0), IsSuccess());
    EXPECT_THAT(print_wrap(model_vector[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(model.populations.array().cast<double>()), 1e-10, 1e-10));

    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSECIRVVS, set_confirmed_cases_data_with_ICU)
{
    auto num_age_groups = 6;
    auto model          = mio::osecirvvs::Model<double>(num_age_groups);
    model.populations.array().setConstant(1);

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
            model.populations[{age_group, mio::osecirvvs::InfectionState::InfectedCriticalNaive}].value();
        ICU_before[(size_t)age_group + 1] =
            model.populations[{age_group, mio::osecirvvs::InfectionState::InfectedCriticalPartialImmunity}].value();
        ICU_before[(size_t)age_group + 2] =
            model.populations[{age_group, mio::osecirvvs::InfectionState::InfectedCriticalImprovedImmunity}].value();
    }

    // get day in mid of the data
    auto mid_day = case_data[(size_t)case_data.size() / 2].date;

    auto model_vector       = std::vector<mio::osecirvvs::Model<double>>{model};
    auto scaling_factor_inf = std::vector<double>(size_t(model.parameters.get_num_groups()), 1.0);
    EXPECT_THAT(
        mio::osecirvvs::details::set_confirmed_cases_data(model_vector, case_data, {1002}, mid_day, scaling_factor_inf),
        IsSuccess());

    // Get new setted ICU compartment
    auto ICU_after = std::vector<double>(size_t(num_age_groups) * 3, 0.0);
    for (auto age_group = mio::AgeGroup(0); age_group < (mio::AgeGroup)num_age_groups; age_group++) {
        ICU_after[(size_t)age_group] =
            model_vector[0].populations[{age_group, mio::osecirvvs::InfectionState::InfectedCriticalNaive}].value();
        ICU_after[(size_t)age_group + 1] =
            model_vector[0]
                .populations[{age_group, mio::osecirvvs::InfectionState::InfectedCriticalPartialImmunity}]
                .value();
        ICU_after[(size_t)age_group + 2] =
            model_vector[0]
                .populations[{age_group, mio::osecirvvs::InfectionState::InfectedCriticalImprovedImmunity}]
                .value();
    }
    // Test if ICU was changed
    EXPECT_NE(ICU_before, ICU_after);
}

TEST(TestOdeSECIRVVS, read_data)
{
    auto num_age_groups = 6; //reading data requires RKI data age groups
    auto model1         = std::vector<mio::osecirvvs::Model<double>>({make_model(num_age_groups)});
    auto model2         = std::vector<mio::osecirvvs::Model<double>>({make_model(num_age_groups)});
    auto model3         = std::vector<mio::osecirvvs::Model<double>>({make_model(num_age_groups)});

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
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecirvvs::InfectionState::Count)) << 8792.15, 175.889,
         3.21484, 0.0633116, 0.221057, 1.42882, 0.0351731, 0.29682, 0, 0, 0, 6.93838, 0.0725173, 0.206715, 0, 0, 0,
         0.0337498, 1.23324e-05, 0.000208293, 0.0292822, 5.8568e-05, 0.000406386, 1340.42, 0, 0, 0, 17067.6, 220.137,
         7.64078, 0.0970237, 0.381933, 4.91193, 0.0779655, 0.741778, 0, 0, 0, 11.7286, 0.0890643, 0.286235, 0, 0, 0,
         0.0434344, 8.40756e-06, 0.000160098, 0.0294125, 3.7932e-05, 0.000296738, 1891.19, 0, 0, 0, 72501, 176.267,
         47.227, 0.113013, 0.490073, 24.4094, 0.0730141, 0.765246, 0, 0, 0, 64.6789, 0.0855947, 0.303032, 0, 0, 0,
         1.23754, 4.5968e-05, 0.000964262, 0.0751837, 1.82724e-05, 0.000157466, 1670.26, 0, 0, 0, 80790.1, 184.645,
         44.5477, 0.100229, 0.50512, 23.6881, 0.0666206, 0.811467, 0, 0, 0, 58.9805, 0.0758111, 0.31192, 0, 0, 0,
         3.75961, 0.000136175, 0.00331973, 0.486628, 0.000111199, 0.00111367, 2022.58, 0, 0, 0, 41581, 177.478, 9.27393,
         0.0389771, 0.216151, 5.77433, 0.030336, 0.4066, 0, 0, 0, 13.3664, 0.0312302, 0.141394, 0, 0, 0, 3.119,
         0.000209444, 0.00561852, 2.60439, 0.00111169, 0.0122515, 2136.6, 0, 0, 0, 13223.8, 216.037, 11.1838, 0.179986,
         0.863926, 3.50537, 0.0705169, 0.818075, 0, 0, 0, 3.52982, 0.0331744, 0.130002, 0, 0, 0, 0.695168, 0.000190699,
         0.00442784, 4.67895, 0.00764769, 0.0729502, 2253.61, 0, 0, 0)
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
                    std::vector<mio::osecirvvs::Model<double>>{model}, tmp_results_dir, {0}, {2020, 12, 01},
                    std::vector<double>(size_t(num_age_groups), 1.0), 1.0, 2,
                    mio::path_join(TEST_DATA_DIR, "county_divi_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "cases_all_county_age_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "county_current_population.json"),
                    mio::path_join(TEST_DATA_DIR, "vacc_county_ageinf_ma7.json")),
                IsSuccess());

    auto data_extrapolated = mio::read_result(mio::path_join(tmp_results_dir, "Results_rki.h5"));
    ASSERT_THAT(data_extrapolated, IsSuccess());

    // Values were generated by the tested function export_input_data_county_timeseries;
    // can only check stability of the implementation, not correctness
    auto expected_results =
        mio::read_result(mio::path_join(TEST_DATA_DIR, "export_time_series_initialization_osecirvvs.h5")).value();

    ASSERT_THAT(print_wrap(data_extrapolated.value()[0].get_groups().matrix()),
                MatrixNear(print_wrap(expected_results[0].get_groups().matrix()), 1e-5, 1e-5));
}

TEST(TestOdeSECIRVVS, export_time_series_init_old_date)
{
    mio::set_log_level(mio::LogLevel::off);
    TempFileRegister temp_file_register;
    auto tmp_results_dir = temp_file_register.get_unique_path();
    ASSERT_THAT(mio::create_directory(tmp_results_dir), IsSuccess());

    auto num_age_groups = 6; // Data to be read requires RKI confirmed cases data age groups
    auto model          = make_model(num_age_groups);

    // set vaccinations to zero
    model.parameters.get<mio::osecirvvs::DailyPartialVaccinations<double>>().array().setConstant(0);
    model.parameters.get<mio::osecirvvs::DailyFullVaccinations<double>>().array().setConstant(0);
    // set all compartments to zero
    model.populations.array().setConstant(0.0);

    // Test exporting time series
    ASSERT_THAT(mio::osecirvvs::export_input_data_county_timeseries(
                    std::vector<mio::osecirvvs::Model<double>>{model}, tmp_results_dir, {0}, {20, 12, 01},
                    std::vector<double>(size_t(num_age_groups), 1.0), 1.0, 0,
                    mio::path_join(TEST_DATA_DIR, "county_divi_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "cases_all_county_age_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "county_current_population.json"),
                    mio::path_join(TEST_DATA_DIR, "vacc_county_ageinf_ma7.json")),
                IsSuccess());

    auto data_extrapolated = mio::read_result(mio::path_join(tmp_results_dir, "Results_rki.h5"));
    ASSERT_THAT(data_extrapolated, IsSuccess());
    auto results_extrapolated = data_extrapolated.value()[0].get_groups().get_value(0);

    // if we enter an old date, the model only should be initialized with the population data.
    // read population data
    std::string path = mio::path_join(TEST_DATA_DIR, "county_current_population.json");
    const std::vector<int> region{0};
    auto population_data = mio::osecirvvs::details::read_population_data(path, region).value();

    // So, the expected values are the population data in the susceptible compartments and zeros in the other compartments.
    for (auto i = 0; i < num_age_groups; i++) {
        EXPECT_NEAR(results_extrapolated(i * Eigen::Index(mio::osecirvvs::InfectionState::Count)),
                    population_data[0][i], 1e-10);
    }
    // sum of all compartments should be equal to the population
    EXPECT_NEAR(results_extrapolated.sum(), std::accumulate(population_data[0].begin(), population_data[0].end(), 0.0),
                1e-5);
    mio::set_log_level(mio::LogLevel::warn);
}

// Model initialization should return same start values as export time series on that day
TEST(TestOdeSECIRVVS, model_initialization)
{
    auto num_age_groups = 6; // Data to be read requires RKI confirmed cases data age groups
    auto model          = make_model(num_age_groups);
    // Vector assignment necessary as read_input_data_county changes model
    auto model_vector = std::vector<mio::osecirvvs::Model<double>>{model};

    ASSERT_THAT(mio::osecirvvs::read_input_data_county(model_vector, {2020, 12, 01}, {0},
                                                       std::vector<double>(size_t(num_age_groups), 1.0), 1.0,
                                                       TEST_DATA_DIR, 2, false),
                IsSuccess());

    // Values from data/export_time_series_init_osecirvvs.h5, for reading in comparison
    // operator for return of mio::read_result and model population needed.
    auto expected_values =
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecirvvs::InfectionState::Count)) << 3.46722e+06, 176.06,
         3.42799, 0.00017139, 0.00059842, 1.52355, 9.52168e-05, 0.000803518, 0, 0, 0, 7.28479, 0.000193296, 0.000551002,
         0, 0, 0, 0.0342843, 3.1805e-08, 5.37183e-07, 0.029746, 1.51045e-07, 1.04806e-06, 1340.35, 0, 0, 0, 7.74734e+06,
         220.4, 7.99917, 0.000224064, 0.000882024, 5.14232, 0.000180051, 0.00171304, 0, 0, 0, 12.1419, 0.000203391,
         0.000653659, 0, 0, 0, 0.0439275, 1.87568e-08, 3.5717e-07, 0.0297464, 8.46241e-08, 6.62006e-07, 1891.15, 0, 0,
         0, 1.92155e+07, 176.538, 47.6768, 0.00043128, 0.00187022, 24.642, 0.000278636, 0.00292033, 0, 0, 0, 65.1411,
         0.000325876, 0.0011537, 0, 0, 0, 1.24042, 1.74173e-07, 3.65358e-06, 0.0753588, 6.9234e-08, 5.96637e-07,
         1671.81, 0, 0, 0, 3.00317e+07, 184.888, 44.9988, 0.000272769, 0.00137466, 23.9279, 0.000181305, 0.00220837, 0,
         0, 0, 59.4274, 0.000205796, 0.000846734, 0, 0, 0, 3.76905, 3.67799e-07, 8.9664e-06, 0.48785, 3.00341e-07,
         3.00797e-06, 2022.51, 0, 0, 0, 1.65123e+07, 177.579, 9.4638, 0.000100211, 0.00055573, 5.89255, 7.79946e-05,
         0.00104538, 0, 0, 0, 13.5709, 7.98864e-05, 0.000361685, 0, 0, 0, 3.13496, 5.30384e-07, 1.4228e-05, 2.61772,
         2.81518e-06, 3.1025e-05, 2136.56, 0, 0, 0, 6.17983e+06, 216.328, 11.9625, 0.000412312, 0.00197908, 3.74944,
         0.00016154, 0.00187405, 0, 0, 0, 3.71387, 7.47535e-05, 0.000292941, 0, 0, 0, 0.707117, 4.15435e-07,
         9.64602e-06, 4.75937, 1.66604e-05, 0.000158922, 2253.59, 0, 0, 0)
            .finished();

    ASSERT_THAT(print_wrap(model_vector[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));
}

TEST(TestOdeSECIRVVS, model_initialization_old_date)
{
    mio::set_log_level(mio::LogLevel::off);
    constexpr auto num_age_groups = 6; // Data to be read requires RKI confirmed cases data age groups
    auto model                    = make_model(num_age_groups);
    // set vaccinations to zero
    model.parameters.get<mio::osecirvvs::DailyPartialVaccinations<double>>().array().setConstant(0);
    model.parameters.get<mio::osecirvvs::DailyFullVaccinations<double>>().array().setConstant(0);
    // set all compartments to zero
    model.populations.array().setConstant(0.0);

    auto model_vector = std::vector<mio::osecirvvs::Model<double>>{model};

    ASSERT_THAT(mio::osecirvvs::read_input_data(model_vector, {100, 12, 01}, {0},
                                                std::vector<double>(size_t(num_age_groups), 1.0), 1.0, TEST_DATA_DIR, 0,
                                                false),
                IsSuccess());

    // if we enter an old date, the model only should be initialized with the population data.
    // read population data
    std::string path = mio::path_join(TEST_DATA_DIR, "county_current_population.json");
    const std::vector<int> region{0};
    auto population_data = mio::osecirvvs::details::read_population_data(path, region).value();

    // So, the expected values are the population data in the susceptible compartments and zeros in the other compartments.
    for (auto i = 0; i < num_age_groups; i++) {
        EXPECT_NEAR(
            model_vector[0].populations.array().cast<double>()(i * Eigen::Index(mio::osecirvvs::InfectionState::Count)),
            population_data[0][i], 1e-5);
    }

    // sum of all compartments should be equal to the population
    EXPECT_NEAR(model_vector[0].populations.array().cast<double>().sum(),
                std::accumulate(population_data[0].begin(), population_data[0].end(), 0.0), 1e-5);
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSECIRVVS, model_initialization_old_date_county)
{
    mio::set_log_level(mio::LogLevel::off);
    constexpr auto num_age_groups = 6; // Data to be read requires RKI confirmed cases data age groups
    auto model                    = make_model(num_age_groups);
    // set vaccinations to zero
    model.parameters.get<mio::osecirvvs::DailyPartialVaccinations<double>>().array().setConstant(0);
    model.parameters.get<mio::osecirvvs::DailyFullVaccinations<double>>().array().setConstant(0);
    // set all compartments to zero
    model.populations.array().setConstant(0.0);

    auto model_vector = std::vector<mio::osecirvvs::Model<double>>{model};

    ASSERT_THAT(mio::osecirvvs::read_input_data_county(model_vector, {100, 12, 01}, {0},
                                                       std::vector<double>(size_t(num_age_groups), 1.0), 1.0,
                                                       TEST_DATA_DIR, 0, false),
                IsSuccess());

    // if we enter an old date, the model only should be initialized with the population data.
    // read population data
    std::string path = mio::path_join(TEST_DATA_DIR, "county_current_population.json");
    const std::vector<int> region{0};
    auto population_data = mio::osecirvvs::details::read_population_data(path, region).value();

    // So, the expected values are the population data in the susceptible compartments and zeros in the other compartments.
    for (auto i = 0; i < num_age_groups; i++) {
        EXPECT_NEAR(
            model_vector[0].populations.array().cast<double>()(i * Eigen::Index(mio::osecirvvs::InfectionState::Count)),
            population_data[0][i], 1e-5);
    }

    // sum of all compartments should be equal to the population
    EXPECT_NEAR(model_vector[0].populations.array().cast<double>().sum(),
                std::accumulate(population_data[0].begin(), population_data[0].end(), 0.0), 1e-5);
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSECIRVVS, set_population_data_overflow_vacc)
{
    auto num_age_groups = 6; // Data to be read requires RKI confirmed cases data age groups
    auto model          = make_model(num_age_groups);
    // set all compartments to zero
    model.populations.array().setConstant(0.0);

    model.parameters
        .template get<mio::osecirvvs::DailyPartialVaccinations<double>>()[{mio::AgeGroup(0), mio::SimulationDay(0)}] =
        1e7 + 1;

    model.parameters
        .template get<mio::osecirvvs::DailyFullVaccinations<double>>()[{mio::AgeGroup(0), mio::SimulationDay(0)}] = 1e7;

    auto model_vector = std::vector<mio::osecirvvs::Model<double>>{model};

    std::string path_pop_data = mio::path_join(TEST_DATA_DIR, "county_current_population.json");
    const std::vector<int> region{0};
    auto population_data = mio::osecirvvs::details::read_population_data(path_pop_data, region).value();

    // we choose the date so that no case data is available
    ASSERT_THAT(mio::osecirvvs::details::set_population_data(
                    model_vector, path_pop_data, mio::path_join(TEST_DATA_DIR, "cases_all_county_age_ma7.json"), {0},
                    {1000, 12, 01}),
                IsSuccess());

    EXPECT_NEAR(
        double(model_vector[0].populations[{mio::AgeGroup(0), mio::osecirvvs::InfectionState::SusceptibleNaive}]), 0.0,
        1e-10);
    EXPECT_NEAR(double(model_vector[0].populations[{mio::AgeGroup(0),
                                                    mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}]),
                1, 1e-10);
    EXPECT_NEAR(double(model_vector[0].populations[{mio::AgeGroup(0),
                                                    mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}]),
                population_data[0][0] - 1, 1e-10);

    EXPECT_NEAR(model_vector[0].populations.get_group_total(mio::AgeGroup(0)), population_data[0][0], 1e-9);
}

TEST(TestOdeSECIRVVS, set_population_data_no_data_avail)
{
    auto num_age_groups = 6; // Data to be read requires RKI confirmed cases data age groups
    auto model          = make_model(num_age_groups);
    // set all compartments to zero
    model.populations.array().setConstant(0.0);

    // if the number of vaccinated individuals is greater than the population, we must limit the number of vaccinated.
    model.parameters
        .template get<mio::osecirvvs::DailyPartialVaccinations<double>>()[{mio::AgeGroup(0), mio::SimulationDay(0)}] =
        0;

    model.parameters
        .template get<mio::osecirvvs::DailyFullVaccinations<double>>()[{mio::AgeGroup(0), mio::SimulationDay(0)}] = 0;

    auto model_vector = std::vector<mio::osecirvvs::Model<double>>{model};

    std::string path_pop_data = mio::path_join(TEST_DATA_DIR, "county_current_population.json");
    const std::vector<int> region{0};
    auto population_data = mio::osecirvvs::details::read_population_data(path_pop_data, region).value();

    // we choose the date so that no case data is available
    ASSERT_THAT(mio::osecirvvs::details::set_population_data(
                    model_vector, path_pop_data, mio::path_join(TEST_DATA_DIR, "cases_all_county_age_ma7.json"), {200},
                    {1000, 12, 01}),
                IsSuccess());

    EXPECT_NEAR(model.populations.get_total(), 0.0, 1e-10);
}

TEST(TestOdeSECIRVVS, run_simulation)
{
    auto num_age_groups = 3;
    auto model          = make_model(num_age_groups);
    auto num_days       = 30;

    auto result = mio::osecirvvs::simulate<double>(0, num_days, 0.1, model);
    result = mio::interpolate_simulation_result(result); // Reduce influence of time steps chosen by the integrator.

    // Load result of a previous run; only tests stability, not correctness.
    auto expected_result =
        mio::read_result(mio::path_join(TEST_DATA_DIR, "results_osecirvvs.h5")).value()[0].get_groups();

    ASSERT_THAT(print_wrap(result.matrix()), MatrixNear(print_wrap(expected_result.matrix()), 1e-5, 1e-5));
}

#endif

TEST(TestOdeSECIRVVS, parameter_percentiles)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    //build small graph
    auto model = make_model(5);
    auto graph = mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>>();
    graph.add_node(0, model);

    //sample a few times
    auto sampled_graphs = std::vector<mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>>>();
    std::generate_n(std::back_inserter(sampled_graphs), 10, [&graph]() {
        return mio::osecirvvs::draw_sample(graph, true);
    });

    //extract nodes from graph
    auto sampled_nodes = std::vector<std::vector<mio::osecirvvs::Model<double>>>();
    std::transform(sampled_graphs.begin(), sampled_graphs.end(), std::back_inserter(sampled_nodes), [](auto&& g) {
        auto models = std::vector<mio::osecirvvs::Model<double>>();
        std::transform(g.nodes().begin(), g.nodes().end(), std::back_inserter(models), [](auto&& n) {
            return n.property;
        });
        return models;
    });

    //compute percentiles
    auto percentile_params = mio::osecirvvs::ensemble_params_percentile(sampled_nodes, 0.6)[0].parameters;

    //spot check parameters
    auto p       = double(percentile_params.get<mio::osecirvvs::ReducTimeInfectedMild<double>>()[mio::AgeGroup(2)]);
    auto samples = std::vector<double>();
    std::transform(
        sampled_nodes.begin(), sampled_nodes.end(), std::back_inserter(samples),
        [](const std::vector<mio::osecirvvs::Model<double>>& nodes) {
            return nodes[0].parameters.get<mio::osecirvvs::ReducTimeInfectedMild<double>>()[mio::AgeGroup(2)];
        });

    std::nth_element(samples.begin(), samples.begin() + 6, samples.end());
    ASSERT_THAT(p, samples[6]);

    p       = double(percentile_params.get<mio::osecirvvs::TimeExposed<double>>()[mio::AgeGroup(2)]);
    samples = std::vector<double>();
    std::transform(sampled_nodes.begin(), sampled_nodes.end(), std::back_inserter(samples),
                   [](const std::vector<mio::osecirvvs::Model<double>>& nodes) {
                       return nodes[0].parameters.get<mio::osecirvvs::TimeExposed<double>>()[mio::AgeGroup(2)];
                   });
    std::nth_element(samples.begin(), samples.begin() + 6, samples.end());
    ASSERT_THAT(p, samples[6]);
}

TEST(TestOdeSECIRVVS, get_infections_relative)
{
    auto model = make_model(2);
    auto sim   = mio::osecirvvs::Simulation<ScalarType>(model);
    auto y     = sim.get_result()[0];

    auto relative_infections = mio::osecirvvs::get_infections_relative<ScalarType>(sim, 0.0, y);

    // see model population init to obtain sum 105=2*(7+7.5+8+9.5+10+10.5)
    ASSERT_DOUBLE_EQ(relative_infections, 105 / model.populations.get_total());
}

TEST(TestOdeSECIRVVS, get_mobility_factors)
{
    auto num_age_groups = 2;
    auto model          = make_model(num_age_groups);
    auto sim            = mio::osecirvvs::Simulation<>(model);
    auto y              = sim.get_result()[0];

    auto mobility_factors = mio::osecirvvs::get_mobility_factors<double>(sim, 0.0, y);

    auto expected_values = (Eigen::VectorXd(Eigen::Index(mio::osecirvvs::InfectionState::Count) * num_age_groups) << 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
                               .finished();
    ASSERT_THAT(print_wrap(mobility_factors), MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));
}

TEST(TestOdeSECIRVVS, test_commuters)
{
    auto model                                      = make_model(2);
    auto mobility_factor                            = 0.1;
    auto non_detection_factor                       = 0.3;
    model.parameters.get_start_commuter_detection() = 0.0;
    model.parameters.get_end_commuter_detection()   = 20.0;
    model.parameters.get_commuter_nondetection()    = non_detection_factor;
    auto sim                                        = mio::osecirvvs::Simulation<>(model);
    auto before_testing                             = sim.get_result().get_last_value().eval();
    auto mobile_population                          = (sim.get_result().get_last_value() * mobility_factor).eval();
    auto mobile_population_tested                   = mobile_population.eval();

    mio::osecirvvs::test_commuters<double>(sim, mobile_population_tested, 0.0);

    ASSERT_NEAR(mobile_population_tested[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsNaive)],
                mobile_population[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsNaive)] *
                    non_detection_factor,
                1e-5);
    ASSERT_NEAR(
        sim.get_result().get_last_value()[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed)],
        before_testing[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed)] +
            mobile_population[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsNaive)] *
                (1 - non_detection_factor),
        1e-5);
    ASSERT_NEAR(mobile_population_tested[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity)],
                mobile_population[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity)] *
                    non_detection_factor,
                1e-5);
    ASSERT_NEAR(
        sim.get_result()
            .get_last_value()[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed)],
        before_testing[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed)] +
            mobile_population[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity)] *
                (1 - non_detection_factor),
        1e-5);
    ASSERT_NEAR(
        mobile_population_tested[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity)],
        mobile_population[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity)] *
            non_detection_factor,
        1e-5);
    ASSERT_NEAR(
        sim.get_result()
            .get_last_value()[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed)],
        before_testing[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed)] +
            mobile_population[Eigen::Index(mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity)] *
                (1 - non_detection_factor),
        1e-5);

    ASSERT_NEAR(mobile_population_tested[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive)],
                mobile_population[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive)] *
                    non_detection_factor,
                1e-5);
    ASSERT_NEAR(sim.get_result()
                    .get_last_value()[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed)],
                before_testing[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed)] +
                    mobile_population[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive)] *
                        (1 - non_detection_factor),
                1e-5);
    ASSERT_NEAR(
        mobile_population_tested[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity)],
        mobile_population[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity)] *
            non_detection_factor,
        1e-5);
    ASSERT_NEAR(
        sim.get_result()
            .get_last_value()[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed)],
        before_testing[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed)] +
            mobile_population[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity)] *
                (1 - non_detection_factor),
        1e-5);
    ASSERT_NEAR(
        mobile_population_tested[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity)],
        mobile_population[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity)] *
            non_detection_factor,
        1e-5);
    ASSERT_NEAR(
        sim.get_result().get_last_value()[Eigen::Index(
            mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed)],
        before_testing[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed)] +
            mobile_population[Eigen::Index(mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity)] *
                (1 - non_detection_factor),
        1e-5);
}

// Test model initialization with total population of 0 and ensure get_flows returns no NaN values
TEST(TestOdeSECIRVVS, population_zero_no_nan)
{
    // initialize simple model with total population 0
    mio::osecirvvs::Model<double> model(1);
    model.populations.set_total(0.0);

    // call the get_flows function
    auto dydt_default = Eigen::VectorXd(45);
    dydt_default.setZero();
    auto y0 = model.get_initial_values();
    model.get_flows(y0, y0, 0, dydt_default);

    // check that there are now NaN values in dydt_default
    for (int i = 0; i < dydt_default.size(); i++) {
        EXPECT_FALSE(std::isnan(dydt_default[i]));
    }
}

TEST(TestOdeSECIRVVS, check_constraints_parameters)
{
    auto model = mio::osecirvvs::Model<double>(1);
    EXPECT_EQ(model.parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::osecirvvs::Seasonality<double>>(-0.2);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::Seasonality<double>>(0.2);
    model.parameters.set<mio::osecirvvs::ICUCapacity<double>>(-2);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ICUCapacity<double>>(2);
    model.parameters.set<mio::osecirvvs::TestAndTraceCapacity<double>>(-1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TestAndTraceCapacity<double>>(1);
    model.parameters.set<mio::osecirvvs::TestAndTraceCapacityMaxRiskNoSymptoms<double>>(-1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TestAndTraceCapacityMaxRiskNoSymptoms<double>>(1);
    model.parameters.set<mio::osecirvvs::TestAndTraceCapacityMaxRiskSymptoms<double>>(-1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TestAndTraceCapacityMaxRiskSymptoms<double>>(1);
    model.parameters.set<mio::osecirvvs::TimeExposed<double>>(-2);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeExposed<double>>(2);
    model.parameters.set<mio::osecirvvs::TimeInfectedNoSymptoms<double>>(0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeInfectedNoSymptoms<double>>(5);
    model.parameters.set<mio::osecirvvs::TimeInfectedSymptoms<double>>(0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeInfectedSymptoms<double>>(2);
    model.parameters.set<mio::osecirvvs::TimeInfectedSevere<double>>(-1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeInfectedSevere<double>>(2);
    model.parameters.set<mio::osecirvvs::TimeInfectedCritical<double>>(0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeInfectedCritical<double>>(2);
    model.parameters.set<mio::osecirvvs::TransmissionProbabilityOnContact<double>>(2.0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TransmissionProbabilityOnContact<double>>(0.5);
    model.parameters.set<mio::osecirvvs::RelativeTransmissionNoSymptoms<double>>(-1.0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::RelativeTransmissionNoSymptoms<double>>(0.5);
    model.parameters.set<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>(3.0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>(0.5);
    model.parameters.set<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>(-0.8);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>(0.5);
    model.parameters.set<mio::osecirvvs::SeverePerInfectedSymptoms<double>>(-0.1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::SeverePerInfectedSymptoms<double>>(0.5);
    model.parameters.set<mio::osecirvvs::CriticalPerSevere<double>>(-1.0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::CriticalPerSevere<double>>(0.5);
    model.parameters.set<mio::osecirvvs::DeathsPerCritical<double>>(1.1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::DeathsPerCritical<double>>(0.5);
    model.parameters.set<mio::osecirvvs::VaccinationGap<double>>(0.2);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::VaccinationGap<double>>(2);
    model.parameters.set<mio::osecirvvs::DaysUntilEffectivePartialImmunity<double>>(-2);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::DaysUntilEffectivePartialImmunity<double>>(30);
    model.parameters.set<mio::osecirvvs::DaysUntilEffectiveImprovedImmunity<double>>(-0.2);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::DaysUntilEffectiveImprovedImmunity<double>>(30);
    model.parameters.set<mio::osecirvvs::ReducExposedPartialImmunity<double>>(0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducExposedPartialImmunity<double>>(0.5);
    model.parameters.set<mio::osecirvvs::ReducExposedImprovedImmunity<double>>(-0.2);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducExposedImprovedImmunity<double>>(0.5);
    model.parameters.set<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>(0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>(0.5);
    model.parameters.set<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>(0.);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>(0.5);
    model.parameters.set<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(-4);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(0.5);
    model.parameters.set<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(-4);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(0.5);
    model.parameters.set<mio::osecirvvs::ReducTimeInfectedMild<double>>(-0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducTimeInfectedMild<double>>(1);
    model.parameters.set<mio::osecirvvs::InfectiousnessNewVariant<double>>(-4);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::InfectiousnessNewVariant<double>>(1);
    EXPECT_EQ(model.parameters.check_constraints(), 0);

    model.parameters.set<mio::osecirvvs::InfectiousnessNewVariant<double>>(1);
    model.parameters.set<mio::osecirvvs::DynamicNPIsImplementationDelay<double>>(-4);
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

    model.parameters.set<mio::osecirvvs::Seasonality<double>>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::Seasonality<double>>(), 0);

    model.parameters.set<mio::osecirvvs::ICUCapacity<double>>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::Seasonality<double>>(), 0);

    model.parameters.set<mio::osecirvvs::TestAndTraceCapacity<double>>(-1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TestAndTraceCapacity<double>>(), 0);

    model.parameters.set<mio::osecirvvs::TestAndTraceCapacityMaxRiskNoSymptoms<double>>(-1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TestAndTraceCapacityMaxRiskNoSymptoms<double>>(), 0);

    model.parameters.set<mio::osecirvvs::TestAndTraceCapacityMaxRiskSymptoms<double>>(-1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TestAndTraceCapacityMaxRiskSymptoms<double>>(), 0);

    model.parameters.set<mio::osecirvvs::TimeExposed<double>>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeExposed<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TimeInfectedNoSymptoms<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedNoSymptoms<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TimeInfectedSymptoms<double>>(1e-10);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedSymptoms<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TimeInfectedSevere<double>>(-1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedSevere<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TimeInfectedCritical<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedCritical<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecirvvs::TransmissionProbabilityOnContact<double>>(2.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR(model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>()[indx_agegroup], 0.0,
                1e-14);

    model.parameters.set<mio::osecirvvs::RelativeTransmissionNoSymptoms<double>>(-1.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::RelativeTransmissionNoSymptoms<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>(3.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>(-0.8);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::SeverePerInfectedSymptoms<double>>(-0.1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::SeverePerInfectedSymptoms<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::CriticalPerSevere<double>>(-1.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::CriticalPerSevere<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::DeathsPerCritical<double>>(1.1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::DeathsPerCritical<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::VaccinationGap<double>>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::VaccinationGap<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::DaysUntilEffectivePartialImmunity<double>>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::DaysUntilEffectivePartialImmunity<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::DaysUntilEffectiveImprovedImmunity<double>>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::DaysUntilEffectiveImprovedImmunity<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::ReducExposedPartialImmunity<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::ReducExposedPartialImmunity<double>>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirvvs::ReducExposedImprovedImmunity<double>>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::ReducExposedImprovedImmunity<double>>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>(0.);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(-4);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(
        model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[indx_agegroup],
        1);

    model.parameters.set<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(-4);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(
        model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[indx_agegroup],
        1);

    model.parameters.set<mio::osecirvvs::ReducTimeInfectedMild<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::ReducTimeInfectedMild<double>>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirvvs::InfectiousnessNewVariant<double>>(-4);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::InfectiousnessNewVariant<double>>()[indx_agegroup], 1);

    model.parameters.set<mio::osecirvvs::DynamicNPIsImplementationDelay<double>>(-4);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::DynamicNPIsImplementationDelay<double>>(), 0);

    EXPECT_EQ(model.parameters.apply_constraints(), 0);
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSECIRVVS, apply_variant_function)
{
    auto model = mio::osecirvvs::Model(1);
    model.parameters.set<mio::osecirvvs::TransmissionProbabilityOnContact<double>>(0.2);

    model.parameters.set<mio::osecirvvs::StartDay>(0);
    model.parameters.set<mio::osecirvvs::StartDayNewVariant>(10);
    model.parameters.set<mio::osecirvvs::InfectiousnessNewVariant<double>>(2.0);
    auto sim = mio::osecirvvs::Simulation<>(model);

    // test that the transmission probability is not changed due to calling the advance function
    sim.advance(0.01);
    EXPECT_NEAR(
        sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)],
        0.2, 1e-10);

    // test if the transmission probability is set to the correct value after applying the variant.
    // The referece values are calculated using equation (36) in doi.org/10.1371/journal.pcbi.1010054
    auto base_infectiousness =
        sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>();

    // however, the parameter should stay unchanged if the new variant is not present in the population.
    sim.apply_variant(0, base_infectiousness);
    EXPECT_NEAR(
        sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)],
        0.2, 1e-10);

    sim.apply_variant(9, base_infectiousness);
    EXPECT_NEAR(
        sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)],
        0.2, 1e-10);

    sim.apply_variant(10, base_infectiousness);
    EXPECT_NEAR(
        sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)],
        0.99 * base_infectiousness[mio::AgeGroup(0)] +
            0.01 * base_infectiousness[mio::AgeGroup(0)] *
                sim.get_model().parameters.get<mio::osecirvvs::InfectiousnessNewVariant<double>>()[mio::AgeGroup(0)],
        1e-10);

    sim.apply_variant(45, base_infectiousness);
    EXPECT_NEAR(
        sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)],
        0.68 * base_infectiousness[mio::AgeGroup(0)] +
            0.32 * base_infectiousness[mio::AgeGroup(0)] *
                sim.get_model().parameters.get<mio::osecirvvs::InfectiousnessNewVariant<double>>()[mio::AgeGroup(0)],
        1e-10);

    sim.apply_variant(1000, base_infectiousness);
    EXPECT_NEAR(
        sim.get_model().parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)],
        0.4, 1e-10);
}
