/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
    model.parameters.get<mio::osecirvvs::DailyFirstVaccination>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirvvs::DailyFirstVaccination>().array().setConstant(0);
    model.parameters.get<mio::osecirvvs::DailyFullVaccination>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirvvs::DailyFullVaccination>().array().setConstant(0);
    mio::TimeSeries<double> result = simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

TEST(TestOdeSECIRVVS, reduceToSecirAndCompareWithPreviousRun)
{
    // double t0   = 0;
    // double tmax = 50;

    mio::osecirvvs::Model model(1);

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

    model.parameters.get<mio::osecirvvs::ICUCapacity>()          = 10000;
    model.parameters.get<mio::osecirvvs::TestAndTraceCapacity>() = 10000;
    model.parameters.get<mio::osecirvvs::DailyFirstVaccination>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirvvs::DailyFirstVaccination>().array().setConstant(0);
    model.parameters.get<mio::osecirvvs::DailyFullVaccination>().resize(mio::SimulationDay(size_t(1000)));
    model.parameters.get<mio::osecirvvs::DailyFullVaccination>().array().setConstant(0);

    auto& contacts       = model.parameters.get<mio::osecirvvs::ContactPatterns>();
    auto& contact_matrix = contacts.get_cont_freq_mat();
    contact_matrix[0]    = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    //times
    model.parameters.get<mio::osecirvvs::IncubationTime>()[mio::AgeGroup(0)]       = 5.2;
    model.parameters.get<mio::osecirvvs::SerialInterval>()[mio::AgeGroup(0)]       = 4.2;
    model.parameters.get<mio::osecirvvs::TimeInfectedSymptoms>()[mio::AgeGroup(0)] = 5;
    model.parameters.get<mio::osecirvvs::TimeInfectedSevere>()[mio::AgeGroup(0)]   = 10;
    model.parameters.get<mio::osecirvvs::TimeInfectedCritical>()[mio::AgeGroup(0)] = 8;

    //probabilities
    model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>()[mio::AgeGroup(0)] = 0.05;
    model.parameters.get<mio::osecirvvs::RelativeTransmissionNoSymptoms>()[mio::AgeGroup(0)]   = 1;
    model.parameters.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic>()[mio::AgeGroup(0)]   = 0.25;
    model.parameters.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms>()[mio::AgeGroup(0)]   = 0.09;
    model.parameters.get<mio::osecirvvs::SeverePerInfectedSymptoms>()[mio::AgeGroup(0)]        = 0.2;
    model.parameters.get<mio::osecirvvs::CriticalPerSevere>()[mio::AgeGroup(0)]                = 0.25;
    model.parameters.get<mio::osecirvvs::DeathsPerCritical>()[mio::AgeGroup(0)]                = 0.3;

    // TODO: Reduction not possible like this, division by zero!
    model.parameters.get<mio::osecirvvs::ReducExposedPartialImmunity>()[mio::AgeGroup(0)]                     = 0;
    model.parameters.get<mio::osecirvvs::ReducExposedImprovedImmunity>()[mio::AgeGroup(0)]                    = 0;
    model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity>()[mio::AgeGroup(0)]            = 0;
    model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity>()[mio::AgeGroup(0)]           = 0;
    model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity>()[mio::AgeGroup(0)]  = 0;
    model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity>()[mio::AgeGroup(0)] = 0;
    model.parameters.get<mio::osecirvvs::ReducTimeInfectedMild>()[mio::AgeGroup(0)]                           = 1;

    model.parameters.get<mio::osecirvvs::Seasonality>() = 0.2;

    model.apply_constraints();

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
    parameters.get<mio::osecirvvs::DailyFirstVaccination>().resize(mio::SimulationDay(size_t(1000)));
    parameters.get<mio::osecirvvs::DailyFirstVaccination>().array().setConstant(5);
    parameters.get<mio::osecirvvs::DailyFullVaccination>().resize(mio::SimulationDay(size_t(1000)));
    parameters.get<mio::osecirvvs::DailyFullVaccination>().array().setConstant(3);
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
    const double incubationTime            = 5.2;
    const double serialIntervalMin         = 0.5 * 2.67 + 0.5 * 5.2;
    const double serialIntervalMax         = 0.5 * 4.00 + 0.5 * 5.2;
    const double timeInfectedSymptomsMin[] = {5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465};
    const double timeInfectedSymptomsMax[] = {8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085};
    const double timeInfectedSevereMin[]   = {3.925, 3.925, 4.85, 6.4, 7.2, 9.};
    const double timeInfectedSevereMax[]   = {6.075, 6.075, 7., 8.7, 9.8, 13.};
    const double timeInfectedCriticalMin[] = {4.95, 4.95, 4.86, 14.14, 14.4, 10.};
    const double timeInfectedCriticalMax[] = {8.95, 8.95, 8.86, 20.58, 19.8, 13.2};

    array_assign_uniform_distribution(params.get<mio::osecirvvs::IncubationTime>(), incubationTime, incubationTime,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::SerialInterval>(), serialIntervalMin,
                                      serialIntervalMax, set_invalid_initial_value);
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

TEST(TestOdeSECIRVVS, draw_sample)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters> graph;

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
    ASSERT_TRUE((parameters0.get<mio::osecirvvs::BaseInfectiousnessB161>().array(),
                 parameters0.get<mio::osecirvvs::TransmissionProbabilityOnContact>().array() * 1.6) //using high variant
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

TEST(TestOdeSECIRVVS, checkPopulationConservation)
{
    auto num_age_groups = 2;
    auto model          = make_model(num_age_groups);
    auto num_days       = 30;

    auto result = mio::osecirvvs::simulate(0, num_days, 0.1, model);

    double num_persons = 0.0;
    for (auto i = 0; i < result.get_last_value().size(); i++) {
        // do not sum up auxiliary compartment
        if ((i + 1) % ((int)mio::osecirvvs::InfectionState::TotalInfections + 1) != 0) {
            EXPECT_GE(result.get_last_value()[i], -1e-3);
            num_persons += result.get_last_value()[i];
        }
    }
    EXPECT_NEAR(num_persons, model.populations.get_total(), 1e-10);
}

#if defined(MEMILIO_HAS_HDF5) && defined(MEMILIO_HAS_JSONCPP)

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
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecirvvs::InfectionState::Count)) << 8792.15, 175.889,
         3.21484, 0.0633116, 0.221057, 1.42882, 0.0351731, 0.29682, 0, 0, 0, 6.93838, 0.0725173, 0.206715, 0, 0, 0,
         0.0337498, 1.23324e-05, 0.000208293, 0.0292822, 5.8568e-05, 0.000406386, 1340.42, 0, 0, 0, 0, 17067.6, 220.137,
         7.64078, 0.0970237, 0.381933, 4.91193, 0.0779655, 0.741778, 0, 0, 0, 11.7286, 0.0890643, 0.286235, 0, 0, 0,
         0.0434344, 8.40756e-06, 0.000160098, 0.0294125, 3.7932e-05, 0.000296738, 1891.19, 0, 0, 0, 0, 72490.2, 176.267,
         47.227, 0.113013, 0.490073, 24.4094, 0.0730141, 0.765246, 0, 0, 0, 64.6789, 0.0855947, 0.303032, 0, 0, 0,
         1.23754, 4.5968e-05, 0.000964262, 0.0751837, 1.82724e-05, 0.000157466, 1681.05, 0, 0, 0, 0, 80790.1, 184.645,
         44.5477, 0.100229, 0.50512, 23.6881, 0.0666206, 0.811467, 0, 0, 0, 58.9805, 0.0758111, 0.31192, 0, 0, 0,
         3.75961, 0.000136175, 0.00331973, 0.486628, 0.000111199, 0.00111367, 2022.58, 0, 0, 0, 0, 41581, 177.478,
         9.27393, 0.0389771, 0.216151, 5.77433, 0.030336, 0.4066, 0, 0, 0, 13.3664, 0.0312302, 0.141394, 0, 0, 0, 3.119,
         0.000209444, 0.00561852, 2.60439, 0.00111169, 0.0122515, 2136.6, 0, 0, 0, 0, 13223.8, 216.037, 11.1838,
         0.179986, 0.863926, 3.50537, 0.0705169, 0.818075, 0, 0, 0, 3.52982, 0.0331744, 0.130002, 0, 0, 0, 0.695168,
         0.000190699, 0.00442784, 4.67895, 0.00764769, 0.0729502, 2253.61, 0, 0, 0, 0)
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
                    std::vector<mio::osecirvvs::Model>{model}, tmp_results_dir, {0}, {2020, 12, 01},
                    std::vector<double>(size_t(num_age_groups), 1.0), 1.0, 2,
                    mio::path_join(TEST_DATA_DIR, "county_divi_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "cases_all_county_age_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "county_current_population.json"), true,
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
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecirvvs::InfectionState::Count)) << 3.46722e+06, 176.06,
         3.42799, 0.00017139, 0.00059842, 1.52355, 9.52168e-05, 0.000803518, 0, 0, 0, 7.28479, 0.000193296, 0.000551002,
         0, 0, 0, 0.0342843, 3.1805e-08, 5.37183e-07, 0.029746, 1.51045e-07, 1.04806e-06, 1340.35, 0, 0, 0, 0,
         7.74734e+06, 220.4, 7.99917, 0.000224064, 0.000882024, 5.14232, 0.000180051, 0.00171304, 0, 0, 0, 12.1419,
         0.000203391, 0.000653659, 0, 0, 0, 0.0439275, 1.87568e-08, 3.5717e-07, 0.0297464, 8.46241e-08, 6.62006e-07,
         1891.15, 0, 0, 0, 0, 1.92155e+07, 176.538, 47.6768, 0.00043128, 0.00187022, 24.642, 0.000278636, 0.00292033, 0,
         0, 0, 65.1411, 0.000325876, 0.0011537, 0, 0, 0, 1.24042, 1.74173e-07, 3.65358e-06, 0.0753588, 6.9234e-08,
         5.96637e-07, 1680.98, 0, 0, 0, 0, 3.00317e+07, 184.888, 44.9988, 0.000272769, 0.00137466, 23.9279, 0.000181305,
         0.00220837, 0, 0, 0, 59.4274, 0.000205796, 0.000846734, 0, 0, 0, 3.76905, 3.67799e-07, 8.9664e-06, 0.48785,
         3.00341e-07, 3.00797e-06, 2022.51, 0, 0, 0, 0, 1.65123e+07, 177.579, 9.4638, 0.000100211, 0.00055573, 5.89255,
         7.79946e-05, 0.00104538, 0, 0, 0, 13.5709, 7.98864e-05, 0.000361685, 0, 0, 0, 3.13496, 5.30384e-07, 1.4228e-05,
         2.61772, 2.81518e-06, 3.1025e-05, 2136.56, 0, 0, 0, 0, 6.17983e+06, 216.328, 11.9625, 0.000412312, 0.00197908,
         3.74944, 0.00016154, 0.00187405, 0, 0, 0, 3.71387, 7.47535e-05, 0.000292941, 0, 0, 0, 0.707117, 4.15435e-07,
         9.64602e-06, 4.75937, 1.66604e-05, 0.000158922, 2253.59, 0, 0, 0, 0)
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
        mio::read_result(mio::path_join(TEST_DATA_DIR, "results_osecirvvs.h5")).value()[0].get_groups();

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
        return mio::osecirvvs::draw_sample(graph, true);
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

    p       = double(percentile_params.get<mio::osecirvvs::SerialInterval>()[mio::AgeGroup(2)]);
    samples = std::vector<double>();
    std::transform(sampled_nodes.begin(), sampled_nodes.end(), std::back_inserter(samples),
                   [](const std::vector<mio::osecirvvs::Model>& nodes) {
                       return nodes[0].parameters.get<mio::osecirvvs::SerialInterval>()[mio::AgeGroup(2)];
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
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
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
    model.parameters.set<mio::osecirvvs::IncubationTime>(-2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::IncubationTime>(2);
    model.parameters.set<mio::osecirvvs::SerialInterval>(1);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::SerialInterval>(5);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::SerialInterval>(1.5);
    model.parameters.set<mio::osecirvvs::TimeInfectedSymptoms>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeInfectedSymptoms>(2);
    model.parameters.set<mio::osecirvvs::TimeInfectedSevere>(-1);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeInfectedSevere>(2);
    model.parameters.set<mio::osecirvvs::TimeInfectedCritical>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::TimeInfectedCritical>(2);
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
    model.parameters.set<mio::osecirvvs::DaysUntilEffectivePartialImmunity>(-2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::DaysUntilEffectivePartialImmunity>(30);
    model.parameters.set<mio::osecirvvs::DaysUntilEffectiveImprovedImmunity>(-0.2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::DaysUntilEffectiveImprovedImmunity>(30);
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
    model.parameters.set<mio::osecirvvs::ReducTimeInfectedMild>(-0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::ReducTimeInfectedMild>(1);
    model.parameters.set<mio::osecirvvs::BaseInfectiousnessB117>(-0.5);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecirvvs::BaseInfectiousnessB117>(0.5);
    model.parameters.set<mio::osecirvvs::BaseInfectiousnessB161>(-4);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSECIRVVS, apply_constraints_parameters)
{
    auto model         = mio::osecirvvs::Model(1);
    auto indx_agegroup = mio::AgeGroup(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    model.parameters.set<mio::osecirvvs::Seasonality>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::Seasonality>(), 0);

    model.parameters.set<mio::osecirvvs::ICUCapacity>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::Seasonality>(), 0);

    model.parameters.set<mio::osecirvvs::IncubationTime>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::IncubationTime>()[indx_agegroup], 2e-4);

    model.parameters.set<mio::osecirvvs::SerialInterval>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR(model.parameters.get<mio::osecirvvs::SerialInterval>()[indx_agegroup], 0.00015, 1e-14);

    model.parameters.set<mio::osecirvvs::SerialInterval>(5);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR(model.parameters.get<mio::osecirvvs::SerialInterval>()[indx_agegroup], 0.00015, 1e-14);

    model.parameters.set<mio::osecirvvs::TimeInfectedSymptoms>(1e-5);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedSymptoms>()[indx_agegroup], 1e-4);

    model.parameters.set<mio::osecirvvs::TimeInfectedSevere>(-1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedSevere>()[indx_agegroup], 1e-4);

    model.parameters.set<mio::osecirvvs::TimeInfectedCritical>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::TimeInfectedCritical>()[indx_agegroup], 1e-4);

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

    model.parameters.set<mio::osecirvvs::DaysUntilEffectivePartialImmunity>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::DaysUntilEffectivePartialImmunity>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::DaysUntilEffectiveImprovedImmunity>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::DaysUntilEffectiveImprovedImmunity>()[indx_agegroup], 0);

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

    model.parameters.set<mio::osecirvvs::BaseInfectiousnessB117>(-0.5);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::BaseInfectiousnessB117>()[indx_agegroup], 0);

    model.parameters.set<mio::osecirvvs::BaseInfectiousnessB161>(-4);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecirvvs::BaseInfectiousnessB161>()[indx_agegroup], 0);

    mio::set_log_level(mio::LogLevel::warn);
}
