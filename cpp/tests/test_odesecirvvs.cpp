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
#include "memilio/mobility/mobility.h"
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

TEST(TestSecir, reduceToSecirAndCompareWithPreviousRun)
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
        EXPECT_GE(result.get_last_value()[i], -1e-3);
        num_persons += result.get_last_value()[i];
    }
    EXPECT_NEAR(num_persons, model.populations.get_total(), 1e-10);
}

#if defined(MEMILIO_HAS_HDF5) && defined(MEMILIO_HAS_JSONCPP)

TEST(TestOdeSECIRVVS, read_data)
{
    auto num_age_groups = 6; //reading data requires RKI data age groups
    auto model          = std::vector<mio::osecirvvs::Model>({make_model(num_age_groups)});

    auto read_result = mio::osecirvvs::read_input_data_county(
        model, {2020, 12, 01}, {1002}, std::vector<double>(size_t(num_age_groups), 1.0), 1.0, TEST_DATA_DIR, 10);
    ASSERT_THAT(read_result, IsSuccess());

    // values were generated by the tested function; can only check stability of the implementation, not correctness
    auto expected_values =
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecirvvs::InfectionState::Count)) << 8795.83, 499.945,
         3.16404, 0.177042, 0.16439, 1.40624, 0.0983565, 0.220732, 0, 0, 0, 6.84143, 0.203161, 0.154011, 0, 0, 0,
         0.0337402, 3.50295e-05, 0.000157341, 0.0292738, 0.000166359, 0.000306977, 1012.73, 0, 0, 0, 17074.9, 482.931,
         7.58957, 0.211334, 0.32517, 4.87901, 0.169822, 0.631534, 0, 0, 0, 11.6593, 0.194153, 0.24389, 0, 0, 0,
         0.0434302, 1.84349e-05, 0.000137211, 0.0294097, 8.31719e-05, 0.000254317, 1621.15, 0, 0, 0, 72497.9, 314.797,
         47.1867, 0.201638, 0.446808, 24.3886, 0.130272, 0.697687, 0, 0, 0, 64.6307, 0.152734, 0.276309, 0, 0, 0,
         1.23753, 8.20856e-05, 0.000879878, 0.0751831, 3.26292e-05, 0.000143686, 1534.92, 0, 0, 0, 80793.3, 375.855,
         44.4994, 0.203793, 0.456367, 23.6624, 0.135457, 0.733146, 0, 0, 0, 58.9247, 0.154166, 0.281854, 0, 0, 0,
         3.75954, 0.000277175, 0.00300253, 0.48662, 0.000226338, 0.00100726, 1828.25, 0, 0, 0, 41583.7, 508.624, 9.2403,
         0.111291, 0.181928, 5.7534, 0.0866177, 0.342224, 0, 0, 0, 13.3241, 0.089212, 0.119063, 0, 0, 0, 3.11879,
         0.000600158, 0.00474586, 2.60421, 0.00318554, 0.0103486, 1802.81, 0, 0, 0, 13230.4, 412.958, 11.1133, 0.341707,
         0.781035, 3.48326, 0.133878, 0.739583, 0, 0, 0, 3.51008, 0.0630277, 0.117614, 0, 0, 0, 0.695112, 0.000364314,
         0.00402809, 4.67857, 0.0146103, 0.0663642, 2050.05, 0, 0, 0)
            .finished();

    ASSERT_THAT(print_wrap(model[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));

    // some more tests which are actually not necessary but can help if the previous tests fails or needs to get changed
    for (mio::AgeGroup i = 0; i < model[0].parameters.get_num_groups(); i++) {
        EXPECT_GE(double(model[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}]), 0);
        EXPECT_GE(double(model[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}]), 0);
        EXPECT_GE(double(model[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(double(model[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}]),
                  0);
        EXPECT_GE(double(model[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}]),
                  0);
        EXPECT_GE(double(model[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}]),
                  0);

        // currently dead and confirmed after commuting compartments are initialized as zero
        EXPECT_EQ(double(model[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}]),
                  0);
        EXPECT_EQ(double(model[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}]), 0);
        EXPECT_EQ(
            double(
                model[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(model[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model[0].populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(
            double(
                model[0].populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]),
            0);
        EXPECT_EQ(double(model[0].populations[{i, mio::osecirvvs::InfectionState::DeadNaive}]), 0);
        EXPECT_EQ(double(model[0].populations[{i, mio::osecirvvs::InfectionState::DeadPartialImmunity}]), 0);
        EXPECT_EQ(double(model[0].populations[{i, mio::osecirvvs::InfectionState::DeadImprovedImmunity}]), 0);
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
                    std::vector<mio::osecirvvs::Model>{model}, TEST_DATA_DIR, tmp_results_dir, {0}, {2020, 12, 01},
                    std::vector<double>(size_t(num_age_groups), 1.0), 1.0, 2, true),
                IsSuccess());

    auto data_extrapolated = mio::read_result(mio::path_join(tmp_results_dir, "Results_rki.h5"));
    ASSERT_THAT(data_extrapolated, IsSuccess());

    // Values were generated by the tested function export_input_data_county_timeseries;
    // can only check stability of the implementation, not correctness
    auto expected_results =
        mio::read_result(mio::path_join(TEST_DATA_DIR, "export_time_series_init_osecirvvs.h5")).value();

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
                                                       TEST_DATA_DIR, 2, true),
                IsSuccess());

    // Values from data/export_time_series_init_osecirvvs.h5, for reading in comparison
    // operator for return of mio::read_result and model population needed.
    auto expected_values =
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecirvvs::InfectionState::Count)) << 3.46722e+06, 500.422,
         3.42784, 0.000487127, 0.000452317, 1.52349, 0.000270626, 0.000607341, 0, 0, 0, 7.28451, 0.000549392,
         0.000416479, 0, 0, 0, 0.0342843, 9.04004e-08, 4.06048e-07, 0.0297459, 4.29321e-07, 7.92212e-07, 1012.64, 0, 0,
         0, 7.74735e+06, 483.505, 7.99905, 0.000491533, 0.000756298, 5.14224, 0.000394982, 0.00146886, 0, 0, 0, 12.1417,
         0.000446185, 0.000560485, 0, 0, 0, 0.0439275, 4.11478e-08, 3.06263e-07, 0.0297464, 1.85645e-07, 5.6765e-07,
         1621.1, 0, 0, 0, 1.92155e+07, 315.279, 47.6767, 0.000770223, 0.00170673, 24.6419, 0.000497616, 0.00266505, 0,
         0, 0, 65.1409, 0.000581982, 0.00105285, 0, 0, 0, 1.24042, 3.11056e-07, 3.33421e-06, 0.0753588, 1.23645e-07,
         5.44483e-07, 1534.84, 0, 0, 0, 3.00317e+07, 376.348, 44.9986, 0.000555233, 0.00124337, 23.9278, 0.000369054,
         0.00199746, 0, 0, 0, 59.4272, 0.000418907, 0.000765866, 0, 0, 0, 3.76905, 7.48674e-07, 8.11008e-06, 0.48785,
         6.11358e-07, 2.7207e-06, 1828.18, 0, 0, 0, 1.65123e+07, 508.914, 9.46371, 0.000287187, 0.000469468, 5.8925,
         0.000223518, 0.000883113, 0, 0, 0, 13.5708, 0.00022894, 0.000305543, 0, 0, 0, 3.13496, 1.51999e-06,
         1.20196e-05, 2.61771, 8.06787e-06, 2.62095e-05, 1802.76, 0, 0, 0, 6.17984e+06, 413.51, 11.9623, 0.000788119,
         0.00180139, 3.74938, 0.000308778, 0.00170579, 0, 0, 0, 3.71382, 0.000142889, 0.00026664, 0, 0, 0, 0.707117,
         7.94101e-07, 8.78009e-06, 4.75937, 3.18463e-05, 0.000144655, 2050.02, 0, 0, 0)
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
                            1, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
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
}