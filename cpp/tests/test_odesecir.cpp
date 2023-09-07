/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#include "load_test_data.h"
#include "memilio/compartments/simulation.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/time_series.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/model.h"
#include "memilio/math/adapt_rk.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/analyze_result.h"
#include "ode_secir/parameters.h"
#include <distributions_helpers.h>
#include <gtest/gtest.h>
#include <iomanip>

#include <fstream>
#include <sstream>
#include <string>

TEST(TestOdeSecir, compareWithPreviousRun)
{
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double cont_freq = 10;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model model(1);

    model.parameters.set<mio::osecir::StartDay>(60);
    model.parameters.set<mio::osecir::Seasonality>(0.2);

    model.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]       = 5.2;
    model.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]       = 4.2;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = 5.8;
    model.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]   = 9.5;
    model.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0] = 7.1;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    model.populations.set_total(nb_total_t0);
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = nb_exp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = nb_car_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]   = nb_inf_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]     = nb_hosp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]   = nb_icu_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]          = nb_rec_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]               = nb_dead_t0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                nb_total_t0);

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0]  = 0.05;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]    = 0.7;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]    = 0.09;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]    = 0.25;
    model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0] = 0.45;
    model.parameters.get<mio::osecir::TestAndTraceCapacity>()                                = 35;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]         = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]                 = 0.25;
    model.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]                 = 0.3;

    model.apply_constraints();

    auto integrator = std::make_shared<mio::RKIntegratorCore>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    mio::TimeSeries<double> secihurd = simulate(t0, tmax, dt, model, integrator);

    auto compare = load_test_data_csv<double>("secihurd-compare.csv");

    ASSERT_EQ(compare.size(), static_cast<size_t>(secihurd.get_num_time_points()));
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size(), static_cast<size_t>(secihurd.get_num_elements()) + 1) << "at row " << i;
        EXPECT_NEAR(secihurd.get_time(i), compare[i][0], 1e-10) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            EXPECT_NEAR(secihurd.get_value(i)[j - 1], compare[i][j], 1e-10) << " at row " << i;
        }
    }
}

TEST(TestOdeSecir, simulateDefault)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.1;

    mio::osecir::Model model(1);
    mio::TimeSeries<double> result = simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

TEST(TestOdeSecir, checkPopulationConservation)
{
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double cont_freq = 10;

    double nb_total_t0 = 10000;

    mio::osecir::Model model(1);

    model.parameters.get<mio::osecir::TestAndTraceCapacity>() = 35;

    model.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]       = 5.2;
    model.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]       = 4.2;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = 5.8;
    model.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]   = 9.5;
    model.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0] = 7.1;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    model.populations.set_total(nb_total_t0);
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]   = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]     = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]   = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]          = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]               = 10;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                nb_total_t0);

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0]  = 0.05;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]    = 1;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]    = 0.09;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]    = 0.25;
    model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0] = 0.45;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]         = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]                 = 0.25;
    model.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]                 = 0.3;

    model.apply_constraints();

    mio::TimeSeries<double> secir = simulate(t0, tmax, dt, model);

    double num_persons = 0.0;
    for (auto i = 0; i < secir.get_last_value().size(); i++) {
        num_persons += secir.get_last_value()[i];
    }
    EXPECT_NEAR(num_persons, nb_total_t0, 1e-10);
}

TEST(TestOdeSecir, testParamConstructors)
{

    double cont_freq = 10;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 54, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 11, nb_dead_t0 = 0;

    double icu_cap   = 4444;
    double start_day = 30, seasonality = 0.3;

    mio::osecir::Model model(1);

    model.parameters.set<mio::osecir::ICUCapacity>(icu_cap);

    model.parameters.set<mio::osecir::StartDay>(start_day);
    model.parameters.set<mio::osecir::Seasonality>(seasonality);

    model.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]       = 5.2;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = 5;
    model.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]       = 4.2;
    model.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]   = 10.;
    model.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0] = 8.;

    model.populations.set_total(nb_total_t0);
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = nb_exp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = nb_car_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]   = nb_inf_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]     = nb_hosp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]   = nb_icu_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]          = nb_rec_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]               = nb_dead_t0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                nb_total_t0);

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0] = 0.05;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]   = 0.67;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]   = 0.09;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = 0.25;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]        = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]                = 0.24;
    model.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]                = 0.3;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    mio::osecir::Model model2{model}; // copy constructor

    EXPECT_EQ(model.parameters.get<mio::osecir::ICUCapacity>(), model2.parameters.get<mio::osecir::ICUCapacity>());
    EXPECT_EQ(model.parameters.get<mio::osecir::StartDay>(), model2.parameters.get<mio::osecir::StartDay>());
    EXPECT_EQ(model.parameters.get<mio::osecir::Seasonality>(), model2.parameters.get<mio::osecir::Seasonality>());

    EXPECT_EQ(model.populations.get_total(), model2.populations.get_total());
    EXPECT_EQ((model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}]),
              (model2.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}]));
    EXPECT_EQ((model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]),
              (model2.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]));
    EXPECT_EQ((model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]),
              (model2.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]));
    EXPECT_EQ((model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]),
              (model2.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]));
    EXPECT_EQ((model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]),
              (model2.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]));
    EXPECT_EQ((model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]),
              (model2.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]));
    EXPECT_EQ((model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]),
              (model2.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]));
    EXPECT_EQ((model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]),
              (model2.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]));

    EXPECT_EQ(model.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::ContactPatterns>().get_cont_freq_mat(),
              model2.parameters.get<mio::osecir::ContactPatterns>().get_cont_freq_mat());

    EXPECT_EQ(model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::ContactPatterns>().get_cont_freq_mat(),
              model2.parameters.get<mio::osecir::ContactPatterns>().get_cont_freq_mat());

    mio::osecir::Model model3 = std::move(model2); // move constructor

    EXPECT_EQ(model.parameters.get<mio::osecir::ICUCapacity>(), model3.parameters.get<mio::osecir::ICUCapacity>());
    EXPECT_EQ(model.parameters.get<mio::osecir::StartDay>(), model3.parameters.get<mio::osecir::StartDay>());
    EXPECT_EQ(model.parameters.get<mio::osecir::Seasonality>(), model3.parameters.get<mio::osecir::Seasonality>());

    EXPECT_EQ(model3.populations.get_total(), model.populations.get_total());
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}]),
              (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]),
              (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]),
              (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]),
              (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]),
              (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]),
              (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]),
              (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]),
              (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]));

    EXPECT_EQ(model3.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0]);

    EXPECT_EQ(model3.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]);

    EXPECT_EQ(model.parameters.get<mio::osecir::ContactPatterns>().get_cont_freq_mat(),
              model3.parameters.get<mio::osecir::ContactPatterns>().get_cont_freq_mat());

    mio::osecir::Model model4 = model3; // copy assignment constructor

    EXPECT_EQ(model4.parameters.get<mio::osecir::ICUCapacity>(), model3.parameters.get<mio::osecir::ICUCapacity>());
    EXPECT_EQ(model4.parameters.get<mio::osecir::StartDay>(), model3.parameters.get<mio::osecir::StartDay>());
    EXPECT_EQ(model4.parameters.get<mio::osecir::Seasonality>(), model3.parameters.get<mio::osecir::Seasonality>());

    EXPECT_EQ(model3.populations.get_total(), model4.populations.get_total());
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}]),
              (model4.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]),
              (model4.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]),
              (model4.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]),
              (model4.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]),
              (model4.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]),
              (model4.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]),
              (model4.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]));
    EXPECT_EQ((model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]),
              (model4.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]));

    EXPECT_EQ(model3.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0]);

    EXPECT_EQ(model3.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]);

    EXPECT_EQ(model4.parameters.get<mio::osecir::ContactPatterns>().get_cont_freq_mat(),
              model3.parameters.get<mio::osecir::ContactPatterns>().get_cont_freq_mat());

    mio::osecir::Model model5 = std::move(model4); // move assignment constructor

    EXPECT_EQ(model5.parameters.get<mio::osecir::ICUCapacity>(), model3.parameters.get<mio::osecir::ICUCapacity>());
    EXPECT_EQ(model5.parameters.get<mio::osecir::StartDay>(), model3.parameters.get<mio::osecir::StartDay>());
    EXPECT_EQ(model5.parameters.get<mio::osecir::Seasonality>(), model3.parameters.get<mio::osecir::Seasonality>());

    EXPECT_EQ(model5.populations.get_total(), model3.populations.get_total());
    EXPECT_EQ((model5.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}]),
              (model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}]));
    EXPECT_EQ((model5.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]),
              (model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]));
    EXPECT_EQ((model5.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]),
              (model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]));
    EXPECT_EQ((model5.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]),
              (model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]));
    EXPECT_EQ((model5.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]),
              (model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]));
    EXPECT_EQ((model5.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]),
              (model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]));
    EXPECT_EQ((model5.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]),
              (model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]));
    EXPECT_EQ((model5.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]),
              (model3.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]));

    EXPECT_EQ(model5.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0]);

    EXPECT_EQ(model5.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]);

    EXPECT_EQ(model5.parameters.get<mio::osecir::ContactPatterns>().get_cont_freq_mat(),
              model3.parameters.get<mio::osecir::ContactPatterns>().get_cont_freq_mat());
}

TEST(TestOdeSecir, testSettersAndGetters)
{
    std::vector<mio::UncertainValue> vec;

    for (int i = 0; i < 22; i++) {
        mio::UncertainValue val = mio::UncertainValue(i);
        val.set_distribution(mio::ParameterDistributionNormal(i, 10 * i, 5 * i, i / 10.0));
        vec.push_back(val);
    }

    mio::osecir::Model model(1);

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    EXPECT_EQ(model.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0].get_distribution().get(), nullptr);

    model.parameters.set<mio::osecir::ICUCapacity>(vec[0]);

    model.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]       = vec[1];
    model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = vec[2];
    model.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]       = vec[3];
    model.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]   = vec[4];
    model.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0] = vec[5];

    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = vec[6];
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = vec[7];
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]   = vec[8];
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]     = vec[9];
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]   = vec[10];
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]          = vec[11];
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]               = vec[12];

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0] = vec[13];
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]   = vec[14];
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]   = vec[15];
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = vec[16];
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]        = vec[17];
    model.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]                = vec[18];
    model.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]                = vec[19];

    EXPECT_NE(model.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0].get_distribution().get(), nullptr);

    check_distribution(*vec[0].get_distribution(),
                       *model.parameters.get<mio::osecir::ICUCapacity>().get_distribution());

    model.parameters.set<mio::osecir::StartDay>(vec[20]);
    model.parameters.set<mio::osecir::Seasonality>(vec[21]);

    EXPECT_NE(model.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0].get_distribution().get(), nullptr);

    check_distribution(*vec[1].get_distribution(),
                       *model.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(*vec[2].get_distribution(),
                       *model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(*vec[3].get_distribution(),
                       *model.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(*vec[4].get_distribution(),
                       *model.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(*vec[5].get_distribution(),
                       *model.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(*vec[6].get_distribution(),
                       *model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}].get_distribution());
    check_distribution(
        *vec[7].get_distribution(),
        *model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}].get_distribution());
    check_distribution(
        *vec[8].get_distribution(),
        *model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}].get_distribution());
    check_distribution(
        *vec[9].get_distribution(),
        *model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}].get_distribution());
    check_distribution(
        *vec[10].get_distribution(),
        *model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}].get_distribution());
    check_distribution(
        *vec[11].get_distribution(),
        *model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}].get_distribution());
    check_distribution(*vec[12].get_distribution(),
                       *model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}].get_distribution());
    check_distribution(
        *vec[13].get_distribution(),
        *model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(
        *vec[14].get_distribution(),
        *model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(
        *vec[15].get_distribution(),
        *model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(
        *vec[16].get_distribution(),
        *model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(
        *vec[17].get_distribution(),
        *model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(*vec[18].get_distribution(),
                       *model.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(*vec[19].get_distribution(),
                       *model.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0].get_distribution());
    // no dist for start day
    check_distribution(*vec[21].get_distribution(),
                       *model.parameters.get<mio::osecir::Seasonality>().get_distribution());

    EXPECT_EQ(vec[0], model.parameters.get<mio::osecir::ICUCapacity>());
    EXPECT_EQ(vec[1], model.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[2], model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[3], model.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[4], model.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[5], model.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[6], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]));
    EXPECT_EQ(vec[7], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]));
    EXPECT_EQ(vec[8], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]));
    EXPECT_EQ(vec[9], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]));
    EXPECT_EQ(vec[10], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]));
    EXPECT_EQ(vec[11], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]));
    EXPECT_EQ(vec[12], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]));
    EXPECT_EQ(vec[13], model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[14], model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[15], model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[16], model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[17], model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[18], model.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[19], model.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[20], model.parameters.get<mio::osecir::StartDay>());
    EXPECT_EQ(vec[21], model.parameters.get<mio::osecir::Seasonality>());
}

TEST(TestOdeSecir, testValueConstraints)
{
    double cont_freq = 10;

    double nb_total_t0 = 10000, nb_exp_t0 = -91, nb_inf_t0 = 39, nb_car_t0 = 36, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 8, nb_dead_t0 = 0;

    mio::osecir::Model model(1);

    model.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]       = 5.1;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = 5.86642;
    model.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]       = 5.08993;
    model.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]   = 11.6138;
    model.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0] = 9.16291;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    model.populations.set_total(nb_total_t0);
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = nb_exp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = nb_car_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]   = nb_inf_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]     = nb_hosp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]   = nb_icu_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]          = nb_rec_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]               = nb_dead_t0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                nb_total_t0);

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0] = 0.064519;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]   = 0.56758;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]   = 2.124921;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = 0.190609;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]        = 0.183693;
    model.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]                = 0.185556;
    model.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]                = 0.245801;

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.check_constraints();
    mio::set_log_level(mio::LogLevel::warn);

    EXPECT_EQ(-91, (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]));
    EXPECT_EQ(2.124921, model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0].value());
    EXPECT_NEAR(5.08993, model.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0], 1e-14);

    model.apply_constraints();

    EXPECT_EQ(0.0, (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]));
    EXPECT_EQ(0.0, model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0].value());
    EXPECT_NEAR(4.6, model.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0], 1e-14);
}

TEST(TestOdeSecir, testModelConstraints)
{
    double t0   = 0;
    double tmax = 57; // after 57 days with cont_freq 10 and winter, the virus would already decline
    double dt   = 0.1;

    double cont_freq = 10;

    double nb_total_t0 = 1000000, nb_exp_t0 = 10000, nb_inf_t0 = 5000, nb_car_t0 = 500, nb_hosp_t0 = 20, nb_icu_t0 = 0,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model model(1);

    model.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]       = 5.2;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = 5;
    model.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]       = 4.2;
    model.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]   = 10.;
    model.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0] = 8.;

    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = nb_exp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = nb_car_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]   = nb_inf_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]     = nb_hosp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]   = nb_icu_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]          = nb_rec_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]               = nb_dead_t0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                nb_total_t0);

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0] = 0.05;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]   = 1;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]   = 0.09;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = 0.25;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]        = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]                = 0.25;
    model.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]                = 0.3;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));

    model.apply_constraints();

    mio::TimeSeries<double> secihurd = simulate(t0, tmax, dt, model);
    double max_icu_cap               = 0;
    for (Eigen::Index i = 0; i < secihurd.get_num_time_points(); i++) {
        if (secihurd.get_value(i)[5] > max_icu_cap) {
            max_icu_cap = secihurd.get_value(i)[5];
        }
    }

    mio::TimeSeries<double> secihurd_interp = mio::interpolate_simulation_result(secihurd);

    model.parameters.set<mio::osecir::StartDay>(100);
    model.parameters.set<mio::osecir::Seasonality>(0.5);

    mio::TimeSeries<double> secihurd_season        = simulate(t0, tmax, dt, model);
    mio::TimeSeries<double> secihurd_season_interp = mio::interpolate_simulation_result(secihurd_season);

    for (Eigen::Index i = 0; i < secihurd_interp.get_num_time_points(); i++) {
        EXPECT_LE(secihurd_season_interp.get_value(i)[3], secihurd_interp.get_value(i)[3]) << " at row " << i;
    }

    model.parameters.set<mio::osecir::StartDay>(280);

    mio::TimeSeries<double> secihurd_season2        = simulate(t0, tmax, dt, model);
    mio::TimeSeries<double> secihurd_season2_interp = mio::interpolate_simulation_result(secihurd_season2);

    for (Eigen::Index i = 0; i < secihurd_interp.get_num_time_points(); i++) {
        EXPECT_GE(secihurd_season2_interp.get_value(i)[3], secihurd_interp.get_value(i)[3]) << " at row " << i;
    }

    // params.set_icu_capacity(max_icu_cap - 3);

    // secihurd = simulate(t0, tmax, dt, params);
    // for (Eigen::Index i = 0; i < secihurd.get_num_time_points(); i++) {
    //     EXPECT_LE(secihurd.get_value(i)[5], max_icu_cap - 2.5) << " at row " << i;
    // }

    // temporary test for random variables
    set_params_distributions_normal(model, t0, tmax, 0.2);

    for (size_t j = 0; j < 10; j++) {
        draw_sample(model);
        model.parameters.set<mio::osecir::ICUCapacity>(8000);
        secihurd = simulate(t0, tmax, dt, model);
        // max_icu_cap = 0;
        // for (Eigen::Index i = 0; i < secihurd.get_num_time_points(); i++) {
        //     if (secihurd.get_value(i)[5] > max_icu_cap) {
        //         max_icu_cap = secihurd.get_value(i)[5];
        //     }
        // }
        // printf("\n max cap: %.4f ", max_icu_cap);
        for (Eigen::Index i = 0; i < secihurd.get_num_time_points(); i++) {
            EXPECT_LE(secihurd.get_value(i)[5], 9000) << " at row " << i;
        }
    }
}

TEST(Secir, testAndTraceCapacity)
{
    double tinc = 5.2, tinf = 6, tserint = 4.2;

    double cont_freq = 10;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50;

    mio::osecir::Model model(1);
    auto& params = model.parameters;

    params.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0]       = tinc;
    params.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = tinf;
    params.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0]       = tserint;

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));

    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = nb_exp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = nb_car_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]   = nb_inf_t0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                nb_total_t0);

    params.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0] = 0.05;
    params.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]   = 1;
    params.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]   = 0.09;
    params.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = 0.25;

    params.apply_constraints();

    auto y = model.populations.get_compartments();

    auto dydt_default = Eigen::VectorXd(Eigen::Index(mio::osecir::InfectionState::Count));
    model.get_derivatives(y, y, 0, dydt_default);

    params.set<mio::osecir::TestAndTraceCapacity>(50);
    params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0] = 0.25 * 3;
    auto dydt_under_capacity = Eigen::VectorXd(Eigen::Index(mio::osecir::InfectionState::Count));
    model.get_derivatives(y, y, 0, dydt_under_capacity);

    params.set<mio::osecir::TestAndTraceCapacity>(10);
    params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0] = 0.25 * 3;
    auto dydt_over_capacity = Eigen::VectorXd(Eigen::Index(mio::osecir::InfectionState::Count));
    model.get_derivatives(y, y, 0, dydt_over_capacity);

    EXPECT_DOUBLE_EQ(dydt_under_capacity[(size_t)mio::osecir::InfectionState::Exposed],
                     dydt_default[(size_t)mio::osecir::InfectionState::Exposed]);
    EXPECT_GT(dydt_over_capacity[(size_t)mio::osecir::InfectionState::Exposed],
              dydt_default[(size_t)mio::osecir::InfectionState::Exposed]);
}

TEST(Secir, getInfectionsRelative)
{
    size_t num_groups = 3;
    mio::osecir::Model model((int)num_groups);
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] = 100.0;
    model.populations.set_difference_from_group_total<mio::AgeGroup>(
        {mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}, 10'000.0);
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedSymptoms}] = 50.0;
    model.populations.set_difference_from_group_total<mio::AgeGroup>(
        {mio::AgeGroup(1), mio::osecir::InfectionState::Susceptible}, 20'000.0);
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedSymptoms}] = 25.0;
    model.populations.set_difference_from_group_total<mio::AgeGroup>(
        {mio::AgeGroup(2), mio::osecir::InfectionState::Susceptible}, 40'000.0);

    mio::osecir::Simulation<> sim(model, 0.0);
    ASSERT_EQ(get_infections_relative(sim, 0.0, sim.get_result().get_last_value()),
              (100. + 50. + 25.) / (10'000 + 20'000 + 40'000));
}

TEST(Secir, get_reproduction_numbers)
{
    size_t num_groups = 3;
    mio::osecir::Model model((int)num_groups);

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(3, 3, 10));

    model.parameters.set<mio::osecir::StartDay>(60);
    model.parameters.set<mio::osecir::Seasonality>(0.2);

    //total population of 10.000
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}]        = 3000;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = 400;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = 50;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]   = 50;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]     = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]   = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]          = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]               = 0;

    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::Susceptible}]        = 4000;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::Exposed}]            = 350;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedNoSymptoms}] = 50;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedSymptoms}]   = 100;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedSevere}]     = 0;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedCritical}]   = 0;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::Recovered}]          = 0;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::Dead}]               = 0;

    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::Susceptible}]        = 1500;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::Exposed}]            = 200;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedNoSymptoms}] = 100;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedSymptoms}]   = 100;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedSevere}]     = 50;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedCritical}]   = 50;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::Recovered}]          = 0;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::Dead}]               = 0;

    for (auto i = mio::AgeGroup(0); i < (mio::AgeGroup)num_groups; i++) {
        model.parameters.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[i] = 5.8;
        model.parameters.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        model.parameters.get<mio::osecir::TimeInfectedSevere>()[i]   = 9.5;
        model.parameters.get<mio::osecir::TimeInfectedCritical>()[i] = 7.1;

        model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[i]  = 0.05;
        model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]    = 0.7;
        model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]    = 0.09;
        model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]    = 0.25;
        model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[i] = 0.45;
        model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[i]         = 0.2;
        model.parameters.get<mio::osecir::CriticalPerSevere>()[i]                 = 0.25;
        model.parameters.get<mio::osecir::DeathsPerCritical>()[i]                 = 0.3;
    }

    mio::TimeSeries<ScalarType> result((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_0((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_1((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_2((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_3((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_4((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_5((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_6((int)mio::osecir::InfectionState::Count * num_groups);

    model.apply_constraints();

    Eigen::VectorXd checkReproductionNumbers(7);
    checkReproductionNumbers << 3.7417747463385571, 3.7146918917498386, 3.6864528889471364, 3.6409596914789653,
        3.5870229607659483, 3.5154768089526982, 3.4647286447599575;

    result_0 << 3000, 400, 50, 50, 0, 0, 0, 0, 4000, 350, 50, 100, 0, 0, 0, 0, 1500, 200, 100, 100, 50, 50, 0, 0;

    result_1 << 2986.2412449968583132, 401.24167503358188469, 59.770932967055152574, 51.623527357069470156,
        0.17417985114868939078, 0.00022729149070066577247, 0.94821218285069908127, 3.1994478027648586745e-07,
        3981.6549933291448724, 357.29614039080206567, 58.339589387963329159, 100.73549052293542161,
        0.34417787494279222793, 0.00045098330804711045904, 1.6291568747836127073, 6.3611987013077957036e-07,
        1493.1206224984291566, 200.62083751679094235, 101.22767332092990955, 102.82956804964119613,
        49.824323129635679663, 49.431134850148950477, 2.7357775149744210097, 0.21006311944956934656;

    result_2 << 2971.5583317377245294, 403.35498239845526314, 69.116677006335692113, 53.650269608953109923,
        0.35280513579635025545, 0.00091208866793335101591, 1.9660194590816646443, 2.5649852579872425173e-06,
        3962.0777756502998272, 365.58091501551371039, 66.510165358127906643, 101.83063681480560092,
        0.68789659464275298983, 0.0017941695505704748168, 3.31081132947206358, 5.0675877399838152961e-06,
        1485.7791658688622647, 201.67749119922763157, 102.42114485557715398, 105.66531768774794386,
        49.660202843865796751, 48.869781815291510441, 5.5091572488792639462, 0.41773848054833506716;

    result_3 << 2942.856400298596327, 409.2171592213586564, 85.012918471965662093, 58.235842719339657947,
        0.68957046457456350197, 0.0033215276324750873937, 3.9847696050448093708, 1.7691487346664204683e-05,
        3923.8085337314628305, 382.81666733073143405, 80.907942542827129273, 104.66793072132175269,
        1.3092281281810831395, 0.0064233169718593535052, 6.4832397042377349905, 3.4524266497013886116e-05,
        1471.4282001492981635, 204.6085796106793282, 104.52885257313833733, 110.78105275523925854,
        49.393601382744165562, 47.878118436204538, 10.595987399938755047, 0.78560769275720121474;

    result_4 << 2907.2133871888117937, 418.67505276798732439, 101.7329866471358315, 64.706265603136884579,
        1.0991226964708620262, 0.0079167447918022521014, 6.5652040345289455203, 6.4317136376874475324e-05,
        3876.2845162517492099, 405.41241335635135101, 96.834733693923240594, 109.13012902817443717,
        2.0196012540506989019, 0.014990460834193003459, 10.30349232671609272, 0.00012362820069712732282,
        1453.6066935944058969, 209.33752638399366219, 106.91548231859164275, 116.55417450338485708,
        49.137205561836246659, 46.790465736158964205, 16.467779593518763193, 1.1906723081098882222;

    result_5 << 2857.0921695713273039, 434.42077730666193247, 121.68015133798982674, 74.740676897795324862,
        1.6806059449234134195, 0.016874713703690334687, 10.368548031693284983, 0.00019619590458404766956,
        3809.4562260951038297, 438.24759156052596154, 117.03004188242431383, 116.6701099172634315,
        2.9553258740962409234, 0.03104121600854261448, 15.609294091390772508, 0.00036936318682803968632,
        1428.5460847856636519, 217.21038865333096624, 110.10383240778212155, 123.94324598867467557,
        48.876141778624635492, 45.44704198659943728, 24.179729493339525703, 1.6935349059846935837;

    result_6 << 2823.7964337008384064, 445.79789336833465541, 133.43010492139592316, 81.843963977092016648,
        2.0791484574835688015, 0.024289771848477689081, 13.027831825199973181, 0.00033397780653605783387,
        3765.0619116011184815, 460.26765919119901582, 129.57755483037871613, 122.32995280769675617,
        3.5572150886264348735, 0.043856585969950019621, 19.161229559313984083, 0.00062033569677932126878,
        1411.8982168504192032, 222.8989466841673277, 112.20485190056430724, 128.5349283083724572, 48.751006448572013596,
        44.640554954647143404, 29.074588258609818325, 1.9969065946475592632;

    result.add_time_point(0, result_0);
    result.add_time_point(0.10000000000000000555, result_1);
    result.add_time_point(0.2000000000000000111, result_2);
    result.add_time_point(0.37998522153999991779, result_3);
    result.add_time_point(0.58252400303872020615, result_4);
    result.add_time_point(0.840598918925461569, result_5);
    result.add_time_point(1, result_6);

    mio::osecir::Simulation<> sim(model, 0.0);
    sim.get_result() = result;

    for (int i = 0; i < sim.get_result().get_num_time_points(); i++) {
        EXPECT_NEAR(checkReproductionNumbers[i], mio::osecir::get_reproduction_numbers(sim)[i], 1e-12);
    }
}

TEST(Secir, get_reproduction_number)
{
    size_t num_groups = 3;
    mio::osecir::Model model((int)num_groups);

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(3, 3, 10));

    model.parameters.set<mio::osecir::StartDay>(60);
    model.parameters.set<mio::osecir::Seasonality>(0.2);

    //total population of 10.000
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}]        = 3000;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = 400;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = 50;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]   = 50;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]     = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]   = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]          = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]               = 0;

    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::Susceptible}]        = 4000;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::Exposed}]            = 350;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedNoSymptoms}] = 50;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedSymptoms}]   = 100;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedSevere}]     = 0;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedCritical}]   = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]          = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]               = 0;

    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::Susceptible}]        = 1500;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::Exposed}]            = 200;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedNoSymptoms}] = 100;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedSymptoms}]   = 100;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedSevere}]     = 50;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedCritical}]   = 50;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]          = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]               = 0;

    for (auto i = mio::AgeGroup(0); i < (mio::AgeGroup)num_groups; i++) {
        model.parameters.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        model.parameters.get<mio::osecir::TimeInfectedSymptoms>()[i] = 5.8;
        model.parameters.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        model.parameters.get<mio::osecir::TimeInfectedSevere>()[i]   = 9.5;
        model.parameters.get<mio::osecir::TimeInfectedCritical>()[i] = 7.1;

        model.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[i]  = 0.05;
        model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]    = 0.7;
        model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]    = 0.09;
        model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]    = 0.25;
        model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[i] = 0.45;
        model.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[i]         = 0.2;
        model.parameters.get<mio::osecir::CriticalPerSevere>()[i]                 = 0.25;
        model.parameters.get<mio::osecir::DeathsPerCritical>()[i]                 = 0.3;
    }

    mio::TimeSeries<ScalarType> result((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_0((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_1((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_2((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_3((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_4((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_5((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_6((int)mio::osecir::InfectionState::Count * num_groups);

    model.apply_constraints();

    result_0 << 3000, 400, 50, 50, 0, 0, 0, 0, 4000, 350, 50, 100, 0, 0, 0, 0, 1500, 200, 100, 100, 50, 50, 0, 0;

    result_1 << 2900, 500, 50, 50, 0, 0, 0, 0, 4000, 350, 50, 100, 0, 0, 0, 0, 1500, 200, 100, 100, 50, 50, 0, 0;

    result_2 << 2850, 550, 50, 50, 0, 0, 0, 0, 4000, 350, 0, 150, 0, 0, 0, 0, 1500, 200, 100, 100, 50, 50, 0, 0;

    result_3 << 2850, 550, 50, 50, 0, 0, 0, 0, 4000, 350, 0, 150, 0, 0, 0, 0, 1300, 400, 100, 100, 50, 50, 0, 0;

    result_4 << 2800, 600, 50, 50, 0, 0, 0, 0, 4000, 300, 0, 200, 0, 0, 0, 0, 1300, 400, 100, 100, 50, 50, 0, 0;

    result_5 << 2800, 600, 50, 50, 0, 0, 0, 0, 4000, 300, 0, 200, 0, 0, 0, 0, 1300, 400, 100, 100, 50, 50, 0, 0;

    result_6 << 2700, 600, 100, 100, 0, 0, 0, 0, 4000, 300, 0, 200, 0, 0, 0, 0, 1300, 400, 100, 100, 0, 100, 0, 0;

    result.add_time_point(0.0, result_0);
    result.add_time_point(0.1000000000000000000, result_1);
    result.add_time_point(0.2000000000000000000, result_2);
    result.add_time_point(0.4000000000000000000, result_3);
    result.add_time_point(0.6000000000000000000, result_4);
    result.add_time_point(0.8000000000000000000, result_5);
    result.add_time_point(1.0, result_6);

    mio::osecir::Simulation<> sim(model, 0.0);
    sim.get_result() = result;

    EXPECT_FALSE(mio::osecir::get_reproduction_number(result.get_time(0) - 0.5, sim)); //Test for indices out of range
    EXPECT_FALSE(mio::osecir::get_reproduction_number(result.get_last_time() + 0.5, sim));
    EXPECT_FALSE(mio::osecir::get_reproduction_number((size_t)result.get_num_time_points(), sim));

    EXPECT_EQ(mio::osecir::get_reproduction_number((size_t)0, sim).value(),
              mio::osecir::get_reproduction_number(0.0, sim).value());

    //Test one function for integer timepoints
    EXPECT_NEAR(mio::osecir::get_reproduction_number((size_t)0, sim).value(), 3.7417747463385571,
                1e-12); //Calculated by hand
    EXPECT_NEAR(mio::osecir::get_reproduction_number((size_t)4, sim).value(), 3.4678495894266841, 1e-12);
    EXPECT_NEAR(mio::osecir::get_reproduction_number((size_t)6, sim).value(), 3.4060279836965339,
                1e-12); //Calculated by hand

    EXPECT_NEAR(mio::osecir::get_reproduction_number(0.05, sim).value(), 3.7153740911442856, 1e-12);
    EXPECT_NEAR(mio::osecir::get_reproduction_number(0.5, sim).value(), 3.4833316698917707, 1e-12);
    EXPECT_NEAR(mio::osecir::get_reproduction_number(0.85, sim).value(), 3.4450376807796337, 1e-12);
}

TEST(Secir, get_migration_factors)
{
    auto beta                                                                              = 0.25;
    auto max_beta                                                                          = 0.5;
    auto model                                                                             = mio::osecir::Model(1);
    model.parameters.get<mio::osecir::IncubationTime>().array()                            = 5.0;
    model.parameters.get<mio::osecir::SerialInterval>().array()                            = 4.0;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>().array()            = 0.1;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>().array()            = beta;
    model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>().array()         = max_beta;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = 100;
    mio::osecir::Simulation<> sim(model, 0.0);
    {
        sim.get_model().parameters.get<mio::osecir::TestAndTraceCapacity>() = 45.;
        auto factors = get_migration_factors(sim, 0.0, sim.get_result().get_last_value());
        auto cmp     = Eigen::VectorXd::Ones(Eigen::Index(mio::osecir::InfectionState::Count)).eval();
        cmp[Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] = beta;
        ASSERT_THAT(print_wrap(factors), MatrixNear(cmp));
    }
    {
        sim.get_model().parameters.get<mio::osecir::TestAndTraceCapacity>() = 45. / 5.;
        auto factors = get_migration_factors(sim, 0.0, sim.get_result().get_last_value());
        auto cmp     = Eigen::VectorXd::Ones(Eigen::Index(mio::osecir::InfectionState::Count)).eval();
        cmp[Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] = max_beta;
        ASSERT_THAT(print_wrap(factors), MatrixNear(cmp));
    }
    {
        sim.get_model().parameters.get<mio::osecir::TestAndTraceCapacity>() = 20.;
        auto factors = get_migration_factors(sim, 0.0, sim.get_result().get_last_value());
        ASSERT_GT(factors[Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)], beta);
        ASSERT_LT(factors[Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)], max_beta);
    }
}

TEST(Secir, check_constraints_parameters)
{
    auto model = mio::osecir::Model(1);
    ASSERT_EQ(model.parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::osecir::Seasonality>(-0.5);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::Seasonality>(0.2);
    model.parameters.set<mio::osecir::ICUCapacity>(-2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::ICUCapacity>(2);
    model.parameters.set<mio::osecir::IncubationTime>(-2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::IncubationTime>(2);
    model.parameters.set<mio::osecir::SerialInterval>(1);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::SerialInterval>(5);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::SerialInterval>(1.5);
    model.parameters.set<mio::osecir::TimeInfectedSymptoms>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::TimeInfectedSymptoms>(2);
    model.parameters.set<mio::osecir::TimeInfectedSevere>(-1);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::TimeInfectedSevere>(2);
    model.parameters.set<mio::osecir::TimeInfectedCritical>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::TimeInfectedCritical>(2);
    model.parameters.set<mio::osecir::TransmissionProbabilityOnContact>(2.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::TransmissionProbabilityOnContact>(0.5);
    model.parameters.set<mio::osecir::RelativeTransmissionNoSymptoms>(-1.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::RelativeTransmissionNoSymptoms>(0.5);
    model.parameters.set<mio::osecir::RecoveredPerInfectedNoSymptoms>(3.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::RecoveredPerInfectedNoSymptoms>(0.5);
    model.parameters.set<mio::osecir::RiskOfInfectionFromSymptomatic>(-2.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::RiskOfInfectionFromSymptomatic>(0.5);
    model.parameters.set<mio::osecir::SeverePerInfectedSymptoms>(-1.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::SeverePerInfectedSymptoms>(0.5);
    model.parameters.set<mio::osecir::CriticalPerSevere>(-1.0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::CriticalPerSevere>(0.5);
    model.parameters.set<mio::osecir::DeathsPerCritical>(1.1);
    ASSERT_EQ(model.parameters.check_constraints(), 1);
    mio::set_log_level(mio::LogLevel::warn);
}
