/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "distributions_helpers.h"
#include "load_test_data.h"
#include "matchers.h"
#include "temp_file_register.h"

#include "ode_secir/model.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/parameters.h"
#include "ode_secir/parameters_io.h"
#include "memilio/data/analyze_result.h"
#include "memilio/math/adapt_rk.h"

#include <gtest/gtest.h>

TEST(TestOdeSecir, compareWithPreviousRun)
{
    /*
    A similar test is implemented in python (without custom integrator) to compare the results of both simulations.
    If this test is change the corresponding python test needs to be changed aswell (also updating the data file).
    */
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double cont_freq = 10;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model<double> model(1);

    model.parameters.set<mio::osecir::StartDay>(60);
    model.parameters.set<mio::osecir::Seasonality<double>>(0.2);

    model.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]            = 3.2;
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0] = 2.0;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0]   = 5.8;
    model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0]     = 9.5;
    model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0]   = 7.1;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    model.populations.set_total(nb_total_t0);
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]                     = nb_exp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]          = nb_car_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]            = nb_inf_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]              = nb_hosp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]            = nb_icu_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]                   = nb_rec_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]                        = nb_dead_t0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                nb_total_t0);

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0]  = 0.05;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]    = 0.7;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]    = 0.09;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]    = 0.25;
    model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0] = 0.45;
    model.parameters.get<mio::osecir::TestAndTraceCapacity<double>>()                                = 35;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0]         = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0]                 = 0.3;
    model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0]                 = 0.3;

    model.apply_constraints();

    auto integrator = std::make_shared<mio::RKIntegratorCore<double>>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-3);
    integrator->set_abs_tolerance(1e-1);
    mio::TimeSeries<double> secihurd =
        mio::simulate<double, mio::osecir::Model<double>>(t0, tmax, dt, model, integrator);

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

    mio::osecir::Model<double> model(1);
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

    mio::osecir::Model<double> model(1);

    model.parameters.get<mio::osecir::TestAndTraceCapacity<double>>() = 35;

    model.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]            = 3.2;
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0] = 2.0;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0]   = 5.8;
    model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0]     = 9.5;
    model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0]   = 7.1;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    model.populations.set_total(nb_total_t0);
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]                     = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]          = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]            = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]              = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]            = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]                   = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]                        = 10;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                nb_total_t0);

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0]  = 0.05;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]    = 1;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]    = 0.09;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]    = 0.25;
    model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0] = 0.45;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0]         = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0]                 = 0.25;
    model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0]                 = 0.3;

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

    mio::osecir::Model<double> model(1);

    model.parameters.set<mio::osecir::ICUCapacity<double>>(icu_cap);

    model.parameters.set<mio::osecir::StartDay>(start_day);
    model.parameters.set<mio::osecir::Seasonality<double>>(seasonality);

    model.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]            = 3.2;
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0] = 2.0;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0]   = 5;
    model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0]     = 10.;
    model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0]   = 8.;

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

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0] = 0.05;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]   = 0.67;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]   = 0.09;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]   = 0.25;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0]        = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0]                = 0.24;
    model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0]                = 0.3;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    mio::osecir::Model<double> model2{model}; // copy constructor

    EXPECT_EQ(model.parameters.get<mio::osecir::ICUCapacity<double>>(),
              model2.parameters.get<mio::osecir::ICUCapacity<double>>());
    EXPECT_EQ(model.parameters.get<mio::osecir::StartDay>(), model2.parameters.get<mio::osecir::StartDay>());
    EXPECT_EQ(model.parameters.get<mio::osecir::Seasonality<double>>(),
              model2.parameters.get<mio::osecir::Seasonality<double>>());

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

    EXPECT_EQ(model.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::ContactPatterns<double>>().get_cont_freq_mat(),
              model2.parameters.get<mio::osecir::ContactPatterns<double>>().get_cont_freq_mat());

    EXPECT_EQ(model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0],
              model2.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<mio::osecir::ContactPatterns<double>>().get_cont_freq_mat(),
              model2.parameters.get<mio::osecir::ContactPatterns<double>>().get_cont_freq_mat());

    mio::osecir::Model<double> model3 = std::move(model2); // move constructor

    EXPECT_EQ(model.parameters.get<mio::osecir::ICUCapacity<double>>(),
              model3.parameters.get<mio::osecir::ICUCapacity<double>>());
    EXPECT_EQ(model.parameters.get<mio::osecir::StartDay>(), model3.parameters.get<mio::osecir::StartDay>());
    EXPECT_EQ(model.parameters.get<mio::osecir::Seasonality<double>>(),
              model3.parameters.get<mio::osecir::Seasonality<double>>());

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

    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0]);

    EXPECT_EQ(model3.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0],
              model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0]);

    EXPECT_EQ(model.parameters.get<mio::osecir::ContactPatterns<double>>().get_cont_freq_mat(),
              model3.parameters.get<mio::osecir::ContactPatterns<double>>().get_cont_freq_mat());

    mio::osecir::Model<double> model4 = model3; // copy assignment constructor

    EXPECT_EQ(model4.parameters.get<mio::osecir::ICUCapacity<double>>(),
              model3.parameters.get<mio::osecir::ICUCapacity<double>>());
    EXPECT_EQ(model4.parameters.get<mio::osecir::StartDay>(), model3.parameters.get<mio::osecir::StartDay>());
    EXPECT_EQ(model4.parameters.get<mio::osecir::Seasonality<double>>(),
              model3.parameters.get<mio::osecir::Seasonality<double>>());

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

    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0]);

    EXPECT_EQ(model3.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0],
              model4.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0]);

    EXPECT_EQ(model4.parameters.get<mio::osecir::ContactPatterns<double>>().get_cont_freq_mat(),
              model3.parameters.get<mio::osecir::ContactPatterns<double>>().get_cont_freq_mat());

    mio::osecir::Model<double> model5 = std::move(model4); // move assignment constructor

    EXPECT_EQ(model5.parameters.get<mio::osecir::ICUCapacity<double>>(),
              model3.parameters.get<mio::osecir::ICUCapacity<double>>());
    EXPECT_EQ(model5.parameters.get<mio::osecir::StartDay>(), model3.parameters.get<mio::osecir::StartDay>());
    EXPECT_EQ(model5.parameters.get<mio::osecir::Seasonality<double>>(),
              model3.parameters.get<mio::osecir::Seasonality<double>>());

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

    EXPECT_EQ(model5.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0]);

    EXPECT_EQ(model5.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0],
              model3.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0]);

    EXPECT_EQ(model5.parameters.get<mio::osecir::ContactPatterns<double>>().get_cont_freq_mat(),
              model3.parameters.get<mio::osecir::ContactPatterns<double>>().get_cont_freq_mat());
}

TEST(TestOdeSecir, testSettersAndGetters)
{
    std::vector<mio::UncertainValue<double>> vec;

    for (int i = 1; i < 23; i++) {
        mio::UncertainValue<double> val = mio::UncertainValue<double>(i);
        val.set_distribution(mio::ParameterDistributionNormal(i, 10 * i, 5 * i, i / 10.0));
        vec.push_back(val);
    }

    mio::osecir::Model<double> model(1);
    EXPECT_EQ(model.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0].get_distribution().get(),
              nullptr);

    model.parameters.set<mio::osecir::ICUCapacity<double>>(vec[0]);

    model.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]            = vec[1];
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0] = vec[2];
    model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0]   = vec[3];
    model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0]     = vec[4];
    model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0]   = vec[5];

    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = vec[6];
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = vec[7];
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]   = vec[8];
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]     = vec[9];
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]   = vec[10];
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]          = vec[11];
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]               = vec[12];

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0] = vec[13];
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]   = vec[14];
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]   = vec[15];
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]   = vec[16];
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0]        = vec[17];
    model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0]                = vec[18];
    model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0]                = vec[19];

    EXPECT_NE(model.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0].get_distribution().get(),
              nullptr);

    check_distribution(*vec[0].get_distribution(),
                       *model.parameters.get<mio::osecir::ICUCapacity<double>>().get_distribution());

    model.parameters.set<mio::osecir::StartDay>(vec[20]);
    model.parameters.set<mio::osecir::Seasonality<double>>(vec[21]);

    check_distribution(*vec[1].get_distribution(),
                       *model.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(
        *vec[2].get_distribution(),
        *model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(
        *vec[3].get_distribution(),
        *model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(
        *vec[4].get_distribution(),
        *model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(
        *vec[5].get_distribution(),
        *model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0].get_distribution());
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
    check_distribution(*vec[13].get_distribution(),
                       *model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0]
                            .get_distribution());
    check_distribution(*vec[14].get_distribution(),
                       *model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]
                            .get_distribution());
    check_distribution(*vec[15].get_distribution(),
                       *model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]
                            .get_distribution());
    check_distribution(*vec[16].get_distribution(),
                       *model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]
                            .get_distribution());
    check_distribution(
        *vec[17].get_distribution(),
        *model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(
        *vec[18].get_distribution(),
        *model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0].get_distribution());
    check_distribution(
        *vec[19].get_distribution(),
        *model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0].get_distribution());
    // no dist for start day
    check_distribution(*vec[21].get_distribution(),
                       *model.parameters.get<mio::osecir::Seasonality<double>>().get_distribution());

    EXPECT_EQ(vec[0], model.parameters.get<mio::osecir::ICUCapacity<double>>());
    EXPECT_EQ(vec[1], model.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[2], model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[3], model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[4], model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[5], model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[6], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]));
    EXPECT_EQ(vec[7], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]));
    EXPECT_EQ(vec[8], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]));
    EXPECT_EQ(vec[9], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]));
    EXPECT_EQ(vec[10], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]));
    EXPECT_EQ(vec[11], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]));
    EXPECT_EQ(vec[12], (model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]));
    EXPECT_EQ(vec[13], model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[14], model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[15], model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[16], model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[17], model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[18], model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[19], model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0]);
    EXPECT_EQ(vec[20], model.parameters.get<mio::osecir::StartDay>());
    EXPECT_EQ(vec[21], model.parameters.get<mio::osecir::Seasonality<double>>());
}

// Test model initialization with total population of 0 and ensure get_flows returns no NaN values
TEST(TestOdeSecir, population_zero_no_nan)
{
    // initialize simple model with total population 0
    mio::osecir::Model<double> model(1);
    model.populations.set_total(0.0);

    // call the get_flows function
    auto dydt_default = Eigen::VectorXd(15);
    dydt_default.setZero();
    auto y0 = model.get_initial_values();
    model.get_flows(y0, y0, 0, dydt_default);

    // check that there are now NaN values in dydt_default
    for (int i = 0; i < dydt_default.size(); i++) {
        EXPECT_FALSE(std::isnan(dydt_default[i]));
    }
}

TEST(TestOdeSecir, testDamping)
{
    // Test functionality of dampings
    // (initially only implemented contact reductions but now also allow contact increases).
    // Contact matrices with dampings are cosine-smoothed in decline/increase along one day to be C1 differentiable.
    // If EulerIntegratorCore with dt=1 is used, we jump across this smoothing so that we can express the relationship
    // between old and new transmission directly, only including damping, contact, and transmission probability values.
    double t0   = 0;
    double dt   = 1;
    double tmax = dt;

    double cont_freq = 10;

    double nb_total_t0 = 1000, nb_inf_t0 = 10;

    auto integrator = std::make_shared<mio::EulerIntegratorCore<ScalarType>>();

    // default model run to be compared against
    mio::osecir::Model model_a(1);
    model_a.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] = nb_inf_t0;
    model_a.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                  nb_total_t0);
    mio::ContactMatrixGroup& contact_matrix_a = model_a.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix_a[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    // set probability of transmission and risk of infection to 1.
    model_a.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>() = 1.0;
    model_a.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()   = 1.0;
    auto result_a = mio::simulate_flows<ScalarType>(t0, tmax, dt, model_a, integrator);

    // reduced transmission
    mio::osecir::Model model_b{model_a};
    model_b.populations.set_total(nb_total_t0);
    model_b.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] = nb_inf_t0;
    model_b.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                  nb_total_t0);
    mio::ContactMatrixGroup& contact_matrix_b = model_b.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix_b[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix_b[0].add_damping(0.5, mio::SimulationTime(0.));
    auto result_b = mio::simulate_flows<ScalarType>(t0, tmax, dt, model_b, integrator);
    EXPECT_EQ(2 * result_b[1].get_last_value()[0], result_a[1].get_last_value()[0]);

    // no transmission
    mio::osecir::Model model_c{model_a};
    model_c.populations.set_total(nb_total_t0);
    model_c.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] = nb_inf_t0;
    model_c.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                  nb_total_t0);
    mio::ContactMatrixGroup& contact_matrix_c = model_c.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix_c[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix_c[0].add_damping(1., mio::SimulationTime(0.));
    auto result_c = mio::simulate_flows<ScalarType>(t0, tmax, dt, model_c, integrator);
    EXPECT_EQ(result_c[1].get_last_value()[0], 0.0);

    // increased transmission to a factor of two (by +1)
    mio::osecir::Model model_d{model_a};
    model_d.populations.set_total(nb_total_t0);
    model_d.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] = nb_inf_t0;
    model_d.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                  nb_total_t0);
    mio::ContactMatrixGroup& contact_matrix_d = model_d.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix_d[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix_d[0].add_damping(-1., mio::SimulationTime(0.));
    auto result_d = mio::simulate_flows<ScalarType>(t0, tmax, dt, model_d, integrator);
    EXPECT_EQ(2 * result_a[1].get_last_value()[0], result_d[1].get_last_value()[0]);
}

TEST(TestOdeSecir, testModelConstraints)
{
    mio::set_log_level(
        mio::LogLevel::err); //as many random things are drawn, warnings are inevitable and cluster output
    double t0   = 0;
    double tmax = 57; // after 57 days with cont_freq 10 and winter, the virus would already decline
    double dt   = 0.1;

    double cont_freq = 10;

    double nb_total_t0 = 1000000, nb_exp_t0 = 10000, nb_inf_t0 = 5000, nb_car_t0 = 500, nb_hosp_t0 = 20, nb_icu_t0 = 0,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model<double> model(1);

    model.parameters.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]            = 2.6;
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0] = 2.6;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0]   = 5;
    model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0]     = 10.;
    model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0]   = 8.;

    model.parameters.get<mio::osecir::Seasonality<double>>()          = 0.0;
    model.parameters.get<mio::osecir::ICUCapacity<double>>()          = 100.0;
    model.parameters.get<mio::osecir::TestAndTraceCapacity<double>>() = 10.0;

    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = nb_exp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = nb_car_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]   = nb_inf_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]     = nb_hosp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]   = nb_icu_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]          = nb_rec_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]               = nb_dead_t0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                nb_total_t0);

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0] = 0.05;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]   = 1;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]   = 0.09;
    model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()                  = 0.85;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]   = 0.25;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[(mio::AgeGroup)0]        = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)0]                = 0.25;
    model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[(mio::AgeGroup)0]                = 0.3;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<double>>();
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

    // Tests that infection numbers are higher in Winter season
    model.parameters.set<mio::osecir::StartDay>(100);
    model.parameters.set<mio::osecir::Seasonality<double>>(0.5);

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

    // temporary test for random variables
    set_params_distributions_normal(model, t0, tmax, 0.2);
    model.parameters.set<mio::osecir::Seasonality<double>>(mio::UncertainValue(0.0));
    model.parameters.set<mio::osecir::ICUCapacity<double>>(mio::UncertainValue(8000));
    for (size_t j = 0; j < 10; j++) {
        draw_sample(model);
        secihurd = simulate(t0, tmax, dt, model);
        for (Eigen::Index i = 0; i < secihurd.get_num_time_points(); i++) {
            EXPECT_LE(secihurd.get_value(i)[5], 9000) << " at row " << i;
        }
    }
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSecir, testAndTraceCapacity)
{
    double cont_freq = 10;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50;

    mio::osecir::Model<double> model(1);
    auto& params = model.parameters;

    params.get<mio::osecir::TimeExposed<double>>()[(mio::AgeGroup)0]            = 3.2;
    params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[(mio::AgeGroup)0] = 2.0;
    params.get<mio::osecir::TimeInfectedSymptoms<double>>()[(mio::AgeGroup)0]   = 6.;

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));

    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]            = nb_exp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = nb_car_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]   = nb_inf_t0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                nb_total_t0);

    params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0] = 0.05;
    params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[(mio::AgeGroup)0]   = 1;
    params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[(mio::AgeGroup)0]   = 0.09;
    params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0]   = 0.25;

    params.apply_constraints();

    auto y = model.populations.get_compartments();

    auto dydt_default = Eigen::VectorXd(Eigen::Index(mio::osecir::InfectionState::Count));
    model.get_derivatives(y.cast<double>(), y.cast<double>(), 0, dydt_default);

    params.set<mio::osecir::TestAndTraceCapacity<double>>(50);
    params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0] = 0.25 * 3;
    auto dydt_under_capacity = Eigen::VectorXd(Eigen::Index(mio::osecir::InfectionState::Count));
    model.get_derivatives(y.cast<double>(), y.cast<double>(), 0, dydt_under_capacity);

    params.set<mio::osecir::TestAndTraceCapacity<double>>(10);
    params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[(mio::AgeGroup)0] = 0.25 * 3;
    auto dydt_over_capacity = Eigen::VectorXd(Eigen::Index(mio::osecir::InfectionState::Count));
    model.get_derivatives(y.cast<double>(), y.cast<double>(), 0, dydt_over_capacity);

    EXPECT_DOUBLE_EQ(dydt_under_capacity[(size_t)mio::osecir::InfectionState::Exposed],
                     dydt_default[(size_t)mio::osecir::InfectionState::Exposed]);
    EXPECT_GT(dydt_over_capacity[(size_t)mio::osecir::InfectionState::Exposed],
              dydt_default[(size_t)mio::osecir::InfectionState::Exposed]);
}

TEST(TestOdeSecir, getInfectionsRelative)
{
    size_t num_groups = 3;
    mio::osecir::Model<double> model((int)num_groups);
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
    ASSERT_EQ(mio::osecir::get_infections_relative<double>(sim, 0.0, sim.get_result().get_last_value()),
              (100. + 50. + 25.) / (10'000 + 20'000 + 40'000));
}

TEST(TestOdeSecir, get_reproduction_number)
{
    const size_t num_groups = 3;
    mio::osecir::Model model((int)num_groups);

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(3, 3, 10));

    model.parameters.set<mio::osecir::StartDay>(60);
    model.parameters.set<mio::osecir::Seasonality<double>>(0.2);

    //total population of 10.000
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}]                 = 3000;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]                     = 400;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]          = 50;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]            = 50;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]              = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]            = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]                   = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]                        = 0;

    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::Susceptible}]                 = 4000;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::Exposed}]                     = 350;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedNoSymptoms}]          = 50;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedSymptoms}]            = 100;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedSevere}]              = 0;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedCritical}]            = 0;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::Recovered}]                   = 0;
    model.populations[{mio::AgeGroup(1), mio::osecir::InfectionState::Dead}]                        = 0;

    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::Susceptible}]                 = 1500;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::Exposed}]                     = 200;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedNoSymptoms}]          = 100;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedSymptoms}]            = 100;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedSevere}]              = 50;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::InfectedCritical}]            = 50;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::Recovered}]                   = 0;
    model.populations[{mio::AgeGroup(2), mio::osecir::InfectionState::Dead}]                        = 0;

    for (auto i = mio::AgeGroup(0); i < (mio::AgeGroup)num_groups; i++) {
        model.parameters.get<mio::osecir::TimeExposed<double>>()[i]            = 3.2;
        model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i] = 2.0;
        model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]   = 5.8;
        model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[i]     = 9.5;
        model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[i]   = 7.1;

        model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i]  = 0.05;
        model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[i]    = 0.7;
        model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i]    = 0.09;
        model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i]    = 0.25;
        model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[i] = 0.45;
        model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i]         = 0.2;
        model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[i]                 = 0.25;
        model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[i]                 = 0.3;
    }
    model.parameters.get<mio::osecir::ICUCapacity<double>>()          = std::numeric_limits<double>::max();
    model.parameters.get<mio::osecir::TestAndTraceCapacity<double>>() = std::numeric_limits<double>::max();

    mio::TimeSeries<ScalarType> time_series1((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_0((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_1((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_2((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_3((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_4((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_5((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_6((int)mio::osecir::InfectionState::Count * num_groups);

    model.apply_constraints();

    result_0 << 3000, 400, 50, 0, 50, 0, 0, 0, 0, 0, 4000, 350, 50, 0, 100, 0, 0, 0, 0, 0, 1500, 200, 100, 0, 100, 0,
        50, 50, 0, 0;

    result_1 << 2900, 500, 50, 0, 50, 0, 0, 0, 0, 0, 4000, 350, 50, 0, 100, 0, 0, 0, 0, 0, 1500, 200, 100, 0, 100, 0,
        50, 50, 0, 0;

    result_2 << 2850, 550, 50, 0, 50, 0, 0, 0, 0, 0, 4000, 350, 0, 0, 150, 0, 0, 0, 0, 0, 1500, 200, 100, 0, 100, 0, 50,
        50, 0, 0;

    result_3 << 2850, 550, 50, 0, 50, 0, 0, 0, 0, 0, 4000, 350, 0, 0, 150, 0, 0, 0, 0, 0, 1300, 400, 100, 0, 100, 0, 50,
        50, 0, 0;

    result_4 << 2800, 600, 50, 0, 50, 0, 0, 0, 0, 0, 4000, 300, 0, 0, 200, 0, 0, 0, 0, 0, 1300, 400, 100, 0, 100, 0, 50,
        50, 0, 0;

    result_5 << 2800, 600, 50, 0, 50, 0, 0, 0, 0, 0, 4000, 300, 0, 0, 200, 0, 0, 0, 0, 0, 1300, 400, 100, 0, 100, 0, 50,
        50, 0, 0;

    result_6 << 2700, 600, 100, 0, 100, 0, 0, 0, 0, 0, 4000, 300, 0, 0, 200, 0, 0, 0, 0, 0, 1300, 400, 100, 0, 100, 0,
        0, 100, 0, 0;

    time_series1.add_time_point(0.0, result_0);
    time_series1.add_time_point(0.1000000000000000000, result_1);
    time_series1.add_time_point(0.2000000000000000000, result_2);
    time_series1.add_time_point(0.4000000000000000000, result_3);
    time_series1.add_time_point(0.6000000000000000000, result_4);
    time_series1.add_time_point(0.8000000000000000000, result_5);
    time_series1.add_time_point(1.0, result_6);

    mio::osecir::Simulation<> sim(model, 0.0);
    sim.get_result() = time_series1;

    EXPECT_FALSE(
        mio::osecir::get_reproduction_number(time_series1.get_time(0) - 0.5, sim)); //Test for indices out of range
    EXPECT_FALSE(mio::osecir::get_reproduction_number(time_series1.get_last_time() + 0.5, sim));
    EXPECT_FALSE(mio::osecir::get_reproduction_number((size_t)time_series1.get_num_time_points(), sim));

    EXPECT_EQ(mio::osecir::get_reproduction_number((size_t)0, sim).value(),
              mio::osecir::get_reproduction_number(0.0, sim).value());

    //Test one function for integer timepoints
    EXPECT_NEAR(mio::osecir::get_reproduction_number((size_t)0, sim).value(), 3.7417747463385571, 1e-12);
    EXPECT_NEAR(mio::osecir::get_reproduction_number((size_t)4, sim).value(), 3.5005445618245297, 1e-12);
    EXPECT_NEAR(mio::osecir::get_reproduction_number((size_t)6, sim).value(), 3.4540372055485653, 1e-12);
    EXPECT_NEAR(mio::osecir::get_reproduction_number(0.05, sim).value(), 3.719862942211813, 1e-12);
    EXPECT_NEAR(mio::osecir::get_reproduction_number(0.5, sim).value(), 3.5121859116705565, 1e-12);
    EXPECT_NEAR(mio::osecir::get_reproduction_number(0.85, sim).value(), 3.4874972585249733, 1e-12);

    //Test handling non-invertibility of V for certain values
    mio::TimeSeries<ScalarType>::Vector result_7((int)mio::osecir::InfectionState::Count * num_groups);
    double icu_occupancy = 0.95 * model.parameters.get<mio::osecir::ICUCapacity<double>>();
    double severe1       = model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[(mio::AgeGroup)0] /
                     (model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[(mio::AgeGroup)0] * 5 *
                      model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[(mio::AgeGroup)1] *
                      3.141592653589793 / (model.parameters.get<mio::osecir::ICUCapacity<double>>()) *
                      std::sin(3.141592653589793 / (0.1 * model.parameters.get<mio::osecir::ICUCapacity<double>>()) *
                               (icu_occupancy - 0.9 * model.parameters.get<mio::osecir::ICUCapacity<double>>())));

    mio::TimeSeries<ScalarType> time_series2((int)mio::osecir::InfectionState::Count * num_groups);
    result_7 << 1000, 0, 0, 0, 0, 0, severe1, 0.95 * model.parameters.get<mio::osecir::ICUCapacity<double>>(), 0, 0,
        1000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    time_series2.add_time_point(0.0, result_7);
    sim.get_result() = time_series2;
    EXPECT_FALSE(mio::osecir::get_reproduction_number((size_t)0, sim));

    //Test in the case of limited test-and-trace capacity:

    //Test for small test and trace
    model.parameters.get<mio::osecir::TestAndTraceCapacity<double>>() = 0;
    mio::osecir::Simulation<> sim2(model, 0.0);
    sim2.get_result() = time_series1;
    EXPECT_NEAR(mio::osecir::get_reproduction_number((size_t)0, sim2).value(), 5.1941804908632792, 1e-12);

    // Test special domain for test-and-trace capacity/requirement:
    model.parameters.get<mio::osecir::TestAndTraceCapacity<double>>() = 1;
    mio::osecir::Simulation<> sim3(model, 0.0);
    mio::TimeSeries<ScalarType> time_series3((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_8((int)mio::osecir::InfectionState::Count * num_groups);
    result_8 << 100, 0, 10, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0;
    time_series3.add_time_point(0.0, result_8);
    sim3.get_result() = time_series3;
    EXPECT_NEAR(mio::osecir::get_reproduction_number((size_t)0, sim3).value(), 1.8462669866786356, 1e-12);

    //Test handling of zero population in at least one agegroup
    mio::TimeSeries<ScalarType> time_series4((int)mio::osecir::InfectionState::Count * num_groups);
    mio::TimeSeries<ScalarType>::Vector result_9((int)mio::osecir::InfectionState::Count * num_groups);
    result_9 << 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    time_series4.add_time_point(0.0, result_9);
    sim.get_result() = time_series4;
    EXPECT_TRUE(mio::osecir::get_reproduction_number((size_t)0, sim));
}

TEST(Secir, get_mobility_factors)
{
    auto beta                                                                              = 0.25;
    auto max_beta                                                                          = 0.5;
    auto model                                                                             = mio::osecir::Model(1);
    model.parameters.get<mio::osecir::TimeExposed<double>>().array()                       = 3.;
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>().array()            = 2.;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>().array()    = 0.1;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>().array()    = beta;
    model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>().array() = max_beta;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = 100;
    mio::osecir::Simulation<> sim(model, 0.0);
    {
        sim.get_model().parameters.get<mio::osecir::TestAndTraceCapacity<double>>() = 45.;
        auto factors = mio::osecir::get_mobility_factors<double>(sim, 0.0, sim.get_result().get_last_value());
        auto cmp     = Eigen::VectorXd::Ones(Eigen::Index(mio::osecir::InfectionState::Count)).eval();
        cmp[Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] = beta;
        ASSERT_THAT(print_wrap(factors), MatrixNear(cmp));
    }
    {
        sim.get_model().parameters.get<mio::osecir::TestAndTraceCapacity<double>>() = 45. / 5.;
        auto factors = mio::osecir::get_mobility_factors<double>(sim, 0.0, sim.get_result().get_last_value());
        auto cmp     = Eigen::VectorXd::Ones(Eigen::Index(mio::osecir::InfectionState::Count)).eval();
        cmp[Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] = max_beta;
        ASSERT_THAT(print_wrap(factors), MatrixNear(cmp));
    }
    {
        sim.get_model().parameters.get<mio::osecir::TestAndTraceCapacity<double>>() = 20.;
        auto factors = mio::osecir::get_mobility_factors<double>(sim, 0.0, sim.get_result().get_last_value());
        ASSERT_GT(factors[Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)], beta);
        ASSERT_LT(factors[Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)], max_beta);
    }
}

TEST(TestOdeSecir, test_commuters)
{
    auto model                                      = mio::osecir::Model(2);
    auto mobility_factor                            = 0.1;
    auto non_detection_factor                       = 0.4;
    model.parameters.get_start_commuter_detection() = 0.0;
    model.parameters.get_end_commuter_detection()   = 20.0;
    model.parameters.get_commuter_nondetection()    = non_detection_factor;
    auto sim                                        = mio::osecir::Simulation<>(model);
    auto before_testing                             = sim.get_result().get_last_value().eval();
    auto mobile_population                          = (sim.get_result().get_last_value() * mobility_factor).eval();
    auto mobile_population_tested                   = mobile_population.eval();

    mio::osecir::test_commuters<double>(sim, mobile_population_tested, 0.0);

    ASSERT_NEAR(mobile_population_tested[Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)],
                mobile_population[Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] * non_detection_factor,
                1e-5);
    ASSERT_NEAR(
        sim.get_result().get_last_value()[Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptomsConfirmed)],
        before_testing[Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptomsConfirmed)] +
            mobile_population[Eigen::Index(mio::osecir::InfectionState::InfectedNoSymptoms)] *
                (1 - non_detection_factor),
        1e-5);
    ASSERT_NEAR(mobile_population_tested[Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)],
                mobile_population[Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] * non_detection_factor,
                1e-5);
    ASSERT_NEAR(sim.get_result().get_last_value()[Eigen::Index(mio::osecir::InfectionState::InfectedSymptomsConfirmed)],
                before_testing[Eigen::Index(mio::osecir::InfectionState::InfectedSymptomsConfirmed)] +
                    mobile_population[Eigen::Index(mio::osecir::InfectionState::InfectedSymptoms)] *
                        (1 - non_detection_factor),
                1e-5);
}

TEST(TestOdeSecir, check_constraints_parameters)
{
    auto model = mio::osecir::Model<double>(1);
    EXPECT_EQ(model.parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::osecir::Seasonality<double>>(-0.5);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::Seasonality<double>>(0.2);
    model.parameters.set<mio::osecir::ICUCapacity<double>>(-2);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::ICUCapacity<double>>(2);
    model.parameters.set<mio::osecir::TestAndTraceCapacity<double>>(-1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::TestAndTraceCapacity<double>>(1);
    model.parameters.set<mio::osecir::TestAndTraceCapacityMaxRisk<double>>(-1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::TestAndTraceCapacityMaxRisk<double>>(1);
    model.parameters.set<mio::osecir::TimeExposed<double>>(-2);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::TimeExposed<double>>(2);
    model.parameters.set<mio::osecir::TimeInfectedNoSymptoms<double>>(-1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::TimeInfectedNoSymptoms<double>>(5);
    model.parameters.set<mio::osecir::TimeInfectedSymptoms<double>>(0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::TimeInfectedSymptoms<double>>(2);
    model.parameters.set<mio::osecir::TimeInfectedSevere<double>>(-1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::TimeInfectedSevere<double>>(2);
    model.parameters.set<mio::osecir::TimeInfectedCritical<double>>(0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::TimeInfectedCritical<double>>(2);
    model.parameters.set<mio::osecir::TransmissionProbabilityOnContact<double>>(2.0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::TransmissionProbabilityOnContact<double>>(0.5);
    model.parameters.set<mio::osecir::RelativeTransmissionNoSymptoms<double>>(-1.0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::RelativeTransmissionNoSymptoms<double>>(0.5);
    model.parameters.set<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>(3.0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>(0.5);
    model.parameters.set<mio::osecir::RiskOfInfectionFromSymptomatic<double>>(-2.0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::RiskOfInfectionFromSymptomatic<double>>(0.5);
    model.parameters.set<mio::osecir::SeverePerInfectedSymptoms<double>>(-1.0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::SeverePerInfectedSymptoms<double>>(0.5);
    model.parameters.set<mio::osecir::CriticalPerSevere<double>>(-1.0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::CriticalPerSevere<double>>(0.5);
    model.parameters.set<mio::osecir::DeathsPerCritical<double>>(1.1);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::DeathsPerCritical<double>>(1.0);
    model.parameters.set<mio::osecir::DynamicNPIsImplementationDelay<double>>(-4);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osecir::DynamicNPIsImplementationDelay<double>>(3);
    EXPECT_EQ(model.parameters.check_constraints(), 0);
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSecir, apply_constraints_parameters)
{
    auto model             = mio::osecir::Model(1);
    auto indx_agegroup     = mio::AgeGroup(0);
    const double tol_times = 1e-1;

    EXPECT_EQ(model.parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::osecir::Seasonality<double>>(-0.5);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::Seasonality<double>>(), 0);

    model.parameters.set<mio::osecir::ICUCapacity<double>>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::ICUCapacity<double>>(), 0);

    model.parameters.set<mio::osecir::TestAndTraceCapacity<double>>(-1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::TestAndTraceCapacity<double>>(), 0);

    model.parameters.set<mio::osecir::TestAndTraceCapacityMaxRisk<double>>(-1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::TestAndTraceCapacityMaxRisk<double>>(), 0);

    model.parameters.set<mio::osecir::TimeExposed<double>>(-2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::TimeExposed<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecir::TimeInfectedNoSymptoms<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecir::TimeInfectedSymptoms<double>>(1e-8);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecir::TimeInfectedSevere<double>>(-1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecir::TimeInfectedCritical<double>>(0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[indx_agegroup], tol_times);

    model.parameters.set<mio::osecir::TransmissionProbabilityOnContact<double>>(2.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR(model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[indx_agegroup], 0.0,
                1e-14);

    model.parameters.set<mio::osecir::RelativeTransmissionNoSymptoms<double>>(-1.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>(3.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecir::RiskOfInfectionFromSymptomatic<double>>(-2.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecir::SeverePerInfectedSymptoms<double>>(-1.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecir::CriticalPerSevere<double>>(-1.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecir::DeathsPerCritical<double>>(1.1);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[indx_agegroup], 0);

    model.parameters.set<mio::osecir::DynamicNPIsImplementationDelay<double>>(-4);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osecir::DynamicNPIsImplementationDelay<double>>(), 0);

    EXPECT_EQ(model.parameters.apply_constraints(), 0);
    mio::set_log_level(mio::LogLevel::warn);
}

#if defined(MEMILIO_HAS_JSONCPP)

TEST(TestOdeSecir, read_population_data_one_age_group)
{
    std::string path = mio::path_join(TEST_DATA_DIR, "county_current_population.json");
    const std::vector<int> region{1001};
    auto result_one_age_group       = mio::osecir::details::read_population_data(path, region, true).value();
    auto result_multiple_age_groups = mio::osecir::details::read_population_data(path, region, false).value();
    EXPECT_EQ(result_one_age_group.size(), 1);
    EXPECT_EQ(result_one_age_group[0].size(), 1);
    EXPECT_EQ(result_one_age_group[0][0], 90163.0);

    EXPECT_EQ(result_multiple_age_groups.size(), 1);
    EXPECT_EQ(result_multiple_age_groups[0].size(), 6);
    EXPECT_EQ(result_multiple_age_groups[0][0], 3433.0);
}
#if defined(MEMILIO_HAS_HDF5)

class ModelTestOdeSecir : public testing::Test
{

public:
    ModelTestOdeSecir()
        : model(6)
        , num_age_groups(6)
    {
    }
    mio::osecir::Model<double> model;
    size_t num_age_groups;

protected:
    void SetUp() override
    {
        model.parameters.set<mio::osecir::StartDay>(60);
        model.parameters.set<mio::osecir::Seasonality<double>>(0.2);

        for (auto i = mio::AgeGroup(0); i < (mio::AgeGroup)num_age_groups; i++) {
            model.parameters.get<mio::osecir::TimeExposed<double>>()[i]            = 3.2;
            model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i] = 2.0;
            model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]   = 5.8;
            model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[i]     = 9.5;
            model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[i]   = 7.1;

            model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i]  = 0.05;
            model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[i]    = 0.7;
            model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i]    = 0.09;
            model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i]    = 0.25;
            model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[i] = 0.45;
            model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i]         = 0.2;
            model.parameters.get<mio::osecir::CriticalPerSevere<double>>()[i]                 = 0.25;
            model.parameters.get<mio::osecir::DeathsPerCritical<double>>()[i]                 = 0.3;
        }
        model.apply_constraints();
    }
};

// Initialize model with data frrom external reported data. Check if the model is initialized correctly.
TEST_F(ModelTestOdeSecir, read_input_data)
{
    auto model1 = std::vector<mio::osecir::Model<double>>{model};
    auto model2 = std::vector<mio::osecir::Model<double>>{model};
    auto model3 = std::vector<mio::osecir::Model<double>>{model};

    auto read_result1 = mio::osecir::read_input_data_county(
        model1, {2020, 12, 01}, {1002}, std::vector<double>(size_t(num_age_groups), 1.0), 1.0, TEST_DATA_DIR, 10);

    auto read_result2 = mio::osecir::read_input_data(
        model2, {2020, 12, 01}, {1002}, std::vector<double>(size_t(num_age_groups), 1.0), 1.0, TEST_DATA_DIR, 10);

    auto read_result_district =
        mio::osecir::read_input_data(model3, {2020, 12, 01}, {1002}, std::vector<double>(size_t(num_age_groups), 1.0),
                                     1.0, mio::path_join(TEST_DATA_DIR, "pydata/District"), 10);

    EXPECT_THAT(read_result1, IsSuccess());
    EXPECT_THAT(read_result2, IsSuccess());
    EXPECT_THAT(read_result_district, IsSuccess());

    // values were generated by the tested function; can only check stability of the implementation, not correctness
    auto expected_values =
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecir::InfectionState::Count)) << 10280, 2.82575, 1.25589,
         0, 5.85714, 0, 1.42857, 1.33333, 28.2857, 0, 19091.3, 6.59341, 4.23862, 0, 11, 0, 2.54286, 1.33333, 88, 0,
         73821.7, 41.9152, 21.6641, 0, 54.5714, 0, 22, 1.33333, 523.857, 0, 82533.6, 39.5604, 21.0361, 0, 51.4286, 0,
         15.1714, 1.33333, 468.571, 0, 43733.1, 8.32025, 5.18053, 0, 11.4286, 0, 4.22857, 1.33333, 158.571, 8, 15642.7,
         10.5181, 3.2967, 0, 3.28571, 0, 0.657143, 1.33333, 49.8571, 7.57143)
            .finished();

    EXPECT_THAT(print_wrap(model1[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));
    EXPECT_THAT(print_wrap(model2[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));
    EXPECT_THAT(print_wrap(model3[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));
}

// Test if the export of the time series is called not failing and results are the same as the initialization.
TEST_F(ModelTestOdeSecir, export_time_series_init)
{
    TempFileRegister temp_file_register;
    auto tmp_results_dir = temp_file_register.get_unique_path();
    EXPECT_THAT(mio::create_directory(tmp_results_dir), IsSuccess());

    // Test exporting time series
    std::vector<mio::osecir::Model<double>> models{model};
    EXPECT_THAT(mio::osecir::export_input_data_county_timeseries(
                    models, tmp_results_dir, {1002}, {2020, 12, 01}, std::vector<double>(size_t(num_age_groups), 1.0),
                    1.0, 2, mio::path_join(TEST_DATA_DIR, "county_divi_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "cases_all_county_age_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "county_current_population.json")),
                IsSuccess());

    auto data_extrapolated = mio::read_result(mio::path_join(tmp_results_dir, "Results_rki.h5"));
    EXPECT_THAT(data_extrapolated, IsSuccess());

    // Values were generated by the tested function export_input_data_county_timeseries;
    // can only check stability of the implementation, not correctness
    auto expected_results =
        mio::read_result(mio::path_join(TEST_DATA_DIR, "export_time_series_initialization_osecir.h5")).value();

    EXPECT_THAT(print_wrap(data_extrapolated.value()[0].get_groups().matrix()),
                MatrixNear(print_wrap(expected_results[0].get_groups().matrix()), 1e-5, 1e-5));
}

// Test the output of the function for a day way in the past. The model should be initialized with the population data since no Case data is available there.
TEST_F(ModelTestOdeSecir, export_time_series_init_old_date)
{
    mio::set_log_level(mio::LogLevel::off);
    TempFileRegister temp_file_register;
    auto tmp_results_dir = temp_file_register.get_unique_path();
    EXPECT_THAT(mio::create_directory(tmp_results_dir), IsSuccess());

    // Test exporting time series
    std::vector<mio::osecir::Model<double>> models{model};
    EXPECT_THAT(mio::osecir::export_input_data_county_timeseries(
                    models, tmp_results_dir, {1002}, {1000, 12, 01}, std::vector<double>(size_t(num_age_groups), 1.0),
                    1.0, 0, mio::path_join(TEST_DATA_DIR, "county_divi_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "cases_all_county_age_ma7.json"),
                    mio::path_join(TEST_DATA_DIR, "county_current_population.json")),
                IsSuccess());

    auto data_extrapolated = mio::read_result(mio::path_join(tmp_results_dir, "Results_rki.h5"));
    EXPECT_THAT(data_extrapolated, IsSuccess());
    auto results_extrapolated = data_extrapolated.value()[0].get_groups().get_value(0);

    // if we enter an old date, the model only should be initialized with the population data.
    // read population data
    std::string path = mio::path_join(TEST_DATA_DIR, "county_current_population.json");
    const std::vector<int> region{1002};
    auto population_data = mio::osecir::details::read_population_data(path, region, false).value();

    // So, the expected values are the population data in the susceptible compartments and zeros in the other compartments.
    for (size_t i = 0; i < num_age_groups; i++) {
        EXPECT_EQ(results_extrapolated(i * Eigen::Index(mio::osecir::InfectionState::Count)), population_data[0][i]);
    }
    // sum of all compartments should be equal to the population
    EXPECT_EQ(results_extrapolated.sum(), std::accumulate(population_data[0].begin(), population_data[0].end(), 0.0));
    mio::set_log_level(mio::LogLevel::warn);
}

// // Model initialization should return same start values as export time series on that day
TEST_F(ModelTestOdeSecir, model_initialization)
{
    // Vector assignment necessary as read_input_data_county changes model
    auto model_vector = std::vector<mio::osecir::Model<double>>{model};

    EXPECT_THAT(mio::osecir::read_input_data_county(model_vector, {2020, 12, 01}, {1002},
                                                    std::vector<double>(size_t(num_age_groups), 1.0), 1.0,
                                                    TEST_DATA_DIR, 2, false),
                IsSuccess());

    // Values from data/export_time_series_init_osecir.h5, for reading in comparison
    // operator for return of mio::read_result and model population needed.
    auto expected_values =
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecir::InfectionState::Count)) << 10280, 2.82575, 1.25589,
         0, 5.85714, 0, 1.42857, 1.33333, 28.2857, 0, 19091.3, 6.59341, 4.23862, 0, 11, 0, 2.54286, 1.33333, 88, 0,
         73821.7, 41.9152, 21.6641, 0, 54.5714, 0, 22, 1.33333, 523.857, 0, 82533.6, 39.5604, 21.0361, 0, 51.4286, 0,
         15.1714, 1.33333, 468.571, 0, 43733.1, 8.32025, 5.18053, 0, 11.4286, 0, 4.22857, 1.33333, 158.571, 8, 15642.7,
         10.5181, 3.2967, 0, 3.28571, 0, 0.657143, 1.33333, 49.8571, 7.57143)
            .finished();

    EXPECT_THAT(print_wrap(model_vector[0].populations.array().cast<double>()),
                MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));
}

// Calling the model initialization with a date way in the past should only initialize the model with the population data.
TEST_F(ModelTestOdeSecir, model_initialization_old_date)
{
    mio::set_log_level(mio::LogLevel::off);
    // Vector assignment necessary as read_input_data_county changes model
    auto model_vector = std::vector<mio::osecir::Model<double>>{model};

    EXPECT_THAT(mio::osecir::read_input_data_county(model_vector, {1000, 12, 01}, {1002},
                                                    std::vector<double>(size_t(num_age_groups), 1.0), 1.0,
                                                    TEST_DATA_DIR, 0, false),
                IsSuccess());

    // if we enter an old date, the model only should be initialized with the population data.
    // read population data
    std::string path = mio::path_join(TEST_DATA_DIR, "county_current_population.json");
    const std::vector<int> region{1002};
    auto population_data = mio::osecir::details::read_population_data(path, region, false).value();

    // So, the expected values are the population data in the susceptible compartments and zeros in the other compartments.
    auto expected_values =
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecir::InfectionState::Count)) << population_data[0][0], 0,
         0, 0, 0, 0, 0, 0, 0, 0, population_data[0][1], 0, 0, 0, 0, 0, 0, 0, 0, 0, population_data[0][2], 0, 0, 0, 0, 0,
         0, 0, 0, 0, population_data[0][3], 0, 0, 0, 0, 0, 0, 0, 0, 0, population_data[0][4], 0, 0, 0, 0, 0, 0, 0, 0, 0,
         population_data[0][5], 0, 0, 0, 0, 0, 0, 0, 0, 0)
            .finished();

    auto results_extrapolated = model_vector[0].populations.array().cast<double>();

    for (size_t i = 0; i < num_age_groups; i++) {
        EXPECT_EQ(results_extrapolated(i * Eigen::Index(mio::osecir::InfectionState::Count)), population_data[0][i]);
    }
    // sum of all compartments should be equal to the population
    EXPECT_EQ(results_extrapolated.sum(), std::accumulate(population_data[0].begin(), population_data[0].end(), 0.0));
    mio::set_log_level(mio::LogLevel::warn);
}

#endif
#endif
