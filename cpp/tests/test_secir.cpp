/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "load_test_data.h"
#include "epidemiology/secir/secir.h"
#include "epidemiology/secir/analyze_result.h"
#include <distributions_helpers.h>
#include <gtest/gtest.h>

TEST(TestSecir, compareWithPreviousRun)
{
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           ticu2death = 5;

    double cont_freq = 10, inf_prob = 0.05, carr_infec = 1, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2,
           theta = 0.25;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    epi::SecirModel model(1);

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    model.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0] = tinc;
    model.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0] = tinfmild;
    model.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0] = tserint;
    model.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0] = thosp2home;
    model.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0] = thome2hosp;
    model.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0] = thosp2icu;
    model.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0] = ticu2home;
    model.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0] = ticu2death;

    epi::ContactMatrixGroup& contact_matrix = model.parameters.get<epi::ContactPatterns>();
    contact_matrix[0]                       = epi::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix[0].add_damping(0.7, epi::SimulationTime(30.));

    model.populations.set_total(nb_total_t0);
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}] = nb_exp_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}] = nb_car_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}] = nb_inf_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}] = nb_hosp_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}] = nb_icu_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}] = nb_rec_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}] = nb_dead_t0;
    model.populations.set_difference_from_total({epi::AgeGroup(0), epi::InfectionState::Susceptible},
                                                nb_total_t0);

    model.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0] = inf_prob;
    model.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0] = carr_infec;
    model.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0] = alpha;
    model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0] = beta;
    model.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0] = rho;
    model.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0] = theta;
    model.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0] = delta;

    model.apply_constraints();

    auto integrator = std::make_shared<epi::RKIntegratorCore>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    epi::TimeSeries<double> secihurd = simulate(t0, tmax, dt, model, integrator);

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

TEST(TestSecir, testParamConstructors)
{
    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           ticu2death = 5;

    double cont_freq = 10, inf_prob = 0.05, carr_infec = 0.67, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2,
           theta = 0.24;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 54, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 11, nb_dead_t0 = 0;

    double icu_cap   = 4444;
    double start_day = 30, seasonality = 0.3;

    epi::SecirModel model(1);

    model.parameters.set<epi::ICUCapacity>(icu_cap);

    model.parameters.set<epi::StartDay>(start_day);
    model.parameters.set<epi::Seasonality>(seasonality);

    model.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0] = tinc;
    model.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0] = tinfmild;
    model.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0] = tserint;
    model.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0] = thosp2home;
    model.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0] = thome2hosp;
    model.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0] = thosp2icu;
    model.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0] = ticu2home;
    model.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0] = ticu2death;

    model.populations.set_total(nb_total_t0);
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}] = nb_exp_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}] = nb_car_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}] = nb_inf_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}] = nb_hosp_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}] = nb_icu_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}] = nb_rec_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}] = nb_dead_t0;
    model.populations.set_difference_from_total({epi::AgeGroup(0), epi::InfectionState::Susceptible},
                                                nb_total_t0);

    model.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0] = inf_prob;
    model.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0] = carr_infec;
    model.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0] = alpha;
    model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0] = beta;
    model.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0] = rho;
    model.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0] = theta;
    model.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0] = delta;

    epi::ContactMatrixGroup& contact_matrix = model.parameters.get<epi::ContactPatterns>();
    contact_matrix[0]                       = epi::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix[0].add_damping(0.7, epi::SimulationTime(30.));

    epi::SecirModel model2{model}; // copy constructor

    EXPECT_EQ(model.parameters.get<epi::ICUCapacity>(), model2.parameters.get<epi::ICUCapacity>());
    EXPECT_EQ(model.parameters.get<epi::StartDay>(), model2.parameters.get<epi::StartDay>());
    EXPECT_EQ(model.parameters.get<epi::Seasonality>(), model2.parameters.get<epi::Seasonality>());

    EXPECT_EQ(model.populations.get_total(), model2.populations.get_total());
    EXPECT_EQ((model.populations[{epi::AgeGroup(0), epi::InfectionState::Susceptible}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Susceptible}]));
    EXPECT_EQ((model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}]));
    EXPECT_EQ((model.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}]));
    EXPECT_EQ((model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}]));
    EXPECT_EQ((model.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}]));
    EXPECT_EQ((model.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}]));
    EXPECT_EQ((model.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}]));
    EXPECT_EQ((model.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}]));

    EXPECT_EQ(model.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0], model2.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0], model2.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0], model2.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::InfectiousTimeAsymptomatic>()[(epi::AgeGroup)0], model2.parameters.get<epi::InfectiousTimeAsymptomatic>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0],
              model2.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0],
              model2.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0],
              model2.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0], model2.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0], model2.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::ContactPatterns>().get_cont_freq_mat(),
              model2.parameters.get<epi::ContactPatterns>().get_cont_freq_mat());

    EXPECT_EQ(model.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0],
              model2.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0],
              model2.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0],
              model2.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0],
              model2.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0],
              model2.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0],
              model2.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0],
              model2.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model.parameters.get<epi::ContactPatterns>().get_cont_freq_mat(),
              model2.parameters.get<epi::ContactPatterns>().get_cont_freq_mat());

    epi::SecirModel model3 = std::move(model2); // move constructor

    EXPECT_EQ(model.parameters.get<epi::ICUCapacity>(), model3.parameters.get<epi::ICUCapacity>());
    EXPECT_EQ(model.parameters.get<epi::StartDay>(), model3.parameters.get<epi::StartDay>());
    EXPECT_EQ(model.parameters.get<epi::Seasonality>(), model3.parameters.get<epi::Seasonality>());

    EXPECT_EQ(model3.populations.get_total(), model.populations.get_total());
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Susceptible}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Susceptible}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}]),
              (model.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}]));

    EXPECT_EQ(model3.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0], model.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0], model.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0], model.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::InfectiousTimeAsymptomatic>()[(epi::AgeGroup)0], model.parameters.get<epi::InfectiousTimeAsymptomatic>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0],
              model.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0],
              model.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0],
              model.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0], model.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0], model.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0]);

    EXPECT_EQ(model3.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0],
              model.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0],
              model.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0],
              model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0],
              model.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0],
              model.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0],
              model.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0],
              model.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0]);

    EXPECT_EQ(model.parameters.get<epi::ContactPatterns>().get_cont_freq_mat(),
              model3.parameters.get<epi::ContactPatterns>().get_cont_freq_mat());

    epi::SecirModel model4 = model3; // copy assignment constructor

    EXPECT_EQ(model4.parameters.get<epi::ICUCapacity>(), model3.parameters.get<epi::ICUCapacity>());
    EXPECT_EQ(model4.parameters.get<epi::StartDay>(), model3.parameters.get<epi::StartDay>());
    EXPECT_EQ(model4.parameters.get<epi::Seasonality>(), model3.parameters.get<epi::Seasonality>());

    EXPECT_EQ(model3.populations.get_total(), model4.populations.get_total());
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Susceptible}]),
              (model4.populations[{epi::AgeGroup(0), epi::InfectionState::Susceptible}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}]),
              (model4.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}]),
              (model4.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}]),
              (model4.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}]),
              (model4.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}]),
              (model4.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}]),
              (model4.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}]));
    EXPECT_EQ((model3.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}]),
              (model4.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}]));

    EXPECT_EQ(model3.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0], model4.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0], model4.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0], model4.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::InfectiousTimeAsymptomatic>()[(epi::AgeGroup)0], model4.parameters.get<epi::InfectiousTimeAsymptomatic>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0],
              model4.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0],
              model4.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0],
              model4.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0], model4.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0], model4.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0]);

    EXPECT_EQ(model3.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0],
              model4.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0],
              model4.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0],
              model4.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0],
              model4.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0],
              model4.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0],
              model4.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model3.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0],
              model4.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0]);

    EXPECT_EQ(model4.parameters.get<epi::ContactPatterns>().get_cont_freq_mat(),
              model3.parameters.get<epi::ContactPatterns>().get_cont_freq_mat());

    epi::SecirModel model5 = std::move(model4); // move assignment constructor

    EXPECT_EQ(model5.parameters.get<epi::ICUCapacity>(), model3.parameters.get<epi::ICUCapacity>());
    EXPECT_EQ(model5.parameters.get<epi::StartDay>(), model3.parameters.get<epi::StartDay>());
    EXPECT_EQ(model5.parameters.get<epi::Seasonality>(), model3.parameters.get<epi::Seasonality>());

    EXPECT_EQ(model5.populations.get_total(), model3.populations.get_total());
    EXPECT_EQ((model5.populations[{epi::AgeGroup(0), epi::InfectionState::Susceptible}]),
              (model3.populations[{epi::AgeGroup(0), epi::InfectionState::Susceptible}]));
    EXPECT_EQ((model5.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}]),
              (model3.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}]));
    EXPECT_EQ((model5.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}]),
              (model3.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}]));
    EXPECT_EQ((model5.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}]),
              (model3.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}]));
    EXPECT_EQ((model5.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}]),
              (model3.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}]));
    EXPECT_EQ((model5.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}]),
              (model3.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}]));
    EXPECT_EQ((model5.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}]),
              (model3.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}]));
    EXPECT_EQ((model5.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}]),
              (model3.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}]));

    EXPECT_EQ(model5.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0], model3.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0], model3.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0], model3.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::InfectiousTimeAsymptomatic>()[(epi::AgeGroup)0], model3.parameters.get<epi::InfectiousTimeAsymptomatic>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0],
              model3.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0],
              model3.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0],
              model3.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0], model3.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0], model3.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0]);

    EXPECT_EQ(model5.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0],
              model3.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0],
              model3.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0],
              model3.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0],
              model3.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0],
              model3.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0],
              model3.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0]);
    EXPECT_EQ(model5.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0],
              model3.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0]);

    EXPECT_EQ(model5.parameters.get<epi::ContactPatterns>().get_cont_freq_mat(),
              model3.parameters.get<epi::ContactPatterns>().get_cont_freq_mat());
}

TEST(TestSecir, testSettersAndGetters)
{
    std::vector<epi::UncertainValue> vec;

    for (int i = 0; i < 26; i++) {
        epi::UncertainValue val = epi::UncertainValue(i);
        val.set_distribution(epi::ParameterDistributionNormal(i, 10 * i, 5 * i, i / 10.0));
        vec.push_back(val);
    }

    epi::SecirModel model(1);

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    EXPECT_EQ(model.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0].get_distribution().get(), nullptr);

    model.parameters.set<epi::ICUCapacity>(vec[0]);

    model.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0] = vec[1];
    model.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0] = vec[2];
    model.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0] = vec[3];
    model.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0] = vec[4];
    model.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0] = vec[5];
    model.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0] = vec[6];
    model.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0] = vec[7];
    model.parameters.get<epi::InfectiousTimeAsymptomatic>()[(epi::AgeGroup)0] = vec[8];
    model.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0] = vec[9];

    model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}] = vec[10];
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}] = vec[11];
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}] = vec[12];
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}] = vec[13];
    model.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}] = vec[14];
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}] = vec[15];
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}] = vec[16];

    model.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0] = vec[17];
    model.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0] = vec[18];
    model.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0] = vec[19];
    model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0] = vec[20];
    model.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0] = vec[21];
    model.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0] = vec[22];
    model.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0] = vec[23];

    EXPECT_NE(model.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0].get_distribution().get(), nullptr);

    check_distribution(*vec[0].get_distribution(), *model.parameters.get<epi::ICUCapacity>().get_distribution());

    model.parameters.set<epi::StartDay>(vec[24]);
    model.parameters.set<epi::Seasonality>(vec[25]);

    EXPECT_NE(model.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0].get_distribution().get(), nullptr);

    check_distribution(*vec[1].get_distribution(), *model.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[2].get_distribution(), *model.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[3].get_distribution(), *model.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[4].get_distribution(),
                       *model.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[5].get_distribution(),
                       *model.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[6].get_distribution(),
                       *model.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[7].get_distribution(), *model.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[8].get_distribution(),
                       *model.parameters.get<epi::InfectiousTimeAsymptomatic>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[9].get_distribution(), *model.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[10].get_distribution(),
                       *model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}].get_distribution());
    check_distribution(*vec[11].get_distribution(),
                       *model.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}].get_distribution());
    check_distribution(*vec[12].get_distribution(),
                       *model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}].get_distribution());
    check_distribution(*vec[13].get_distribution(),
                       *model.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}].get_distribution());
    check_distribution(*vec[14].get_distribution(),
                       *model.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}].get_distribution());
    check_distribution(*vec[15].get_distribution(),
                       *model.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}].get_distribution());
    check_distribution(*vec[16].get_distribution(),
                       *model.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}].get_distribution());
    check_distribution(*vec[17].get_distribution(),
                       *model.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[18].get_distribution(),
                       *model.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[19].get_distribution(),
                       *model.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[20].get_distribution(),
                       *model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[21].get_distribution(),
                       *model.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[22].get_distribution(),
                       *model.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0].get_distribution());
    check_distribution(*vec[23].get_distribution(),
                       *model.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0].get_distribution());
    // no dist for start day
    check_distribution(*vec[25].get_distribution(), *model.parameters.get<epi::Seasonality>().get_distribution());

    EXPECT_EQ(vec[0], model.parameters.get<epi::ICUCapacity>());
    EXPECT_EQ(vec[1], model.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[2], model.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[3], model.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[4], model.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[5], model.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[6], model.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[7], model.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[8], model.parameters.get<epi::InfectiousTimeAsymptomatic>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[9], model.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[10], (model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}]));
    EXPECT_EQ(vec[11], (model.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}]));
    EXPECT_EQ(vec[12], (model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}]));
    EXPECT_EQ(vec[13], (model.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}]));
    EXPECT_EQ(vec[14], (model.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}]));
    EXPECT_EQ(vec[15], (model.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}]));
    EXPECT_EQ(vec[16], (model.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}]));
    EXPECT_EQ(vec[17], model.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[18], model.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[19], model.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[20], model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[21], model.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[22], model.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[23], model.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0]);
    EXPECT_EQ(vec[24], model.parameters.get<epi::StartDay>());
    EXPECT_EQ(vec[25], model.parameters.get<epi::Seasonality>());
}

TEST(TestSecir, testValueConstraints)
{
    double tinc    = 5.1, // R_2^(-1)+R_3^(-1)
        tinfmild   = 5.86642, // 4-14  (=R4^(-1))
        tserint    = 5.08993, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        thosp2home = 11.6138, // 7-16 (=R5^(-1))
        thome2hosp = 4.45361, // 2.5-7 (=R6^(-1))
        thosp2icu  = 2.15791, // 1-3.5 (=R7^(-1))
        ticu2home  = 9.16291, // 5-16 (=R8^(-1))
        ticu2death = 5.90264; // 3.5-7 (=R5^(-1))

    double cont_freq = 10, // 0.2-0.75
        inf_prob = 0.064519, carr_infec = 0.56758,
           alpha = 2.124921, // 0.01-0.16
        beta     = 0.190609, // 0.05-0.5
        delta    = 0.245801, // 0.15-0.77
        rho      = 0.183693, // 0.1-0.35
        theta    = 0.185556; // 0.15-0.4

    double nb_total_t0 = 10000, nb_exp_t0 = -91, nb_inf_t0 = 39, nb_car_t0 = 36, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 8, nb_dead_t0 = 0;

    epi::SecirModel model(1);

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    model.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0] = tinc;
    model.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0] = tinfmild;
    model.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0] = tserint;
    model.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0] = thosp2home;
    model.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0] = thome2hosp;
    model.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0] = thosp2icu;
    model.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0] = ticu2home;
    model.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0] = ticu2death;

    epi::ContactMatrixGroup& contact_matrix = model.parameters.get<epi::ContactPatterns>();
    contact_matrix[0]                       = epi::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix[0].add_damping(0.7, epi::SimulationTime(30.));

    model.populations.set_total(nb_total_t0);
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}] = nb_exp_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}] = nb_car_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}] = nb_inf_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}] = nb_hosp_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}] = nb_icu_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}] = nb_rec_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}] = nb_dead_t0;
    model.populations.set_difference_from_total({epi::AgeGroup(0), epi::InfectionState::Susceptible},
                                                nb_total_t0);

    model.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0] = inf_prob;
    model.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0] = carr_infec;
    model.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0] = alpha;
    model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0] = beta;
    model.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0] = rho;
    model.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0] = theta;
    model.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0] = delta;

    epi::set_log_level(epi::LogLevel::off);
    model.parameters.check_constraints();
    epi::set_log_level(epi::LogLevel::warn);

    EXPECT_EQ(-91, (model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}]));
    EXPECT_EQ(2.124921, model.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0].value());
    EXPECT_NEAR(5.08993, model.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0], 1e-14);

    model.apply_constraints();

    EXPECT_EQ(0.0, (model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}]));
    EXPECT_EQ(0.0, model.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0].value());
    EXPECT_NEAR(4.6, model.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0], 1e-14);
}

TEST(TestSecir, testModelConstraints)
{
    double t0   = 0;
    double tmax = 57; // after 57 days with cont_freq 10 and winter, the virus would already decline
    double dt   = 0.1;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           ticu2death = 5;

    double cont_freq = 10, inf_prob = 0.05, carr_infec = 1, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2,
           theta = 0.25;

    double nb_total_t0 = 1000000, nb_exp_t0 = 10000, nb_inf_t0 = 5000, nb_car_t0 = 500, nb_hosp_t0 = 20, nb_icu_t0 = 0,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    epi::SecirModel model(1);

    model.parameters.get<epi::IncubationTime>()[(epi::AgeGroup)0] = tinc;
    model.parameters.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0] = tinfmild;
    model.parameters.get<epi::SerialInterval>()[(epi::AgeGroup)0] = tserint;
    model.parameters.get<epi::HospitalizedToHomeTime>()[(epi::AgeGroup)0] = thosp2home;
    model.parameters.get<epi::HomeToHospitalizedTime>()[(epi::AgeGroup)0] = thome2hosp;
    model.parameters.get<epi::HospitalizedToICUTime>()[(epi::AgeGroup)0] = thosp2icu;
    model.parameters.get<epi::ICUToHomeTime>()[(epi::AgeGroup)0] = ticu2home;
    model.parameters.get<epi::ICUToDeathTime>()[(epi::AgeGroup)0] = ticu2death;

    model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}] = nb_exp_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}] = nb_car_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}] = nb_inf_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}] = nb_hosp_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}] = nb_icu_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}] = nb_rec_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}] = nb_dead_t0;
    model.populations.set_difference_from_total({epi::AgeGroup(0), epi::InfectionState::Susceptible},
                                                nb_total_t0);

    model.parameters.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0] = inf_prob;
    model.parameters.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0] = carr_infec;
    model.parameters.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0] = alpha;
    model.parameters.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0] = beta;
    model.parameters.get<epi::HospitalizedCasesPerInfectious>()[(epi::AgeGroup)0] = rho;
    model.parameters.get<epi::ICUCasesPerHospitalized>()[(epi::AgeGroup)0] = theta;
    model.parameters.get<epi::DeathsPerHospitalized>()[(epi::AgeGroup)0] = delta;

    epi::ContactMatrixGroup& contact_matrix = model.parameters.get<epi::ContactPatterns>();
    contact_matrix[0]                           = epi::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));

    model.apply_constraints();

    epi::TimeSeries<double> secihurd = simulate(t0, tmax, dt, model);
    double max_icu_cap               = 0;
    for (Eigen::Index i = 0; i < secihurd.get_num_time_points(); i++) {
        if (secihurd.get_value(i)[5] > max_icu_cap) {
            max_icu_cap = secihurd.get_value(i)[5];
        }
    }

    epi::TimeSeries<double> secihurd_interp = epi::interpolate_simulation_result(secihurd);

    model.parameters.set<epi::StartDay>(100);
    model.parameters.set<epi::Seasonality>(0.5);

    epi::TimeSeries<double> secihurd_season        = simulate(t0, tmax, dt, model);
    epi::TimeSeries<double> secihurd_season_interp = epi::interpolate_simulation_result(secihurd_season);

    for (Eigen::Index i = 0; i < secihurd_interp.get_num_time_points(); i++) {
        EXPECT_LE(secihurd_season_interp.get_value(i)[3], secihurd_interp.get_value(i)[3]) << " at row " << i;
    }

    model.parameters.set<epi::StartDay>(280);

    epi::TimeSeries<double> secihurd_season2        = simulate(t0, tmax, dt, model);
    epi::TimeSeries<double> secihurd_season2_interp = epi::interpolate_simulation_result(secihurd_season2);

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
        model.parameters.set<epi::ICUCapacity>(8000);
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
    double tinc = 5.2, tinfmild = 6, tserint = 4.2;

    double cont_freq = 10, inf_prob = 0.05, carr_infec = 1, alpha = 0.09, beta = 0.25;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50;

    epi::SecirModel model(1);
    auto& params = model.parameters;

    params.get<epi::IncubationTime>()[(epi::AgeGroup)0] = tinc;
    params.get<epi::InfectiousTimeMild>()[(epi::AgeGroup)0] = tinfmild;
    params.get<epi::SerialInterval>()[(epi::AgeGroup)0] = tserint;

    epi::ContactMatrixGroup& contact_matrix = params.get<epi::ContactPatterns>();
    contact_matrix[0]                       = epi::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));

    model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}] = nb_exp_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}] = nb_car_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}] = nb_inf_t0;
    model.populations.set_difference_from_total({epi::AgeGroup(0), epi::InfectionState::Susceptible},
                                                nb_total_t0);

    params.get<epi::InfectionProbabilityFromContact>()[(epi::AgeGroup)0] = inf_prob;
    params.get<epi::RelativeCarrierInfectability>()[(epi::AgeGroup)0] = carr_infec;
    params.get<epi::AsymptoticCasesPerInfectious>()[(epi::AgeGroup)0] = alpha;
    params.get<epi::RiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0] = beta;

    params.apply_constraints();

    auto y = model.populations.get_compartments();

    auto dydt_default = Eigen::VectorXd(Eigen::Index(epi::InfectionState::Count));
    model.get_derivatives(y, y, 0, dydt_default);

    params.set<epi::TestAndTraceCapacity>(50);
    params.get<epi::MaxRiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0] = beta * 3;
    auto dydt_under_capacity = Eigen::VectorXd(Eigen::Index(epi::InfectionState::Count));
    model.get_derivatives(y, y, 0, dydt_under_capacity);

    params.set<epi::TestAndTraceCapacity>(10);
    params.get<epi::MaxRiskOfInfectionFromSympomatic>()[(epi::AgeGroup)0] = beta * 3;
    auto dydt_over_capacity = Eigen::VectorXd(Eigen::Index(epi::InfectionState::Count));
    model.get_derivatives(y, y, 0, dydt_over_capacity);

    EXPECT_DOUBLE_EQ(dydt_under_capacity[(size_t)epi::InfectionState::Exposed], dydt_default[(size_t)epi::InfectionState::Exposed]);
    EXPECT_GT(dydt_over_capacity[(size_t)epi::InfectionState::Exposed], dydt_default[(size_t)epi::InfectionState::Exposed]);
}

TEST(Secir, getInfectionsRelative)
{
    size_t num_groups = 3;
    epi::SecirModel model((int)num_groups);
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}] = 100.0;
    model.populations.set_difference_from_group_total<epi::AgeGroup>({epi::AgeGroup(0), epi::InfectionState::Susceptible}, 10'000.0);
    model.populations[{epi::AgeGroup(1), epi::InfectionState::Infected}] = 50.0;
    model.populations.set_difference_from_group_total<epi::AgeGroup>({epi::AgeGroup(1), epi::InfectionState::Susceptible}, 20'000.0);
    model.populations[{epi::AgeGroup(2), epi::InfectionState::Infected}] = 25.0;
    model.populations.set_difference_from_group_total<epi::AgeGroup>({epi::AgeGroup(2), epi::InfectionState::Susceptible}, 40'000.0);

    epi::SecirSimulation<> sim(model, 0.0);
    ASSERT_EQ(epi::get_infections_relative(sim, 0.0, sim.get_result().get_last_value()), (100. + 50. + 25.) / (10'000 + 20'000 + 40'000));
}
