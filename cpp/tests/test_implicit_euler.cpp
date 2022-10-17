/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Alexander Ruettgers, Martin J. Kuehn
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
#include "secir/implicit_euler.h"
#include <gtest/gtest.h>

TEST(TestImplicitEuler, compareOneTimeStep)
{
    double t0 = 0;
    double dt = 0.1;

    // working_params
    double tinc    = 5.2, // R_2^(-1)+R_3^(-1)
        tinfmild   = 6, // 4-14  (=R4^(-1))
        tserint    = 4.2, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        thosp2home = 12, // 7-16 (=R5^(-1))
        thome2hosp = 5, // 2.5-7 (=R6^(-1))
        thosp2icu  = 2, // 1-3.5 (=R7^(-1))
        ticu2home  = 8, // 5-16 (=R8^(-1))
        tinfasy    = 6.2, // (=R9^(-1)=R_3^(-1)+0.5*R_4^(-1))
        ticu2death = 5; // 3.5-7 (=R5^(-1))

    double cont_freq = 0.5, // 0.2-0.75
        alpha        = 0.09, // 0.01-0.16
        beta         = 0.25, // 0.05-0.5
        delta        = 0.3, // 0.15-0.77
        rho          = 0.2, // 0.1-0.35
        theta        = 0.25; // 0.15-0.4

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::SecirModel model(1);
    auto& params = model.parameters;

    params.get<mio::IncubationTime>()[(mio::AgeGroup)0]             = tinc;
    params.get<mio::InfectiousTimeMild>()[(mio::AgeGroup)0]         = tinfmild;
    params.get<mio::SerialInterval>()[(mio::AgeGroup)0]             = tserint;
    params.get<mio::HospitalizedToHomeTime>()[(mio::AgeGroup)0]     = thosp2home;
    params.get<mio::HomeToHospitalizedTime>()[(mio::AgeGroup)0]     = thome2hosp;
    params.get<mio::HospitalizedToICUTime>()[(mio::AgeGroup)0]      = thosp2icu;
    params.get<mio::ICUToHomeTime>()[(mio::AgeGroup)0]              = ticu2home;
    params.get<mio::InfectiousTimeAsymptomatic>()[(mio::AgeGroup)0] = tinfasy;
    params.get<mio::ICUToDeathTime>()[(mio::AgeGroup)0]             = ticu2death;

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix.add_damping(0.7, mio::SimulationTime(30.));

    model.populations.set_total(nb_total_t0);
    model.populations[{mio::AgeGroup(0), mio::InfectionState::Exposed}]      = nb_exp_t0;
    model.populations[{mio::AgeGroup(0), mio::InfectionState::Carrier}]      = nb_car_t0;
    model.populations[{mio::AgeGroup(0), mio::InfectionState::Infected}]     = nb_inf_t0;
    model.populations[{mio::AgeGroup(0), mio::InfectionState::Hospitalized}] = nb_hosp_t0;
    model.populations[{mio::AgeGroup(0), mio::InfectionState::ICU}]          = nb_icu_t0;
    model.populations[{mio::AgeGroup(0), mio::InfectionState::Recovered}]    = nb_rec_t0;
    model.populations[{mio::AgeGroup(0), mio::InfectionState::Dead}]         = nb_dead_t0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::InfectionState::Susceptible}, nb_total_t0);

    params.get<mio::InfectionProbabilityFromContact>()[(mio::AgeGroup)0] = 1.;
    params.get<mio::AsymptomaticCasesPerInfectious>()[(mio::AgeGroup)0]    = alpha;
    params.get<mio::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = beta;
    params.get<mio::HospitalizedCasesPerInfectious>()[(mio::AgeGroup)0]  = rho;
    params.get<mio::ICUCasesPerHospitalized>()[(mio::AgeGroup)0]         = theta;
    params.get<mio::DeathsPerICU>()[(mio::AgeGroup)0]                    = delta;

    mio::DerivFunction dummy_f; // only required for explicit time integrators

    auto y0 = model.populations.get_compartments();
    Eigen::VectorXd y1(y0.size()); // solution at time t=\Delta t=0.1

    mio::ImplicitEulerIntegratorCore(model).step(dummy_f, y0, t0, dt,
                                                 y1); // just one iteration of implicit Euler scheme

    EXPECT_NEAR(y1[0], 9756.897564185243, 1e-10);
    EXPECT_NEAR(y1[1], 99.97811957794602, 1e-10);
    EXPECT_NEAR(y1[2], 50.74190209182218, 1e-10);
    EXPECT_NEAR(y1[3], 51.41751953982101, 1e-10);
    EXPECT_NEAR(y1[4], 19.833786579788253, 1e-10);
    EXPECT_NEAR(y1[5], 10.098962633404634, 1e-10);
    EXPECT_NEAR(y1[6], 10.97155161617429, 1e-10);
    EXPECT_NEAR(y1[7], 0.0605937758004278, 1e-10);
}
