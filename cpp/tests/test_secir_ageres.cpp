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
#include "epidemiology/model/simulation.h"
#include <gtest/gtest.h>

TEST(TestSecir, compareAgeResWithSingleRun)
{
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           ticu2death = 5;

    double cont_freq = 0.5, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    epi::SecirModel model(3);
    epi::AgeGroup nb_groups = model.parameters.get_num_groups();
    double fact      = 1.0 / (double)(size_t)nb_groups;

    auto& params = model.parameters;
    for (auto i = epi::AgeGroup(0); i < nb_groups; i++) {
        params.get<epi::IncubationTime>()[i] = tinc;
        params.get<epi::InfectiousTimeMild>()[i] = tinfmild;
        params.get<epi::SerialInterval>()[i] = tserint;
        params.get<epi::HospitalizedToHomeTime>()[i] = thosp2home;
        params.get<epi::HomeToHospitalizedTime>()[i] = thome2hosp;
        params.get<epi::HospitalizedToICUTime>()[i] = thosp2icu;
        params.get<epi::ICUToHomeTime>()[i] = ticu2home;
        params.get<epi::ICUToDeathTime>()[i] = ticu2death;

        model.populations[{i, epi::InfectionState::Exposed}] = fact * nb_exp_t0;
        model.populations[{i, epi::InfectionState::Carrier}] = fact * nb_car_t0;
        model.populations[{i, epi::InfectionState::Infected}] = fact * nb_inf_t0;
        model.populations[{i, epi::InfectionState::Hospitalized}] = fact * nb_hosp_t0;
        model.populations[{i, epi::InfectionState::ICU}] = fact * nb_icu_t0;
        model.populations[{i, epi::InfectionState::Recovered}] = fact * nb_rec_t0;
        model.populations[{i, epi::InfectionState::Dead}] = fact * nb_dead_t0;
        model.populations.set_difference_from_group_total<epi::AgeGroup>({i, epi::InfectionState::Susceptible},
                                                                         fact * nb_total_t0);

        params.get<epi::InfectionProbabilityFromContact>()[i] = 1.;
        params.get<epi::RelativeCarrierInfectability>()[i] = 1.;
        params.get<epi::AsymptoticCasesPerInfectious>()[i] = alpha;
        params.get<epi::RiskOfInfectionFromSympomatic>()[i] = beta;
        params.get<epi::HospitalizedCasesPerInfectious>()[i] = rho;
        params.get<epi::ICUCasesPerHospitalized>()[i] = theta;
        params.get<epi::DeathsPerHospitalized>()[i] = delta;
    }

    params.apply_constraints();

    epi::ContactMatrixGroup& contact_matrix = params.get<epi::ContactPatterns>();
    contact_matrix[0] = epi::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    contact_matrix[0].add_damping(0.7, epi::SimulationTime(30.));

    auto integrator = std::make_shared<epi::RKIntegratorCore>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    epi::TimeSeries<double> secihurd = simulate(t0, tmax, dt, model, integrator);

    // char vars[] = {'S', 'E', 'C', 'I', 'H', 'U', 'R', 'D'};
    // printf("People in\n");
    // for (size_t k = 0; k < epi::InfectionState::SecirCount; k++) {
    //     double dummy = 0;

    //     for (size_t i = 0; i < params.get_num_groups(); i++) {
    //         printf("\t %c[%d]: %.0f", vars[k], (int)i, secir[secir.size() - 1][k]);
    //         dummy += secir[secir.size() - 1][k];
    //     }

    //     printf("\t %c_otal: %.0f\n", vars[k], dummy);
    // }

    auto compare = load_test_data_csv<double>("secihurd-compare.csv");

    ASSERT_EQ(compare.size(), static_cast<size_t>(secihurd.get_num_time_points()));
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size() - 1, secihurd.get_num_elements() / (size_t)nb_groups) << "at row " << i;
        ASSERT_NEAR(secihurd.get_time(i), compare[i][0], 1e-10) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            double dummy = 0;
            for (size_t k = 0; k < (size_t)nb_groups; k++) {
                dummy += secihurd.get_value(i)[j - 1 + k * (size_t)epi::InfectionState::Count];
            }
            EXPECT_NEAR(dummy, compare[i][j], 1e-10) << " at row " << i;
        }
    }
}
