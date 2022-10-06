/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Wadim Koslow
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
#include "memilio/compartments/simulation.h"
#include "secir/secir.h"
#include "memilio/utils/time_series.h"
#include "memilio/io/result_io.h"
#include "temp_file_register.h"
#include <gtest/gtest.h>

TEST(TestSaveResult, compareResultWithH5)
{
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double cont_freq = 10, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::SecirModel model(1);
    auto& params            = model.parameters;
    mio::AgeGroup nb_groups = params.get_num_groups();
    ;

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::IncubationTime>()[i]             = 5.2;
        params.get<mio::TimeInfectedSymptoms>()[i]         = 5.;
        params.get<mio::SerialInterval>()[i]             = 4.2;
        params.get<mio::TimeInfectedSevere>()[i]     = 10.;
        params.get<mio::TimeInfectedCritical>()[i]              = 8.;

        model.populations[{i, mio::InfectionState::Exposed}]      = nb_exp_t0;
        model.populations[{i, mio::InfectionState::Carrier}]      = nb_car_t0;
        model.populations[{i, mio::InfectionState::Infected}]     = nb_inf_t0;
        model.populations[{i, mio::InfectionState::Hospitalized}] = nb_hosp_t0;
        model.populations[{i, mio::InfectionState::ICU}]          = nb_icu_t0;
        model.populations[{i, mio::InfectionState::Recovered}]    = nb_rec_t0;
        model.populations[{i, mio::InfectionState::Dead}]         = nb_dead_t0;
        model.populations.set_difference_from_total({i, mio::InfectionState::Susceptible}, nb_total_t0);

        params.get<mio::TransmissionProbabilityOnContact>()[i] = 0.06;
        params.get<mio::RelativeTransmissionNoSymptoms>()[i]    = 0.67;
        params.get<mio::RecoveredPerInfectedNoSymptoms>()[i]    = alpha;
        params.get<mio::RiskOfInfectionFromSymptomatic>()[i]   = beta;
        params.get<mio::SeverePerInfectedSymptoms>()[i]  = rho;
        params.get<mio::CriticalPerSevere>()[i]         = theta;
        params.get<mio::DeathsPerCritical>()[i]                    = delta;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::ContactPatterns>();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    auto result_from_sim                                  = simulate(t0, tmax, dt, model);
    std::vector<mio::TimeSeries<double>> results_from_sim = {result_from_sim, result_from_sim};
    std::vector<int> ids                                  = {1, 2};

    TempFileRegister file_register;
    auto results_file_path  = file_register.get_unique_path("test_result-%%%%-%%%%.h5");
    auto save_result_status = mio::save_result(results_from_sim, ids, (int)(size_t)nb_groups, results_file_path);
    ASSERT_TRUE(save_result_status);

    auto results_from_file = mio::read_result(results_file_path);
    ASSERT_TRUE(results_from_file);
    auto result_from_file = results_from_file.value()[0];

    ASSERT_EQ(result_from_file.get_groups().get_num_time_points(), result_from_sim.get_num_time_points());
    ASSERT_EQ(result_from_file.get_totals().get_num_time_points(), result_from_sim.get_num_time_points());
    for (Eigen::Index i = 0; i < result_from_sim.get_num_time_points(); i++) {
        ASSERT_EQ(result_from_file.get_groups().get_num_elements(), result_from_sim.get_num_elements())
            << "at row " << i;
        ASSERT_EQ(result_from_file.get_totals().get_num_elements(),
                  result_from_sim.get_num_elements() / static_cast<Eigen::Index>((size_t)nb_groups))
            << "at row " << i;
        ASSERT_NEAR(result_from_sim.get_time(i), result_from_file.get_groups().get_time(i), 1e-10) << "at row " << i;
        ASSERT_NEAR(result_from_sim.get_time(i), result_from_file.get_totals().get_time(i), 1e-10) << "at row " << i;
        for (Eigen::Index l = 0; l < result_from_file.get_totals().get_num_elements(); l++) {
            double total = 0.0;
            for (Eigen::Index j = 0; j < Eigen::Index((size_t)nb_groups); j++) {
                total += result_from_sim[i][j * (size_t)mio::InfectionState::Count + l];
                EXPECT_NEAR(result_from_file.get_groups()[i][j * (size_t)mio::InfectionState::Count + l],
                            result_from_sim[i][j * (size_t)mio::InfectionState::Count + l], 1e-10)
                    << " at row " << i << " at row " << l << " at Group " << j;
            }
            EXPECT_NEAR(result_from_file.get_totals()[i][l], total, 1e-10) << " at row " << i << " at row " << l;
        }
    }
}
