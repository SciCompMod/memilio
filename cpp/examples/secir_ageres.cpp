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
#include <epidemiology/secir/secir.h>
#include <epidemiology/utils/time_series.h>
#include <epidemiology/utils/logging.h>
#include <epidemiology/model/simulation.h>

int main()
{

    epi::set_log_level(epi::LogLevel::debug);

    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    epi::log_info("Simulating SECIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    double tinc    = 5.2, // R_2^(-1)+R_3^(-1)
        tinfmild   = 6, // 4-14  (=R4^(-1))
        tserint    = 4.2, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        thosp2home = 12, // 7-16 (=R5^(-1))
        thome2hosp = 5, // 2.5-7 (=R6^(-1))
        thosp2icu  = 2, // 1-3.5 (=R7^(-1))
        ticu2home  = 8, // 5-16 (=R8^(-1))
        // tinfasy    = 6.2, // (=R9^(-1)=R_3^(-1)+0.5*R_4^(-1))
        ticu2death = 5; // 3.5-7 (=R5^(-1))

    double cont_freq = 10, // see Polymod study
        inf_prob = 0.05, carr_infec = 0.67,
           alpha = 0.09, // 0.01-0.16
        beta     = 0.25, // 0.05-0.5
        delta    = 0.3, // 0.15-0.77
        rho      = 0.2, // 0.1-0.35
        theta    = 0.25; // 0.15-0.4

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    epi::SecirModel model(3);
    auto nb_groups = model.parameters.get_num_groups();
    double fact    = 1.0 / (double)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<epi::ICUCapacity>(std::numeric_limits<double>::max());
    params.set<epi::StartDay>(0);
    params.set<epi::Seasonality>(0);

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

        params.get<epi::InfectionProbabilityFromContact>()[i] = inf_prob;
        params.get<epi::RelativeCarrierInfectability>()[i] = carr_infec;
        params.get<epi::AsymptoticCasesPerInfectious>()[i] = alpha;
        params.get<epi::RiskOfInfectionFromSympomatic>()[i] = beta;
        params.get<epi::HospitalizedCasesPerInfectious>()[i] = rho;
        params.get<epi::ICUCasesPerHospitalized>()[i] = theta;
        params.get<epi::DeathsPerICU>()[i] = delta;
    }

    epi::ContactMatrixGroup& contact_matrix = params.get<epi::ContactPatterns>();
    contact_matrix[0] = epi::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7), epi::SimulationTime(30.));

    model.apply_constraints();

    epi::TimeSeries<double> secir = simulate(t0, tmax, dt, model);

    char vars[] = {'S', 'E', 'C', 'I', 'H', 'U', 'R', 'D'};
    printf("Number of time points :%d\n", static_cast<int>(secir.get_num_time_points()));
    printf("People in\n");

    for (size_t k = 0; k < (size_t)epi::InfectionState::Count; k++) {
        double dummy = 0;

        for (size_t i = 0; i < (size_t)params.get_num_groups(); i++) {
            printf("\t %c[%d]: %.0f", vars[k], (int)i,
                   secir.get_last_value()[k + (size_t)epi::InfectionState::Count * (int)i]);
            dummy += secir.get_last_value()[k + (size_t)epi::InfectionState::Count * (int)i];
        }

        printf("\t %c_otal: %.0f\n", vars[k], dummy);
    }
}
