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
#include "secir/secir.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    mio::log_info("Simulating SECIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    double cont_freq = 10; // see Polymod study

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::SecirModel model(1);

    model.parameters.set<mio::StartDay>(60);
    model.parameters.set<mio::Seasonality>(0.2);

    model.parameters.get<mio::IncubationTime>()       = 5.2;
    model.parameters.get<mio::TimeInfectedSymptoms>() = 5.8;
    model.parameters.get<mio::SerialInterval>()       = 4.2;
    model.parameters.get<mio::TimeInfectedSevere>()   = 9.5;
    model.parameters.get<mio::TimeInfectedCritical>() = 7.1;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    model.populations.set_total(nb_total_t0);
    model.populations[{mio::AgeGroup(0), mio::InfectionState::Exposed}]            = nb_exp_t0;
    model.populations[{mio::AgeGroup(0), mio::InfectionState::InfectedNoSymptoms}] = nb_car_t0;
    model.populations[{mio::AgeGroup(0), mio::InfectionState::InfectedSymptoms}]   = nb_inf_t0;
    model.populations[{mio::AgeGroup(0), mio::InfectionState::InfectedSevere}]     = nb_hosp_t0;
    model.populations[{mio::AgeGroup(0), mio::InfectionState::InfectedCritical}]   = nb_icu_t0;
    model.populations[{mio::AgeGroup(0), mio::InfectionState::Recovered}]          = nb_rec_t0;
    model.populations[{mio::AgeGroup(0), mio::InfectionState::Dead}]               = nb_dead_t0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::InfectionState::Susceptible}, nb_total_t0);

    model.parameters.get<mio::TransmissionProbabilityOnContact>()  = 0.05;
    model.parameters.get<mio::RelativeTransmissionNoSymptoms>()    = 0.7;
    model.parameters.get<mio::RecoveredPerInfectedNoSymptoms>()    = 0.09;
    model.parameters.get<mio::RiskOfInfectionFromSymptomatic>()    = 0.25;
    model.parameters.get<mio::MaxRiskOfInfectionFromSymptomatic>() = 0.45;
    model.parameters.get<mio::TestAndTraceCapacity>()              = 35;
    model.parameters.get<mio::SeverePerInfectedSymptoms>()         = 0.2;
    model.parameters.get<mio::CriticalPerSevere>()                 = 0.25;
    model.parameters.get<mio::DeathsPerCritical>()                 = 0.3;

    model.apply_constraints();

    auto integrator = std::make_shared<mio::RKIntegratorCore>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    mio::TimeSeries<double> secir = simulate(t0, tmax, dt, model, integrator);

    bool print_to_terminal = true;

    if (print_to_terminal) {
        char vars[] = {'S', 'E', 'C', 'I', 'H', 'U', 'R', 'D'};
        printf("\n # t");
        for (size_t k = 0; k < (size_t)mio::InfectionState::Count; k++) {
            printf(" %c", vars[k]);
        }
        auto num_points = static_cast<size_t>(secir.get_num_time_points());
        for (size_t i = 0; i < num_points; i++) {
            printf("\n%.14f ", secir.get_time(i));
            Eigen::VectorXd res_j = secir.get_value(i);
            for (size_t j = 0; j < (size_t)mio::InfectionState::Count; j++) {
                printf(" %.14f", res_j[j]);
            }
        }

        Eigen::VectorXd res_j = secir.get_last_value();
        printf("number total: %f",
               res_j[0] + res_j[1] + res_j[2] + res_j[3] + res_j[4] + res_j[5] + res_j[6] + res_j[7]);
    }
}
