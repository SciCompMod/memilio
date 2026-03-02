/*
* Copyright (C) 2020-2026 MEmilio
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
#include "ode_secir/model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0;
    ScalarType tmax = 50;
    ScalarType dt   = 0.1;

    mio::log_info("Simulating SECIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    ScalarType cont_freq = 10; // see Polymod study

    ScalarType nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
               nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model<ScalarType> model(1);

    model.parameters.template set<mio::osecir::StartDay<ScalarType>>(60);
    model.parameters.set<mio::osecir::Seasonality<ScalarType>>(0.2);

    model.parameters.get<mio::osecir::TimeExposed<ScalarType>>()            = 3.2;
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>() = 2.0;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()   = 5.8;
    model.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()     = 9.5;
    model.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()   = 7.1;

    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, cont_freq));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime<ScalarType>(30.));

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

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()  = 0.05;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()    = 0.7;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()    = 0.09;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()    = 0.25;
    model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<ScalarType>>() = 0.45;
    model.parameters.get<mio::osecir::TestAndTraceCapacity<ScalarType>>()              = 35;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()         = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()                 = 0.25;
    model.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()                 = 0.3;
    // The function apply_constraints() ensures that all parameters are within their defined bounds.
    // Note that negative values are set to zero instead of stopping the simulation.
    model.apply_constraints();

    // Using default Integrator
    mio::TimeSeries<ScalarType> secir = mio::osecir::simulate<ScalarType>(t0, tmax, dt, model);

    /*
    Example of using a different integrator
   All available integrators are listed in cpp/memilio/math/README.md

    auto integrator = std::make_unique<mio::RKIntegratorCore>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    mio::TimeSeries<ScalarType> secir = simulate<ScalarType>(t0, tmax, dt, model, std::move(integrator));
    */

    bool print_to_terminal = true;

    if (print_to_terminal) {
        std::vector<std::string> vars = {"S", "E", "C", "C_confirmed", "I", "I_confirmed", "H", "U", "R", "D"};
        printf("\n # t");
        for (size_t k = 0; k < (size_t)mio::osecir::InfectionState::Count; k++) {
            printf(" %s", vars[k].c_str());
        }

        auto num_points = static_cast<size_t>(secir.get_num_time_points());
        for (size_t i = 0; i < num_points; i++) {
            printf("\n%.14f ", secir.get_time(i));
            Eigen::VectorX<ScalarType> res_j = secir.get_value(i);
            for (size_t j = 0; j < (size_t)mio::osecir::InfectionState::Count; j++) {
                printf(" %.14f", res_j[j]);
            }
        }

        Eigen::VectorX<ScalarType> res_j = secir.get_last_value();
        printf("number total: %f",
               res_j[0] + res_j[1] + res_j[2] + res_j[3] + res_j[4] + res_j[5] + res_j[6] + res_j[7]);
    }
}
