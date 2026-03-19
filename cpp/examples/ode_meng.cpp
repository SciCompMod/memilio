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

    mio::log_info("Simulating meningitis model; t={} ... {} with dt = {}.", t0, tmax, dt);

    ScalarType cont_freq = 10; // set in line with the transmission risk below

    ScalarType nb_total_t0 = 10000;

    mio::omeng::Model<ScalarType> model(1);

    model.parameters.get<mio::omeng::RateCarrierToInfected<ScalarType>>()           = 0.00022; // sigma
    model.parameters.get<mio::omeng::RateCarrierToRecovered<ScalarType>>()          = 0.8; // eta_2
    model.parameters.get<mio::omeng::RateInfectedToRecovered<ScalarType>>()         = 0.43; // eta_1
    model.parameters.get<mio::omeng::RateInfectedToDead<ScalarType>>()              = 0.495; // d
    model.parameters.get<mio::omeng::RateNaturalDeath<ScalarType>>()                = 0.0152207; // mu
    model.parameters.get<mio::omeng::RateImmunityLoss<ScalarType>>()                = 0.851; // xi
    model.parameters.get<mio::omeng::ProbabilityImmunityLossSusLow<ScalarType>>()   = 0.6997; //  theta
    model.parameters.get<mio::omeng::ModificationRate<ScalarType>>()                = 0.23; //  a
    model.parameters.get<mio::omeng::RiskOfInfectionFromFromCarrier<ScalarType>>()  = 0.742; // omega
    model.parameters.get<mio::omeng::RiskOfInfectionFromFromInfected<ScalarType>>() = 0.425; // 0 to 0.85; omega1
    model.parameters.get<mio::omeng::TransmissionProbabilityOnContact<ScalarType>>() =
        0.05; // set in line with contact frequency = 10
    model.parameters.get<mio::omeng::IncomeFractionSusLow<ScalarType>>() = 0.585; //  Delta
    model.parameters.get<mio::omeng::IncomeRate<ScalarType>>()           = 19787; // Pi

    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.get<mio::omeng::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, cont_freq));
    // contact_matrix[0].add_damping(0.7, mio::SimulationTime<ScalarType>(30.));

    model.populations.set_total(nb_total_t0);
    model.populations[{mio::AgeGroup(0), mio::omeng::InfectionState::SusceptibleHigh}] = 0;
    model.populations[{mio::AgeGroup(0), mio::omeng::InfectionState::Carrier}]         = 0;
    model.populations[{mio::AgeGroup(0), mio::omeng::InfectionState::Infected}]        = 0;
    model.populations[{mio::AgeGroup(0), mio::omeng::InfectionState::Recovered}]       = 0;
    model.populations[{mio::AgeGroup(0), mio::omeng::InfectionState::Dead}]            = 0;
    model.populations[{mio::AgeGroup(0), mio::omeng::InfectionState::DeadNatural}]     = 0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::omeng::InfectionState::SusceptibleLow},
                                                nb_total_t0);

    model.check_constraints();

    // Using default Integrator
    mio::TimeSeries<ScalarType> secir = mio::omeng::simulate<ScalarType>(t0, tmax, dt, model);

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
        std::vector<std::string> vars = {"S_H", "S_L", "C", "I", "R", "D_D", "D_N"};
        printf("\n # t");
        for (size_t k = 0; k < (size_t)mio::omeng::InfectionState::Count; k++) {
            printf(" %s", vars[k].c_str());
        }

        auto num_points = static_cast<size_t>(secir.get_num_time_points());
        for (size_t i = 0; i < num_points; i++) {
            printf("\n%.14f ", secir.get_time(i));
            Eigen::VectorX<ScalarType> res_j = secir.get_value(i);
            for (size_t j = 0; j < (size_t)mio::omeng::InfectionState::Count; j++) {
                printf(" %.14f", res_j[j]);
            }
        }

        Eigen::VectorX<ScalarType> res_j = secir.get_last_value();
        printf("number total: %f", res_j[0] + res_j[1] + res_j[2] + res_j[3] + res_j[4] + res_j[5] + res_j[6]);
    }
}
