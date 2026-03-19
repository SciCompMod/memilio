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
#include "ode_meningitis/model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/data/analyze_result.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0;
    ScalarType tmax = 10;
    ScalarType dt   = 0.1;

    mio::log_info("Simulating meningitis model; t={} ... {} with dt = {}.", t0, tmax, dt);

    ScalarType cont_freq = 10; // set in line with the transmission risk below

    ScalarType nb_total_t0 = 10000;
    ScalarType nb_inf_t0   = 100;

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
    model.populations[{mio::AgeGroup(0), mio::omeng::InfectionState::Incoming}]        = 0;
    model.populations[{mio::AgeGroup(0), mio::omeng::InfectionState::SusceptibleHigh}] = 0;
    model.populations[{mio::AgeGroup(0), mio::omeng::InfectionState::Carrier}]         = 0;
    model.populations[{mio::AgeGroup(0), mio::omeng::InfectionState::Infected}]        = nb_inf_t0;
    model.populations[{mio::AgeGroup(0), mio::omeng::InfectionState::Recovered}]       = 0;
    model.populations[{mio::AgeGroup(0), mio::omeng::InfectionState::Dead}]            = 0;
    model.populations[{mio::AgeGroup(0), mio::omeng::InfectionState::DeadNatural}]     = 0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::omeng::InfectionState::SusceptibleLow},
                                                nb_total_t0 - nb_inf_t0);

    model.check_constraints();

    // Using default Integrator
    auto result = mio::omeng::simulate_flows<ScalarType>(t0, tmax, dt, model);

    /*
    Example of using a different integrator
    All available integrators are listed in cpp/memilio/math/README.md

    auto integrator = std::make_unique<mio::RKIntegratorCore>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    auto result = mio::omeng::simulate_flows<ScalarType>(t0, tmax, dt, model, std::move(integrator));
    */

    // Note: The "Incoming" compartment (Inc) acts as a pure source node for recruitment.
    // Since the FlowModel framework only supports conservative flows between compartments, there is
    // no external inflow into Incoming. Its value therefore becomes increasingly negative over time.
    // This is expected model behaviour: Inc accumulates the total number of recruits with a negative sign
    // and does NOT represent a real population. The living population is S_L + S_H + C + I + R.

    // simulate_flows additionally records all inter-compartment flows at each time step.
    // The return value is a pair: [0] = compartment TimeSeries, [1] = flow TimeSeries.
    printf("Compartments:\n");
    auto result_interpolated = mio::interpolate_simulation_result(result[0]);
    result_interpolated.print_table({"Inc", "S_L", "S_H", "C", "I", "R", "D_D", "D_N"});

    // auto result_flows_interp = mio::interpolate_simulation_result(result[1]);
    // printf("\nFlows:\n");
    // result_flows_interp.print_table({"Inc->S_H", "Inc->S_L", "S_H->C", "S_L->C", "C->I", "C->R", "I->R", "I->D_D",
    //                                  "R->S_L", "R->S_H", "S_H->D_N", "S_L->D_N", "C->D_N", "I->D_N", "R->D_N"});
}
