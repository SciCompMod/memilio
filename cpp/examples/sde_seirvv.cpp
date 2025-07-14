/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Nils Wassmuth, Rene Schmieding, Martin J. Kuehn
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
#include "memilio/utils/logging.h"
#include "memilio/utils/uncertain_value.h"
#include "sde_seirvv/model.h"
#include "sde_seirvv/simulation.h"

#include <vector>

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0.;
    ScalarType tmid = 1.;
    ScalarType tmax = 3.;
    ScalarType dt   = 0.1;

    mio::log_info("Simulating SEIRVV; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::sseirvv::Model<ScalarType> model;

    ScalarType total_population = 180000;

    model.populations[{mio::sseirvv::InfectionState::ExposedV1}]     = 0;
    model.populations[{mio::sseirvv::InfectionState::ExposedV2}]     = 0;
    model.populations[{mio::sseirvv::InfectionState::InfectedV1}]    = 7200;
    model.populations[{mio::sseirvv::InfectionState::InfectedV2}]    = 0;
    model.populations[{mio::sseirvv::InfectionState::RecoveredV1}]   = 0;
    model.populations[{mio::sseirvv::InfectionState::RecoveredV2}]   = 0;
    model.populations[{mio::sseirvv::InfectionState::ExposedV1V2}]   = 0;
    model.populations[{mio::sseirvv::InfectionState::InfectedV1V2}]  = 0;
    model.populations[{mio::sseirvv::InfectionState::RecoveredV1V2}] = 0;
    model.populations[{mio::sseirvv::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::sseirvv::InfectionState::ExposedV1}] -
        model.populations[{mio::sseirvv::InfectionState::ExposedV2}] -
        model.populations[{mio::sseirvv::InfectionState::InfectedV1}] -
        model.populations[{mio::sseirvv::InfectionState::InfectedV2}] -
        model.populations[{mio::sseirvv::InfectionState::RecoveredV1}] -
        model.populations[{mio::sseirvv::InfectionState::RecoveredV2}] -
        model.populations[{mio::sseirvv::InfectionState::ExposedV1V2}] -
        model.populations[{mio::sseirvv::InfectionState::InfectedV1V2}] -
        model.populations[{mio::sseirvv::InfectionState::RecoveredV1V2}];

    // It is assumed that both variants have the same transmission probability
    // on contact and the same time exposed. The time infected is scaled by
    // 1.35 for the second variant.
    model.parameters.get<mio::sseirvv::ContactPatterns<ScalarType>>().get_baseline()(0, 0) = 1;
    model.parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV1<ScalarType>>(0.076);
    model.parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV2<ScalarType>>(0.076);
    model.parameters.set<mio::sseirvv::TimeExposedV1<ScalarType>>(5.33);
    model.parameters.set<mio::sseirvv::TimeExposedV2<ScalarType>>(5.33);
    model.parameters.set<mio::sseirvv::TimeInfectedV1<ScalarType>>(17.2);
    model.parameters.set<mio::sseirvv::TimeInfectedV2<ScalarType>>(17.2 * 1.35);

    model.check_constraints();

    // Simulate the model up until tmid, with only the first variant.
    auto sseirv = mio::sseirvv::simulate<ScalarType>(t0, tmid, dt, model);
    // Set the model population to the simulation result, so it is used as initial value for the second simulation.
    model.populations.array() = sseirv.get_last_value().cast<mio::UncertainValue<ScalarType>>();
    // The second variant enters with 100 individuals. This increases the model population to total_population + 100.
    model.populations[{mio::sseirvv::InfectionState::InfectedV2}] = 100;
    // Simulate the model from tmid to tmax, now with both variants.
    auto sseirv2 = mio::sseirvv::simulate<ScalarType>(tmid, tmax, dt, model);
    sseirv.print_table({"Susceptible", "ExposedV1", "InfectedV1", "RecoveredV1", "ExposedV2", "InfectedV2",
                        "RecoveredV2", "ExposedV1V2", "InfectedV1V2", "RecoveredV1V2"});
    sseirv2.print_table({"Susceptible", "ExposedV1", "InfectedV1", "RecoveredV1", "ExposedV2", "InfectedV2",
                         "RecoveredV2", "ExposedV1V2", "InfectedV1V2", "RecoveredV1V2"});
}
