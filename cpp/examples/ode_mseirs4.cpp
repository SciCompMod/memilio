/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#include "ode_mseirs4/model.h"
#include "memilio/compartments/simulation.h"
#include <iostream>

int main()
{
    using FP = double;
    mio::omseirs4::Model<FP> model;
    auto& params = model.parameters;

    // Example parameter values for day-based time unit (t in days)
    params.get<mio::omseirs4::BaseTransmissionRate<FP>>()     = 0.4; // b0 per day
    params.get<mio::omseirs4::SeasonalAmplitude<FP>>()        = 0.15; // b1
    params.get<mio::omseirs4::SeasonalPhase<FP>>()            = 0.0; // phi; for phase shift use 2*pi*offsetDays/365
    params.get<mio::omseirs4::NaturalBirthDeathRate<FP>>()    = 1.0 / (70.0 * 365.0); // mu per day
    params.get<mio::omseirs4::LossMaternalImmunityRate<FP>>() = 1.0 / 90.0; // xi per day (~3 months)
    params.get<mio::omseirs4::ProgressionRate<FP>>()          = 1.0 / 7.0; // sigma per day (≈7 days latent)
    params.get<mio::omseirs4::RecoveryRate<FP>>()             = 1.0 / 14.0; // nu per day (≈14 days infectious)
    params.get<mio::omseirs4::ImmunityWaningRate<FP>>()       = 1.0 / (5.0 * 365.0); // gamma per day (5 years)
    // factors default already set (0.5, 0.35, 0.25) as in the paper https://doi.org/10.1016/S0025-5564(01)00066-9.

    // Initial population (absolute counts), set all compartments explicitly.
    double N = 1e6;

    // Initial populations.
    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::MaternalImmune)}] =
        5000.0;

    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::E1)}] = 300.0;
    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::E2)}] = 150.0;
    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::E3)}] = 80.0;
    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::E4)}] = 70.0;

    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::I1)}] = 200.0;
    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::I2)}] = 100.0;
    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::I3)}] = 50.0;
    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::I4)}] = 50.0;

    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::R1)}] = 40000.0;
    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::R2)}] = 30000.0;
    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::R3)}] = 20000.0;
    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::R4)}] = 10000.0;

    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::S2)}] = 100000.0;
    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::S3)}] = 50000.0;
    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::S4)}] = 50000.0;

    // Compute S1 as residual to match N
    double assigned = 5000.0 + (300.0 + 150.0 + 80.0 + 70.0) + (200.0 + 100.0 + 50.0 + 50.0) +
                      (40000.0 + 30000.0 + 20000.0 + 10000.0) + (100000.0 + 50000.0 + 50000.0);
    double S1 = N - assigned;
    if (S1 < 0)
        S1 = 0;
    model.populations[{mio::Index<mio::omseirs4::InfectionState>(mio::omseirs4::InfectionState::S1)}] = S1;

    model.check_constraints();

    // simulate
    double t0   = 0.0;
    double tmax = 20.0; // days
    double dt   = 1.0; // daily output
    auto result = mio::simulate(t0, tmax, dt, model);

    // print header
    std::cout << "t M S1 S2 S3 S4 E1 E2 E3 E4 I1 I2 I3 I4 R1 R2 R3 R4\n";
    for (size_t i = 0; i < (size_t)result.get_num_time_points(); ++i) {
        std::cout << result.get_time(i);
        const auto& y = result.get_value(i);
        for (size_t k = 0; k < (size_t)mio::omseirs4::InfectionState::Count; ++k) {
            std::cout << ' ' << y[(Eigen::Index)k];
        }
        std::cout << '\n';
    }
    return 0;
}
