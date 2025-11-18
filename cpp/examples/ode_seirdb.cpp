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

#include "ode_seirdb/model.h"
#include "ode_seirdb/infection_state.h"
#include "ode_seirdb/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0;
    ScalarType tmax = 30.;
    ScalarType dt   = 1.0;

    mio::log_info("Simulating ODE SEIRDB; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::oseirdb::Model<ScalarType> model(1);

    ScalarType total_population                                                    = 10000;
    model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Exposed}]   = 100;
    model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Infected}]  = 100;
    model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Recovered}] = 100;
    model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Dead}]      = 100;
    model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Buried}]    = 100;
    model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Exposed}] -
        model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Recovered}] -
        model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Dead}] -
        model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Buried}];

    model.parameters.set<mio::oseirdb::TimeExposed<ScalarType>>(5.2);
    model.parameters.set<mio::oseirdb::TimeInfected<ScalarType>>(6.0);
    model.parameters.set<mio::oseirdb::TimeToBurial<ScalarType>>(5.5);
    model.parameters.set<mio::oseirdb::TransmissionProbabilityOnContact<ScalarType>>(0.1);
    model.parameters.set<mio::oseirdb::TransmissionProbabilityFromDead<ScalarType>>(0.01);
    model.parameters.set<mio::oseirdb::ProbabilityToRecover<ScalarType>>(0.75);

    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.get<mio::oseirdb::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(2.7);
    contact_matrix[0].add_damping(0.7, mio::SimulationTime<ScalarType>(10.));

    model.check_constraints();

    auto result = mio::simulate<ScalarType>(t0, tmax, dt, model);

    result.print_table({"S", "E", "I", "R", "D", "B"});
    std::cout << "\nnumber total: " << result.get_last_value().sum() << "\n";
}
