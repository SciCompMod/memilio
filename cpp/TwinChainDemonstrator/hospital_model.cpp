/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Julia Bicker
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

#include <iostream>
#include "abm/model.h"
#include "abm/simulation.h"
#include "abm/common_abm_loggers.h"

extern "C" {
#include "julia_init.h"
}
#include "JuliaTestPkg.h"

int main()
{
    auto model = mio::abm::Model(1);
    auto home = model.add_location(mio::abm::LocationType::Home);
    auto pid = model.add_person(home, mio::AgeGroup(0), mio::abm::ActivityType::Home);
    auto& person = model.get_person(pid);
    person.set_assigned_location(mio::abm::ActivityType::Home, home, model.get_id());
    auto sim = mio::abm::Simulation(mio::abm::TimePoint(0), std::move(model));
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)};
    sim.advance(mio::abm::TimePoint(1), historyTimeSeries);

    // init_julia(0, nullptr);

    // std::cout << "Test\n";
    // printint(39);
    // std::cout << "mean = " << mean_of_two(2, 7) << "\n";

    // shutdown_julia(0);
    return 0;
}
