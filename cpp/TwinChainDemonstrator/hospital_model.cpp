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
    // Set severe to crittical to zero so no one goes to ICU.
    model.parameters.get<mio::abm::CriticalPerInfectedSevere>() = 0.0;
    // let all agents go to work and nobody goes to school
    model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({mio::AgeGroup(0)}, true);
    model.parameters.get<mio::abm::AgeGroupGotoSchool>().set_multiple({mio::AgeGroup(0)}, false);
    auto home1    = model.add_location(mio::abm::LocationType::Home);
    auto home2    = model.add_location(mio::abm::LocationType::Home);
    auto hospital = model.add_location(mio::abm::LocationType::Hospital);
    auto office   = model.add_location(mio::abm::LocationType::Office);
    // Add 10 persons and assign them to homes
    for (int i = 0; i < 10; ++i) {
        // Assign have of the persons to home1 and the other half to home2
        auto pid     = model.add_person(i < 5 ? home1 : home2, mio::AgeGroup(0), mio::abm::ActivityType::Home);
        auto& person = model.get_person(pid);
        person.set_assigned_location(mio::abm::ActivityType::Home, i < 5 ? home1 : home2, model.get_id());
        // Assign all persons the hospital as their hospital location
        person.set_assigned_location(mio::abm::ActivityType::Hospital, hospital, model.get_id());
        // Assign every second person to the office as work place and the others to the hospital as work place
        person.set_assigned_location(mio::abm::ActivityType::Work, i % 2 == 0 ? office : hospital, model.get_id());
    }
    auto sim = mio::abm::Simulation(mio::abm::TimePoint(0), std::move(model));
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)};
#ifdef MEMILIO_USE_JULIA_LIB
    init_julia(0, nullptr); // once, before the simulation
#endif

    sim.advance(mio::abm::TimePoint(0) + mio::abm::days(5), historyTimeSeries);
#ifdef MEMILIO_USE_JULIA_LIB
    shutdown_julia(0); // once, after the simulation
#endif

    return 0;
}
