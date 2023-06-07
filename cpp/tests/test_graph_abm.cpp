/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Julia Bicker, Martin J. Kuehn
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
#include "models/graph_abm/graph_world.h"
#include "models/abm/world.h"
#include "models/abm/person.h"
#include "abm_helpers.h"
#include <gtest/gtest.h>

TEST(TestGraphWorld, test_find_location)
{
    mio::abm::GlobalInfectionParameters infection_params;
    //auto world  = mio::abm::World(infection_params, 0);
    auto world1 = mio::graph_abm::GraphWorld(infection_params, 0);
    auto world2 = mio::graph_abm::GraphWorld(infection_params, 1);
    auto loc1   = world1.add_location(mio::abm::LocationType::Home);
    auto loc2   = world2.add_location(mio::abm::LocationType::Work);
    auto& loc   = world1.get_individualized_location(loc1);
    mio::unused(loc);
    // auto person = mio::abm::Person(loc, mio::abm::AgeGroup::Age35to59, 0, 0);
    // person.set_assigned_location(loc1);
    // person.set_assigned_location(loc2);
    // auto& found_loc1 = world1.find_location(mio::abm::LocationType::Home, person);
    // auto& found_loc2 = world1.find_location(mio::abm::LocationType::Work, person);
    // ASSERT_EQ(found_loc1, loc1);
    // ASSERT_EQ(found_loc2, loc2);
}