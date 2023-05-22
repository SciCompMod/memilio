/*
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Khoa Nguyen
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
#include "abm/abm.h"
#include "graph_abm/graph_world.h"
#include "abm/household.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/meta_mobility_instant.h"
#include <cstdio>

int main()
{
    // Set global infection parameters (similar to infection parameters in SECIR model) and initialize the world
    mio::abm::GlobalInfectionParameters infection_params;

    // Set same infection parameter for all age groups. For example, the incubation period is 4 days.
    infection_params.get<mio::abm::IncubationPeriod>() = 4.;

    // Create the worlds with infection parameters.
    auto world1 = mio::graph_abm::GraphWorld(infection_params);
    auto world2 = mio::graph_abm::GraphWorld(infection_params);

    //add households for world 1
    auto child = mio::abm::HouseholdMember(); // A child is 50/50% 0-4 or 5-14.
    child.set_age_weight(mio::abm::AgeGroup::Age0to4, 1);
    child.set_age_weight(mio::abm::AgeGroup::Age5to14, 1);

    auto parent = mio::abm::HouseholdMember(); // A parent is 50/50% 15-34 or 35-59.
    parent.set_age_weight(mio::abm::AgeGroup::Age15to34, 1);
    parent.set_age_weight(mio::abm::AgeGroup::Age35to59, 1);

    auto two_person_household_group1 = mio::abm::HouseholdGroup();
    auto two_person_household_group2 = mio::abm::HouseholdGroup();
    auto two_person_family           = mio::abm::Household();
    two_person_family.add_members(child, 1);
    two_person_family.add_members(parent, 1);
    auto two_person_adults = mio::abm::Household();
    two_person_adults.add_members(parent, 2);
    two_person_household_group1.add_households(two_person_family, 2);
    two_person_household_group1.add_households(two_person_adults, 1);
    two_person_household_group2.add_households(two_person_family, 1);
    two_person_household_group2.add_households(two_person_adults, 3);

    auto three_person_household_group1 = mio::abm::HouseholdGroup();
    auto three_person_household_group2 = mio::abm::HouseholdGroup();
    auto three_person_family           = mio::abm::Household();
    three_person_family.add_members(child, 1);
    three_person_family.add_members(parent, 2);
    three_person_household_group1.add_households(three_person_family, 2);
    three_person_household_group2.add_households(three_person_family, 1);

    //add households to world 1
    add_household_group_to_world(world1, two_person_household_group1);
    add_household_group_to_world(world1, three_person_household_group1);
    //add households to world 2
    add_household_group_to_world(world2, two_person_household_group2);
    add_household_group_to_world(world2, three_person_household_group2);

    //add internal locations for both worlds
    //auto event =

    auto t0 = mio::abm::TimePoint(0);

    mio::Graph<mio::SimulationNode<mio::abm::Simulation>, mio::MigrationEdge> g;
    g.add_node(1, t0, std::move(world1));
    g.add_node(2, t0, std::move(world2));

    // The results are saved in a table with 9 rows.
    // The first row is t = time, the others correspond to the number of people with a certain infection state at this time:
    // S = Susceptible, E = Exposed, C = Carrier, I = Infected, I_s = Infected_Severe,
    // I_c = Infected_Critical, R_C = Recovered_Carrier, R_I = Recovered_Infected, D = Dead
    // auto f_abm = fopen("abm_minimal.txt", "w");
    // fprintf(f_abm, "# t S E C I I_s I_c R_C R_I D\n");
    // for (auto i = 0; i < sim.get_result().get_num_time_points(); ++i) {
    //     fprintf(f_abm, "%f ", sim.get_result().get_time(i));
    //     auto v = sim.get_result().get_value(i);
    //     for (auto j = 0; j < v.size(); ++j) {
    //         fprintf(f_abm, "%f", v[j]);
    //         if (j < v.size() - 1) {
    //             fprintf(f_abm, " ");
    //         }
    //     }
    //     if (i < sim.get_result().get_num_time_points() - 1) {
    //         fprintf(f_abm, "\n");
    //     }
    // }
    // fclose(f_abm);
}
