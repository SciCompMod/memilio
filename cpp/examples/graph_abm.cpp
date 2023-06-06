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
#include "models/abm/abm.h"
#include "graph_abm/graph_world.h"
#include "models/abm/household.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include <cstdio>

int main()
{
    // Set global infection parameters (similar to infection parameters in SECIR model) and initialize the world
    mio::abm::GlobalInfectionParameters infection_params;

    // Set same infection parameter for all age groups. For example, the incubation period is 4 days.
    infection_params.get<mio::abm::IncubationPeriod>() = 4.;

    // Create the worlds with infection parameters.
    auto world1 = mio::graph_abm::GraphWorld(infection_params, 0);
    auto world2 = mio::graph_abm::GraphWorld(infection_params, 1);

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
    //world1
    std::vector<mio::abm::LocationId> events_w1;
    std::vector<mio::abm::LocationId> works_w1;
    //add two social events
    auto event1_w1 = world1.add_location(mio::abm::LocationType::SocialEvent);
    world1.get_individualized_location(event1_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto event2_w1 = world1.add_location(mio::abm::LocationType::SocialEvent);
    world1.get_individualized_location(event2_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    events_w1.push_back(event1_w1);
    events_w1.push_back(event2_w1);
    // Add hospital and ICU with 5 maximum contacs.
    auto hospital = world1.add_location(mio::abm::LocationType::Hospital);
    world1.get_individualized_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto icu = world1.add_location(mio::abm::LocationType::ICU);
    world1.get_individualized_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    // Add one supermarket, maximum constacts are assumed to be 20.
    auto shop_w1 = world1.add_location(mio::abm::LocationType::BasicsShop);
    world1.get_individualized_location(shop_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    // At every school, the maximum contacts are 20.
    auto school = world1.add_location(mio::abm::LocationType::School);
    world1.get_individualized_location(school).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    // At every workplace, maximum contacts are 10.
    auto work1_w1 = world1.add_location(mio::abm::LocationType::Work);
    world1.get_individualized_location(work1_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    auto work2_w1 = world1.add_location(mio::abm::LocationType::Work);
    world1.get_individualized_location(work2_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    works_w1.push_back(work1_w1);
    works_w1.push_back(work2_w1);

    //world2
    std::vector<mio::abm::LocationId> events_w2;
    std::vector<mio::abm::LocationId> works_w2;
    //add one social event
    auto event_w2 = world2.add_location(mio::abm::LocationType::SocialEvent);
    world2.get_individualized_location(event_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    events_w2.push_back(event_w2);
    // Add one supermarket, maximum constacts are assumed to be 20.
    auto shop_w2 = world2.add_location(mio::abm::LocationType::BasicsShop);
    world2.get_individualized_location(shop_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    // At every workplace, maximum contacts are 10.
    auto work_w2 = world2.add_location(mio::abm::LocationType::Work);
    world2.get_individualized_location(work_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    works_w2.push_back(work_w2);

    //add external locations
    //world 1
    works_w1.push_back(work_w2);

    //world 2
    //TODO: MaximumContacts for external locations setzen
    events_w2.push_back(event2_w1);
    works_w2.push_back(work2_w1);

    auto start_date = mio::abm::TimePoint(0);

    //assign infection states
    //world 1
    auto persons_w1 = world1.get_persons();
    for (auto& person : persons_w1) {
        mio::abm::InfectionState infection_state =
            (mio::abm::InfectionState)(rand() % ((uint32_t)mio::abm::InfectionState::Count - 1));
        if (infection_state != mio::abm::InfectionState::Susceptible)
            person.add_new_infection(mio::abm::Infection(mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         world1.get_global_infection_parameters(), start_date,
                                                         infection_state));
    }
    //world 2
    auto persons_w2 = world2.get_persons();
    for (auto& person : persons_w2) {
        mio::abm::InfectionState infection_state =
            (mio::abm::InfectionState)(rand() % ((uint32_t)mio::abm::InfectionState::Count - 1));
        if (infection_state != mio::abm::InfectionState::Susceptible)
            person.add_new_infection(mio::abm::Infection(mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         world2.get_global_infection_parameters(), start_date,
                                                         infection_state));
    }

    //assign locations
    //word 1
    for (auto& person : persons_w1) {
        //choose one of the events with same probability and assign to person
        auto event =
            events_w1[mio::DiscreteDistribution<size_t>::get_instance()(std::vector<double>(events_w1.size(), 1.0))];
        person.set_assigned_location(event);
        //assign shop
        person.set_assigned_location(shop_w1);
        //assign hospital and ICU
        person.set_assigned_location(hospital);
        person.set_assigned_location(icu);
        //assign work and school dependent on person's age
        if (person.get_age() == mio::abm::AgeGroup::Age5to14) {
            person.set_assigned_location(school);
        }
        if (person.get_age() == mio::abm::AgeGroup::Age15to34 || person.get_age() == mio::abm::AgeGroup::Age35to59) {
            auto work =
                works_w1[mio::DiscreteDistribution<size_t>::get_instance()(std::vector<double>(works_w1.size(), 1.0))];
            person.set_assigned_location(work);
        }
    }
    //world 2
    for (auto& person : persons_w2) {
        //choose one of the events with same probability and assign to person
        auto event =
            events_w2[mio::DiscreteDistribution<size_t>::get_instance()(std::vector<double>(events_w2.size(), 1.0))];
        person.set_assigned_location(event);
        //assign shop
        person.set_assigned_location(shop_w2);
        //assign hospital and ICU
        person.set_assigned_location(hospital);
        person.set_assigned_location(icu);
        //assign work and school dependent on person's age
        if (person.get_age() == mio::abm::AgeGroup::Age5to14) {
            person.set_assigned_location(school);
        }
        if (person.get_age() == mio::abm::AgeGroup::Age15to34 || person.get_age() == mio::abm::AgeGroup::Age35to59) {
            auto work =
                works_w2[mio::DiscreteDistribution<size_t>::get_instance()(std::vector<double>(works_w2.size(), 1.0))];
            person.set_assigned_location(work);
        }
    }

    auto t0   = mio::abm::TimePoint(0);
    auto dt   = mio::abm::hours(12);
    auto tmax = mio::abm::TimePoint(0) + mio::abm::days(10);

    mio::Graph<mio::SimulationNode<mio::graph_abm::GraphSimulation>, mio::MigrationEdgeABM> g;
    g.add_node(0, t0, std::move(world1));
    g.add_node(1, t0, std::move(world2));
    g.add_edge(0, 1);
    g.add_edge(1, 0);

    auto sim = mio::make_migration_sim(t0, dt, std::move(g));

    sim.advance(tmax, true);
    for (auto& n : sim.get_graph().nodes()) {
        std::cout << "node " << n.id << "\n";
        std::cout << "\n";

        std::cout << "# t S E C I I_s I_c R_C R_I D\n";
        for (auto i = 0; i < n.property.get_result().get_num_time_points(); ++i) {
            std::cout << n.property.get_result().get_time(i) << " ";
            auto v = n.property.get_result().get_value(i);
            for (auto j = 0; j < v.size(); ++j) {
                std::cout << v[j] << " ";
                if (j < v.size() - 1) {
                    std::cout << " ";
                }
            }
            if (i < n.property.get_result().get_num_time_points() - 1) {
                std::cout << "\n";
            }
        }
        std::cout << "\n";
    }
}
