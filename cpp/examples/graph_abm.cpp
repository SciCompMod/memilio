/* 
* Copyright (C) 2020-2024 MEmilio
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

#include "abm/household.h"
#include "abm/model.h"
#include "abm/infection_state.h"
#include "abm/location_type.h"
#include "abm/time.h"
#include "abm/person_id.h"
#include "graph_abm/graph_abm_mobility.h"
#include "graph_abm/graph_abmodel.h"
#include "memilio/io/history.h"
#include "memilio/mobility/graph.h"
#include <cstddef>
#include <cstdint>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

//Logger
struct Logger : mio::LogAlways {
    /**
    * A vector of tuples with the Location information i.e. each tuple contains the following information:
    * - The LocationId (including the model id)
    * - The total number of Persons at the location
    * - A map containing the number of Persons per InfectionState at the location
    */
    using Type = std::vector<std::tuple<int, mio::abm::LocationType, mio::abm::LocationId, size_t,
                                        std::map<mio::abm::InfectionState, size_t>>>;
    static Type log(const mio::abm::Simulation<mio::GraphABModel>& sim)
    {
        Type location_information{};
        location_information.reserve(size_t(mio::abm::LocationType::Count));
        auto t = sim.get_time();
        for (auto&& loc : sim.get_model().get_locations()) {
            std::map<mio::abm::InfectionState, size_t> persons_per_infection_state;
            for (size_t i = 0; i < static_cast<size_t>(mio::abm::InfectionState::Count); ++i) {
                auto inf_state = mio::abm::InfectionState(i);
                persons_per_infection_state.insert(
                    {inf_state, sim.get_model().get_subpopulation(loc.get_id(), t, inf_state)});
            }
            location_information.push_back(std::make_tuple(loc.get_model_id(), loc.get_type(), loc.get_id(),
                                                           sim.get_model().get_number_persons(loc.get_id()),
                                                           persons_per_infection_state));
        }
        return location_information;
    }
};

int main()
{
    // This is an example with three age groups representing children, adults and seniors.
    size_t num_age_groups         = 3;
    const auto age_group_children = mio::AgeGroup(0);
    const auto age_group_adults   = mio::AgeGroup(1);
    const auto age_group_seniors  = mio::AgeGroup(2);

    auto model1 = mio::GraphABModel(num_age_groups, 0);

    //Set infection parameters
    model1.parameters.get<mio::abm::IncubationPeriod>()              = 4.;
    model1.parameters.get<mio::abm::InfectedNoSymptomsToSymptoms>()  = 2.;
    model1.parameters.get<mio::abm::InfectedNoSymptomsToRecovered>() = 4.;
    model1.parameters.get<mio::abm::InfectedSymptomsToRecovered>()   = 5.;
    model1.parameters.get<mio::abm::InfectedSymptomsToSevere>()      = 6.;
    model1.parameters.get<mio::abm::SevereToRecovered>()             = 8.;
    model1.parameters.get<mio::abm::SevereToCritical>()              = 7.;
    model1.parameters.get<mio::abm::CriticalToRecovered>()           = 10.;
    model1.parameters.get<mio::abm::CriticalToDead>()                = 11.;

    //Age group 0 goes to school and age group 1 goes to work
    model1.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_children] = true;
    model1.parameters.get<mio::abm::AgeGroupGotoWork>()[age_group_adults]     = true;

    //Household members can be child, parent or senior
    auto child = mio::abm::HouseholdMember(num_age_groups);
    child.set_age_weight(age_group_children, 1);
    auto parent = mio::abm::HouseholdMember(num_age_groups);
    parent.set_age_weight(age_group_adults, 1);
    auto adult = mio::abm::HouseholdMember(num_age_groups);
    adult.set_age_weight(age_group_adults, 1);
    adult.set_age_weight(age_group_seniors, 1);

    //Single-Person households
    auto single_hh = mio::abm::Household();
    single_hh.add_members(adult, 1);

    //Two-Adult household
    auto two_adult_hh = mio::abm::Household();
    two_adult_hh.add_members(adult, 2);

    //Single-Parent household
    auto single_parent_hh = mio::abm::Household();
    single_parent_hh.add_members(child, 1);
    single_parent_hh.add_members(parent, 1);

    //Family household
    auto family_hh = mio::abm::Household();
    family_hh.add_members(child, 1);
    family_hh.add_members(parent, 2);

    //Household groups for model 1
    auto single_hh_group_m1 = mio::abm::HouseholdGroup();
    single_hh_group_m1.add_households(single_hh, 5);
    auto two_adult_hh_group_m1 = mio::abm::HouseholdGroup();
    two_adult_hh_group_m1.add_households(two_adult_hh, 3);
    auto single_parent_hh_group_m1 = mio::abm::HouseholdGroup();
    single_parent_hh_group_m1.add_households(single_parent_hh, 5);
    auto family_hh_group_m1 = mio::abm::HouseholdGroup();
    family_hh_group_m1.add_households(family_hh, 10);
    add_household_group_to_model(model1, single_hh_group_m1);
    add_household_group_to_model(model1, two_adult_hh_group_m1);
    add_household_group_to_model(model1, single_hh_group_m1);
    add_household_group_to_model(model1, family_hh_group_m1);

    auto model2 = mio::GraphABModel(num_age_groups, 1);

    //Set infection parameters
    model2.parameters.get<mio::abm::IncubationPeriod>()              = 4.;
    model2.parameters.get<mio::abm::InfectedNoSymptomsToSymptoms>()  = 2.;
    model2.parameters.get<mio::abm::InfectedNoSymptomsToRecovered>() = 4.;
    model2.parameters.get<mio::abm::InfectedSymptomsToRecovered>()   = 5.;
    model2.parameters.get<mio::abm::InfectedSymptomsToSevere>()      = 6.;
    model2.parameters.get<mio::abm::SevereToRecovered>()             = 8.;
    model2.parameters.get<mio::abm::SevereToCritical>()              = 7.;
    model2.parameters.get<mio::abm::CriticalToRecovered>()           = 10.;
    model2.parameters.get<mio::abm::CriticalToDead>()                = 11.;

    //Age group 0 goes to school and age group 1 goes to work
    model2.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_children] = true;
    model2.parameters.get<mio::abm::AgeGroupGotoWork>()[age_group_adults]     = true;

    //Household groups for model 2
    auto single_hh_group_m2 = mio::abm::HouseholdGroup();
    single_hh_group_m2.add_households(single_hh, 6);
    auto two_adult_hh_group_m2 = mio::abm::HouseholdGroup();
    two_adult_hh_group_m2.add_households(two_adult_hh, 2);
    auto single_parent_hh_group_m2 = mio::abm::HouseholdGroup();
    single_parent_hh_group_m2.add_households(single_parent_hh, 10);
    auto family_hh_group_m2 = mio::abm::HouseholdGroup();
    family_hh_group_m2.add_households(family_hh, 11);
    add_household_group_to_model(model2, single_hh_group_m2);
    add_household_group_to_model(model2, two_adult_hh_group_m2);
    add_household_group_to_model(model2, single_hh_group_m2);
    add_household_group_to_model(model2, family_hh_group_m2);

    //Create locations for both models
    //model 1
    auto event_m1 = model1.add_location(mio::abm::LocationType::SocialEvent);
    model1.get_location(event_m1).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    auto hospital_m1 = model1.add_location(mio::abm::LocationType::Hospital);
    model1.get_location(hospital_m1).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    auto icu_m1 = model1.add_location(mio::abm::LocationType::ICU);
    model1.get_location(icu_m1).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto shop_m1 = model1.add_location(mio::abm::LocationType::BasicsShop);
    model1.get_location(shop_m1).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto school_m1 = model1.add_location(mio::abm::LocationType::School);
    model1.get_location(school_m1).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto work_m1 = model1.add_location(mio::abm::LocationType::Work);
    model1.get_location(work_m1).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    //model 2
    auto event_m2 = model2.add_location(mio::abm::LocationType::SocialEvent);
    model2.get_location(event_m2).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    auto hospital_m2 = model2.add_location(mio::abm::LocationType::Hospital);
    model2.get_location(hospital_m2).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    auto icu_m2 = model2.add_location(mio::abm::LocationType::ICU);
    model2.get_location(icu_m2).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto shop_m2 = model2.add_location(mio::abm::LocationType::BasicsShop);
    model2.get_location(shop_m2).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto school_m2 = model2.add_location(mio::abm::LocationType::School);
    model2.get_location(school_m2).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto work_m2 = model2.add_location(mio::abm::LocationType::Work);
    model2.get_location(work_m2).get_infection_parameters().set<mio::abm::MaximumContacts>(10);

    auto start_date = mio::abm::TimePoint(0);
    auto end_date   = mio::abm::TimePoint(0) + mio::abm::days(30);

    //Assign infection states and locations to persons from model 1
    std::vector<double> infection_distribution_m1{0.5, 0.3, 0.05, 0.05, 0.05, 0.05, 0.0, 0.0};
    for (auto& person : model1.get_persons()) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(model1.get_rng(), infection_distribution_m1));
        auto rng = mio::abm::PersonalRandomNumberGenerator(model1.get_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model1.parameters, start_date, infection_state));
        }
        person.set_assigned_location(mio::abm::LocationType::SocialEvent, event_m1, model1.get_id());
        person.set_assigned_location(mio::abm::LocationType::BasicsShop, shop_m1, model1.get_id());
        person.set_assigned_location(mio::abm::LocationType::Hospital, hospital_m1, model1.get_id());
        person.set_assigned_location(mio::abm::LocationType::ICU, icu_m1, model1.get_id());
        if (person.get_age() == age_group_children) {
            person.set_assigned_location(mio::abm::LocationType::School, school_m1, model1.get_id());
        }
        if (person.get_age() == age_group_adults) {
            //10% of adults in model 1 work in model 2
            size_t work_model = mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(),
                                                                                  std::vector<double>{0.9, 0.1});
            if (work_model == 1) { //person works in other model
                person.set_assigned_location(mio::abm::LocationType::Work, work_m2, model2.get_id());
            }
            else { //person works in same model
                person.set_assigned_location(mio::abm::LocationType::Work, work_m1, model1.get_id());
            }
        }
    }

    //Assign infection states and locations to persons from model 2
    std::vector<double> infection_distribution_m2{0.7, 0.1, 0.05, 0.05, 0.1, 0.0, 0.0, 0.0};
    for (auto& person : model2.get_persons()) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(model2.get_rng(), infection_distribution_m2));
        auto rng = mio::abm::PersonalRandomNumberGenerator(model2.get_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model2.parameters, start_date, infection_state));
        }
        person.set_assigned_location(mio::abm::LocationType::SocialEvent, event_m2, model2.get_id());
        person.set_assigned_location(mio::abm::LocationType::BasicsShop, shop_m2, model2.get_id());
        person.set_assigned_location(mio::abm::LocationType::Hospital, hospital_m2, model2.get_id());
        person.set_assigned_location(mio::abm::LocationType::ICU, icu_m2, model2.get_id());
        if (person.get_age() == age_group_children) {
            person.set_assigned_location(mio::abm::LocationType::School, school_m2, model2.get_id());
        }
        if (person.get_age() == age_group_adults) {
            //20% of adults in model 2 work in model 1
            size_t work_model = mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(),
                                                                                  std::vector<double>{0.2, 0.8});
            if (work_model == 1) { //person works in same model
                person.set_assigned_location(mio::abm::LocationType::Work, work_m2, model2.get_id());
            }
            else { //person works in other model
                person.set_assigned_location(mio::abm::LocationType::Work, work_m1, model1.get_id());
            }
        }
    }

    using HistoryType = mio::History<mio::DataWriterToMemory, Logger>;
    mio::Graph<mio::ABMSimulationNode<HistoryType>, mio::ABMMobilityEdge<HistoryType>> graph;
    graph.add_node(model1.get_id(), HistoryType{}, start_date, std::move(model1));
    graph.add_node(model2.get_id(), HistoryType{}, start_date, std::move(model2));
    graph.add_edge(model1.get_id(), model2.get_id());
    graph.add_edge(model2.get_id(), model1.get_id());

    auto exchange_time_span = mio::abm::hours(12);
    auto sim                = mio::make_abm_graph_sim<HistoryType>(start_date, exchange_time_span, std::move(graph));
    sim.advance(end_date);

    return 0;
}
