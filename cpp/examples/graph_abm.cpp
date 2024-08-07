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
#include "graph_abm/model_wrapper.h"
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
    static Type log(const mio::abm::Simulation& sim)
    {
        Type location_information{};
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

    auto model1 = mio::abm::WrapperModel(num_age_groups, 0);

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

    // Vector holding all persons for the graph simulation. This vector is copied to all models at the end.
    std::vector<mio::abm::Person> persons;

    //Household groups for model 1
    auto single_hh_group_w1 = mio::abm::HouseholdGroup();
    single_hh_group_w1.add_households(single_hh, 5);
    auto two_adult_hh_group_w1 = mio::abm::HouseholdGroup();
    two_adult_hh_group_w1.add_households(two_adult_hh, 3);
    auto single_parent_hh_group_w1 = mio::abm::HouseholdGroup();
    single_parent_hh_group_w1.add_households(single_parent_hh, 5);
    auto family_hh_group_w1 = mio::abm::HouseholdGroup();
    family_hh_group_w1.add_households(family_hh, 10);
    add_household_group_to_model(model1, single_hh_group_w1);
    add_household_group_to_model(model1, two_adult_hh_group_w1);
    add_household_group_to_model(model1, single_hh_group_w1);
    add_household_group_to_model(model1, family_hh_group_w1);

    //add persons from model 0 to vector
    for (auto& person : model1.get_persons()) {
        mio::abm::PersonId new_id{static_cast<uint32_t>(persons.size())};
        persons.emplace_back(person, new_id);
    }

    auto model2 = mio::abm::Model(num_age_groups, 1);

    //Household groups for model 2
    auto single_hh_group_w2 = mio::abm::HouseholdGroup();
    single_hh_group_w2.add_households(single_hh, 6);
    auto two_adult_hh_group_w2 = mio::abm::HouseholdGroup();
    two_adult_hh_group_w2.add_households(two_adult_hh, 2);
    auto single_parent_hh_group_w2 = mio::abm::HouseholdGroup();
    single_parent_hh_group_w2.add_households(single_parent_hh, 10);
    auto family_hh_group_w2 = mio::abm::HouseholdGroup();
    family_hh_group_w2.add_households(family_hh, 11);
    add_household_group_to_model(model2, single_hh_group_w2);
    add_household_group_to_model(model2, two_adult_hh_group_w2);
    add_household_group_to_model(model2, single_hh_group_w2);
    add_household_group_to_model(model2, family_hh_group_w2);

    //add persons from model 1 to vector
    for (auto& person : model2.get_persons()) {
        mio::abm::PersonId new_id{static_cast<uint32_t>(persons.size())};
        persons.emplace_back(person, new_id);
    }

    //Create locations for both models
    //model 0
    auto event_w1 = model1.add_location(mio::abm::LocationType::SocialEvent);
    model1.get_location(event_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    auto hospital_w1 = model1.add_location(mio::abm::LocationType::Hospital);
    model1.get_location(hospital_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    auto icu_w1 = model1.add_location(mio::abm::LocationType::ICU);
    model1.get_location(icu_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto shop_w1 = model1.add_location(mio::abm::LocationType::BasicsShop);
    model1.get_location(shop_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto school_w1 = model1.add_location(mio::abm::LocationType::School);
    model1.get_location(school_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto work_w1 = model1.add_location(mio::abm::LocationType::Work);
    model1.get_location(work_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    //model 1
    auto event_w2 = model2.add_location(mio::abm::LocationType::SocialEvent);
    model2.get_location(event_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    auto hospital_w2 = model2.add_location(mio::abm::LocationType::Hospital);
    model2.get_location(hospital_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    auto icu_w2 = model2.add_location(mio::abm::LocationType::ICU);
    model2.get_location(icu_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto shop_w2 = model2.add_location(mio::abm::LocationType::BasicsShop);
    model2.get_location(shop_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto school_w2 = model2.add_location(mio::abm::LocationType::School);
    model2.get_location(school_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto work_w2 = model2.add_location(mio::abm::LocationType::Work);
    model2.get_location(work_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(10);

    auto start_date = mio::abm::TimePoint(0);
    auto end_date   = mio::abm::TimePoint(0) + mio::abm::days(30);
    std::vector<uint32_t> params_e1;
    std::vector<uint32_t> params_e2;

    //Assign infection states and locations
    std::vector<double> infection_distribution{0.5, 0.3, 0.05, 0.05, 0.05, 0.05, 0.0, 0.0};
    for (auto& person : persons) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), infection_distribution));
        auto rng = mio::abm::PersonalRandomNumberGenerator(mio::thread_local_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model1.parameters, start_date, infection_state));
        }
        // Assign locations to persons from model 1
        if (person.get_assigned_location_model_id(mio::abm::LocationType::Home) == model1.get_id()) {
            person.set_assigned_location(mio::abm::LocationType::SocialEvent, event_w1, model1.get_id());
            person.set_assigned_location(mio::abm::LocationType::BasicsShop, shop_w1, model1.get_id());
            person.set_assigned_location(mio::abm::LocationType::Hospital, hospital_w1, model1.get_id());
            person.set_assigned_location(mio::abm::LocationType::ICU, icu_w1, model1.get_id());
            if (person.get_age() == age_group_children) {
                person.set_assigned_location(mio::abm::LocationType::School, school_w1, model1.get_id());
            }
            if (person.get_age() == age_group_adults) {
                //10% of adults in model 0 work in model 1
                size_t work_model = mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(),
                                                                                      std::vector<double>{0.9, 0.1});
                if (work_model == 1) { //person works in other model
                    person.set_assigned_location(mio::abm::LocationType::Work, work_w2, model2.get_id());
                    //add person to edge parameters
                    params_e1.push_back(person.get_id().get());
                }
                else { //person works in same model
                    person.set_assigned_location(mio::abm::LocationType::Work, work_w1, model1.get_id());
                }
            }
        }
        else {
            person.set_assigned_location(mio::abm::LocationType::SocialEvent, event_w2, model2.get_id());
            person.set_assigned_location(mio::abm::LocationType::BasicsShop, shop_w2, model2.get_id());
            person.set_assigned_location(mio::abm::LocationType::Hospital, hospital_w2, model2.get_id());
            person.set_assigned_location(mio::abm::LocationType::ICU, icu_w2, model2.get_id());
            if (person.get_age() == age_group_children) {
                person.set_assigned_location(mio::abm::LocationType::School, school_w2, model2.get_id());
            }
            if (person.get_age() == age_group_adults) {
                //20% of adults in model 1 work in model 0
                size_t work_model = mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(),
                                                                                      std::vector<double>{0.2, 0.8});
                if (work_model == 0) { //person works in other model
                    person.set_assigned_location(mio::abm::LocationType::Work, work_w1, model1.get_id());
                    //add person to edge parameters
                    params_e2.push_back(person.get_id().get());
                }
                else { //person works in same model
                    person.set_assigned_location(mio::abm::LocationType::Work, work_w2, model2.get_id());
                }
            }
        }
    }

    //copy persons to both models
    model1.set_persons(persons);
    model2.set_persons(persons);

    using HistoryType = mio::History<mio::DataWriterToMemory, Logger>;
    mio::Graph<mio::ABMSimulationNode<HistoryType>, mio::ABMMobilityEdge<HistoryType>> graph;
    graph.add_node(model1.get_id(), HistoryType{}, start_date, std::move(model1));
    graph.add_node(model2.get_id(), HistoryType{}, start_date, std::move(model2));
    graph.add_edge(0, 1, params_e1,
                   std::vector<mio::ABMMobilityEdge<HistoryType>::MobilityRuleType>{&mio::apply_commuting});
    graph.add_edge(1, 0, params_e2,
                   std::vector<mio::ABMMobilityEdge<HistoryType>::MobilityRuleType>{&mio::apply_commuting});

    auto sim = mio::make_abm_graph_sim<HistoryType>(start_date, mio::abm::hours(12), std::move(graph));
    sim.advance(end_date);

    auto& log_n1 = std::get<0>(sim.get_graph().nodes()[0].property.get_history()).get_log();
    auto& log_n2 = std::get<0>(sim.get_graph().nodes()[1].property.get_history()).get_log();

    return 0;
}
