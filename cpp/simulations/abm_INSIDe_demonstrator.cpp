/*
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Khoa Nguyen
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
#include "abm/state.h"
#include "abm/location.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/io/io.h"
#include "memilio/utils/random_number_generator.h"
#include "boost/filesystem.hpp"
#include <fstream>
#include <sstream>

namespace fs = boost::filesystem;

struct LocationMapping {
    std::string inputId;
    std::vector<std::string> modelId{};
};

void read_txt(std::vector<std::pair<std::string, std::string>>& areas, std::vector<int>& inhabitants,
              const fs::path& input_dir)
{
    std::ifstream file;
    file.open(input_dir.string());
    std::string line = "";
    while (std::getline(file, line)) {
        int inhabitant;
        std::string Id;
        std::string area;
        std::string tmpString;
        std::stringstream inputString(line);
        std::getline(inputString, Id, ',');
        std::getline(inputString, tmpString, ',');
        inhabitant = atoi(tmpString.c_str());
        inhabitants.push_back(inhabitant);
        std::getline(inputString, area, ',');
        areas.push_back(std::make_pair(Id, area));
        line = "";
    }
}

mio::abm::LocationId add_location(mio::abm::World& world, mio::abm::LocationType type, double max_contacts, int persons,
                                  int volume)
{
    //add a social event to the world
    auto location = world.add_location(type);
    //set maximum number of contacts and capacity //TODO
    world.get_individualized_location(location).get_infection_parameters().set<mio::abm::MaximumContacts>(max_contacts);
    world.get_individualized_location(location).set_capacity(persons, volume);

    return location;
}

void insert_locations_to_map(std::vector<LocationMapping>& locationIds, std::string& inputId,
                             std::vector<mio::abm::LocationId>& locations)
{
    LocationMapping map;
    map.inputId = inputId;
    std::string locationType;
    std::string locationIndex;
    for (auto& location : locations) {
        locationType  = std::to_string(int(location.type));
        locationIndex = std::to_string(location.index);
        if (int(location.type) < 10) {
            locationType = "0" + locationType;
        }
        if (location.index < 10) {
            locationIndex = "0" + locationIndex;
        }
        map.modelId.push_back((locationType + locationIndex));
    }

    locationIds.push_back(map);
}

mio::abm::HouseholdGroup make_one_person_households(const mio::abm::HouseholdMember& member, int household_size,
                                                    int number_of_households)
{
    auto household_group = mio::abm::HouseholdGroup();
    for (int hh = 0; hh < number_of_households; ++hh) {
        auto household = mio::abm::Household();
        household.add_members(member, household_size);
        household_group.add_households(household, 1);
    }
    return household_group;
}

mio::abm::HouseholdGroup make_multiple_person_households(const mio::abm::HouseholdMember& child,
                                                         const mio::abm::HouseholdMember& parent,
                                                         const mio::abm::HouseholdMember& other, int household_size,
                                                         int number_of_two_parent_households,
                                                         int number_of_one_parent_households,
                                                         int number_of_other_households)
{
    auto household_group = mio::abm::HouseholdGroup();

    //Add two parent households
    auto hh_two_parents = mio::abm::Household();
    hh_two_parents.add_members(child, household_size - 2);
    hh_two_parents.add_members(parent, 2);
    household_group.add_households(hh_two_parents, number_of_two_parent_households);

    //Add one parent households
    auto hh_one_parent = mio::abm::Household();
    hh_one_parent.add_members(child, household_size - 1);
    hh_one_parent.add_members(parent, 1);
    household_group.add_households(hh_one_parent, number_of_one_parent_households);

    //add other households
    auto hh_other = mio::abm::Household();
    hh_other.add_members(other, household_size);
    household_group.add_households(hh_other, number_of_other_households);

    return household_group;
}

std::vector<mio::abm::LocationId> add_households(mio::abm::World& world, std::vector<double>& distribution,
                                                 int num_inhabitants)
{
    std::vector<int> households(distribution.size());
    size_t household_size;
    std::vector<mio::abm::LocationId> locations;
    int new_index = (int)world.get_locations()[(uint32_t)mio::abm::LocationType::Home].size();
    while (num_inhabitants > 0) {
        household_size = mio::DiscreteDistribution<size_t>::get_instance()(distribution);
        households[household_size] += 1;
        num_inhabitants -= (int)(household_size + 1);
    }

    //One-Person Households
    auto one_person_household_member = mio::abm::HouseholdMember();
    one_person_household_member.set_age_weight(mio::abm::AgeGroup::Age15to34, 5);
    one_person_household_member.set_age_weight(mio::abm::AgeGroup::Age35to59, 6);
    one_person_household_member.set_age_weight(mio::abm::AgeGroup::Age60to79, 4);
    one_person_household_member.set_age_weight(mio::abm::AgeGroup::Age80plus, 2);
    int number_of_one_person_households = households[0];

    auto one_person_household_group =
        make_one_person_households(one_person_household_member, 1, number_of_one_person_households);
    add_household_group_to_world(world, one_person_household_group);

    //Members for multiple person households
    auto child = mio::abm::HouseholdMember();
    child.set_age_weight(mio::abm::AgeGroup::Age0to4, 1);
    child.set_age_weight(mio::abm::AgeGroup::Age5to14, 1);

    auto parent = mio::abm::HouseholdMember();
    parent.set_age_weight(mio::abm::AgeGroup::Age15to34, 2);
    parent.set_age_weight(mio::abm::AgeGroup::Age35to59, 2);
    parent.set_age_weight(mio::abm::AgeGroup::Age60to79, 1);

    auto other = mio::abm::HouseholdMember();
    other.set_age_weight(mio::abm::AgeGroup::Age0to4, 1);
    other.set_age_weight(mio::abm::AgeGroup::Age5to14, 2);
    other.set_age_weight(mio::abm::AgeGroup::Age15to34, 3);
    other.set_age_weight(mio::abm::AgeGroup::Age35to59, 3);
    other.set_age_weight(mio::abm::AgeGroup::Age60to79, 2);
    other.set_age_weight(mio::abm::AgeGroup::Age80plus, 2);

    //Two-Person Households
    int two_person_two_parents = (int)(0.4 * households[1]);
    int two_person_one_parent  = (int)(0.4 * households[1]);
    int two_person_others      = households[1] - two_person_two_parents - two_person_one_parent;

    auto two_person_household_group = make_multiple_person_households(child, parent, other, 2, two_person_two_parents,
                                                                      two_person_one_parent, two_person_others);
    add_household_group_to_world(world, two_person_household_group);

    //Three-Person Households
    int three_person_two_parents = (int)(0.4 * households[2]);
    int three_person_one_parent  = (int)(0.4 * households[2]);
    int three_person_others      = households[2] - three_person_two_parents - three_person_one_parent;

    auto three_person_household_group = make_multiple_person_households(
        child, parent, other, 3, three_person_two_parents, three_person_one_parent, three_person_others);
    add_household_group_to_world(world, three_person_household_group);

    //Four-Person Households
    int four_person_two_parents = (int)(0.5 * households[3]);
    int four_person_one_parent  = (int)(0.2 * households[3]);
    int four_person_others      = households[3] - four_person_two_parents - four_person_one_parent;

    auto four_person_household_group = make_multiple_person_households(child, parent, other, 4, four_person_two_parents,
                                                                       four_person_one_parent, four_person_others);
    add_household_group_to_world(world, four_person_household_group);

    //Five-Person Households
    int five_person_two_parents = (int)(0.6 * households[4]);
    int five_person_one_parent  = (int)(0.1 * households[4]);
    int five_person_others      = households[4] - five_person_two_parents - five_person_one_parent;

    auto five_person_household_group = make_multiple_person_households(child, parent, other, 5, five_person_two_parents,
                                                                       five_person_one_parent, five_person_others);
    add_household_group_to_world(world, five_person_household_group);

    //fill location id vector with new locations for mapping
    mio::abm::LocationId id;
    id.type = mio::abm::LocationType::Home;
    //add LocationIds for Mapping
    for (int hh = 0; hh < std::accumulate(households.begin(), households.end(), 0); ++hh) {
        id.index = new_index;
        locations.emplace_back(id);
        ++new_index;
    }

    return locations;
}

void create_locations_from_input(std::vector<std::pair<std::string, std::string>>& areas, std::vector<int>& inhabitants,
                                 mio::abm::World& world, std::vector<LocationMapping>& locationIds,
                                 ScalarType one_person_hh, ScalarType two_person_hh, ScalarType three_person_hh,
                                 ScalarType four_person_hh, ScalarType five_person_hh)
{
    assert(areas.size() == inhabitants.size());
    std::string res                                = "residential";
    std::vector<ScalarType> household_distribution = {one_person_hh, two_person_hh, three_person_hh, four_person_hh,
                                                      five_person_hh};
    bool has_school                                = false;
    for (auto loc = 0; loc < areas.size(); ++loc) {
        std::vector<mio::abm::LocationId> locations;
        if (std::search((areas[loc].second).begin(), (areas[loc].second).end(), res.begin(), res.end()) !=
            (areas[loc].second).end()) {
            //Home
            locations = add_households(world, household_distribution, inhabitants[loc]);
        }
        else if (areas[loc].second == "mixed") {
            size_t location_type = mio::DiscreteDistribution<size_t>::get_instance()(std::vector<double>{1, 1});
            if (location_type) {
                locations.emplace_back(add_location(world, mio::abm::LocationType::Work, 40., 100, 2000));
            }
            else {
                //Home
                locations = add_households(world, household_distribution, inhabitants[loc]);
            }
        }
        else {
            if (areas[loc].second == "recreational") {
                //Social Event
                locations.emplace_back(add_location(world, mio::abm::LocationType::SocialEvent, 30., 30, 40));
            }
            else if (areas[loc].second == "shopping_business") {

                if (!has_school) {
                    locations.emplace_back(add_location(world, mio::abm::LocationType::School, 40., 500, 2000));
                    has_school = true;
                }
                else {
                    //uniformly draw whether shopping_business should be work or BasicsShop
                    size_t location_type = mio::DiscreteDistribution<size_t>::get_instance()(std::vector<double>{1, 1});
                    if (location_type) {
                        locations.emplace_back(add_location(world, mio::abm::LocationType::BasicsShop, 20., 100, 1000));
                    }
                    else {
                        locations.emplace_back(add_location(world, mio::abm::LocationType::Work, 40., 300, 2000));
                    }
                }
            }
            else if (areas[loc].second == "university") {
                locations.emplace_back(add_location(world, mio::abm::LocationType::Work, 50., 200, 4000));
            }
            else {
                mio::log_error("Area input type does not match to abm location type.");
            }
        }

        insert_locations_to_map(locationIds, areas[loc].first, locations);
    }
}

mio::abm::InfectionState get_infection_state(ScalarType exposed, ScalarType infected_no_symptoms,
                                             ScalarType infected_symptoms, ScalarType infected_severe,
                                             ScalarType infected_critical, ScalarType recovered_infected_no_symptoms,
                                             ScalarType recovered_infected)
{
    ScalarType susceptible = 1 - exposed - infected_no_symptoms - infected_symptoms - infected_severe -
                             infected_critical - recovered_infected_no_symptoms - recovered_infected;
    std::vector<ScalarType> weights = {
        susceptible,     exposed,           infected_no_symptoms,           infected_symptoms,
        infected_severe, infected_critical, recovered_infected_no_symptoms, recovered_infected};
    if (weights.size() != (size_t)mio::abm::InfectionState::Count - 1) {
        mio::log_error("Initialization in ABM wrong, please correct vector length.");
    }
    size_t state = mio::DiscreteDistribution<size_t>::get_instance()(weights);
    return (mio::abm::InfectionState)state;
}

void assign_infection_states(mio::abm::World& world, ScalarType exposed, ScalarType infected_no_symptoms,
                             ScalarType infected_symptoms, ScalarType infected_severe, ScalarType infected_critical,
                             ScalarType recovered_infected_no_symptoms, ScalarType recovered_infected)
{
    auto persons = world.get_persons();
    for (auto& person : persons) {
        world.set_infection_state(person, get_infection_state(exposed, infected_no_symptoms, infected_symptoms,
                                                              infected_severe, infected_critical,
                                                              recovered_infected_no_symptoms, recovered_infected));
    }
}

void assign_locations(mio::abm::World& world)
{
    //TODO: add hospital and ICU
    std::vector<mio::abm::Location> schools = world.get_locations()[(uint32_t)mio::abm::LocationType::School];
    std::vector<ScalarType> school_weights(schools.size(), 1);
    std::vector<mio::abm::Location> workplaces = world.get_locations()[(uint32_t)mio::abm::LocationType::Work];
    std::vector<ScalarType> workplaces_weights(workplaces.size(), 1);
    std::vector<mio::abm::Location> basic_shops = world.get_locations()[(uint32_t)mio::abm::LocationType::BasicsShop];
    std::vector<ScalarType> basic_shops_weights(basic_shops.size(), 1);
    std::vector<mio::abm::Location> social_events =
        world.get_locations()[(uint32_t)mio::abm::LocationType::SocialEvent];
    std::vector<ScalarType> social_event_weights(social_events.size(), 1);

    auto persons = world.get_persons();
    for (auto& person : persons) {
        //assign shop
        size_t shop = mio::DiscreteDistribution<size_t>::get_instance()(basic_shops_weights);
        person.set_assigned_location(basic_shops[shop]);
        //assign event
        size_t event = mio::DiscreteDistribution<size_t>::get_instance()(social_event_weights);
        person.set_assigned_location(social_events[event]);
        //assign work and school
        if (person.get_age() == mio::abm::AgeGroup::Age5to14) {
            size_t school = mio::DiscreteDistribution<size_t>::get_instance()(school_weights);
            person.set_assigned_location(schools[school]);
        }
        if (person.get_age() == mio::abm::AgeGroup::Age15to34 || person.get_age() == mio::abm::AgeGroup::Age35to59) {
            size_t work = mio::DiscreteDistribution<size_t>::get_instance()(workplaces_weights);
            person.set_assigned_location(workplaces[work]);
        }
    }
}

mio::abm::Simulation create_sampled_simulation(const mio::abm::TimePoint& t0, const fs::path& input_dir,
                                               std::vector<LocationMapping>& LocationIds)
{
    std::vector<std::pair<std::string, std::string>> areas;
    std::vector<int> inhabitants;
    // Read input file containing the locations and the number of inhabitants per location
    read_txt(areas, inhabitants, input_dir);

    //Assumed percentage for 1-Person, 2-Person, 3-Person, 4-Person and 5-Persons households
    ScalarType one_person_hh_pct = 0.37, two_person_hh_pct = 0.33, three_person_hh_pct = 0.15, four_person_hh_pct = 0.1,
               five_person_hh_pct = 0.05;
    // Assumed percentage of infection state at the beginning of the simulation.
    ScalarType exposed_pct = 0.005, infected_no_symptoms_pct = 0.001, infected_symptoms_pct = 0.001,
               infected_severe_pct = 0.0001, infected_critical_pct = 0.0, recovered_infected_no_symptoms_pct = 0.0,
               recovered_infected_symptoms_pct = 0.0;
    //Set global infection parameters
    mio::abm::GlobalInfectionParameters infection_params;

    //TODO: set infection parameters

    // Create the world with infection parameters.
    auto world = mio::abm::World(infection_params);

    //Transform the input location to correspondin abm location
    create_locations_from_input(areas, inhabitants, world, LocationIds, one_person_hh_pct, two_person_hh_pct,
                                three_person_hh_pct, four_person_hh_pct, five_person_hh_pct);

    //Assign an infection state to every person
    assign_infection_states(world, exposed_pct, infected_no_symptoms_pct, infected_symptoms_pct, infected_severe_pct,
                            infected_critical_pct, recovered_infected_no_symptoms_pct, recovered_infected_symptoms_pct);
    //Assign the locations to persons
    assign_locations(world);

    auto sim = mio::abm::Simulation(t0, std::move(world));
    return sim;
}

mio::IOResult<void> run(const fs::path& input_dir)
{
    auto t0      = mio::abm::TimePoint(0); // Start time per simulation
    auto tmax    = mio::abm::TimePoint(0) + mio::abm::days(14); // End time per simulation
    auto results = std::vector<mio::TimeSeries<ScalarType>>{}; // Vector of collected results

    std::vector<LocationMapping> LocationIds;

    auto sim = create_sampled_simulation(t0, input_dir, LocationIds);
    sim.advance(tmax);
    return mio::success();
}

template <class T>
void print(T& data)
{
    for (auto item : data) {
        std::cout << item << " ";
    }
    std::cout << std::endl;
}

int main()
{
    const fs::path input_dir =
        "C:/Users/bick_ju/Documents/INSIDe/Demonstrator/INSIDeDemonstrator/INSIDe_Demonstrator_AreaList.csv";
    auto result = run(input_dir);
    return 0;
}