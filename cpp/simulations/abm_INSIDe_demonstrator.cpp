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
#include "abm/infection_state.h"
#include "abm/location.h"
#include "abm/simulation.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/io/io.h"
#include "memilio/io/history.h"
#include "memilio/utils/random_number_generator.h"
#include "boost/filesystem.hpp"
#include <cstddef>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = boost::filesystem;

// Assign the name to general age group.
size_t num_age_groups         = 6;
const auto age_group_0_to_4   = mio::AgeGroup(0);
const auto age_group_5_to_14  = mio::AgeGroup(1);
const auto age_group_15_to_34 = mio::AgeGroup(2);
const auto age_group_35_to_59 = mio::AgeGroup(3);
const auto age_group_60_to_79 = mio::AgeGroup(4);
const auto age_group_80_plus  = mio::AgeGroup(5);

struct LocationMapping {
    std::string inputId;
    std::vector<std::string> modelId{};
};

std::string convert_loc_id_to_string(std::tuple<mio::abm::LocationType, uint32_t> tuple_id)
{
    std::string locationType  = std::to_string(static_cast<std::uint32_t>(std::get<0>(tuple_id)));
    std::string locationIndex = std::to_string(std::get<1>(tuple_id));
    if (static_cast<std::uint32_t>(std::get<0>(tuple_id)) < 10) {
        locationType = "0" + locationType;
    }
    if (std::get<1>(tuple_id) < 10) {
        locationIndex = "0" + locationIndex;
    }

    return "I" + locationType + locationIndex;
}

std::vector<std::tuple<uint32_t, mio::abm::TimeSpan>> get_agents_per_location(
    std::tuple<mio::abm::LocationType, uint32_t> loc_id,
    std::vector<std::tuple<mio::abm::LocationId, uint32_t, mio::abm::TimeSpan, mio::abm::InfectionState>>& log)
{
    std::vector<std::tuple<uint32_t, mio::abm::TimeSpan>> agents_per_location;
    for (auto& log_tuple : log) {
        if (std::get<0>(log_tuple).type == std::get<0>(loc_id) && std::get<0>(log_tuple).index == std::get<1>(loc_id)) {
            agents_per_location.push_back(std::make_tuple(std::get<1>(log_tuple), std::get<2>(log_tuple)));
        }
    }
    return agents_per_location;
}

/**
 * read input test file and save input areas and inhabitants per area
 * @param[in, out] areas vector that is filled with area types and area ids from input file
 * @param[in, out] inhabitants vector that is filled with the inhabitants per area given in the input file
 * @param[in] input_dir path to input txt file
*/
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

/**
 * Add location to world
 * @param[in] world world object the location is added to
 * @param[in] type location type
 * @param[in] max_contacts maximum number of contacts at the location
 * @param[in] persons maximum number of persons for location capacity
 * @param[in] volume maximum volume for location capacity
 * @return added location
*/
mio::abm::LocationId add_location(mio::abm::World& world, mio::abm::LocationType type, double max_contacts, int persons,
                                  int volume)
{
    //add a location of given type to the world
    auto location = world.add_location(type);
    //set maximum number of contacts and capacity //TODO
    world.get_individualized_location(location).get_infection_parameters().set<mio::abm::MaximumContacts>(max_contacts);
    world.get_individualized_location(location).set_capacity(persons, volume);

    return location;
}

/**
 * Insert the abm location Ids to mapping for the corresponding input area id.
 * @param[in, out] locationIds vector that maps every input area id to the corresponding abm location ids
 * @param[in] inputId input area id that is added to mapping vector
 * @param[in] locations abm locations that correspond to inputId
 * @return mapping vector
*/
void insert_locations_to_map(std::vector<LocationMapping>& locationIds, std::string& inputId,
                             std::vector<mio::abm::LocationId>& locations)
{
    LocationMapping map;
    map.inputId = inputId;
    std::string locationType;
    std::string locationIndex;
    //An abm locationId consists of a type and an index.
    //For the mapping the locationId is stored as a string of the form xxyy where xx specifies
    //the location type and yy the location index
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

/**
 * Make a one-person household group
 * @param[in] member household member
 * @param[in] number_of_households number of one-person households in the household group
 * @return household_group The one-person household group
 *
*/
mio::abm::HouseholdGroup make_one_person_households(const mio::abm::HouseholdMember& member, int number_of_households)
{
    auto household_group = mio::abm::HouseholdGroup();
    for (int hh = 0; hh < number_of_households; ++hh) {
        auto household = mio::abm::Household();
        household.add_members(member, 1);
        //add one-person household to household group
        household_group.add_households(household, 1);
    }
    return household_group;
}

/**
 * Make a multiple-person household group
 * @param[in] child household member representing a child
 * @param[in] parent household member representing a parent
 * @param[in] other random household member
 * @param[in] household_size household size e.g 2-person household, 3-person household etc.
 * @param[in] number_of_two_parent_households number of households with two parents and (household_size - 2) children
 * @param[in] number_of_one_parent_households number of households with one parent and (household_size - 1) children
 * @param[in] number_of_other_households number of households with random members
 * @return household_group multiple-person household group
*/
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

/**
 * Add households to the world for a given number of inhabitants.
 * @param[in, out] world
 * @param[in] distribution vector containing the percentages of 1-person, 2-person, ... households
 * @param[in] num_inhabitants number of inhabitants that should be distributed to households
 * @return locations vector with location ids of the added households
*/
std::vector<mio::abm::LocationId> add_households(mio::abm::World& world, std::vector<double>& distribution,
                                                 int num_inhabitants)
{
    //vector that saves the number of households for every household size
    std::vector<int> households(distribution.size());
    size_t household_size;
    std::vector<mio::abm::LocationId> locations;
    //index of the first new household
    size_t new_index = world.get_locations().size();
    while (num_inhabitants > 0) {
        //draw household size from the given distribution
        household_size = mio::DiscreteDistribution<size_t>::get_instance()(world.get_rng(), distribution);
        //increase the number of households of the drawn household size
        households[household_size] += 1;
        num_inhabitants -= (int)(household_size + 1);
    }

    //One-Person Households
    auto one_person_household_member = mio::abm::HouseholdMember(num_age_groups);
    one_person_household_member.set_age_weight(age_group_15_to_34, 5);
    one_person_household_member.set_age_weight(age_group_35_to_59, 6);
    one_person_household_member.set_age_weight(age_group_60_to_79, 4);
    one_person_household_member.set_age_weight(age_group_80_plus, 2);
    int number_of_one_person_households = households[0];

    auto one_person_household_group =
        make_one_person_households(one_person_household_member, number_of_one_person_households);
    add_household_group_to_world(world, one_person_household_group);

    //Members for multiple person households
    auto child = mio::abm::HouseholdMember(num_age_groups);
    child.set_age_weight(age_group_0_to_4, 1);
    child.set_age_weight(age_group_5_to_14, 1);

    auto parent = mio::abm::HouseholdMember(num_age_groups);
    parent.set_age_weight(age_group_15_to_34, 2);
    parent.set_age_weight(age_group_35_to_59, 2);
    parent.set_age_weight(age_group_60_to_79, 1);

    auto other = mio::abm::HouseholdMember(num_age_groups);
    other.set_age_weight(age_group_0_to_4, 1);
    other.set_age_weight(age_group_5_to_14, 2);
    other.set_age_weight(age_group_15_to_34, 3);
    other.set_age_weight(age_group_35_to_59, 3);
    other.set_age_weight(age_group_60_to_79, 2);
    other.set_age_weight(age_group_80_plus, 2);

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
        id.index = (int)new_index;
        locations.emplace_back(id);
        ++new_index;
    }

    return locations;
}

/**
 * Creates abm locations from input areas
 * @param[in] areas input areas consisting of a type and an id
 * @param[in] inhabitants number of inhabitants per input area
 * @param[in, out] world
 * @param[in, out] locationIds mapping of abm location ids to corresponding input area ids
 * @param[in] one_person_hh percentage of one-person households
 * @param[in] two_person_hh percentage of two-person households
 * @param[in] three_person_hh percentage of three-person households
 * @param[in] four_person_hh percentage of four-person households
 * @param[in] five_person_hh percentage of five-person households
*/
void create_locations_from_input(std::vector<std::pair<std::string, std::string>>& areas, std::vector<int>& inhabitants,
                                 mio::abm::World& world, std::vector<LocationMapping>& locationIds,
                                 ScalarType one_person_hh, ScalarType two_person_hh, ScalarType three_person_hh,
                                 ScalarType four_person_hh, ScalarType five_person_hh)
{
    assert(areas.size() == inhabitants.size());
    std::string residential                        = "residential";
    std::vector<ScalarType> household_distribution = {one_person_hh, two_person_hh, three_person_hh, four_person_hh,
                                                      five_person_hh};
    //school and hospital is needed for migration rules
    bool has_school   = false;
    bool has_hospital = false;
    for (size_t loc = 0; loc < areas.size(); ++loc) {
        std::vector<mio::abm::LocationId> locations;
        if (std::search((areas[loc].second).begin(), (areas[loc].second).end(), residential.begin(),
                        residential.end()) != (areas[loc].second).end()) {
            //Home
            locations = add_households(world, household_distribution, inhabitants[loc]);
        }
        else if (areas[loc].second == "mixed" || areas[loc].second == "mixed\r") {
            //areas of type "mixed" are equally distributed to  of location of type Home amd type Work
            size_t location_type =
                mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), std::vector<double>{1., 1.});
            if (location_type) {
                locations.emplace_back(add_location(world, mio::abm::LocationType::Work, 40., 100, 2000));
            }
            else {
                //Home
                locations = add_households(world, household_distribution, inhabitants[loc]);
            }
        }
        else {
            if (areas[loc].second == "recreational" || areas[loc].second == "recreational\r") {
                //Social Event
                locations.emplace_back(add_location(world, mio::abm::LocationType::SocialEvent, 30., 30, 40));
            }
            else if (areas[loc].second == "shopping_business" || areas[loc].second == "shopping_business\r") {
                //school, hospital and icu are added first
                if (!has_school) {
                    locations.emplace_back(add_location(world, mio::abm::LocationType::School, 40., 500, 2000));
                    has_school = true;
                }
                else if (!has_hospital) {
                    locations.emplace_back(add_location(world, mio::abm::LocationType::Hospital, 5., 300, 10000));
                    locations.emplace_back(add_location(world, mio::abm::LocationType::ICU, 5., 30, 1000));
                    has_hospital = true;
                }
                else {
                    //the other areas of type 'shopping_business' are equally distributed to locations of type BasicsShop and type Work
                    size_t location_type = mio::DiscreteDistribution<size_t>::get_instance()(
                        mio::thread_local_rng(), std::vector<double>{1., 1.});
                    if (location_type) {
                        locations.emplace_back(add_location(world, mio::abm::LocationType::BasicsShop, 20., 100, 1000));
                    }
                    else {
                        locations.emplace_back(add_location(world, mio::abm::LocationType::Work, 40., 300, 2000));
                    }
                }
            }
            else if (areas[loc].second == "university" || areas[loc].second == "university\r") {
                //area of type 'university' is converted to location of type Work
                locations.emplace_back(add_location(world, mio::abm::LocationType::Work, 50., 200, 4000));
            }
            else {
                mio::log_error("Area input type does not match to abm location type.");
            }
        }
        //insert locations to input area mapping
        insert_locations_to_map(locationIds, areas[loc].first, locations);
    }
}

/**
* Returns parameters for LogNormalDistributiom given desired expected value and standard deviation.
* @param[in] mean desired expected value
* @param[in] std desired standard deviation
* @return pair with parameters for LogNormalDistribtuion
*/
std::pair<double, double> get_my_and_sigma(double mean, double std)
{
    double my    = log(mean * mean / sqrt(mean * mean + std * std));
    double sigma = sqrt(log(1 + std * std / (mean * mean)));
    return {my, sigma};
}

/**
 * Set infection parameters
 * @param[in, out] infection_params infection parameters
*/
void set_infection_parameters(mio::abm::Parameters& infection_params)
{
    //set parameters for every agegroup
    auto incubation_period_params                      = get_my_and_sigma(3, 1.2);
    infection_params.get<mio::abm::IncubationPeriod>() = {incubation_period_params.first,
                                                          incubation_period_params.second};

    //0-4
    auto TimeInfectedNoSymptomsToSymptoms = get_my_and_sigma(2.2, 0.5);
    infection_params
        .get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = {
        TimeInfectedNoSymptomsToSymptoms.first, TimeInfectedNoSymptomsToSymptoms.second};

    auto TimeInfectedNoSymptomsToRecovered = get_my_and_sigma(9.2, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = {
        TimeInfectedNoSymptomsToRecovered.first, TimeInfectedNoSymptomsToRecovered.second};

    auto TimeInfectedSymptomsToSevere = get_my_and_sigma(10.5, 1.1);
    infection_params
        .get<mio::abm::TimeInfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = {
        TimeInfectedSymptomsToSevere.first, TimeInfectedSymptomsToSevere.second};

    auto TimeInfectedSymptomsToRecovered = get_my_and_sigma(7.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = {
        TimeInfectedSymptomsToRecovered.first, TimeInfectedSymptomsToRecovered.second};

    auto TimeInfectedSevereToRecovered = get_my_and_sigma(5.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = {
        TimeInfectedSevereToRecovered.first, TimeInfectedSevereToRecovered.second};

    auto TimeInfectedSevereToCritical = get_my_and_sigma(5.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = {
        TimeInfectedSevereToCritical.first, TimeInfectedSevereToCritical.second};

    auto TimeInfectedCriticalToRecovered = get_my_and_sigma(7.0, 3.0);
    infection_params
        .get<mio::abm::TimeInfectedCriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = {
        TimeInfectedCriticalToRecovered.first, TimeInfectedCriticalToRecovered.second};

    auto TimeInfectedCriticalToDead = get_my_and_sigma(6.0, 2.0);
    infection_params.get<mio::abm::TimeInfectedCriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        {TimeInfectedCriticalToDead.first, TimeInfectedCriticalToDead.second};

    //5-14
    TimeInfectedNoSymptomsToSymptoms = get_my_and_sigma(2.2, 0.5);
    infection_params
        .get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = {
        TimeInfectedNoSymptomsToSymptoms.first, TimeInfectedNoSymptomsToSymptoms.second};

    TimeInfectedNoSymptomsToRecovered = get_my_and_sigma(9.2, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = {
        TimeInfectedNoSymptomsToRecovered.first, TimeInfectedNoSymptomsToRecovered.second};

    TimeInfectedSymptomsToSevere = get_my_and_sigma(10.5, 1.1);
    infection_params
        .get<mio::abm::TimeInfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = {
        TimeInfectedSymptomsToSevere.first, TimeInfectedSymptomsToSevere.second};

    TimeInfectedSymptomsToRecovered = get_my_and_sigma(7.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = {
        TimeInfectedSymptomsToRecovered.first, TimeInfectedSymptomsToRecovered.second};

    TimeInfectedSevereToRecovered = get_my_and_sigma(5.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = {
        TimeInfectedSevereToRecovered.first, TimeInfectedSevereToRecovered.second};

    TimeInfectedSevereToCritical = get_my_and_sigma(5.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = {
        TimeInfectedSevereToCritical.first, TimeInfectedSevereToCritical.second};

    TimeInfectedCriticalToRecovered = get_my_and_sigma(7.0, 3.0);
    infection_params
        .get<mio::abm::TimeInfectedCriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = {
        TimeInfectedCriticalToRecovered.first, TimeInfectedCriticalToRecovered.second};

    TimeInfectedCriticalToDead = get_my_and_sigma(6.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedCriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = {
        TimeInfectedCriticalToDead.first, TimeInfectedCriticalToDead.second};

    //15-34
    TimeInfectedNoSymptomsToSymptoms = get_my_and_sigma(2.2, 0.5);
    infection_params
        .get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = {
        TimeInfectedNoSymptomsToSymptoms.first, TimeInfectedNoSymptomsToSymptoms.second};

    TimeInfectedNoSymptomsToRecovered = get_my_and_sigma(9.2, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = {
        TimeInfectedNoSymptomsToRecovered.first, TimeInfectedNoSymptomsToRecovered.second};

    TimeInfectedSymptomsToSevere = get_my_and_sigma(10.5, 1.1);
    infection_params
        .get<mio::abm::TimeInfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = {
        TimeInfectedSymptomsToSevere.first, TimeInfectedSymptomsToSevere.second};

    TimeInfectedSymptomsToRecovered = get_my_and_sigma(7.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = {
        TimeInfectedSymptomsToRecovered.first, TimeInfectedSymptomsToRecovered.second};

    TimeInfectedSevereToRecovered = get_my_and_sigma(6.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = {
        TimeInfectedSevereToRecovered.first, TimeInfectedSevereToRecovered.second};

    TimeInfectedSevereToCritical = get_my_and_sigma(5.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = {
        TimeInfectedSevereToCritical.first, TimeInfectedSevereToCritical.second};

    TimeInfectedCriticalToRecovered = get_my_and_sigma(7.0, 3.0);
    infection_params
        .get<mio::abm::TimeInfectedCriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = {
        TimeInfectedCriticalToRecovered.first, TimeInfectedCriticalToRecovered.second};

    TimeInfectedCriticalToDead = get_my_and_sigma(6.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedCriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = {
        TimeInfectedCriticalToDead.first, TimeInfectedCriticalToDead.second};

    //35-59
    TimeInfectedNoSymptomsToSymptoms = get_my_and_sigma(2.2, 0.5);
    infection_params
        .get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = {
        TimeInfectedNoSymptomsToSymptoms.first, TimeInfectedNoSymptomsToSymptoms.second};

    TimeInfectedNoSymptomsToRecovered = get_my_and_sigma(9.2, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = {
        TimeInfectedNoSymptomsToRecovered.first, TimeInfectedNoSymptomsToRecovered.second};

    TimeInfectedSymptomsToSevere = get_my_and_sigma(6.0, 1.1);
    infection_params
        .get<mio::abm::TimeInfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = {
        TimeInfectedSymptomsToSevere.first, TimeInfectedSymptomsToSevere.second};

    TimeInfectedSymptomsToRecovered = get_my_and_sigma(7.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = {
        TimeInfectedSymptomsToRecovered.first, TimeInfectedSymptomsToRecovered.second};

    TimeInfectedSevereToRecovered = get_my_and_sigma(8.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = {
        TimeInfectedSevereToRecovered.first, TimeInfectedSevereToRecovered.second};

    TimeInfectedSevereToCritical = get_my_and_sigma(5.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = {
        TimeInfectedSevereToCritical.first, TimeInfectedSevereToCritical.second};

    TimeInfectedCriticalToRecovered = get_my_and_sigma(17.5, 3.0);
    infection_params
        .get<mio::abm::TimeInfectedCriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = {
        TimeInfectedCriticalToRecovered.first, TimeInfectedCriticalToRecovered.second};

    TimeInfectedCriticalToDead = get_my_and_sigma(16.5, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedCriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = {
        TimeInfectedCriticalToDead.first, TimeInfectedCriticalToDead.second};

    //60-79
    TimeInfectedNoSymptomsToSymptoms = get_my_and_sigma(2.2, 0.5);
    infection_params
        .get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = {
        TimeInfectedNoSymptomsToSymptoms.first, TimeInfectedNoSymptomsToSymptoms.second};

    TimeInfectedNoSymptomsToRecovered = get_my_and_sigma(9.2, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = {
        TimeInfectedNoSymptomsToRecovered.first, TimeInfectedNoSymptomsToRecovered.second};

    TimeInfectedSymptomsToSevere = get_my_and_sigma(6.0, 1.1);
    infection_params
        .get<mio::abm::TimeInfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = {
        TimeInfectedSymptomsToSevere.first, TimeInfectedSymptomsToSevere.second};

    TimeInfectedSymptomsToRecovered = get_my_and_sigma(7.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = {
        TimeInfectedSymptomsToRecovered.first, TimeInfectedSymptomsToRecovered.second};

    TimeInfectedSevereToRecovered = get_my_and_sigma(10.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = {
        TimeInfectedSevereToRecovered.first, TimeInfectedSevereToRecovered.second};

    TimeInfectedSevereToCritical = get_my_and_sigma(5.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = {
        TimeInfectedSevereToCritical.first, TimeInfectedSevereToCritical.second};

    TimeInfectedCriticalToRecovered = get_my_and_sigma(17.5, 3.0);
    infection_params
        .get<mio::abm::TimeInfectedCriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = {
        TimeInfectedCriticalToRecovered.first, TimeInfectedCriticalToRecovered.second};

    TimeInfectedCriticalToDead = get_my_and_sigma(16.5, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedCriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = {
        TimeInfectedCriticalToDead.first, TimeInfectedCriticalToDead.second};

    //80+
    TimeInfectedNoSymptomsToSymptoms = get_my_and_sigma(2.2, 0.5);
    infection_params
        .get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] = {
        TimeInfectedNoSymptomsToSymptoms.first, TimeInfectedNoSymptomsToSymptoms.second};

    TimeInfectedNoSymptomsToRecovered = get_my_and_sigma(9.2, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] = {
        TimeInfectedNoSymptomsToRecovered.first, TimeInfectedNoSymptomsToRecovered.second};

    TimeInfectedSymptomsToSevere = get_my_and_sigma(6.0, 1.1);
    infection_params
        .get<mio::abm::TimeInfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] = {
        TimeInfectedSymptomsToSevere.first, TimeInfectedSymptomsToSevere.second};

    TimeInfectedSymptomsToRecovered = get_my_and_sigma(7.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] = {
        TimeInfectedSymptomsToRecovered.first, TimeInfectedSymptomsToRecovered.second};

    TimeInfectedSevereToRecovered = get_my_and_sigma(15.0, 3.0);
    infection_params
        .get<mio::abm::TimeInfectedSevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] = {
        TimeInfectedSevereToRecovered.first, TimeInfectedSevereToRecovered.second};

    TimeInfectedSevereToCritical = get_my_and_sigma(5.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedSevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] = {
        TimeInfectedSevereToCritical.first, TimeInfectedSevereToCritical.second};

    TimeInfectedCriticalToRecovered = get_my_and_sigma(12.5, 3.0);
    infection_params
        .get<mio::abm::TimeInfectedCriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] = {
        TimeInfectedCriticalToRecovered.first, TimeInfectedCriticalToRecovered.second};

    TimeInfectedCriticalToDead = get_my_and_sigma(11.0, 2.0);
    infection_params
        .get<mio::abm::TimeInfectedCriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] = {
        TimeInfectedCriticalToDead.first, TimeInfectedCriticalToDead.second};

    // Set percentage parameters
    infection_params.get<mio::abm::SymptomsPerInfectedNoSymptoms>() = 0.79;
    infection_params.get<mio::abm::SeverePerInfectedSymptoms>()     = 0.08;
    infection_params.get<mio::abm::CriticalPerInfectedSevere>()     = 0.18;
    infection_params.get<mio::abm::DeathsPerInfectedCritical>()     = 0.1;
}

/**
 * Get a person's infection state according to a given distribution
 * @param[in] exposed percentage of infection state 'exposed'
 * @param[in] infected_no_symptoms percentage of infection state 'infected no symptoms'
 * @param[in] infected_symptoms percentage of of infection state 'infected symptoms'
 * @param[in] infected_severe percentage of infection state 'infected severe'
 * @param[in] infected_critical percentage of infection state 'infected critical'
 * @param[in] recovered_infected_no_symptoms percentage of infection state 'recovered of state infected no symptoms'
 * @param[in] recovered_infected percentage of infection state 'recovered of other infection state'
 * @return state drawn infection state
*/
mio::abm::InfectionState get_infection_state(mio::abm::Person::RandomNumberGenerator& rng, ScalarType exposed,
                                             ScalarType infected_no_symptoms, ScalarType infected_symptoms,
                                             ScalarType infected_severe, ScalarType infected_critical,
                                             ScalarType recovered)
{
    ScalarType susceptible =
        1 - exposed - infected_no_symptoms - infected_symptoms - infected_severe - infected_critical - recovered;
    std::vector<ScalarType> weights = {
        susceptible, exposed, infected_no_symptoms, infected_symptoms, infected_severe, infected_critical, recovered};
    if (weights.size() != (size_t)mio::abm::InfectionState::Count - 1) {
        mio::log_error("Initialization in ABM wrong, please correct vector length.");
    }
    size_t state = mio::DiscreteDistribution<size_t>::get_instance()(rng, weights);
    return (mio::abm::InfectionState)state;
}

/**
 * Assign an infection state to each person.
 * @param[in, out] world
 * @param[in] t0 starting time point
 * @param[in] exposed percentage of infection state 'exposed'
 * @param[in] infected_no_symptoms percentage of infection state 'infected no symptoms'
 * @param[in] infected_symptoms percentage of of infection state 'infected symptoms'
 * @param[in] infected_severe percentage of infection state 'infected severe'
 * @param[in] infected_critical percentage of infection state 'infected critical'
 * @param[in] recovered percentage of infection state 'recovered'
*/
void assign_infection_states(mio::abm::World& world, mio::abm::TimePoint t0, ScalarType exposed,
                             ScalarType infected_no_symptoms, ScalarType infected_symptoms, ScalarType infected_severe,
                             ScalarType infected_critical, ScalarType recovered)
{
    auto persons = world.get_persons();
    for (auto& person : persons) {
        auto rng             = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);
        auto infection_state = get_infection_state(rng, exposed, infected_no_symptoms, infected_symptoms,
                                                   infected_severe, infected_critical, recovered);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         world.parameters, t0, infection_state,
                                                         person.get_latest_protection(), false),
                                     t0);
        }
    }
}

std::vector<mio::abm::LocationId> find_all_locations_of_type(mio::abm::World& world, mio::abm::LocationType type)
{
    std::vector<mio::abm::LocationId> locations;
    for (auto& loc : world.get_locations()) {
        if (loc.get_type() == type) {
            locations.push_back(mio::abm::LocationId{loc.get_index(), type});
        }
    }
    return locations;
}

/**
 * Assign locations to persons.
 * @param[in, out] world
*/
void assign_locations(mio::abm::World& world)
{
    //get locations from world
    //schools
    std::vector<mio::abm::LocationId> schools = find_all_locations_of_type(world, mio::abm::LocationType::School);
    std::vector<ScalarType> school_weights(schools.size(), 1);
    //hispitals
    std::vector<mio::abm::LocationId> hospitals = find_all_locations_of_type(world, mio::abm::LocationType::Hospital);
    std::vector<ScalarType> hospital_weights(hospitals.size(), 1);
    //icu
    std::vector<mio::abm::LocationId> icus = find_all_locations_of_type(world, mio::abm::LocationType::ICU);
    std::vector<ScalarType> icu_weights(icus.size(), 1);
    //workplaces
    std::vector<mio::abm::LocationId> workplaces = find_all_locations_of_type(world, mio::abm::LocationType::Work);
    std::vector<ScalarType> workplaces_weights(workplaces.size(), 1);
    //shops
    std::vector<mio::abm::LocationId> basic_shops =
        find_all_locations_of_type(world, mio::abm::LocationType::BasicsShop);
    std::vector<ScalarType> basic_shops_weights(basic_shops.size(), 1);
    //social events
    std::vector<mio::abm::LocationId> social_events =
        find_all_locations_of_type(world, mio::abm::LocationType::SocialEvent);
    std::vector<ScalarType> social_event_weights(social_events.size(), 1);

    auto persons = world.get_persons();
    for (auto& person : persons) {
        //assign shop
        size_t shop = mio::DiscreteDistribution<size_t>::get_instance()(world.get_rng(), basic_shops_weights);
        person.set_assigned_location(basic_shops[shop]);
        //assign hospital
        size_t hospital = mio::DiscreteDistribution<size_t>::get_instance()(world.get_rng(), hospital_weights);
        person.set_assigned_location(hospitals[hospital]);
        //assign icu
        size_t icu = mio::DiscreteDistribution<size_t>::get_instance()(world.get_rng(), icu_weights);
        person.set_assigned_location(icus[icu]);
        //assign event
        size_t event = mio::DiscreteDistribution<size_t>::get_instance()(world.get_rng(), social_event_weights);
        person.set_assigned_location(social_events[event]);
        //assign work and school
        if (person.get_age() == age_group_5_to_14) {
            size_t school = mio::DiscreteDistribution<size_t>::get_instance()(world.get_rng(), school_weights);
            person.set_assigned_location(schools[school]);
        }
        if (person.get_age() == age_group_15_to_34 || person.get_age() == age_group_35_to_59) {
            size_t work = mio::DiscreteDistribution<size_t>::get_instance()(world.get_rng(), workplaces_weights);
            person.set_assigned_location(workplaces[work]);
        }
    }
}

/**
 * Create sampled simulation with start time t0.
 * @param[in] t0 start time of the simulation
 * @param[in] input_dir text file with the locations and the inhabitants per residential area
 * @param[in, out] locationIds mapping of input area to abm locations
 * @return sim abm simulation
*/
mio::abm::Simulation create_sampled_simulation(const mio::abm::TimePoint& t0, const fs::path& input_dir,
                                               std::vector<LocationMapping>& locationIds)
{
    std::vector<std::pair<std::string, std::string>> areas;
    std::vector<int> inhabitants;
    mio::unused(locationIds);
    mio::unused(input_dir);
    // Read input file containing the locations and the number of inhabitants per location
    read_txt(areas, inhabitants, input_dir);

    //Assumed percentage for 1-Person, 2-Person, 3-Person, 4-Person and 5-Persons households
    ScalarType one_person_hh_pct = 0.37, two_person_hh_pct = 0.33, three_person_hh_pct = 0.15, four_person_hh_pct = 0.1,
               five_person_hh_pct = 0.05;
    // Assumed percentage of infection state at the beginning of the simulation.
    ScalarType exposed_pct = 0.005, infected_no_symptoms_pct = 0.001, infected_symptoms_pct = 0.001,
               infected_severe_pct = 0.0001, infected_critical_pct = 0.0, recovered_pct = 0.0;
    //Set global infection parameters
    auto world = mio::abm::World(num_age_groups);
    set_infection_parameters(world.parameters);

    //Transform the input location to correspondin abm location
    create_locations_from_input(areas, inhabitants, world, locationIds, one_person_hh_pct, two_person_hh_pct,
                                three_person_hh_pct, four_person_hh_pct, five_person_hh_pct);

    //Assign an infection state to every person
    assign_infection_states(world, t0, exposed_pct, infected_no_symptoms_pct, infected_symptoms_pct,
                            infected_severe_pct, infected_critical_pct, recovered_pct);
    // //Assign the locations to persons
    assign_locations(world);

    auto sim = mio::abm::Simulation(t0, std::move(world));
    return sim;
}

//Loggers used for output object

//time point logger
struct LogTimePoint : mio::LogAlways {
    using Type = double;
    static Type log(const mio::abm::Simulation& sim)
    {
        return sim.get_time().hours();
    }
};

//LocationId logger
struct LogLocationIds : mio::LogOnce {
    using Type = std::vector<std::tuple<mio::abm::LocationType, uint32_t>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        std::vector<std::tuple<mio::abm::LocationType, uint32_t>> location_ids{};
        for (auto&& location : sim.get_world().get_locations()) {
            location_ids.push_back(std::make_tuple(location.get_type(), location.get_index()));
        }
        return location_ids;
    }
};

//agent logger
struct LogPersonsPerLocationAndInfectionTime : mio::LogAlways {
    using Type = std::vector<std::tuple<mio::abm::LocationId, uint32_t, mio::abm::TimeSpan, mio::abm::InfectionState>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        std::vector<std::tuple<mio::abm::LocationId, uint32_t, mio::abm::TimeSpan, mio::abm::InfectionState>>
            location_ids_person{};
        for (auto&& person : sim.get_world().get_persons()) {
            location_ids_person.push_back(std::make_tuple(person.get_location().get_id(), person.get_person_id(),
                                                          person.get_time_since_transmission(),
                                                          person.get_infection_state(sim.get_time())));
        }
        return location_ids_person;
    }
};

void write_results_to_file(std::string path,
                           mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                        LogPersonsPerLocationAndInfectionTime>::WriteWrapper::Data& logg)
{
    auto location_ids = std::get<1>(logg);
    auto agents       = std::get<2>(logg);
    auto time_points  = std::get<0>(logg);

    std::string input;
    std::ofstream myfile(path);
    for (size_t loc_id_index = 0; loc_id_index < location_ids[0].size(); ++loc_id_index) {
        input = convert_loc_id_to_string(location_ids[0][loc_id_index]) + " " + std::to_string(time_points.size());
        for (size_t t = 0; t < time_points.size(); ++t) {
            auto a_per_loc = get_agents_per_location(location_ids[0][loc_id_index], agents[t]);
            input += " " + std::to_string(time_points[t]) + " " + std::to_string(a_per_loc.size());
            for (auto& agent : a_per_loc) {
                double time_since_transmission;
                if (std::get<1>(agent) > mio::abm::TimeSpan(std::numeric_limits<int>::max() / 4)) {
                    time_since_transmission = -1;
                }
                else {
                    time_since_transmission = std::get<1>(agent).hours();
                }
                input += " " + std::to_string(std::get<0>(agent)) + " " + std::to_string(time_since_transmission);
            }
        }
        myfile << input << "\n";
    }
    myfile.close();
}

void write_location_mapping_to_file(std::string path, std::vector<LocationMapping>& LocationIds)
{
    std::string input;
    std::ofstream myfile(path);
    for (auto& id : LocationIds) {
        input = id.inputId + " ";
        for (auto& model_id : id.modelId) {
            input += model_id + " ";
        }
        myfile << input << "\n";
    }

    myfile.close();
}

std::map<int, std::vector<std::string>> initialize_model(mio::abm::Model& model, std::string person_file)
{

    std::map<int, std::vector<std::string>> loc_area_mapping;
    std::map<int, mio::abm::LocationId> locations;

    const fs::path p = filename;
    if (!fs::exists(p)) {
        mio::log_error("Cannot read in data. File does not exist.");
    }
    // File pointer
    std::fstream fin;

    // Open an existing file
    fin.open(filename, std::ios::in);
    std::vector<int32_t> row;
    std::vector<std::string> row_string;
    std::string line;

    // Read the Titles from the Data file
    std::getline(fin, line);
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    std::vector<std::string> titles;
    boost::split(titles, line, boost::is_any_of(","));
    uint32_t count_of_titles              = 0;
    std::map<std::string, uint32_t> index = {};
    for (auto const& title : titles) {
        index.insert({title, count_of_titles});
        row_string.push_back(title);
        count_of_titles++;
    }

    while (std::getline(fin, line)) {
        row.clear();

        // read columns in this row
        split_line(line, &row);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        uint32_t age = row[index["age"]];

        int home_id   = row[index["home_id"]];
        int home_zone = row[index["home_zone"]];

        auto iter_home = locations.find(home_id);
        if (iter_home == locations.end()) {
            home = model.add_location(mio::abm::LocationType::Home);
            locations.insert({home_id, home});
            std::string loc = "0" + std::to_string(mio::abm::LocationType::Home) + std::to_string(home.get());
            auto zone_iter  = loc_area_mapping.find(home_zone);
            if (zone_iter == loc_area_mapping.end()) {
                loc_area_mapping.insert({home_zone, {loc}});
            }
            else {
                loc_area_mapping[home_zone].push_back(loc);
            }
        }
        else {
            home = locations[home_id];
        }
        auto pid     = model.add_person(home, determine_age_group(age));
        auto& person = model.get_person(pid);
        person.set_assigned_location(mio::abm::LocationType::Home, home);

        int shop_id   = row[index["shop_id"]];
        int shop_zone = row[index["shop_zone"]];

        auto iter_shop = locations.find(shop_id);
        if (iter_shop == locations.end()) {
            shop = model.add_location(mio::abm::LocationType::BasicsShop);
            if (shop_id = !-1) {
                locations.insert({shop_id, shop});
            }
            std::string loc = "0" + std::to_string(mio::abm::LocationType::BasicsShop) + std::to_string(shop.get());
            auto zone_iter  = loc_area_mapping.find(shop_zone);
            if (zone_iter == loc_area_mapping.end()) {
                loc_area_mapping.insert({shop_zone, {loc}});
            }
            else {
                loc_area_mapping[shop_zone].push_back(loc);
            }
        }
        else {
            shop = locations[shop_id];
        }
        person.set_assigned_location(mio::abm::LocationType::BasicsShop, shop);

        int event_id   = row[index["event_id"]];
        int event_zone = row[index["event_zone"]];

        auto iter_event = locations.find(event_id);
        if (iter_event == locations.end()) {
            event = model.add_location(mio::abm::LocationType::SocialEvent);
            if (event_id != -1) {
                locations.insert({event_id, event});
            }
            std::string loc = "0" + std::to_string(mio::abm::LocationType::SocialEvent) + std::to_string(event.get());
            auto zone_iter  = loc_area_mapping.find(event_zone);
            if (zone_iter == loc_area_mapping.end()) {
                loc_area_mapping.insert({event_zone, {loc}});
            }
            else {
                loc_area_mapping[event_zone].push_back(loc);
            }
        }
        else {
            event = locations[event_id];
        }
        person.set_assigned_location(mio::abm::LocationType::SocialEvent, event);

        if (person.get_age() == mio::AgeGroup(1)) {
            int school_id   = row[index["school_id"]];
            int school_zone = row[index["school_zone"]];

            auto iter_school = locations.find(school_id);
            if (iter_school == locations.end()) {
                school = model.add_location(mio::abm::LocationType::School);
                if (school_id != -1) {
                    locations.insert({school_id, school});
                }
                std::string loc = "0" + std::to_string(mio::abm::LocationType::School) + std::to_string(school.get());
                auto zone_iter  = loc_area_mapping.find(school_zone);
                if (zone_iter == loc_area_mapping.end()) {
                    loc_area_mapping.insert({school_zone, {loc}});
                }
                else {
                    loc_area_mapping[school_zone].push_back(loc);
                }
            }
            else {
                school = locations[school_id];
            }
            person.set_assigned_location(mio::abm::LocationType::School, school);
        }

        if (person.get_age() == mio::AgeGroup(2) || person.get_age() == mio::AgeGroup(3)) {
            int work_id   = row[index["work_id"]];
            int work_zone = row[index["work_zone"]];

            auto iter_work = locations.find(work_id);
            if (iter_work == locations.end()) {
                work = model.add_location(mio::abm::LocationType::Work);
                if (work_id != -1) {
                    locations.insert({work_id, work});
                }
                std::string loc = "0" + std::to_string(mio::abm::LocationType::Work) + std::to_string(work.get());
                auto zone_iter  = loc_area_mapping.find(work_zone);
                if (zone_iter == loc_area_mapping.end()) {
                    loc_area_mapping.insert({work_zone, {loc}});
                }
                else {
                    loc_area_mapping[work_zone].push_back(loc);
                }
            }
            else {
                work = locations[work_id];
            }
            person.set_assigned_location(mio::abm::LocationType::Work, work);
        }
    }

    return loc_area_mapping;
}

mio::IOResult<void> run(const fs::path& input_dir)
{
    mio::set_log_level(mio::LogLevel::warn);
    // auto t0   = mio::abm::TimePoint(0); // Start time per simulation
    // auto tmax = mio::abm::TimePoint(0) + mio::abm::days(14); // End time per simulation

    auto model = mio::abm::Model(size_t(6));
    auto dict  = initialize_model(model, "../../pycode/examples/simulation/ABM Demonstrator/input/persons.csv");

    mio::unused(model);
    mio::unused(dict);
    // //mapping of input areas to abm locations
    // std::vector<LocationMapping> LocationIds;
    // //create sampled simulation
    // auto sim = create_sampled_simulation(t0, input_dir, LocationIds);

    // //output object
    // mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds, LogPersonsPerLocationAndInfectionTime> history;

    // //advance until tmax
    // sim.advance(tmax, history);

    // //output
    // auto logg = history.get_log();
    // write_results_to_file("output_abm_demonstrator.txt", logg);
    // write_location_mapping_to_file("location_mapping.txt", LocationIds);

    // std::cout << "# t S E C I I_s I_c R_C R_I D\n";
    // for (auto i = 0; i < sim.get_result().get_num_time_points(); ++i) {
    //     std::cout << sim.get_result().get_time(i) << " ";
    //     auto v = sim.get_result().get_value(i);
    //     for (auto j = 0; j < v.size(); ++j) {
    //         std::cout << v[j] << " ";
    //         if (j < v.size() - 1) {
    //             std::cout << " ";
    //         }
    //     }
    //     if (i < sim.get_result().get_num_time_points() - 1) {
    //         std::cout << "\n";
    //     }
    // }

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
