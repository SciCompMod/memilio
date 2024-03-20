/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Sascha Korf, Carlotta Gerstein
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

#include <fstream>
#include <vector>
#include <iostream>
#include <chrono>
#include "abm/abm.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/uncertain_value.h"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "abm/vaccine.h"
#include "abm/common_abm_loggers.h"
#include "generate_graph_from_data.cpp"
#include "memilio/utils/miompi.h"
#include "memilio/io/binary_serializer.h"

namespace fs = boost::filesystem;

// Assign the name to general age group.
size_t num_age_groups         = 6;
const auto age_group_0_to_4   = mio::AgeGroup(0);
const auto age_group_5_to_14  = mio::AgeGroup(1);
const auto age_group_15_to_34 = mio::AgeGroup(2);
const auto age_group_35_to_59 = mio::AgeGroup(3);
const auto age_group_60_to_79 = mio::AgeGroup(4);
const auto age_group_80_plus  = mio::AgeGroup(5);

const std::map<mio::osecir::InfectionState, mio::abm::InfectionState> infection_state_map{
    {mio::osecir::InfectionState::Susceptible, mio::abm::InfectionState::Susceptible},
    {mio::osecir::InfectionState::Exposed, mio::abm::InfectionState::Exposed},
    {mio::osecir::InfectionState::InfectedNoSymptoms, mio::abm::InfectionState::InfectedNoSymptoms},
    {mio::osecir::InfectionState::InfectedNoSymptomsConfirmed, mio::abm::InfectionState::InfectedNoSymptoms},
    {mio::osecir::InfectionState::InfectedSymptoms, mio::abm::InfectionState::InfectedSymptoms},
    {mio::osecir::InfectionState::InfectedSymptomsConfirmed, mio::abm::InfectionState::InfectedSymptoms},
    {mio::osecir::InfectionState::InfectedSevere, mio::abm::InfectionState::InfectedSevere},
    {mio::osecir::InfectionState::InfectedCritical, mio::abm::InfectionState::InfectedCritical},
    {mio::osecir::InfectionState::Recovered, mio::abm::InfectionState::Recovered},
    {mio::osecir::InfectionState::Dead, mio::abm::InfectionState::Dead}};

mio::CustomIndexArray<double, mio::AgeGroup, mio::osecir::InfectionState> initial_infection_distribution{
    {mio::AgeGroup(num_age_groups), mio::osecir::InfectionState::Count}, 0.02};

/**
 * Determine initial distribution of infection states.
*/
void determine_initial_infection_states_world(const fs::path& input_dir, const mio::Date date)
{
    // estimate intial population by ODE compartiments
    auto initial_graph                     = get_graph(date, 1, input_dir);
    size_t braunschweig_id                 = 16; // Braunschweig has ID 16
    auto braunschweig_node                 = initial_graph.value()[braunschweig_id];
    initial_infection_distribution.array() = braunschweig_node.populations.array().cast<double>();
}

/**
 * Assign an infection state to each person according to real world data read in through the ODE secir model.
 */
void assign_infection_state(mio::abm::World& world, mio::abm::TimePoint t)
{
    // convert initial population to ABM initial infections
    for (auto& person : world.get_persons()) {
        auto rng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);

        auto infection_state = mio::osecir::InfectionState(mio::DiscreteDistribution<size_t>::get_instance()(
            rng, initial_infection_distribution.slice(person.get_age()).as_array().array()));

        //bool detected = false;
        if (infection_state != mio::osecir::InfectionState::Susceptible) {
            //if (infection_state == mio::osecir::InfectionState::InfectedNoSymptomsConfirmed ||
            //    infection_state == mio::osecir::InfectionState::InfectedSymptomsConfirmed) {
            //    detected = true;
            //}
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         world.parameters, t, infection_state_map.at(infection_state)));
        }
    }
}

int stringToMinutes(const std::string& input)
{
    size_t colonPos = input.find(":");
    if (colonPos == std::string::npos) {
        // Handle invalid input (no colon found)
        return -1; // You can choose a suitable error code here.
    }

    std::string xStr = input.substr(0, colonPos);
    std::string yStr = input.substr(colonPos + 1);

    int x = std::stoi(xStr);
    int y = std::stoi(yStr);
    return x * 60 + y;
}

int longLatToInt(const std::string& input)
{
    double y = std::stod(input) * 1e+5; //we want the 5 numbers after digit
    return (int)y;
}
void split_line(std::string string, std::vector<int32_t>* row)
{
    std::vector<std::string> strings;

    std::string x = ",,", y = ",-1,";
    size_t pos;
    while ((pos = string.find(x)) != std::string::npos) {
        string.replace(pos, 2, y);
    } // Temporary fix to handle empty cells.
    boost::split(strings, string, boost::is_any_of(","));
    std::transform(strings.begin(), strings.end(), std::back_inserter(*row), [&](std::string s) {
        if (s.find(":") != std::string::npos) {
            return stringToMinutes(s);
        }
        else if (s.find(".") != std::string::npos) {
            return longLatToInt(s);
        }
        else if (s == "null") {
            return 10; // This shouldnt be too often, just assume a short time after 12:00 o'clock for now
        }
        else {
            return std::stoi(s);
        }
    });
}

mio::abm::LocationType get_location_type(uint32_t acitivity_end)
{
    mio::abm::LocationType type;
    switch (acitivity_end) {
    case 1:
        type = mio::abm::LocationType::Work;
        break;
    case 2:
        type = mio::abm::LocationType::School;
        break;
    case 3:
        type = mio::abm::LocationType::BasicsShop;
        break;
    case 4:
        type = mio::abm::LocationType::SocialEvent; // Freizeit
        break;
    case 5:
        type = mio::abm::LocationType::BasicsShop; // Private Erledigung
        break;
    case 6:
        type = mio::abm::LocationType::SocialEvent; // Sonstiges
        break;
    default:
        type = mio::abm::LocationType::Home;
        break;
    }
    return type;
}

mio::AgeGroup determine_age_group(uint32_t age)
{
    if (age <= 4) {
        return age_group_0_to_4;
    }
    else if (age <= 14) {
        return age_group_5_to_14;
    }
    else if (age <= 34) {
        return age_group_15_to_34;
    }
    else if (age <= 59) {
        return age_group_35_to_59;
    }
    else if (age <= 79) {
        return age_group_60_to_79;
    }
    else {
        return age_group_80_plus;
    }
}

void create_world_from_data(mio::abm::World& world, const std::string& filename, const mio::abm::TimePoint t0,
                            int max_number_persons)
{
    // Open File
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

    std::map<uint32_t, mio::abm::LocationId> locations = {};
    std::map<uint32_t, mio::abm::Person&> persons      = {};
    std::vector<std::tuple<int, uint32_t, uint32_t, uint32_t>> schools =
        {}; //schools need to save the id, the zone it is in and the amount of persons who will be assigned to it, we need it double when we actually assign persons
    std::map<uint32_t, uint32_t> person_ids = {};
    std::map<uint32_t, std::pair<uint32_t, int>> locations_before;
    std::map<uint32_t, std::pair<uint32_t, int>> locations_after;

    uint32_t max_number_of_persons_in_school = 600;

    // For the world we need: Hospitals, ICUs (for both we just create one for now), Homes for each unique householdID, One Person for each person_id with respective age and home_id.

    // We assume that no person goes to an hospital, altough e.g. "Sonstiges" could be a hospital
    auto hospital = world.add_location(mio::abm::LocationType::Hospital);
    world.get_individualized_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    world.get_individualized_location(hospital).set_capacity(std::numeric_limits<uint32_t>::max(),
                                                             std::numeric_limits<uint32_t>::max());
    auto icu = world.add_location(mio::abm::LocationType::ICU);
    world.get_individualized_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    world.get_individualized_location(icu).set_capacity(std::numeric_limits<uint32_t>::max(),
                                                        std::numeric_limits<uint32_t>::max());

    // First we determine the persons number and their starting locations
    int number_of_persons = 0;

    while (std::getline(fin, line)) {
        row.clear();

        // read columns in this row
        split_line(line, &row);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        uint32_t person_id = row[index["puid"]];
        auto it_person_id  = person_ids.find(person_id);
        if (it_person_id == person_ids.end()) {
            if (number_of_persons >= max_number_persons)
                break; //This is okay because the data is sorted by person_id
            person_ids.insert({person_id, number_of_persons});
            number_of_persons++;
        }

        // The starting location of a person is the end location of the last trip he made, either on the same day or on
        // the day before
        uint32_t target_location_id = std::abs(row[index["loc_id_end"]]);
        int trip_start              = row[index["start_time"]];
        if (trip_start < t0.hour_of_day()) {
            auto it_person = locations_before.find(person_id);
            if (it_person == locations_before.end()) {
                locations_before.insert({person_id, std::make_pair(target_location_id, trip_start)});
            }
            else {
                if (it_person->second.second <= trip_start) {
                    it_person->second.first  = target_location_id;
                    it_person->second.second = trip_start;
                }
            }
        }
        else {
            auto it_person = locations_after.find(person_id);
            if (it_person == locations_after.end()) {
                locations_after.insert({person_id, std::make_pair(target_location_id, trip_start)});
            }
            else {
                if (it_person->second.second <= trip_start) {
                    it_person->second.first  = target_location_id;
                    it_person->second.second = trip_start;
                }
            }
        }
    }

    fin.clear();
    fin.seekg(0);
    std::getline(fin, line); // Skip header row

    // Add all locations to the world
    int school_id_counter   = -100;
    uint32_t last_person_id = 0;
    while (std::getline(fin, line)) {
        row.clear();

        // read columns in this row
        split_line(line, &row);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        uint32_t person_id = row[index["puid"]];
        if (person_ids.find(person_id) == person_ids.end())
            break;

        uint32_t home_id                                 = row[index["huid"]];
        int target_location_id                           = row[index["loc_id_end"]];
        uint32_t end_zone                                = row[index["end_zone"]];
        uint32_t activity_end                            = row[index["activity_end"]];
        mio::abm::GeographicalLocation location_long_lat = {(double)row[index["lon_end"]] / 1e+5,
                                                            (double)row[index["lat_end"]] / 1e+5};
        mio::abm::LocationId home;
        auto it_home = locations.find(home_id);
        if (it_home == locations.end()) {
            home = world.add_location(mio::abm::LocationType::Home, 1);
            locations.insert({home_id, home});
        }
        else {
            home = it_home->second;
        }

        mio::abm::LocationId location;

        switch (target_location_id) {
        case -1: //home has to be anonymized
            target_location_id = home_id;
        case -2: //schools have to be anonymized we place them into a general school and every 600 persons into a new school
        {
            if (last_person_id != person_id) {
                bool found_school = false;
                for (auto& school : schools) {
                    if (std::get<1>(school) == end_zone && std::get<2>(school) < max_number_of_persons_in_school) {
                        target_location_id = std::get<0>(school);
                        std::get<2>(school)++;
                        last_person_id = person_id;
                        found_school   = true;
                        break;
                    }
                }
                if (!found_school) {
                    target_location_id = school_id_counter--;
                    schools.push_back(std::make_tuple(target_location_id, end_zone, 1, 0));
                }
            }
        }
        default:
            break;
        }

        auto it_location = locations.find(
            target_location_id); // Check if location already exists also for home which have the same id (home_id = target_location_id)
        if (it_location == locations.end()) {
            location = world.add_location(
                get_location_type(activity_end),
                1); // Assume one place has one activity, this may be untrue but not important for now(?)
            locations.insert({target_location_id, location});
            world.get_individualized_location(location).set_geographical_location(location_long_lat);
        }
    }
    fin.clear();
    fin.seekg(0);
    std::getline(fin, line); // Skip header row
    std::unordered_set<uint32_t> ids_in_bs;

    // Add the persons and trips
    last_person_id = 0;
    while (std::getline(fin, line)) {
        row.clear();

        // read columns in this row
        split_line(line, &row);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        uint32_t person_id = row[index["puid"]];
        if (person_ids.find(person_id) == person_ids.end())
            break;

        uint32_t age            = row[index["age"]];
        uint32_t home_id        = row[index["huid"]];
        int target_location_id  = row[index["loc_id_end"]];
        int start_location_id   = row[index["loc_id_start"]];
        uint32_t end_zone       = row[index["end_zone"]];
        uint32_t trip_start     = row[index["start_time"]];
        uint32_t transport_mode = row[index["travel_mode"]];
        uint32_t acticity_end   = row[index["activity_end"]];
        bool home_in_bs         = true;

        switch (target_location_id) {
        case -1: //home has to be anonymized
            target_location_id = home_id;
        case -2: //schools have to be anonymized we place them into a general school and every 600 persons into a new school
            for (auto& school : schools) {
                if (std::get<1>(school) == end_zone && std::get<3>(school) < max_number_of_persons_in_school) {
                    if (last_person_id != person_id) {
                        last_person_id = person_id;
                        std::get<3>(school)++;
                    }
                    target_location_id = std::get<0>(school);
                    break;
                }
            }
        default:
            break;
        }

        // Add the trip to the trip list person and location must exist at this point
        auto target_location = locations.find(target_location_id)->second;
        auto start_location  = locations.find(start_location_id)->second;

        auto it_person = persons.find(person_id);

        if (it_person == persons.end()) {
            auto it_first_location_id = locations_before.find(person_id);
            if (it_first_location_id == locations_before.end()) {
                it_first_location_id = locations_after.find(person_id);
            }
            auto first_location_id = it_first_location_id->second.first;
            auto first_location    = locations.find(first_location_id)->second;
            auto& person           = world.add_person(first_location, determine_age_group(age));
            if (home_in_bs) {
                ids_in_bs.insert(person_id);
            }
            auto home = locations.find(home_id)->second;
            person.set_assigned_location(home);
            person.set_assigned_location(hospital);
            person.set_assigned_location(icu);
            persons.insert({person_id, person});
            it_person = persons.find(person_id);
        }

        it_person->second.set_assigned_location(
            target_location); //This assumes that we only have in each tripchain only one location type for each person
        if (locations.find(start_location_id) == locations.end()) {
            // For trips where the start location is not known use Home instead
            start_location = {it_person->second.get_assigned_location_index(mio::abm::LocationType::Home),
                              mio::abm::LocationType::Home};
        }
        world.get_trip_list().add_trip(mio::abm::Trip(
            it_person->second.get_person_id(), mio::abm::TimePoint(0) + mio::abm::minutes(trip_start), target_location,
            start_location, mio::abm::TransportMode(transport_mode), mio::abm::ActivityType(acticity_end)));
    }
    world.get_trip_list().use_weekday_trips_on_weekend();
    world.parameters.get<mio::abm::LogAgentIds>() = ids_in_bs;

    // //Some Data about the world:
    // //write how many persons are in each home
    // std::map<uint32_t, uint32_t> persons_in_home;
    // for (auto& person : world.get_persons()) {
    //     auto home = person.get_assigned_location_index(mio::abm::LocationType::Home);
    //     if (persons_in_home.find(home) == persons_in_home.end()) {
    //         persons_in_home.insert({home, 1});
    //     }
    //     else {
    //         persons_in_home[home]++;
    //     }
    // }
    // for (auto& home : persons_in_home) {
    //     mio::log_info("Home ", home.first, " has ", home.second, " persons.");
    // }
}

std::pair<double, double> get_my_and_sigma(std::pair<double, double> mean_and_std)
{
    auto mean    = mean_and_std.first;
    auto stddev  = mean_and_std.second;
    double my    = log(mean * mean / sqrt(mean * mean + stddev * stddev));
    double sigma = sqrt(log(1 + stddev * stddev / (mean * mean)));
    return {my, sigma};
}

void set_parameters(mio::abm::Parameters& params)
{
    mio::RandomNumberGenerator rng;

    // Set the Time parameters for the infection same for every age group for now

    auto incubation_period_my_sigma          = get_my_and_sigma({4.5, 1.5});
    params.get<mio::abm::IncubationPeriod>() = {incubation_period_my_sigma.first, incubation_period_my_sigma.second};

    auto InfectedNoSymptoms_to_symptoms_my_sigma             = get_my_and_sigma({1.1, 0.9});
    params.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>() = {InfectedNoSymptoms_to_symptoms_my_sigma.first,
                                                                InfectedNoSymptoms_to_symptoms_my_sigma.second};

    auto TimeInfectedNoSymptomsToRecovered_my_sigma           = get_my_and_sigma({8.0, 2.0});
    params.get<mio::abm::TimeInfectedNoSymptomsToRecovered>() = {TimeInfectedNoSymptomsToRecovered_my_sigma.first,
                                                                 TimeInfectedNoSymptomsToRecovered_my_sigma.second};

    auto TimeInfectedSymptomsToSevere_my_sigma           = get_my_and_sigma({6.6, 4.9});
    params.get<mio::abm::TimeInfectedSymptomsToSevere>() = {TimeInfectedSymptomsToSevere_my_sigma.first,
                                                            TimeInfectedSymptomsToSevere_my_sigma.second};

    auto TimeInfectedSymptomsToRecovered_my_sigma           = get_my_and_sigma({8.0, 2.0});
    params.get<mio::abm::TimeInfectedSymptomsToRecovered>() = {TimeInfectedSymptomsToRecovered_my_sigma.first,
                                                               TimeInfectedSymptomsToRecovered_my_sigma.second};

    auto TimeInfectedSevereToCritical_my_sigma           = get_my_and_sigma({1.5, 2.0});
    params.get<mio::abm::TimeInfectedSevereToCritical>() = {TimeInfectedSevereToCritical_my_sigma.first,
                                                            TimeInfectedSevereToCritical_my_sigma.second};

    auto TimeInfectedSevereToRecovered_my_sigma           = get_my_and_sigma({18.1, 6.3});
    params.get<mio::abm::TimeInfectedSevereToRecovered>() = {TimeInfectedSevereToRecovered_my_sigma.first,
                                                             TimeInfectedSevereToRecovered_my_sigma.second};

    auto TimeInfectedCriticalToDead_my_sigma           = get_my_and_sigma({10.7, 4.8});
    params.get<mio::abm::TimeInfectedCriticalToDead>() = {TimeInfectedCriticalToDead_my_sigma.first,
                                                          TimeInfectedCriticalToDead_my_sigma.second};

    auto TimeInfectedCriticalToRecovered_my_sigma           = get_my_and_sigma({18.1, 6.3});
    params.get<mio::abm::TimeInfectedCriticalToRecovered>() = {TimeInfectedCriticalToRecovered_my_sigma.first,
                                                               TimeInfectedCriticalToRecovered_my_sigma.second};

    // Set percentage parameters
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>() = 0.5;
    params.get<mio::abm::SeverePerInfectedSymptoms>()     = 0.1;
    params.get<mio::abm::CriticalPerInfectedSevere>()     = 0.05;
    params.get<mio::abm::DeathsPerInfectedCritical>()     = 0.002;

    // Set infection parameters
    // Set protection level from high viral load. Information based on: https://doi.org/10.1093/cid/ciaa886
    params.get<mio::abm::HighViralLoadProtectionFactor>() = [](ScalarType days) -> ScalarType {
        return mio::linear_interpolation_of_data_set<ScalarType, ScalarType>(
            {{0, 0.863}, {1, 0.969}, {7, 0.029}, {10, 0.002}, {14, 0.0014}, {21, 0}}, days);
    };

    // Set protection level against an severe infection. Information based on: https://doi.org/10.1093/cid/ciaa886
    params.get<mio::abm::SeverityProtectionFactor>() = [](ScalarType days) -> ScalarType {
        return mio::linear_interpolation_of_data_set<ScalarType, ScalarType>({{0, 0.967},
                                                                              {30, 0.975},
                                                                              {60, 0.977},
                                                                              {90, 0.974},
                                                                              {120, 0.963},
                                                                              {150, 0.947},
                                                                              {180, 0.93},
                                                                              {210, 0.929},
                                                                              {240, 0.923},
                                                                              {270, 0.908},
                                                                              {300, 0.893},
                                                                              {330, 0.887},
                                                                              {360, 0.887},
                                                                              {360, 0.5}},
                                                                             days);
    };

    //Set other parameters
    params.get<mio::abm::MaskProtection>()           = 0.66;
    params.get<mio::abm::AerosolTransmissionRates>() = 0.0;
}

// set location specific parameters
void set_local_parameters(mio::abm::World& world)
{
    const int n_age_groups = world.parameters.get_num_groups();

    // setting this up in matrix-form would be much nicer,
    // but we somehow can't construct Eigen object with initializer lists
    /* baseline_home
        0.4413 0.4504 1.2383 0.8033 0.0494 0.0017
        0.0485 0.7616 0.6532 1.1614 0.0256 0.0013
        0.1800 0.1795 0.8806 0.6413 0.0429 0.0032
        0.0495 0.2639 0.5189 0.8277 0.0679 0.0014
        0.0087 0.0394 0.1417 0.3834 0.7064 0.0447
        0.0292 0.0648 0.1248 0.4179 0.3497 0.1544
    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_home(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_home[{age_group_0_to_4, age_group_0_to_4}]     = 0.4413;
    contacts_home[{age_group_0_to_4, age_group_5_to_14}]    = 0.0504;
    contacts_home[{age_group_0_to_4, age_group_15_to_34}]   = 1.2383;
    contacts_home[{age_group_0_to_4, age_group_35_to_59}]   = 0.8033;
    contacts_home[{age_group_0_to_4, age_group_60_to_79}]   = 0.0494;
    contacts_home[{age_group_0_to_4, age_group_80_plus}]    = 0.0017;
    contacts_home[{age_group_5_to_14, age_group_0_to_4}]    = 0.0485;
    contacts_home[{age_group_5_to_14, age_group_5_to_14}]   = 0.7616;
    contacts_home[{age_group_5_to_14, age_group_15_to_34}]  = 0.6532;
    contacts_home[{age_group_5_to_14, age_group_35_to_59}]  = 1.1614;
    contacts_home[{age_group_5_to_14, age_group_60_to_79}]  = 0.0256;
    contacts_home[{age_group_5_to_14, age_group_80_plus}]   = 0.0013;
    contacts_home[{age_group_15_to_34, age_group_0_to_4}]   = 0.1800;
    contacts_home[{age_group_15_to_34, age_group_5_to_14}]  = 0.1795;
    contacts_home[{age_group_15_to_34, age_group_15_to_34}] = 0.8806;
    contacts_home[{age_group_15_to_34, age_group_35_to_59}] = 0.6413;
    contacts_home[{age_group_15_to_34, age_group_60_to_79}] = 0.0429;
    contacts_home[{age_group_15_to_34, age_group_80_plus}]  = 0.0032;
    contacts_home[{age_group_35_to_59, age_group_0_to_4}]   = 0.0495;
    contacts_home[{age_group_35_to_59, age_group_5_to_14}]  = 0.2639;
    contacts_home[{age_group_35_to_59, age_group_15_to_34}] = 0.5189;
    contacts_home[{age_group_35_to_59, age_group_35_to_59}] = 0.8277;
    contacts_home[{age_group_35_to_59, age_group_60_to_79}] = 0.0679;
    contacts_home[{age_group_35_to_59, age_group_80_plus}]  = 0.0014;
    contacts_home[{age_group_60_to_79, age_group_0_to_4}]   = 0.0087;
    contacts_home[{age_group_60_to_79, age_group_5_to_14}]  = 0.0394;
    contacts_home[{age_group_60_to_79, age_group_15_to_34}] = 0.1417;
    contacts_home[{age_group_60_to_79, age_group_35_to_59}] = 0.3834;
    contacts_home[{age_group_60_to_79, age_group_60_to_79}] = 0.7064;
    contacts_home[{age_group_60_to_79, age_group_80_plus}]  = 0.0447;
    contacts_home[{age_group_80_plus, age_group_0_to_4}]    = 0.0292;
    contacts_home[{age_group_80_plus, age_group_5_to_14}]   = 0.0648;
    contacts_home[{age_group_80_plus, age_group_15_to_34}]  = 0.1248;
    contacts_home[{age_group_80_plus, age_group_35_to_59}]  = 0.4179;
    contacts_home[{age_group_80_plus, age_group_60_to_79}]  = 0.3497;
    contacts_home[{age_group_80_plus, age_group_80_plus}]   = 0.1544;

    /* baseline_school
        1.1165 0.2741 0.2235 0.1028 0.0007 0.0000
        0.1627 1.9412 0.2431 0.1780 0.0130 0.0000
        0.0148 0.1646 1.1266 0.0923 0.0074 0.0000
        0.0367 0.1843 0.3265 0.0502 0.0021 0.0005
        0.0004 0.0370 0.0115 0.0014 0.0039 0.0000
        0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_school(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_school[{age_group_0_to_4, age_group_0_to_4}]     = 1.1165;
    contacts_school[{age_group_0_to_4, age_group_5_to_14}]    = 0.2741;
    contacts_school[{age_group_0_to_4, age_group_15_to_34}]   = 0.2235;
    contacts_school[{age_group_0_to_4, age_group_35_to_59}]   = 0.1028;
    contacts_school[{age_group_0_to_4, age_group_60_to_79}]   = 0.0007;
    contacts_school[{age_group_0_to_4, age_group_80_plus}]    = 0.0000;
    contacts_school[{age_group_5_to_14, age_group_0_to_4}]    = 0.1627;
    contacts_school[{age_group_5_to_14, age_group_5_to_14}]   = 1.9412;
    contacts_school[{age_group_5_to_14, age_group_15_to_34}]  = 0.2431;
    contacts_school[{age_group_5_to_14, age_group_35_to_59}]  = 0.1780;
    contacts_school[{age_group_5_to_14, age_group_60_to_79}]  = 0.0130;
    contacts_school[{age_group_5_to_14, age_group_80_plus}]   = 0.0000;
    contacts_school[{age_group_15_to_34, age_group_0_to_4}]   = 0.0148;
    contacts_school[{age_group_15_to_34, age_group_5_to_14}]  = 0.1646;
    contacts_school[{age_group_15_to_34, age_group_15_to_34}] = 1.1266;
    contacts_school[{age_group_15_to_34, age_group_35_to_59}] = 0.0923;
    contacts_school[{age_group_15_to_34, age_group_60_to_79}] = 0.0074;
    contacts_school[{age_group_15_to_34, age_group_80_plus}]  = 0.0000;
    contacts_school[{age_group_35_to_59, age_group_0_to_4}]   = 0.0367;
    contacts_school[{age_group_35_to_59, age_group_5_to_14}]  = 0.1843;
    contacts_school[{age_group_35_to_59, age_group_15_to_34}] = 0.3265;
    contacts_school[{age_group_35_to_59, age_group_35_to_59}] = 0.0502;
    contacts_school[{age_group_35_to_59, age_group_60_to_79}] = 0.0021;
    contacts_school[{age_group_35_to_59, age_group_80_plus}]  = 0.0005;
    contacts_school[{age_group_60_to_79, age_group_0_to_4}]   = 0.0004;
    contacts_school[{age_group_60_to_79, age_group_5_to_14}]  = 0.0370;
    contacts_school[{age_group_60_to_79, age_group_15_to_34}] = 0.0115;
    contacts_school[{age_group_60_to_79, age_group_35_to_59}] = 0.0014;
    contacts_school[{age_group_60_to_79, age_group_60_to_79}] = 0.0039;
    contacts_school[{age_group_60_to_79, age_group_80_plus}]  = 0.0000;
    contacts_school[{age_group_80_plus, age_group_0_to_4}]    = 0.0000;
    contacts_school[{age_group_80_plus, age_group_5_to_14}]   = 0.0000;
    contacts_school[{age_group_80_plus, age_group_15_to_34}]  = 0.0000;
    contacts_school[{age_group_80_plus, age_group_35_to_59}]  = 0.0000;
    contacts_school[{age_group_80_plus, age_group_60_to_79}]  = 0.0000;
    contacts_school[{age_group_80_plus, age_group_80_plus}]   = 0.0000;

    /* baseline_work
        0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
        0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
        0.0000 0.0127 1.7570 1.6050 0.0133 0.0000
        0.0000 0.0020 1.0311 2.3166 0.0098 0.0000
        0.0000 0.0002 0.0194 0.0325 0.0003 0.0000
        0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_work(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_work[{age_group_0_to_4, age_group_0_to_4}]     = 0.0000;
    contacts_work[{age_group_0_to_4, age_group_5_to_14}]    = 0.0000;
    contacts_work[{age_group_0_to_4, age_group_15_to_34}]   = 0.0000;
    contacts_work[{age_group_0_to_4, age_group_35_to_59}]   = 0.0000;
    contacts_work[{age_group_0_to_4, age_group_60_to_79}]   = 0.0000;
    contacts_work[{age_group_0_to_4, age_group_80_plus}]    = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_0_to_4}]    = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_5_to_14}]   = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_15_to_34}]  = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_35_to_59}]  = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_60_to_79}]  = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_80_plus}]   = 0.0000;
    contacts_work[{age_group_15_to_34, age_group_0_to_4}]   = 0.0000;
    contacts_work[{age_group_15_to_34, age_group_5_to_14}]  = 0.0127;
    contacts_work[{age_group_15_to_34, age_group_15_to_34}] = 1.7570;
    contacts_work[{age_group_15_to_34, age_group_35_to_59}] = 1.6050;
    contacts_work[{age_group_15_to_34, age_group_60_to_79}] = 0.0133;
    contacts_work[{age_group_15_to_34, age_group_80_plus}]  = 0.0000;
    contacts_work[{age_group_35_to_59, age_group_0_to_4}]   = 0.0000;
    contacts_work[{age_group_35_to_59, age_group_5_to_14}]  = 0.0020;
    contacts_work[{age_group_35_to_59, age_group_15_to_34}] = 1.0311;
    contacts_work[{age_group_35_to_59, age_group_35_to_59}] = 2.3166;
    contacts_work[{age_group_35_to_59, age_group_60_to_79}] = 0.0098;
    contacts_work[{age_group_35_to_59, age_group_80_plus}]  = 0.0000;
    contacts_work[{age_group_60_to_79, age_group_0_to_4}]   = 0.0000;
    contacts_work[{age_group_60_to_79, age_group_5_to_14}]  = 0.0002;
    contacts_work[{age_group_60_to_79, age_group_15_to_34}] = 0.0194;
    contacts_work[{age_group_60_to_79, age_group_35_to_59}] = 0.0325;
    contacts_work[{age_group_60_to_79, age_group_60_to_79}] = 0.0003;
    contacts_work[{age_group_60_to_79, age_group_80_plus}]  = 0.0000;
    contacts_work[{age_group_80_plus, age_group_0_to_4}]    = 0.0000;
    contacts_work[{age_group_80_plus, age_group_5_to_14}]   = 0.0000;
    contacts_work[{age_group_80_plus, age_group_15_to_34}]  = 0.0000;
    contacts_work[{age_group_80_plus, age_group_35_to_59}]  = 0.0000;
    contacts_work[{age_group_80_plus, age_group_60_to_79}]  = 0.0000;
    contacts_work[{age_group_80_plus, age_group_80_plus}]   = 0.0000;

    /* baseline_other
        0.5170 0.3997 0.7957 0.9958 0.3239 0.0428
        0.0632 0.9121 0.3254 0.4731 0.2355 0.0148
        0.0336 0.1604 1.7529 0.8622 0.1440 0.0077
        0.0204 0.1444 0.5738 1.2127 0.3433 0.0178
        0.0371 0.0393 0.4171 0.9666 0.7495 0.0257
        0.0791 0.0800 0.3480 0.5588 0.2769 0.0180

    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_other(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_other[{age_group_0_to_4, age_group_0_to_4}]     = 0.5170;
    contacts_other[{age_group_0_to_4, age_group_5_to_14}]    = 0.3997;
    contacts_other[{age_group_0_to_4, age_group_15_to_34}]   = 0.7957;
    contacts_other[{age_group_0_to_4, age_group_35_to_59}]   = 0.9958;
    contacts_other[{age_group_0_to_4, age_group_60_to_79}]   = 0.3239;
    contacts_other[{age_group_0_to_4, age_group_80_plus}]    = 0.0428;
    contacts_other[{age_group_5_to_14, age_group_0_to_4}]    = 0.0632;
    contacts_other[{age_group_5_to_14, age_group_5_to_14}]   = 0.9121;
    contacts_other[{age_group_5_to_14, age_group_15_to_34}]  = 0.3254;
    contacts_other[{age_group_5_to_14, age_group_35_to_59}]  = 0.4731;
    contacts_other[{age_group_5_to_14, age_group_60_to_79}]  = 0.2355;
    contacts_other[{age_group_5_to_14, age_group_80_plus}]   = 0.0148;
    contacts_other[{age_group_15_to_34, age_group_0_to_4}]   = 0.0336;
    contacts_other[{age_group_15_to_34, age_group_5_to_14}]  = 0.1604;
    contacts_other[{age_group_15_to_34, age_group_15_to_34}] = 1.7529;
    contacts_other[{age_group_15_to_34, age_group_35_to_59}] = 0.8622;
    contacts_other[{age_group_15_to_34, age_group_60_to_79}] = 0.1440;
    contacts_other[{age_group_15_to_34, age_group_80_plus}]  = 0.0077;
    contacts_other[{age_group_35_to_59, age_group_0_to_4}]   = 0.0204;
    contacts_other[{age_group_35_to_59, age_group_5_to_14}]  = 0.1444;
    contacts_other[{age_group_35_to_59, age_group_15_to_34}] = 0.5738;
    contacts_other[{age_group_35_to_59, age_group_35_to_59}] = 1.2127;
    contacts_other[{age_group_35_to_59, age_group_60_to_79}] = 0.3433;
    contacts_other[{age_group_35_to_59, age_group_80_plus}]  = 0.0178;
    contacts_other[{age_group_60_to_79, age_group_0_to_4}]   = 0.0371;
    contacts_other[{age_group_60_to_79, age_group_5_to_14}]  = 0.0393;
    contacts_other[{age_group_60_to_79, age_group_15_to_34}] = 0.4171;
    contacts_other[{age_group_60_to_79, age_group_35_to_59}] = 0.9666;
    contacts_other[{age_group_60_to_79, age_group_60_to_79}] = 0.7495;
    contacts_other[{age_group_60_to_79, age_group_80_plus}]  = 0.0257;
    contacts_other[{age_group_80_plus, age_group_0_to_4}]    = 0.0791;
    contacts_other[{age_group_80_plus, age_group_5_to_14}]   = 0.0800;
    contacts_other[{age_group_80_plus, age_group_15_to_34}]  = 0.3480;
    contacts_other[{age_group_80_plus, age_group_35_to_59}]  = 0.5588;
    contacts_other[{age_group_80_plus, age_group_60_to_79}]  = 0.2769;
    contacts_other[{age_group_80_plus, age_group_80_plus}]   = 0.0180;

    for (auto& loc : world.get_locations()) {
        switch (loc.get_type()) {
        case mio::abm::LocationType::Home:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_home;
            break;
        case mio::abm::LocationType::School:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_work;
            break;
        case mio::abm::LocationType::Work:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_work;
            break;
        default:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            break;
        }
    }
}

/**
 * @brief Add testing strategies to the world.
*/
void add_testing_strategies(mio::abm::World& world, bool school, bool work, bool symptomatic)
{
    // Tests in schools
    auto testing_criteria_school = mio::abm::TestingCriteria();

    auto testing_min_time = mio::abm::days(7);
    auto start_date       = mio::abm::TimePoint(0);
    auto end_date         = mio::abm::TimePoint(0) + mio::abm::days(60);
    auto test_type        = mio::abm::AntigenTest();
    auto probability      = mio::UncertainValue();
    assign_uniform_distribution(probability, 0.5, 0.5);

    auto testing_scheme_school = mio::abm::TestingScheme(testing_criteria_school, testing_min_time, start_date,
                                                         end_date, test_type, probability.draw_sample());
    if (school)
        world.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::School, testing_scheme_school);

    // Tests in work places
    auto testing_criteria_work = mio::abm::TestingCriteria();

    assign_uniform_distribution(probability, 0.5, 0.5);
    auto testing_scheme_work = mio::abm::TestingScheme(testing_criteria_work, testing_min_time, start_date, end_date,
                                                       test_type, probability.draw_sample());
    if (work)
        world.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme_work);

    // Test when symptomatic
    auto testing_criteria_symptomatic = mio::abm::TestingCriteria({}, {mio::abm::InfectionState::InfectedSymptoms});
    auto testing_scheme_symptomatic =
        mio::abm::TestingScheme(testing_criteria_symptomatic, testing_min_time, start_date, end_date, test_type, 0.7);

    if (symptomatic)
        world.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Home, testing_scheme_symptomatic);
}

/**
 * Create a sampled simulation with start time t0.
 * @param t0 The start time of the Simulation.
 */
void create_sampled_world(mio::abm::World& world, const fs::path& input_dir, const mio::abm::TimePoint& t0,
                          int max_num_persons)
{
    //Set global infection parameters (similar to infection parameters in SECIR model) and initialize the world

    set_parameters(world.parameters);
    set_local_parameters(world);

    // Create the world object from statistical data.
    create_world_from_data(world, (input_dir / "mobility/braunschweig_result.csv").generic_string(), t0,
                           max_num_persons);
    world.use_migration_rules(false);

    // Assign an infection state to each person.
    assign_infection_state(world, t0);

    //auto t_lockdown = mio::abm::TimePoint(0) + mio::abm::days(20);

    // During the lockdown, 25% of people work from home and schools are closed for 90% of students.
    // Social events are very rare.
    //mio::abm::set_home_office(t_lockdown, 0.25, world.parameters);
    //mio::abm::set_school_closure(t_lockdown, 0.9, world.parameters);
    //mio::abm::close_social_events(t_lockdown, 0.9, world.parameters);

    // Verschiedene Fälle:
    //1. ohne testing scheme
    //2. testing scheme: Test bei Symptomen (unabh. von Location, 1x am Tag, InfectedSymptomatic, 70% der Bevölkerung) ab Tag 0
    //3. testen in schulen und Arbeitsplätzen (unabh. von Alter, 1x am Tag, unabh. von InfectionState) ab Tag 0
    //4. 2.+3.

    add_testing_strategies(world, false, false, false);
}

template <typename T>
void write_log_to_file_person_and_location_data(const T& history)
{
    auto logg     = history.get_log();
    auto loc_id   = std::get<0>(logg)[0];
    auto agent_id = std::get<1>(logg)[0];
    // Write lo to a text file.
    std::ofstream myfile("locations_lookup.txt");
    myfile << "location_id, location_type, latitude, longitude\n";
    for (uint32_t loc_id_index = 0; loc_id_index < loc_id.size(); ++loc_id_index) {
        auto id            = std::get<0>(loc_id[loc_id_index]);
        auto location_type = (int)std::get<1>(loc_id[loc_id_index]);
        auto id_longitute  = std::get<2>(loc_id[loc_id_index]).longitude;
        auto id_latitude   = std::get<2>(loc_id[loc_id_index]).latitude;
        myfile << id << ", " << location_type << ", " << id_longitute << ", " << id_latitude << "\n";
    }
    myfile.close();

    std::ofstream myfile2("agents_lookup.txt");
    myfile2 << "agent_id, home_id, age\n";
    for (uint32_t agent_id_index = 0; agent_id_index < agent_id.size(); ++agent_id_index) {
        auto id      = std::get<0>(agent_id[agent_id_index]);
        auto home_id = std::get<1>(agent_id[agent_id_index]);
        auto age     = std::get<2>(agent_id[agent_id_index]);
        myfile2 << id << ", " << home_id << ", " << age << "\n";
    }
    myfile2.close();
}

template <typename T>
void write_log_to_file_trip_data(const T& history)
{

    auto movement_data = std::get<0>(history.get_log());
    std::ofstream myfile3("movement_data.txt");
    myfile3 << "agent_id, trip_id, start_location, end_location, start_time, end_time, transport_mode, activity, "
               "infection_state \n";
    int trips_id = 0;
    for (uint32_t movement_data_index = 2; movement_data_index < movement_data.size(); ++movement_data_index) {
        myfile3 << "timestep Nr.: " << movement_data_index - 1 << "\n";
        for (uint32_t trip_index = 0; trip_index < movement_data[movement_data_index].size(); trip_index++) {
            auto agent_id = (int)std::get<0>(movement_data[movement_data_index][trip_index]);

            int start_index = movement_data_index - 1;
            using Type      = std::tuple<uint32_t, uint32_t, mio::abm::TimePoint, mio::abm::TransportMode,
                                    mio::abm::ActivityType, mio::abm::InfectionState>;
            while (!std::binary_search(std::begin(movement_data[start_index]), std::end(movement_data[start_index]),
                                       movement_data[movement_data_index][trip_index],
                                       [](const Type& v1, const Type& v2) {
                                           return std::get<0>(v1) < std::get<0>(v2);
                                       })) {
                start_index--;
            }
            auto start_location_pointer =
                std::lower_bound(std::begin(movement_data[start_index]), std::end(movement_data[start_index]),
                                 movement_data[movement_data_index][trip_index], [](const Type& v1, const Type& v2) {
                                     return std::get<0>(v1) < std::get<0>(v2);
                                 });
            int start_location = (int)std::get<1>(*start_location_pointer);

            auto end_location = (int)std::get<1>(movement_data[movement_data_index][trip_index]);

            auto start_time = (int)std::get<2>(movement_data[movement_data_index][trip_index]).seconds();
            auto end_time   = (int)std::get<2>(movement_data[movement_data_index][trip_index]).seconds();

            auto transport_mode  = (int)std::get<3>(movement_data[movement_data_index][trip_index]);
            auto activity        = (int)std::get<4>(movement_data[movement_data_index][trip_index]);
            auto infection_state = (int)std::get<5>(movement_data[movement_data_index][trip_index]);
            myfile3 << agent_id << ", " << trips_id << ", " << start_location << " , " << end_location << " , "
                    << start_time << " , " << end_time << " , " << transport_mode << " , " << activity << " , "
                    << infection_state << "\n";
            trips_id++;
        }
    }
    myfile3.close();
}

void write_txt_file_for_graphical_compartment_output(std::vector<std::vector<mio::TimeSeries<ScalarType>>> input_file)
{
    // mio::unused(input_file);
    // In the input file is a h5 file there are multiple runs with each having the amount of people in each compartment for each timestep.
    // The output folder should have the following format:
    // for each run there is a file with the name "run_1.txt" and so on.
    // in each file the the rows represent the timesteps and the columns represent the compartments.
    // The first row is the header with the compartment names.
    // Time = Time in days, S = Susceptible, E = Exposed, I_NS = InfectedNoSymptoms, I_Sy = InfectedSymptoms, I_Sev = InfectedSevere,
    // I_Crit = InfectedCritical, R = Recovered, D = Dead

    // Output folder name:
    std::string folderName = "folder_run_bs";
    // Create folder
    fs::create_directory(folderName);
    // Loop over all runs
    for (int run = 0; run < (int)input_file.size(); run++) {
        // Create file name
        std::string fileName = folderName + "/run_" + std::to_string(run) + ".txt";
        // Create file
        std::ofstream myfile;
        myfile.open(fileName);
        // Write header
        input_file.at(run).at(0).print_table({"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4, myfile);
    }
}

template <typename T>
void write_log_to_file_infection_per_location_type(const T& history, const fs::path& result_dir)
{
    std::ofstream myfile4((result_dir / "infection_per_location_type.txt").string());
    const std::vector<std::string>& labels = {
        "Home",     "Work", "School", "SocialEvent",     "BasicsShop",
        "Hospital", "ICU",  "Car",    "PublicTransport", "TransportWithoutContact",
        "Cemetery"};
    std::get<0>(history.get_log()).print_table(labels, 12, 4, myfile4);
}

template <typename T>
void write_log_to_file_infection_per_age_group(const T& history, const fs::path& result_dir)
{
    std::ofstream myfile5((result_dir / "infection_per_age_group.txt").string());
    const std::vector<std::string>& labels = {"0_to_4", "5_to_14", "15_to_34", "35_to_59", "60_to_79", "80_plus"};
    std::get<0>(history.get_log()).print_table(labels, 7, 4, myfile5);
}

struct LogInfectionStatePerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    /** 
     * @brief Log the TimeSeries of the number of Person%s in an #InfectionState.
     * @param[in] sim The simulation of the abm.
     * @return A pair of the TimePoint and the TimeSeries of the number of Person%s in an #InfectionState.
     */
    static Type log(const mio::abm::Simulation& sim)
    {

        Eigen::VectorXd sum = Eigen::VectorXd::Zero(
            Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups()));
        auto curr_time     = sim.get_time();
        const auto persons = sim.get_world().get_persons();

        // PRAGMA_OMP(parallel for)
        for (auto i = size_t(0); i < persons.size(); ++i) {
            auto& p    = persons[i];
            auto index = (((size_t)(mio::abm::InfectionState::Count)) * ((uint32_t)p.get_age().get())) +
                         ((uint32_t)p.get_infection_state(curr_time));
            // PRAGMA_OMP(atomic)
            sum[index] += 1;
        }
        return std::make_pair(curr_time, sum);
    }
};

template <typename T>
T gather_results(int rank, int num_procs, int num_runs, T ensemble_vec)
{
    auto gathered_ensemble_vec = T{};

    if (rank == 0) {
        gathered_ensemble_vec.reserve(num_runs);
        std::copy(ensemble_vec.begin(), ensemble_vec.end(), std::back_inserter(gathered_ensemble_vec));
        for (int src_rank = 1; src_rank < num_procs; ++src_rank) {
            int bytes_size;
            MPI_Recv(&bytes_size, 1, MPI_INT, src_rank, 0, mio::mpi::get_world(), MPI_STATUS_IGNORE);
            mio::ByteStream bytes(bytes_size);
            MPI_Recv(bytes.data(), bytes.data_size(), MPI_BYTE, src_rank, 0, mio::mpi::get_world(), MPI_STATUS_IGNORE);

            auto src_ensemble_results = mio::deserialize_binary(bytes, mio::Tag<decltype(ensemble_vec)>{});
            if (!src_ensemble_results) {
                mio::log_error("Error receiving ensemble results from rank {}.", src_rank);
            }
            std::copy(src_ensemble_results.value().begin(), src_ensemble_results.value().end(),
                      std::back_inserter(gathered_ensemble_vec));
        }
    }
    else {
        auto bytes      = mio::serialize_binary(ensemble_vec);
        auto bytes_size = int(bytes.data_size());
        MPI_Send(&bytes_size, 1, MPI_INT, 0, 0, mio::mpi::get_world());
        MPI_Send(bytes.data(), bytes.data_size(), MPI_BYTE, 0, 0, mio::mpi::get_world());
    }
    return gathered_ensemble_vec;
}

std::vector<size_t> distribute_runs(size_t num_runs, int num_procs)
{
    //evenly distribute runs
    //lower processes do one more run if runs are not evenly distributable
    auto num_runs_local = num_runs / num_procs; //integer division!
    auto remainder      = num_runs % num_procs;

    std::vector<size_t> run_distribution(num_procs);
    std::fill(run_distribution.begin(), run_distribution.begin() + remainder, num_runs_local + 1);
    std::fill(run_distribution.begin() + remainder, run_distribution.end(), num_runs_local);

    return run_distribution;
}

mio::IOResult<void> run(const fs::path& input_dir, const fs::path& result_dir, size_t num_runs,
                        bool save_single_runs = true)
{
    int num_procs, rank;
#ifdef MEMILIO_ENABLE_MPI
    MPI_Comm_size(mio::mpi::get_world(), &num_procs);
    MPI_Comm_rank(mio::mpi::get_world(), &rank);
#else
    num_procs = 1;
    rank      = 0;
#endif

    auto run_distribution = distribute_runs(num_runs, num_procs);
    auto start_run_idx = std::accumulate(run_distribution.begin(), run_distribution.begin() + size_t(rank), size_t(0));
    auto end_run_idx   = start_run_idx + run_distribution[size_t(rank)];

    mio::Date start_date{2021, 3, 1};
    auto t0              = mio::abm::TimePoint(0); // Start time per simulation
    auto tmax            = mio::abm::TimePoint(0) + mio::abm::days(120); // End time per simulation
    auto max_num_persons = 50000;
    auto ensemble_infection_per_loc_type =
        std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of infection per location type results
    // ensemble_infection_per_loc_type.reserve(size_t(num_runs));
    // auto ensemble_infection_per_age_group =
    //     std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of infection per age group results
    // ensemble_infection_per_age_group.reserve(size_t(num_runs));

    auto ensemble_infection_state_per_age_group =
        std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of infection state per age group results
    ensemble_infection_state_per_age_group.reserve(size_t(num_runs));
    auto ensemble_params = std::vector<std::vector<mio::abm::World>>{}; // Vector of all worlds
    ensemble_params.reserve(size_t(num_runs));
    // auto run_idx            = size_t(1); // The run index
    // auto save_result_result = mio::IOResult<void>(mio::success()); // Variable informing over successful IO operations

    int tid = -1;
#pragma omp parallel private(tid) // Start of parallel region: forks threads
    {
        tid = omp_get_thread_num(); // default is number of CPUs on machine
        printf("Hello World from thread = %d and rank = %d\n", tid, rank);
        if (tid == 0) {
            printf("Number of threads = %d\n", omp_get_num_threads());
        }
    } // ** end of the the parallel: joins threads

    // Determine inital infection state distribution
    //Time this
    auto start0 = std::chrono::high_resolution_clock::now();
    // determine_initial_infection_states_world(input_dir, start_date);
    auto stop0     = std::chrono::high_resolution_clock::now();
    auto duration0 = std::chrono::duration<double>(stop0 - start0);
    std::cout << "Time taken by determine_initial_infection_states_world: " << duration0.count() << " seconds"
              << std::endl;

    // Create one world for all simulations that will be copied
    // auto world = mio::abm::World(num_age_groups);
    // create_sampled_world(world, input_dir, t0, max_num_persons);

    // Loop over a number of runs
    for (size_t run_idx = start_run_idx; run_idx < end_run_idx; run_idx++) {

        // Start the clock before create_sampled_world
        auto start1 = std::chrono::high_resolution_clock::now();
        // Create the sampled simulation with start time t0.
        auto world = mio::abm::World(num_age_groups);
        create_sampled_world(world, input_dir, t0, max_num_persons);
        // Stop the clock after create_sampled_world and calculate the duration
        auto stop1     = std::chrono::high_resolution_clock::now();
        auto duration1 = std::chrono::duration<double>(stop1 - start1);
        std::cout << "Time taken by create_sampled_world: " << duration1.count() << " seconds" << std::endl;

        // auto world_copy = world; // COPY CONSTRUCTOR DOESN'T WORK. LOCATIONS AREN'T ASSIGNED!
        auto sim = mio::abm::Simulation(t0, std::move(world));

        //output object
        // mio::History<mio::DataWriterToMemory, mio::abm::LogLocationInformation, mio::abm::LogPersonInformation,
        //              mio::abm::LogDataForMovement>
        //     historyPersonInf;
        // mio::History<mio::abm::DataWriterToMemoryDelta, mio::abm::LogDataForMovement> historyPersonInfDelta;
        mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionPerLocationType> historyInfectionPerLocationType{
            Eigen::Index(mio::abm::LocationType::Count)};
        // mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionPerAgeGroup> historyInfectionPerAgeGroup{
        //     Eigen::Index(sim.get_world().parameters.get_num_groups())};
        mio::History<mio::abm::TimeSeriesWriter, LogInfectionStatePerAgeGroup> historyInfectionStatePerAgeGroup{
            Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups())};

        // Collect the id of location in world.
        // std::vector<int> loc_ids;
        // for (auto& location : sim.get_world().get_locations()) {
        //     loc_ids.push_back(location.get_index());
        // }

        // Start the clock before sim.advance
        auto start2 = std::chrono::high_resolution_clock::now();
        // Advance the world to tmax
        // sim.advance(tmax, historyPersonInf, historyInfectionPerLocationType, historyInfectionPerAgeGroup,
        //             historyPersonInfDelta, historyInfectionStatePerAgeGroup);
        sim.advance(tmax, historyInfectionStatePerAgeGroup, historyInfectionPerLocationType);
        // Stop the clock after sim.advance and calculate the duration
        auto stop2     = std::chrono::high_resolution_clock::now();
        auto duration2 = std::chrono::duration<double>(stop2 - start2);
        std::cout << "Time taken by sim.advance: " << duration2.count() << " seconds" << std::endl;

        // TODO: update result of the simulation to be a vector of location result.
        auto temp_sim_infection_per_loc_tpye =
            std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionPerLocationType.get_log())};
        // auto temp_sim_infection_per_age_group =
        //     std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionPerAgeGroup.get_log())};
        auto temp_sim_infection_state_per_age_group =
            std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionStatePerAgeGroup.get_log())};

        // Push result of the simulation back to the result vector
        ensemble_infection_per_loc_type.emplace_back(temp_sim_infection_per_loc_tpye);
        // ensemble_infection_per_age_group.push_back(temp_sim_infection_per_age_group);
        // ensemble_params.push_back(std::vector<mio::abm::World>{sim.get_world()});
        ensemble_infection_state_per_age_group.emplace_back(temp_sim_infection_state_per_age_group);
        // Option to save the current run result to file
        /* if (save_result_result && save_single_runs) {
            auto result_dir_run = result_dir / ("abm_result_run_" + std::to_string(run_idx) + ".h5");
            BOOST_OUTCOME_TRY(save_result(ensemble_results.back(), {0}, 1, result_dir_run.string()));
            result_dir_run = result_dir / ("abm_result_infections_per_loc_run_" + std::to_string(run_idx) + ".h5");
            BOOST_OUTCOME_TRY(
                save_result(ensemble_results_infections_per_loc.back(), loc_type_ids, 1, result_dir_run.string()));
            result_dir_run = result_dir / ("abm_result_infections_per_age_run_" + std::to_string(run_idx) + ".h5");
            BOOST_OUTCOME_TRY(
                save_result(ensemble_results_infections_per_age.back(), age_group_ids, 1, result_dir_run.string()));
        } */
        // write_log_to_file_person_and_location_data(historyPersonInf);
        // write_log_to_file_trip_data(historyPersonInfDelta);
        // write_log_to_file_infection_per_age_group(historyInfectionPerAgeGroup, result_dir);
        write_log_to_file_infection_per_location_type(historyInfectionPerLocationType, result_dir);

        std::cout << "Run " << run_idx + 1 << " of " << num_runs << " finished." << std::endl;

        //HACK since // gather_results(rank, num_procs, num_runs, ensemble_params);
        //for now this doesnt work, but we can still save the results of the last world since the
        //parameters are the same for each run
        if (rank == 0 && run_idx == end_run_idx - 1) {
            for (size_t i = 0; i < num_runs; i++) {
                ensemble_params.emplace_back(std::vector<mio::abm::World>{sim.get_world()});
            }
        }
    }

    printf("Saving results ... ");

#ifdef MEMILIO_ENABLE_MPI
    //gather results
    auto final_ensemble_infection_state_per_age_group =
        gather_results(rank, num_procs, num_runs, ensemble_infection_state_per_age_group);
    auto final_ensemble_infection_per_loc_type =
        gather_results(rank, num_procs, num_runs, ensemble_infection_per_loc_type);
    if (rank == 0) {
        BOOST_OUTCOME_TRY(save_results(final_ensemble_infection_state_per_age_group, ensemble_params, {0},
                                       result_dir / "infection_state_per_age_group/", save_single_runs));
        BOOST_OUTCOME_TRY(save_results(final_ensemble_infection_per_loc_type, ensemble_params, {0},
                                       result_dir / "infection_per_location_type/", save_single_runs));
    }
#else

    BOOST_OUTCOME_TRY(save_results(ensemble_infection_state_per_age_group, ensemble_params, {0},
                                   result_dir / "infection_state_per_age_group/", save_single_runs));
    BOOST_OUTCOME_TRY(save_results(ensemble_infection_per_loc_type, ensemble_params, {0},
                                   result_dir / "infection_per_location_type/", save_single_runs));

#endif

    // BOOST_OUTCOME_TRY(save_results(ensemble_infection_per_loc_type, ensemble_params, {0},
    //                                result_dir / "infection_per_location_type/", save_single_runs));
    // BOOST_OUTCOME_TRY(save_results(ensemble_infection_per_age_group, ensemble_params, {0},
    //                                result_dir / "infection_per_age_group/", save_single_runs));

    printf("done.\n");
    //write_txt_file_for_graphical_compartment_output(ensemble_infection_state_per_age_group);
    // BOOST_OUTCOME_TRY(save_result_result);
    return mio::success();
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::err);
#ifdef MEMILIO_ENABLE_MPI
    mio::mpi::init();
#endif

    std::string input_dir  = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data";
    std::string result_dir = input_dir + "/results";
    size_t num_runs;
    bool save_single_runs = true;

    if (argc == 2) {
        num_runs = atoi(argv[1]);
        printf("Number of run is %s.\n", argv[1]);
        printf("Saving results to the current directory.\n");
    }

    else if (argc == 3) {
        num_runs   = atoi(argv[1]);
        result_dir = argv[2];
        printf("Number of run is %s.\n", argv[1]);
        printf("Saving results to \"%s\".\n", result_dir.c_str());
    }
    else {
        printf("Usage:\n");
        printf("abm_example <num_runs>\n");
        printf("\tRun the simulation for <num_runs> time(s).\n");
        printf("\tStore the results in the current directory.\n");
        printf("abm_braunschweig <num_runs> <result_dir>\n");
        printf("\tRun the simulation for <num_runs> time(s).\n");
        printf("\tStore the results in <result_dir>.\n");

        num_runs = 1;
        printf("Running with number of runs = %d.\n", (int)num_runs);
    }

    // mio::thread_local_rng().seed({...}); //set seeds, e.g., for debugging
    // printf("Seeds: ");
    // for (auto s : mio::thread_local_rng().get_seeds()) {
    //     printf("%u, ", s);
    // }
    // printf("\n");

    auto result = run(input_dir, result_dir, num_runs, save_single_runs);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        mio::mpi::finalize();
        return -1;
    }

    mio::mpi::finalize();
    return 0;
}
