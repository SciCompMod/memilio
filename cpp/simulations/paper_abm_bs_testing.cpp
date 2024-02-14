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
#include "abm/abm.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/uncertain_value.h"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "abm/vaccine.h"
#include "abm/common_abm_loggers.h"

namespace fs = boost::filesystem;

// Assign the name to general age group.
size_t num_age_groups         = 6;
const auto age_group_0_to_4   = mio::AgeGroup(0);
const auto age_group_5_to_14  = mio::AgeGroup(1);
const auto age_group_15_to_34 = mio::AgeGroup(2);
const auto age_group_35_to_59 = mio::AgeGroup(3);
const auto age_group_60_to_79 = mio::AgeGroup(4);
const auto age_group_80_plus  = mio::AgeGroup(5);

/**
 * Determine the infection state of a person at the beginning of the simulation.
 * The infection states are chosen randomly. They are distributed according to the probabilites set in the example.
 * @return random infection state
 */
mio::abm::InfectionState determine_infection_state(mio::abm::Person::RandomNumberGenerator& rng, ScalarType exposed,
                                                   ScalarType infected_no_symptoms, ScalarType infected_symptoms,
                                                   ScalarType recovered)
{
    ScalarType susceptible          = 1 - exposed - infected_no_symptoms - infected_symptoms - recovered;
    std::vector<ScalarType> weights = {
        susceptible,           exposed,  infected_no_symptoms, infected_symptoms / 3, infected_symptoms / 3,
        infected_symptoms / 3, recovered};
    if (weights.size() != (size_t)mio::abm::InfectionState::Count - 1) {
        mio::log_error("Initialization in ABM wrong, please correct vector length.");
    }
    auto state = mio::DiscreteDistribution<size_t>::get_instance()(rng, weights);
    return (mio::abm::InfectionState)state;
}

/**
 * Assign an infection state to each person.
 */
void assign_infection_state(mio::abm::World& world, mio::abm::TimePoint t, double exposed_prob,
                            double infected_no_symptoms_prob, double infected_symptoms_prob, double recovered_prob)
{
    auto persons = world.get_persons();
    for (auto& person : persons) {
        auto rng             = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);
        auto infection_state = determine_infection_state(rng, exposed_prob, infected_no_symptoms_prob,
                                                         infected_symptoms_prob, recovered_prob);
        if (infection_state != mio::abm::InfectionState::Susceptible)
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         world.parameters, t, infection_state));
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
    std::map<uint32_t, uint32_t> person_ids            = {};
    std::map<uint32_t, std::pair<uint32_t, int>> locations_before;
    std::map<uint32_t, std::pair<uint32_t, int>> locations_after;

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
    while (std::getline(fin, line)) {
        row.clear();

        // read columns in this row
        split_line(line, &row);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        uint32_t person_id = row[index["puid"]];
        if (person_ids.find(person_id) == person_ids.end())
            break;

        uint32_t home_id                                 = row[index["huid"]];
        uint32_t target_location_id                      = std::abs(row[index["loc_id_end"]]);
        uint32_t activity_end                            = row[index["activity_end"]];
        mio::abm::GeographicalLocation location_long_lat = {(double)row[index["lon_end"]] / 1e+5,
                                                            (double)row[index["lat_end"]] / 1e+5};
        mio::abm::LocationId home;
        auto it_home = locations.find(home_id);
        if (it_home == locations.end()) {
            home = world.add_location(mio::abm::LocationType::Home, 1);
            locations.insert({home_id, home});
            mio::abm::GeographicalLocation location_long_lat_home = {(double)row[index["lon_start"]] / 1e+5,
                                                                     (double)row[index["lat_start"]] / 1e+5};
            world.get_individualized_location(home).set_geographical_location(location_long_lat_home);
        }
        else {
            home = it_home->second;
        }

        mio::abm::LocationId location;
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

    // Add the persons and trips
    while (std::getline(fin, line)) {
        row.clear();

        // read columns in this row
        split_line(line, &row);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        uint32_t person_id = row[index["puid"]];
        if (person_ids.find(person_id) == person_ids.end())
            break;

        uint32_t age                = row[index["age"]];
        uint32_t home_id            = row[index["huid"]];
        uint32_t target_location_id = std::abs(row[index["loc_id_end"]]);
        uint32_t start_location_id  = std::abs(row[index["loc_id_start"]]);
        uint32_t trip_start         = row[index["start_time"]];
        uint32_t transport_mode     = row[index["travel_mode"]];
        uint32_t acticity_end       = row[index["activity_end"]];

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
            auto home              = locations.find(home_id)->second;
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
}

std::pair<double, double> get_my_and_sigma(std::pair<double, double> mean_and_std)
{
    auto mean    = mean_and_std.first;
    auto stddev  = mean_and_std.second;
    double my    = log(mean * mean / sqrt(mean * mean + stddev * stddev));
    double sigma = sqrt(log(1 + stddev * stddev / (mean * mean)));
    return {my, sigma};
}

void set_parameters(mio::abm::Parameters params)
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

    // Set protection level from low viral load. Information based on: https://doi.org/10.1093/cid/ciaa886
    // params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ExposureType::NaturalInfection, age_group_60_to_79,
    //                                                   mio::abm::VirusVariant::Wildtype}] =
    //     [](ScalarType days) -> ScalarType {
    //     return mio::linear_interpolation_of_data_set<ScalarType, ScalarType>({{0, 0.967},
    //                                                                           {30, 0.975},
    //                                                                           {60, 0.977},
    //                                                                           {90, 0.974},
    //                                                                           {120, 0.963},
    //                                                                           {150, 0.947},
    //                                                                           {180, 0.93},
    //                                                                           {210, 0.929},
    //                                                                           {240, 0.923},
    //                                                                           {270, 0.908},
    //                                                                           {300, 0.893},
    //                                                                           {330, 0.887},
    //                                                                           {360, 0.887},
    //                                                                           {360, 0.5}},
    //                                                                          days);
    // };

    //Set other parameters
    params.get<mio::abm::MaskProtection>()           = 0.5;
    params.get<mio::abm::AerosolTransmissionRates>() = 1.0;
}

/**
 * Create a sampled simulation with start time t0.
 * @param t0 The start time of the Simulation.
 */
mio::abm::Simulation create_sampled_simulation(const std::string& input_file, const mio::abm::TimePoint& t0,
                                               int max_num_persons)
{
    // Assumed percentage of infection state at the beginning of the simulation.
    ScalarType exposed_prob = 0.05, infected_no_symptoms_prob = 0.001, infected_symptoms_prob = 0.001,
               recovered_prob = 0.0;

    //Set global infection parameters (similar to infection parameters in SECIR model) and initialize the world
    auto world = mio::abm::World(num_age_groups);

    set_parameters(world.parameters);

    // Create the world object from statistical data.
    create_world_from_data(world, input_file, t0, max_num_persons);
    world.use_migration_rules(false);

    // Assign an infection state to each person.
    assign_infection_state(world, t0, exposed_prob, infected_no_symptoms_prob, infected_symptoms_prob, recovered_prob);

    auto t_lockdown = mio::abm::TimePoint(0) + mio::abm::days(1);

    // During the lockdown, 25% of people work from home and schools are closed for 90% of students.
    // Social events are very rare.
    mio::abm::set_home_office(t_lockdown, 0.25, world.parameters);
    mio::abm::set_school_closure(t_lockdown, 0.9, world.parameters);
    mio::abm::close_social_events(t_lockdown, 0.9, world.parameters);

    auto sim = mio::abm::Simulation(t0, std::move(world));
    return sim;
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

mio::IOResult<void> run(const std::string& input_file, const fs::path& result_dir, size_t num_runs,
                        bool save_single_runs = true)
{
    auto t0               = mio::abm::TimePoint(0); // Start time per simulation
    auto tmax             = mio::abm::TimePoint(0) + mio::abm::days(60); // End time per simulation
    auto ensemble_results = std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of collected results
    ensemble_results.reserve(size_t(num_runs));
    auto run_idx            = size_t(1); // The run index
    auto save_result_result = mio::IOResult<void>(mio::success()); // Variable informing over successful IO operations
    auto max_num_persons    = 10000;

    // Loop over a number of runs
    while (run_idx <= num_runs) {

        // Create the sampled simulation with start time t0.
        auto sim = create_sampled_simulation(input_file, t0, max_num_persons);
        //output object
        mio::History<mio::DataWriterToMemory, mio::abm::LogLocationInformation, mio::abm::LogPersonInformation,
                     mio::abm::LogDataForMovement>
            historyPersonInf;
        mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
            Eigen::Index(mio::abm::InfectionState::Count)};
        mio::History<mio::abm::DataWriterToMemoryDelta, mio::abm::LogDataForMovement> historyPersonInfDelta;
        // Collect the id of location in world.
        std::vector<int> loc_ids;
        for (auto& location : sim.get_world().get_locations()) {
            loc_ids.push_back(location.get_index());
        }
        // Advance the world to tmax
        sim.advance(tmax, historyPersonInf, historyTimeSeries, historyPersonInfDelta);
        // TODO: update result of the simulation to be a vector of location result.
        auto temp_sim_result = std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyTimeSeries.get_log())};
        // Push result of the simulation back to the result vector
        ensemble_results.push_back(temp_sim_result);
        // Option to save the current run result to file
        if (save_result_result && save_single_runs) {
            auto result_dir_run = result_dir / ("abm_result_run_" + std::to_string(run_idx) + ".h5");
            BOOST_OUTCOME_TRY(save_result(ensemble_results.back(), loc_ids, 1, result_dir_run.string()));
        }
        write_log_to_file_person_and_location_data(historyPersonInf);
        write_log_to_file_trip_data(historyPersonInfDelta);

        std::cout << "Run " << run_idx << " of " << num_runs << " finished." << std::endl;
        ++run_idx;
    }
    write_txt_file_for_graphical_compartment_output(ensemble_results);
    BOOST_OUTCOME_TRY(save_result_result);
    return mio::success();
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::warn);

    std::string result_dir = ".";
    std::string input_file =
        "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/cpp/simulations/bs_und_umgebung.csv";
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
        printf("Running with number of runs = 1.\n");
        num_runs = 1;
    }

    // mio::thread_local_rng().seed({...}); //set seeds, e.g., for debugging
    //printf("Seeds: ");
    //for (auto s : mio::thread_local_rng().get_seeds()) {
    //    printf("%u, ", s);
    //}
    //printf("\n");

    auto result = run(input_file, result_dir, num_runs, save_single_runs);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }
    return 0;
}
