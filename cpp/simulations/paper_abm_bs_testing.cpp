/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Sascha Korf, David Kerkmann, Khoa Nguyen, Carlotta Gerstein
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
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"

#include <chrono>
#define TIME_TYPE std::chrono::steady_clock::time_point
#define TIME_NOW std::chrono::steady_clock::now()
#define PRINTABLE_TIME(_time) (std::chrono::duration_cast<std::chrono::duration<double>>(_time)).count()

#define restart_timer(timer, description)                                                                              \
    {                                                                                                                  \
        TIME_TYPE new_time = TIME_NOW;                                                                                 \
        std::cout << "\r" << description << " :: " << PRINTABLE_TIME(new_time - timer) << std::endl << std::flush;     \
        timer = new_time;                                                                                              \
    }
#define DEBUG(cout_args) std::cerr << cout_args << std::endl << std::flush;

namespace fs = boost::filesystem;

TIME_TYPE timer;

// Assign the name to general age group.
size_t num_age_groupss        = 6;
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

/**
 * Determine initial distribution of infection states.
*/
mio::CustomIndexArray<double, mio::AgeGroup, mio::osecir::InfectionState>
determine_initial_infection_states_world(const fs::path& input_dir, const mio::Date date, double sclaling_infected)
{
    // estimate intial population by ODE compartiments
    auto initial_graph           = get_graph(date, 1, input_dir, sclaling_infected);
    const size_t braunschweig_id = 16; // Braunschweig has ID 16
    auto braunschweig_node       = initial_graph.value()[braunschweig_id];
    mio::CustomIndexArray<double, mio::AgeGroup, mio::osecir::InfectionState> initial_infection_distribution{
        {mio::AgeGroup(num_age_groupss), mio::osecir::InfectionState::Count}, 0.5};

    initial_infection_distribution.array() = braunschweig_node.populations.array().cast<double>();
    return initial_infection_distribution;
}

/**
 * Assign an infection state to each person according to real world data read in through the ODE secir model.
 * Infections are set with the rounded values in the rows in initial_infection_distribution.
 * Only works if enough persons in the all age groups exist.
 * The number of agents in the model should fit to the sum of the rows in initial_infection_distribution,
 * otherwise many agents will be susceptible.
 */
void assign_infection_state(
    mio::abm::World& world, mio::abm::TimePoint t,
    mio::CustomIndexArray<double, mio::AgeGroup, mio::osecir::InfectionState> initial_infection_distribution)
{
    // save all persons with age groups
    std::vector<std::vector<uint32_t>> persons_by_age(num_age_groupss);

    for (auto& person : world.get_persons()) {
        if (person.get_should_be_logged()) {
            persons_by_age[person.get_age().get()].push_back(person.get_person_id());
        }
    }

    for (size_t age = 0; age < num_age_groupss; ++age) {
        auto age_slice = initial_infection_distribution.slice(mio::AgeGroup(age)).as_array().array();
        auto age_grp   = mio::AgeGroup(age).get();
        // Check that the world has enough persons in each age group to initialize infections.
        // (All persons minus the susceptibles.)
        // For lower population sizes use the same method with _prob at the end.
        assert(age_slice.sum() - age_slice[0] <= persons_by_age[age_grp].size() &&
               "Not enough persons to initialize with given amount of infections.");

        // Iterate over all InfectionStates except the susceptibles.
        for (auto i = 1; i < age_slice.size(); ++i) {
            for (auto j = 0; j < std::floor(age_slice[i]); ++j) {
                // select random person and assign Infection
                uint32_t id_rnd          = persons_by_age[age_grp][mio::UniformIntDistribution<size_t>::get_instance()(
                    world.get_rng(), 0U, persons_by_age[age_grp].size() - 1)];
                mio::abm::Person& person = world.get_person(id_rnd);
                auto rng                 = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);
                person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Alpha, person.get_age(),
                                                             world.parameters, t,
                                                             infection_state_map.at(mio::osecir::InfectionState(i))));
                persons_by_age[age_grp].erase(
                    std::remove(persons_by_age[age_grp].begin(), persons_by_age[age_grp].end(), id_rnd),
                    persons_by_age[age_grp].end());
            }
        }
    }
}

size_t determine_age_group_from_rki(mio::AgeGroup age)
{
    if (age == mio::AgeGroup(0)) {
        return 0;
    }
    else if (age == mio::AgeGroup(1)) {
        return 1;
    }
    else if (age == mio::AgeGroup(2)) {
        return 2;
    }
    else if (age == mio::AgeGroup(3)) {
        return 3;
    }
    else if (age == mio::AgeGroup(4)) {
        return 4;
    }
    else if (age == mio::AgeGroup(5)) {
        return 5;
    }
    else {
        return 2;
    }
}

std::map<mio::Date, std::vector<std::pair<uint32_t, uint32_t>>> prepare_vaccination_state(mio::Date simulation_end,
                                                                                          const std::string& filename)
{
    // for saving previous day of vaccination
    std::vector<std::pair<uint32_t, uint32_t>> vacc_vector_prev(num_age_groupss);
    //inizialize the vector with 0
    for (size_t i = 0; i < num_age_groupss; ++i) {
        vacc_vector_prev[i] = std::make_pair(0, 0);
    }
    std::map<mio::Date, std::vector<std::pair<uint32_t, uint32_t>>> vacc_map{};

    //Read in file with vaccination data
    auto vacc_data = mio::read_vaccination_data(filename).value();
    for (auto& vacc_entry : vacc_data) {
        // we need ot filter out braunschweig with zip code 3101
        if (vacc_entry.county_id.value() == mio::regions::CountyId(3101)) {
            //we need the vaccination from the beginning till the end of the simulaiton (2021-05-30)
            if (vacc_entry.date <= simulation_end && vacc_entry.date >= mio::Date(2020, 12, 01)) {
                // if the date isn't in the map we need to add a vector of size num_age_groupss
                if (vacc_map.find(vacc_entry.date) == vacc_map.end()) {
                    vacc_map[vacc_entry.date] = std::vector<std::pair<uint32_t, uint32_t>>(num_age_groupss);
                }
                // we need to add the number of persons to the vector of the date, but these are cumulative so we need to substract the day before
                vacc_map[vacc_entry.date][determine_age_group_from_rki(vacc_entry.age_group)].first =
                    (int)vacc_entry.num_vaccinations_partially -
                    (int)vacc_vector_prev[determine_age_group_from_rki(vacc_entry.age_group)].first;
                vacc_map[vacc_entry.date][determine_age_group_from_rki(vacc_entry.age_group)].second =
                    (int)vacc_entry.num_vaccinations_completed -
                    (int)vacc_vector_prev[determine_age_group_from_rki(vacc_entry.age_group)].second;

                //update the vector for the next iteration
                vacc_vector_prev[determine_age_group_from_rki(vacc_entry.age_group)].first =
                    (int)vacc_entry.num_vaccinations_partially;
                vacc_vector_prev[determine_age_group_from_rki(vacc_entry.age_group)].second =
                    (int)vacc_entry.num_vaccinations_completed;
            }
        }
    }
    return vacc_map;
}

/**
 * @brief assign an vaccination state to each person according to real world data read in through the ODE secir model.
 * 
 * @param input 
 * @return int 
 */
void assign_vaccination_state(mio::abm::World& world, mio::Date simulation_beginning,
                              std::map<mio::Date, std::vector<std::pair<uint32_t, uint32_t>>> vacc_map)
{
    // we check if we even have enough people to vaccinate in each respective age group
    std::vector<size_t> num_persons_by_age(num_age_groupss);
    for (auto& person : world.get_persons()) {
        if (person.get_should_be_logged()) {
            num_persons_by_age[determine_age_group_from_rki(person.get_age())]++;
        }
    }
    //sum over all dates in the vacc_map to check if we have enough persons to vaccinate
    std::vector<size_t> num_persons_by_age_vaccinate(num_age_groupss);
    for (auto& vacc_entry : vacc_map) {
        for (size_t age = 0; age < vacc_entry.second.size(); ++age) {
            num_persons_by_age_vaccinate[age] += vacc_entry.second[age].first;
        }
    }
    //check
    for (size_t age = 0; age < num_persons_by_age.size(); ++age) {
        if (num_persons_by_age[age] < num_persons_by_age_vaccinate[age]) {
            // mio::log_error(
            //     "Not enough persons to vaccinate in age group we dont vaccinate if an age group is fully vaccinated! ",
            //     age);
        }
    }

    // save all persons with age groups
    std::vector<std::vector<uint32_t>> persons_by_age(num_age_groupss);

    for (auto& person : world.get_persons()) {
        persons_by_age[person.get_age().get()].push_back(person.get_person_id());
    }

    // vaccinate random persons according to the data and we also beforehand need to keep a list of persons which are already vaccinated and their age to vaccinate them with the second dose

    // first we need a vector with a list of ids of already vaccinated persons for each age group
    std::vector<std::vector<uint32_t>> vaccinated_persons(num_age_groupss);
    for (auto& vacc_entry : vacc_map) {
        for (size_t age = 0; age < vacc_entry.second.size(); ++age) {
            for (uint32_t i = 0; i < vacc_entry.second[age].first; ++i) {
                if (persons_by_age[age].size() == 0) {
                    // mio::log_error("Not enough to vacc people 1st time");
                }
                else {
                    // select random person and assign Vaccination
                    uint32_t id_rnd          = persons_by_age[age][mio::UniformIntDistribution<size_t>::get_instance()(
                        world.get_rng(), 0U, persons_by_age[age].size() - 1)];
                    mio::abm::Person& person = world.get_person(id_rnd);
                    auto timePoint           = mio::abm::TimePoint(
                        mio::get_offset_in_days(vacc_entry.first, simulation_beginning) * 24 * 60 * 60);
                    person.add_new_vaccination(
                        mio::abm::Vaccination(mio::abm::ExposureType::GenericVaccine, timePoint));
                    persons_by_age[age].erase(
                        std::remove(persons_by_age[age].begin(), persons_by_age[age].end(), id_rnd),
                        persons_by_age[age].end());
                    vaccinated_persons[age].push_back(id_rnd);
                }
            }
            // for (uint32_t i = 0; i < vacc_entry.second[age].second; ++i) {
            //     if (vaccinated_persons[age].size() == 0) {
            //         // mio::log_error("Not enough vaccinated people to vacc 2nd time! ");
            //     }
            //     else {
            //         // select random already vaccinated person and assign Vaccination
            //         uint32_t id_rnd = vaccinated_persons[age][mio::UniformIntDistribution<size_t>::get_instance()(
            //             world.get_rng(), 0U, vaccinated_persons[age].size() - 1)];
            //         mio::abm::Person& person = world.get_person(id_rnd);
            //         auto timePoint           = mio::abm::TimePoint(
            //             mio::get_offset_in_days(vacc_entry.first, simulation_beginning) * 24 * 60 * 60);
            //         person.add_new_vaccination(
            //             mio::abm::Vaccination(mio::abm::ExposureType::GenericVaccine, timePoint));
            //         vaccinated_persons[age].erase(
            //             std::remove(vaccinated_persons[age].begin(), vaccinated_persons[age].end(), id_rnd),
            //             vaccinated_persons[age].end());
            //     }
            // }
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
    boost::split(strings, string, boost::is_any_of(","));
    std::transform(strings.begin(), strings.end(), std::back_inserter(*row), [&](std::string s) {
        if (s.find(":") != std::string::npos) {
            return stringToMinutes(s);
        }
        else if (s.find(".") != std::string::npos) {
            return longLatToInt(s);
        }
        else {
            return std::stoi(s);
        }
    });
}

mio::abm::LocationType get_location_type(const int location_type)
{
    mio::abm::LocationType type;
    switch (location_type) {
    case 0:
        type = mio::abm::LocationType::Home;
        break;
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
        type = mio::abm::LocationType::SocialEvent;
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
    else if (age > 79) {
        return age_group_80_plus;
    }
    else {
        return age_group_0_to_4;
    }
}

void create_world_from_data(mio::abm::World& world, const std::string& filename, const int max_number_persons,
                            const int run_number)
{

    // Open File; we use the cleaned up version of https://zenodo.org/records/13318436
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

    std::map<uint32_t, mio::abm::LocationId> locations = {}; //uint is the location id from the data file
    std::map<uint32_t, mio::abm::Person&> persons      = {}; //uint is the person id from the data file

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
    // Add all home locations to the world and count persons
    uint32_t last_person_id       = 0;
    int number_of_persons_read_in = 0;
    while (std::getline(fin, line)) {
        row.clear();

        // read columns in this row
        split_line(line, &row);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        uint32_t person_id = row[index["puid"]];
        if (person_id != last_person_id) {
            number_of_persons_read_in++;
            last_person_id = person_id;
            if (number_of_persons_read_in > max_number_persons) {
                break;
            }
        }

        uint32_t home_id = row[index["huid"]];

        mio::abm::LocationId home;
        auto it_home = locations.find(home_id);
        if (it_home == locations.end()) {
            home = world.add_location(mio::abm::LocationType::Home);
            locations.insert({home_id, home});
        }
    }

    fin.clear();
    fin.seekg(0);
    std::getline(fin, line); // Skip header row
    last_person_id            = 0;
    number_of_persons_read_in = 0;
    while (std::getline(fin, line)) {
        row.clear();

        // read columns in this row
        split_line(line, &row);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        uint32_t person_id = row[index["puid"]];
        if (person_id != last_person_id) {
            number_of_persons_read_in++;
            last_person_id = person_id;
            if (number_of_persons_read_in > max_number_persons) {
                break;
            }
        }

        mio::abm::LocationId location;
        int target_location_id                           = row[index["loc_id_end"]];
        uint32_t location_type                           = row[index["location_type"]];
        mio::abm::GeographicalLocation location_long_lat = {(double)row[index["lon_end"]] / 1e+5,
                                                            (double)row[index["lat_end"]] / 1e+5};
        auto it_location                                 = locations.find(
            target_location_id); // Check if location already exists also for home which have the same id (home_id = target_location_id)
        if (it_location == locations.end()) {
            location = world.add_location(get_location_type(location_type));
            locations.insert({target_location_id, location});
            world.get_individualized_location(location).set_geographical_location(location_long_lat);
        }
    }

    fin.clear();
    fin.seekg(0);
    std::getline(fin, line); // Skip header row
    last_person_id            = 0;
    number_of_persons_read_in = 0;
    std::vector<mio::abm::Trip> trip_vec;
    while (std::getline(fin, line)) {
        row.clear();

        // read columns in this row
        split_line(line, &row);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        uint32_t person_id = row[index["puid"]];
        bool has_home_trip = true;
        if (person_id != last_person_id) {
            number_of_persons_read_in++;
            last_person_id = person_id;
            if (number_of_persons_read_in > max_number_persons) {
                break;
            }
            has_home_trip = row[index["has_home_trip"]];
        }

        uint32_t age                = row[index["age"]];
        uint32_t home_id            = row[index["huid"]];
        uint32_t target_location_id = row[index["loc_id_end"]];
        uint32_t trip_start         = row[index["start_time"]];
        bool home_in_bs             = row[index["home_in_bs"]];

        // Add the trip to the trip list, person and location must exist at this point
        auto target_location = locations.find(target_location_id)->second;
        auto it_person       = persons.find(person_id);
        if (it_person == persons.end()) {
            auto home    = locations.find(home_id)->second;
            auto& person = world.add_person(home, determine_age_group(age), run_number * 1000000);
            person.set_mask_preferences({0.0, -0.1, -0.1, -0.3, -0.2, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, -0.0});
            person.set_assigned_location(home);
            person.set_assigned_location(hospital);
            person.set_assigned_location(icu);
            persons.insert({person_id, person});
            home_in_bs ? person.set_should_be_logged(true) : person.set_should_be_logged(false);
            it_person = persons.find(person_id);
        }
        if (target_location.type != mio::abm::LocationType::Home) {
            it_person->second.set_assigned_location(target_location);
        }
        trip_vec.push_back(mio::abm::Trip(it_person->second.get_person_id(),
                                          mio::abm::TimePoint(0) + mio::abm::minutes(trip_start), target_location));
        if (!has_home_trip) {
            trip_vec.push_back(mio::abm::Trip(it_person->second.get_person_id(),
                                              mio::abm::TimePoint(0) + mio::abm::hours(17),
                                              locations.find(home_id)->second));
        }
    }
    world.get_trip_list().add_several_trips(trip_vec);
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

    //Set testing parameters
    auto pcr_test_values     = mio::abm::TestParameters{0.9, 0.995};
    auto antigen_test_values = mio::abm::TestParameters{0.71, 0.996}; //https://doi.org/10.1016/j.eclinm.2021.100954
    auto generic_test_values = mio::abm::TestParameters{0.7, 0.95};
    params.get<mio::abm::TestData>()[mio::abm::TestType::PCR]     = pcr_test_values;
    params.get<mio::abm::TestData>()[mio::abm::TestType::Antigen] = antigen_test_values;
    params.get<mio::abm::TestData>()[mio::abm::TestType::Generic] = generic_test_values;

    // Set percentage parameters
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_0_to_4}]   = 0.50;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_5_to_14}]  = 0.55;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_15_to_34}] = 0.60;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_35_to_59}] = 0.70;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_60_to_79}] = 0.83;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_80_plus}]  = 0.90;

    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_0_to_4}]   = 0.02;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_5_to_14}]  = 0.03;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_15_to_34}] = 0.04;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_35_to_59}] = 0.07;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_60_to_79}] = 0.17;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_80_plus}]  = 0.24;

    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_0_to_4}]   = 0.1;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_5_to_14}]  = 0.11;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_15_to_34}] = 0.12;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_35_to_59}] = 0.14;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_60_to_79}] = 0.33;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_80_plus}]  = 0.62;

    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_0_to_4}]   = 0.12;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_5_to_14}]  = 0.13;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_15_to_34}] = 0.15;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_35_to_59}] = 0.26;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_60_to_79}] = 0.40;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_80_plus}]  = 0.48;

    // Set infection parameters
    // Set protection level against an severe infection. Information based on: https://doi.org/10.1016/j.ebiom.2023.104734, https://www.sciencedirect.com/science/article/pii/S2590113322000062
    params.get<mio::abm::SeverityProtectionFactor>() = [](ScalarType days) -> ScalarType {
        return mio::linear_interpolation_of_data_set<ScalarType, ScalarType>({{0, 0.8}, {150, 0.8}}, days);
    };

    // Set protection level against an infection. Information based on: 10.1016/j.vaccine.2023.03.069
    params.get<mio::abm::InfectionProtectionFactor>() = [](ScalarType days) -> ScalarType {
        return mio::linear_interpolation_of_data_set<ScalarType, ScalarType>({{0, 0.0}, {150, 0.0}}, days);
    };

    //Set other parameters
    params.get<mio::abm::AerosolTransmissionRates>() = 0.0;
}

// set location specific parameters
void set_local_parameters(mio::abm::World& world)
{
    const int n_age_groups = (int)world.parameters.get_num_groups();

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

    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_random(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 1.0);

    for (auto& loc : world.get_locations()) {
        switch (loc.get_type()) {
        case mio::abm::LocationType::Home:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_home;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1.6*1; //15 hours
            break;
        case mio::abm::LocationType::School:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_school;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 12.0*0; //2 hours
            break;
        case mio::abm::LocationType::Work:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_work;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 8.0*0; // 3 hours
            break;
        case mio::abm::LocationType::SocialEvent:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1.2;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 8.0*0; // 3 hours
            break;
        case mio::abm::LocationType::BasicsShop:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.8;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 12.0*0; // 2 hours
            break;
        default:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_random;
            break;
        }
    }
}

std::vector<int> read_in_icu(std::vector<mio::DiviEntry> divi_data, mio::Date start_date, int max_num_days)
{
    std::vector<int> icu_data;
    for (auto& entry : divi_data) {
        if (entry.county_id->get() == 3101) {
            if (entry.date >= start_date && entry.date < mio::offset_date_by_days(start_date, max_num_days)) {
                icu_data.push_back(entry.num_icu);
            }
        }
    }
    return icu_data;
}

std::vector<int> read_in_deaths(std::vector<mio::ConfirmedCasesDataEntry> rki_data, mio::Date start_date,
                                int max_num_days)
{
    std::vector<std::vector<int>> death_data_age{num_age_groupss};
    auto date_need = mio::offset_date_by_days(start_date, -18);
    for (auto& entry : rki_data) {
        if (entry.county_id->get() == 3101) {
            if (entry.date >= date_need && entry.date < mio::offset_date_by_days(date_need, max_num_days)) {
                int age_group = entry.age_group.get();
                death_data_age.at(age_group).push_back(entry.num_deaths);
            }
        }
    }
    std::vector<int> death_data;
    for (int i = 0; i < max_num_days; i++) {
        int sum = 0;
        for (size_t j = 0; j < death_data_age.size(); j++) {
            sum += death_data_age[j][i];
        }
        death_data.push_back(sum);
    }

    return death_data;
}

std::vector<int> read_in_detected(std::vector<mio::ConfirmedCasesDataEntry> rki_data, mio::Date start_date,
                                  int max_num_days)
{
    std::vector<std::vector<double>> conf_data_age{num_age_groupss};
    for (auto& entry : rki_data) {
        if (entry.county_id->get() == 3101) {
            if (entry.date >= start_date && entry.date < mio::offset_date_by_days(start_date, max_num_days)) {
                auto age_group = entry.age_group.get();
                conf_data_age.at(age_group).push_back(entry.num_confirmed);
            }
        }
    }
    double sum_day_minus_one = 0;
    for (size_t j = 0; j < conf_data_age.size(); j++) {
        sum_day_minus_one += conf_data_age[j][0];
    }

    std::vector<int> conf_data;
    for (int i = 0; i < max_num_days; i++) {
        double sum = 0;
        for (size_t j = 0; j < conf_data_age.size(); j++) {
            sum += conf_data_age[j][i];
        }
        conf_data.push_back((int)(sum - sum_day_minus_one));
    }

    return conf_data;
}

double calculate_rmse_from_results(const fs::path& data_dir, mio::TimeSeries<ScalarType> sim_inf_states,
                                   mio::TimeSeries<ScalarType> sim_det, int max_num_days, mio::Date start_date)
{
    // We need to read in the results from the results directory
    auto real_data_dead_path =
        mio::path_join((data_dir / "pydata" / "Germany").string(), "cases_all_county_age_ma1.json");
    auto real_data_icu_path = mio::path_join((data_dir / "pydata" / "Germany").string(), "county_divi.json");
    auto real_detected_path =
        mio::path_join((data_dir / "pydata" / "Germany").string(), "cases_all_county_age_repdate_ma1.json");
    auto divi_data     = mio::read_divi_data(real_data_icu_path);
    auto death_data    = mio::read_confirmed_cases_data(real_data_dead_path);
    auto detected_data = mio::read_confirmed_cases_data(real_detected_path);
    // read in the real data
    auto real_data_dead_vec = read_in_deaths(death_data.value(), start_date, max_num_days);
    auto real_data_icu_vec  = read_in_icu(divi_data.value(), start_date, max_num_days);
    auto real_data_conf_vec = read_in_detected(detected_data.value(), start_date, max_num_days);

    // Simulated data
    std::vector<int> sim_data_vec_icu;
    std::vector<int> sim_data_vec_dead;
    std::vector<int> sim_data_vec_conf;
    for (int i = 0; i < max_num_days; i++) {
        int number_of_persons_in_icu = 0;
        int number_of_persons_dead   = 0;
        int number_of_det            = 0;
        for (size_t j = 0; j < (size_t)num_age_groupss; j++) {
            auto index_icu = (((size_t)(mio::abm::InfectionState::Count)) * (j)) +
                             ((uint32_t)(mio::abm::InfectionState::InfectedCritical));
            auto index_dead =
                (((size_t)(mio::abm::InfectionState::Count)) * (j)) + ((uint32_t)(mio::abm::InfectionState::Dead));
            number_of_persons_in_icu += sim_inf_states[i * 24][index_icu];
            number_of_persons_dead += sim_inf_states[i * 24][index_dead];
            number_of_det += sim_det[i * 24][j];
        }
        sim_data_vec_conf.push_back(number_of_det);
        sim_data_vec_icu.push_back(number_of_persons_in_icu);
        sim_data_vec_dead.push_back(number_of_persons_dead);
    }

    // for debugging we print the simulated and real data
    // std::cout << "Real data dead: ";
    // for (auto& entry : real_data_dead_vec) {
    //     std::cout << entry << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "Simulated data dead: ";
    // for (auto& entry : sim_data_vec_dead) {
    //     std::cout << entry << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "Real data icu: ";
    // for (auto& entry : real_data_icu_vec) {
    //     std::cout << entry << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "Simulated data icu: ";
    // for (auto& entry : sim_data_vec_icu) {
    //     std::cout << entry << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "Real data conf: ";
    // for (auto& entry : real_data_conf_vec) {
    //     std::cout << entry << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "Simulated data conf: ";
    // for (auto& entry : sim_data_vec_conf) {
    //     std::cout << entry << " ";
    // }
    // std::cout << std::endl;
    // now we calculate the RMSE
    double rmse_dead = 0;
    double rmse_icu  = 0;
    double rmse_conf = 0;
    for (size_t i = 0; i < real_data_dead_vec.size(); i++) {
        rmse_dead += pow(real_data_dead_vec[i] - sim_data_vec_dead[i], 2);
        rmse_icu += pow((int)(sim_data_vec_icu[i] * 0.55) - real_data_icu_vec[i], 2);
        rmse_conf += pow(real_data_conf_vec[i] - sim_data_vec_conf[i], 2);
    }
    rmse_dead = rmse_dead / real_data_dead_vec.size();
    rmse_icu  = rmse_icu / real_data_icu_vec.size();
    rmse_conf = rmse_conf / real_data_conf_vec.size();

    //write  to terminal

    std::cout << "RMSE dead: " << rmse_dead << std::endl;
    std::cout << "RMSE icu: " << rmse_icu << std::endl;
    std::cout << "RMSE conf: " << rmse_conf << std::endl;

    return (1.00 * rmse_dead) + (0.1 * rmse_icu) + (0.01 * 0.01 * 3.0 * rmse_conf);
}

/**
 * @brief Calculate a grid search for a given set of parameters.
 * @input std::vector where size is the amount of parameters, and the first entry is the min and the second is the max value
 */
std::vector<std::vector<double>> grid_points(std::vector<std::pair<double, double>> parameter_boundaries,
                                             std::vector<int> number_of_points)
{
    std::vector<std::vector<double>> grid;
    for (size_t i = 0; i < parameter_boundaries.size(); i++) {
        std::vector<double> temp;
        double step = (parameter_boundaries[i].second - parameter_boundaries[i].first) / (number_of_points.at(i) - 1);
        if (number_of_points.at(i) == 1) {
            step = 0;
        }
        for (int j = 0; j < number_of_points.at(i); j++) {
            temp.push_back(parameter_boundaries[i].first + j * step);
        }
        grid.push_back(temp);
    }
    return grid;
}

/**
 * @brief Calculate a grid search for a given set of parameters.
 * @input std::vector where size is the amount of parameters, and the first entry is the min and the second is the max value
 */
std::vector<std::vector<double>> grid_points(const std::vector<double>& parameter_points,
                                             const std::vector<int>& number_of_points)
{
    std::vector<std::vector<double>> grid;
    for (size_t i = 0; i < parameter_points.size(); i++) {
        std::vector<double> temp;
        double min_value = parameter_points[i] * 0.95;
        double max_value = parameter_points[i] * 1.05;
        double step      = (max_value - min_value) / (number_of_points.at(i) - 1);
        if (number_of_points.at(i) == 1) {
            step = 0;
        }
        for (int j = 0; j < number_of_points.at(i); j++) {
            temp.push_back(min_value + j * step);
        }
        grid.push_back(temp);
    }
    return grid;
}

std::vector<std::vector<double>> every_combination_of_parameters(std::vector<std::vector<double>> grid)
{
    // we get n parameters and several values for each parameter. We want to get all possible combinations of the parameters
    std::vector<std::vector<double>> parameter_combinations;
    std::vector<int> counter_per_dimension(grid.size(), 0);
    int number_of_points = 1;
    for (size_t i = 0; i < grid.size(); i++) {
        number_of_points = number_of_points * grid[i].size();
    }
    for (int i = 0; i < number_of_points; i++) {
        std::vector<double> temp;
        for (size_t j = 0; j < grid.size(); j++) {
            temp.push_back(grid[j][counter_per_dimension[j]]);
        }
        parameter_combinations.push_back(temp);
        // we increase the counter for the last dimension
        counter_per_dimension.back()++;
        // we increase the counter for the other dimensions if the last dimension has reached the end
        for (int k = (int)grid.size(); k > 0; k--) {
            if (counter_per_dimension[k] == (int)grid[k].size()) {
                counter_per_dimension[k] = 0;
                counter_per_dimension[k - 1]++;
            }
        }
    }
    return parameter_combinations;
}

/**
 * @brief Distribute the grid search over the MPI ranks.
 */
std::vector<std::vector<double>> distribute_grid_search(int rank, int num_procs, std::vector<std::vector<double>> grid)
{
    //Calculate how many grid points there are, assuming that each parameter has the same amount of points
    int number_of_points = 1;
    for (size_t i = 0; i < grid.size(); i++) {
        number_of_points = number_of_points * grid[i].size();
    }
    //Calculate how many points each rank should calculate
    int points_per_rank = number_of_points / num_procs;
    //leftover points goes to the last rank
    if (rank == num_procs - 1) {
        points_per_rank = number_of_points - points_per_rank * (num_procs - 1);
    }
    // we calculate every possible combination of the grid, independently of the rank
    std::vector<std::vector<double>> grid_search{};
    std::vector<size_t> counter_per_dimension(grid.size(), 0);
    for (int i = 0; i < number_of_points; i++) {
        std::vector<double> temp;
        for (size_t j = 0; j < grid.size(); j++) {
            temp.push_back(grid[j][counter_per_dimension[j]]);
        }
        grid_search.push_back(temp);
        // we increase the counter for the last dimension
        counter_per_dimension[grid.size() - 1]++;
        // we increase the counter for the other dimensions if the last dimension has reached the end
        for (size_t k = grid.size() - 1; k > 0; k--) {
            if (counter_per_dimension[k] == grid[k].size()) {
                counter_per_dimension[k] = 0;
                counter_per_dimension[k - 1]++;
            }
        }
    }
    // we calculate the grid search for the rank
    std::vector<std::vector<double>> grid_search_ranks;
    for (int i = 0; i < points_per_rank; i++) {
        grid_search_ranks.push_back(grid_search[i + rank * points_per_rank]);
    }

    return grid_search_ranks;
}

void create_easter_social_event(mio::abm::World& world, double perc_easter_event)
{
    // int number_of_persons_per_easter_event = 6;
    // int actual_number_of_persons           = 0;
    // auto event                             = world.add_location(mio::abm::LocationType::Event);

    // for (auto& p : world.get_persons()) {
    //     if (actual_number_of_persons >= number_of_persons_per_easter_event) {
    //         actual_number_of_persons = 0;
    //         event                    = world.add_location(mio::abm::LocationType::Event);
    //     }
    //     p.set_goes_to_easter(true);
    //     p.set_assigned_location(event);
    //     actual_number_of_persons++;
    // }
    // we take 50% of each age group and move them to the social event location in a random fashion with 1 person per age group and 2 of the age group 15-34
    std::vector<std::vector<uint32_t>> person_per_age_group;
    person_per_age_group.resize(num_age_groupss);

    for (size_t i = 0; i < world.get_persons().size(); i++) {
        auto& p = world.get_persons()[i];
        person_per_age_group[determine_age_group_from_rki(p.get_age())].push_back(i);
    }

    for (size_t i = 0; i < person_per_age_group.size(); i++) {
        std::shuffle(person_per_age_group[i].begin(), person_per_age_group[i].end(), world.get_rng());
    }
    int number_of_persons       = (int)(perc_easter_event * world.get_persons().size());
    int number_of_social_events = number_of_persons / (num_age_groupss + 1);

    for (auto i = 0; i < number_of_social_events; i++) {
        auto easter_event = world.add_location(mio::abm::LocationType::Event);
        for (size_t j = 0; j < num_age_groupss; j++) {
            if (person_per_age_group[j].empty()) {
                continue;
            }
            world.get_persons()[person_per_age_group[j].back()].set_assigned_location(easter_event);
            world.get_persons()[person_per_age_group[j].back()].set_goes_to_easter(true);
            person_per_age_group[j].pop_back();

            if (j == 3) {
                if (person_per_age_group[j].empty()) {
                    continue;
                }
                world.get_persons()[person_per_age_group[j].back()].set_assigned_location(easter_event);
                world.get_persons()[person_per_age_group[j].back()].set_goes_to_easter(true);
                person_per_age_group[j].pop_back();
            }
            if (j == 4) {
                if (person_per_age_group[j].empty()) {
                    continue;
                }
                world.get_persons()[person_per_age_group[j].back()].set_assigned_location(easter_event);
                world.get_persons()[person_per_age_group[j].back()].set_goes_to_easter(true);
                person_per_age_group[j].pop_back();
            }
            if (j == 5) {
                if (person_per_age_group[j].empty()) {
                    continue;
                }
                world.get_persons()[person_per_age_group[j].back()].set_assigned_location(easter_event);
                world.get_persons()[person_per_age_group[j].back()].set_goes_to_easter(true);
                person_per_age_group[j].pop_back();
            }
        }
    }
}

/**
 * Create a sampled simulation with start time t0.
 * @param t0 The start time of the Simulation.
 */
void create_sampled_world(
    mio::abm::World& world, const fs::path& input_dir, const mio::abm::TimePoint& t0, int max_num_persons,
    mio::Date start_date_sim, double perc_easter_event,
    mio::CustomIndexArray<double, mio::AgeGroup, mio::osecir::InfectionState> initial_infection_distribution,
    std::map<mio::Date, std::vector<std::pair<uint32_t, uint32_t>>> vacc_map, uint32_t run_number)
{
    mio::unused(start_date_sim);
    //Set global infection parameters (similar to infection parameters in SECIR model) and initialize the world

    set_parameters(world.parameters);

    restart_timer(timer, "time taken for setting up parameters and local parameters");

    // Create the world object from statistical data.
    create_world_from_data(world, (input_dir / "mobility/braunschweig_result_ffa8_modified2.csv").generic_string(),
                           max_num_persons, run_number);
    world.use_migration_rules(false);
    restart_timer(timer, "time taken for braunschweig trip input");

    create_easter_social_event(world, perc_easter_event);
    restart_timer(timer, "time taken for setting up easter social ebent");

    // Assign an infection state to each person.
    assign_infection_state(world, t0, initial_infection_distribution);
    restart_timer(timer, "time taken for assigning infection state");

    // Assign vaccination status to each person.
    assign_vaccination_state(world, start_date_sim, vacc_map);
    restart_timer(timer, "time taken for assigning vaccination state");
    set_local_parameters(world);
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
        const auto curr_time = sim.get_time();
        const auto persons   = sim.get_world().get_persons();

        // PRAGMA_OMP(parallel for)
        for (auto i = size_t(0); i < persons.size(); ++i) {
            auto& p = persons[i];
            if (p.get_should_be_logged()) {
                auto index = (((size_t)(mio::abm::InfectionState::Count)) * ((uint32_t)p.get_age().get())) +
                             ((uint32_t)p.get_infection_state(curr_time));
                // PRAGMA_OMP(atomic)
                sum[index] += 1;
            }
        }
        return std::make_pair(curr_time, sum);
    }
};

struct LogInfectionPerLocationTypePerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    /** 
     * @brief Log the TimeSeries of the number of Person%s in an #InfectionState.
     * @param[in] sim The simulation of the abm.
     * @return A pair of the TimePoint and the TimeSeries of the number of Person%s in an #InfectionState.
     */
    static Type log(const mio::abm::Simulation& sim)
    {

        Eigen::VectorXd sum = Eigen::VectorXd::Zero(
            Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups()));
        auto curr_time     = sim.get_time();
        auto prev_time     = sim.get_prev_time();
        const auto persons = sim.get_world().get_persons();

        // PRAGMA_OMP(parallel for)
        for (auto i = size_t(0); i < persons.size(); ++i) {
            auto& p = persons[i];
            if (p.get_should_be_logged()) {
                // PRAGMA_OMP(atomic)
                if ((p.get_infection_state(prev_time) != mio::abm::InfectionState::Exposed) &&
                    (p.get_infection_state(curr_time) == mio::abm::InfectionState::Exposed)) {
                    auto index = (((size_t)(mio::abm::LocationType::Count)) * ((uint32_t)p.get_age().get())) +
                                 ((uint32_t)p.get_location().get_type());
                    sum[index] += 1;
                }
            }
        }
        return std::make_pair(curr_time, sum);
    }
};

struct LogTestPerLocationTypePerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    /** 
     * @brief Log the TimeSeries of the number of Person%s that have been tested since the last time step.
     * @param[in] sim The simulation of the abm.
     * @return A pair of the TimePoint and the TimeSeries of the number of Person%s that have been tested.
     */
    static Type log(const mio::abm::Simulation& sim)
    {

        Eigen::VectorXd sum = Eigen::VectorXd::Zero(
            Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups()));
        const auto curr_time = sim.get_time();
        const auto prev_time = sim.get_prev_time();
        const auto persons   = sim.get_world().get_persons();

        // PRAGMA_OMP(parallel for)
        for (auto i = size_t(0); i < persons.size(); ++i) {
            auto& p = persons[i];
            if (p.get_should_be_logged()) {
                // PRAGMA_OMP(atomic)
                if ((p.get_time_of_last_test() == prev_time)) {
                    auto index = (((size_t)(mio::abm::LocationType::Count)) * ((uint32_t)p.get_age().get())) +
                                 ((uint32_t)p.get_location().get_type());
                    sum[index] += 1;
                }
            }
        }
        return std::make_pair(curr_time, sum);
    }
};

struct LogPositiveTestPerLocationTypePerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    /** 
     * @brief Log the TimeSeries of the number of Person%s that have been tested positive since the last time step.
     * @param[in] sim The simulation of the abm.
     * @return A pair of the TimePoint and the TimeSeries of the number of Person%s that have been tested positive.
     */
    static Type log(const mio::abm::Simulation& sim)
    {

        Eigen::VectorXd sum = Eigen::VectorXd::Zero(
            Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups()));
        const auto curr_time = sim.get_time();
        const auto prev_time = sim.get_prev_time();
        const auto persons   = sim.get_world().get_persons();

        // PRAGMA_OMP(parallel for)
        for (auto i = size_t(0); i < persons.size(); ++i) {
            auto& p = persons[i];
            if (p.get_should_be_logged()) {
                // PRAGMA_OMP(atomic)
                // careful: this check assumes a persons goes into quarantine when a test is positive and that the duration of quarantine is >0
                if (p.get_time_of_last_test() == prev_time &&
                    p.is_in_quarantine(curr_time, sim.get_world().parameters)) {
                    auto index = (((size_t)(mio::abm::LocationType::Count)) * ((uint32_t)p.get_age().get())) +
                                 ((uint32_t)p.get_location().get_type());
                    sum[index] += 1;
                }
            }
        }
        return std::make_pair(curr_time, sum);
    }
};

struct LogCumulativeDetectedInfectionsPerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    /** 
     * @brief Log the TimeSeries of the number of Person%s in an #InfectionState.
     * @param[in] sim The simulation of the abm.
     * @return A pair of the TimePoint and the TimeSeries of the number of Person%s in an #InfectionState.
     */
    static Type log(const mio::abm::Simulation& sim)
    {

        Eigen::VectorXd sum  = Eigen::VectorXd::Zero(Eigen::Index(sim.get_world().parameters.get_num_groups()));
        const auto curr_time = sim.get_time();
        const auto persons   = sim.get_world().get_persons();

        // PRAGMA_OMP(parallel for)
        for (auto i = size_t(0); i < persons.size(); ++i) {
            auto& p = persons[i];
            if (p.get_should_be_logged()) {
                // PRAGMA_OMP(atomic)
                if (p.was_person_tested_positive()) {
                    uint32_t index = (uint32_t)p.get_age().get();
                    sum[index] += 1;
                }
            }
        }
        return std::make_pair(curr_time, sum);
    }
};

struct LogEstimatedReproductionNumber : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    /** 
     * @brief Log the TimeSeries of the estimated reproduction number.
     * We use the formula from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3816335/
     * @param[in] sim The simulation of the abm.
     * @return A pair of the TimePoint and the TimeSeries of the estimated reproduction number.
     */
    static Type log(const mio::abm::Simulation& sim)
    {
        Eigen::VectorXd estimation = Eigen::VectorXd::Zero(Eigen::Index(1));

        // time period to take into account for estimating the reproduction number
        // longer periods lead to more averaged results
        mio::abm::TimeSpan time_frame = sim.get_dt();
        const auto t                  = sim.get_time();
        const auto persons            = sim.get_world().get_persons();

        // PRAGMA_OMP(parallel for)
        int number_newly_infected  = 0;
        double infection_incidence = 0;
        for (auto i = size_t(0); i < persons.size(); ++i) {
            auto& p = persons[i];
            if (p.get_should_be_logged()) {
                // PRAGMA_OMP(atomic)
                if (p.is_infected(t) && !p.is_infected(t - time_frame)) {
                    number_newly_infected += 1;
                }
                // sum the total infection incidence from all infected people in the desired time frame
                if (p.is_infected(t) || p.is_infected(t - time_frame)) {
                    // this is the integral of the viral shed between t-time_frame and t
                    // divided by the total integral to normalize it to a density function
                    infection_incidence += p.get_infection().get_viral_shed_integral(t - time_frame, t) /
                                           p.get_infection().get_viral_shed_integral();
                }
            }
        }
        // if the infection dies out, the reproduction number has no meaningful value
        if (infection_incidence == 0) {
            return std::make_pair(t, estimation);
        }

        estimation[0] = number_newly_infected / infection_incidence;
        return std::make_pair(t, estimation);
    }
};

mio::TimeSeries<ScalarType> get_new_detected_infections(mio::TimeSeries<ScalarType> cumulative_infections)
{
    mio::TimeSeries<ScalarType> new_detected_infection(mio::TimeSeries<ScalarType>::zero(
        cumulative_infections.get_num_time_points() / 24, cumulative_infections.get_num_elements()));
    for (Eigen::Index i = 1; i < cumulative_infections.get_num_time_points() / 24; i++) {
        new_detected_infection[i] = cumulative_infections[i * 24] - cumulative_infections[(i - 1) * 24];
    }

    return new_detected_infection;
}

#ifdef MEMILIO_ENABLE_MPI
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

void get_grid_search_results_and_write_them_to_file(
    int rank, int num_procs, const fs::path& result_dir,
    std::vector<std::pair<std::vector<double>, double>> grid_my_rank_with_rmse)
{
    // we just send every grid point and the corresponding rmse to rank 0 who writes it to a file
    auto gathered_grid_vector_with_rmse = std::vector<std::pair<std::vector<double>, double>>{};

    if (rank == 0) {
        std::copy(grid_my_rank_with_rmse.begin(), grid_my_rank_with_rmse.end(),
                  std::back_inserter(gathered_grid_vector_with_rmse));

        for (int src_rank = 1; src_rank < num_procs; ++src_rank) {
            int bytes_size;
            MPI_Recv(&bytes_size, 1, MPI_INT, src_rank, 0, mio::mpi::get_world(), MPI_STATUS_IGNORE);
            mio::ByteStream bytes(bytes_size);
            MPI_Recv(bytes.data(), bytes.data_size(), MPI_BYTE, src_rank, 0, mio::mpi::get_world(), MPI_STATUS_IGNORE);

            auto src_grid_results = mio::deserialize_binary(bytes, mio::Tag<decltype(grid_my_rank_with_rmse)>{});
            std::copy(src_grid_results.value().begin(), src_grid_results.value().end(),
                      std::back_inserter(gathered_grid_vector_with_rmse));
        }
    }
    else {
        auto bytes_grid = mio::serialize_binary(grid_my_rank_with_rmse);
        auto bytes_size = int(bytes_grid.data_size());
        MPI_Send(&bytes_size, 1, MPI_INT, 0, 0, mio::mpi::get_world());
        MPI_Send(bytes_grid.data(), bytes_grid.data_size(), MPI_BYTE, 0, 0, mio::mpi::get_world());
    }

    // write the gathered grid search results to a file
    if (rank == 0) {
        std::ofstream file((result_dir / "/grid_search/grid_search_results.txt").string());
        if (!file.is_open()) {
            mio::log_error("Error opening file for writing grid search results.");
        }
        for (size_t i = 0; i < gathered_grid_vector_with_rmse.size(); i++) {
            file << "Grid point: ";
            for (size_t j = 0; j < gathered_grid_vector_with_rmse[i].first.size(); j++) {
                file << gathered_grid_vector_with_rmse[i].first.at(j) << " ";
            }
            file << "RMSE: " << gathered_grid_vector_with_rmse[i].second << std::endl;
        }

        // Sort the grid search results by RMSE in ascending order
        std::sort(gathered_grid_vector_with_rmse.begin(), gathered_grid_vector_with_rmse.end(),
                  [](const auto& a, const auto& b) {
                      return a.second < b.second;
                  });

        // Write the 5 best RMSE at the end of the file
        file << std::endl << "Top 20 RMSE:" << std::endl;
        for (size_t i = 0; i < std::min<size_t>(20, gathered_grid_vector_with_rmse.size()); i++) {
            file << "Grid point: ";
            for (size_t j = 0; j < gathered_grid_vector_with_rmse[i].first.size(); j++) {
                file << gathered_grid_vector_with_rmse[i].first.at(j) << " ";
            }
            file << "RMSE: " << gathered_grid_vector_with_rmse[i].second << std::endl;
        }

        file.close();
    }
}

#endif

void add_npi_testing_strategies_to_world(mio::abm::Simulation& sim, mio::abm::TimePoint tmax,
                                         double testing_probability_sympt, double ratio_asympt_to_sympt,
                                         double lockdown_prob, double after_lockdown_prob)
{
    double testing_probability_asympt = testing_probability_sympt / ratio_asympt_to_sympt;
    auto start_date_test              = mio::abm::TimePoint(mio::abm::days(0).seconds());
    auto lockdown_start_date          = mio::abm::TimePoint(mio::abm::days(29).seconds());
    auto easter_end_date              = mio::abm::TimePoint(mio::abm::days(37).seconds());
    auto lockdown_end_date            = mio::abm::TimePoint(mio::abm::days(60).seconds());
    auto end_date_test                = tmax;

    auto antigen_test = mio::abm::TestType::Antigen;
    auto antigen_test_parameters =
        sim.get_world().parameters.get<mio::abm::TestData>()[antigen_test]; // Test parameters
    auto testing_min_time     = mio::abm::days(3);
    auto vector_sympt_states  = std::vector<mio::abm::InfectionState>{mio::abm::InfectionState::InfectedSymptoms,
                                                                      mio::abm::InfectionState::InfectedSevere,
                                                                      mio::abm::InfectionState::InfectedCritical};
    auto vector_asympt_states = std::vector<mio::abm::InfectionState>{
        mio::abm::InfectionState::InfectedNoSymptoms, mio::abm::InfectionState::Exposed,
        mio::abm::InfectionState::Susceptible, mio::abm::InfectionState::Recovered};

    auto testing_criteria_asympt = mio::abm::TestingCriteria({}, vector_asympt_states);
    auto testing_criteria_sympt  = mio::abm::TestingCriteria({}, vector_sympt_states);

    // We go through each location type and add the testing strategies

    auto testing_scheme_asympt =
        mio::abm::TestingScheme(testing_criteria_asympt, testing_min_time, start_date_test, lockdown_start_date,
                                antigen_test_parameters, testing_probability_asympt);
    auto testing_scheme_sympt =
        mio::abm::TestingScheme(testing_criteria_sympt, testing_min_time, start_date_test, lockdown_start_date,
                                antigen_test_parameters, testing_probability_sympt);
    auto testing_scheme_asympt_easter =
        mio::abm::TestingScheme(testing_criteria_asympt, testing_min_time, lockdown_start_date, easter_end_date,
                                antigen_test_parameters, 0.66 * testing_probability_asympt);
    auto testing_scheme_sympt_easter =
        mio::abm::TestingScheme(testing_criteria_sympt, testing_min_time, lockdown_start_date, easter_end_date,
                                antigen_test_parameters, 0.66 * testing_probability_sympt);
    auto testing_scheme_asympt_wl =
        mio::abm::TestingScheme(testing_criteria_asympt, testing_min_time, easter_end_date, lockdown_end_date,
                                antigen_test_parameters, lockdown_prob * testing_probability_asympt);
    auto testing_scheme_sympt_wl =
        mio::abm::TestingScheme(testing_criteria_sympt, testing_min_time, easter_end_date, lockdown_end_date,
                                antigen_test_parameters, lockdown_prob * testing_probability_sympt);
    auto testing_scheme_asympt_al =
        mio::abm::TestingScheme(testing_criteria_asympt, testing_min_time, lockdown_end_date, end_date_test,
                                antigen_test_parameters, after_lockdown_prob * testing_probability_asympt);
    auto testing_scheme_sympt_al =
        mio::abm::TestingScheme(testing_criteria_sympt, testing_min_time, lockdown_end_date, end_date_test,
                                antigen_test_parameters, after_lockdown_prob * testing_probability_sympt);

    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::School, testing_scheme_asympt);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::School, testing_scheme_sympt);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::School,
                                                              testing_scheme_asympt_easter);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::School,
                                                              testing_scheme_sympt_easter);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::School, testing_scheme_asympt_wl);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::School, testing_scheme_sympt_wl);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::School, testing_scheme_asympt_al);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::School, testing_scheme_sympt_al);

    // Work
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme_asympt);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme_sympt);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work,
                                                              testing_scheme_asympt_easter);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work,
                                                              testing_scheme_sympt_easter);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme_asympt_wl);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme_sympt_wl);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme_asympt_al);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme_sympt_al);

    //basic shops:
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::BasicsShop,
                                                              testing_scheme_asympt);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::BasicsShop, testing_scheme_sympt);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::BasicsShop,
                                                              testing_scheme_asympt_easter);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::BasicsShop,
                                                              testing_scheme_sympt_easter);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::BasicsShop,
                                                              testing_scheme_asympt_wl);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::BasicsShop,
                                                              testing_scheme_sympt_wl);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::BasicsShop,
                                                              testing_scheme_asympt_al);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::BasicsShop,
                                                              testing_scheme_sympt_al);

    // social events
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::SocialEvent,
                                                              testing_scheme_asympt);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::SocialEvent,
                                                              testing_scheme_sympt);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::SocialEvent,
                                                              testing_scheme_asympt_easter);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::SocialEvent,
                                                              testing_scheme_sympt_easter);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::SocialEvent,
                                                              testing_scheme_asympt_wl);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::SocialEvent,
                                                              testing_scheme_sympt_wl);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::SocialEvent,
                                                              testing_scheme_asympt_al);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::SocialEvent,
                                                              testing_scheme_sympt_al);

    // easter event
    // sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Event, testing_scheme_asympt_wl);
    // sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Event, testing_scheme_sympt_wl);
}

void write_grid_search_prematurely_to_file(int rank, const fs::path& result_dir,
                                           std::pair<std::vector<double>, double> rmse)
{
    //make file if there is no file for my rank

    //file name  is grid_search_results_"my_rank".txt
    std::string filename = "grid_search_results_" + std::to_string(rank) + ".txt";
    // look if file exists in folder
    bool file_exists = boost::filesystem::exists(mio::path_join((result_dir / "grid_search").string(), filename));
    if (file_exists) {
        // we open it and attach the new rmse double to the end
        std::ofstream curr_file((result_dir / "/grid_search/" / filename).string(), std::ios::app);
        // write the new rmse with parameters
        curr_file << "Grid point: ";
        for (size_t j = 0; j < rmse.first.size(); j++) {
            curr_file << rmse.first.at(j) << " ";
        }
        curr_file << "RMSE: " << rmse.second << std::endl;
    }
    else {
        // we create a new file and write the rmse double to it
        std::ofstream curr_file((result_dir / "/grid_search/" / filename).string());
        curr_file << "Grid point: ";
        for (size_t j = 0; j < rmse.first.size(); j++) {
            curr_file << rmse.first.at(j) << " ";
        }
        curr_file << "RMSE: " << rmse.second << std::endl;
    }
}

mio::IOResult<void> run_with_grid_search(const fs::path& input_dir, const fs::path& result_dir, int num_runs,
                                         std::vector<std::vector<double>> grid_points, mio::RandomNumberGenerator rng)
{
    int num_procs, rank;
#ifdef MEMILIO_ENABLE_MPI
    MPI_Comm_size(mio::mpi::get_world(), &num_procs);
    MPI_Comm_rank(mio::mpi::get_world(), &rank);
#else
    num_procs = 1;
    rank      = 0;
#endif

    // define parameters for grid search
    // Distribute the grid search over the MPI ranks
    auto grid_search_rank = distribute_grid_search(rank, num_procs, grid_points);
    // short debug print to see if everything worked. Printing rank and amount of grid points as well as first point
    std::cout << "Rank: " << rank << " has " << grid_search_rank.size() << " grid points" << std::endl;
    std::cout << "First grid point: ";
    for (size_t j = 0; j < grid_search_rank[0].size(); j++) {
        std::cout << grid_search_rank[0].at(j) << " ";
    }
    std::cout << std::endl;

    // Run the grid search
    std::vector<double> rmse_results_per_grid_point;
    rmse_results_per_grid_point.resize(grid_search_rank.size());

#pragma omp parallel for num_threads(96) firstprivate(rng)
    for (size_t i = 0; i < grid_search_rank.size(); i++) {
        auto params = grid_search_rank[i];

        printf("I am Thread %d\n", omp_get_thread_num());
        // print my parameters
        std::cout << "Rank: " << rank << " Thread: " << omp_get_thread_num() << "place bound: " << omp_get_place_num()
                  << " Parameters: ";
        for (size_t j = 0; j < params.size(); j++) {
            std::cout << params.at(j) << " ";
        }
        std::cout << std::endl;

        const double viral_shedding_rate        = params[0];
        const double dark_figure                = params[1];
        const double contact_red_lockdown       = params[2];
        const double damping_community_lockdown = 0.5;
        // const double testing_probability_sympt  = 0.037;
        // const double ratio_asympt_to_sympt      = 20.0;
        const double testing_probability_sympt = params[3];
        const double ratio_asympt_to_sympt     = params[4];

        const double lockdown_test_prob     = 1.2;
        const auto after_lockdown_test_prob = 1.0;

        const auto seasonality_april = 0.95;
        const auto seasonality_may   = 0.85;

        const double masks                            = 0.25;
        const double after_lockdown_contact_reduction = 0.50;

        const double perc_easter_event        = 0.2;
        const auto quarantine_duration        = mio::abm::days(10);
        const double quarantine_effectiveness = 0.5;

        mio::Date start_date{2021, 3, 1};
        int date_of_lockdown     = 29;
        int end_date_of_lockdown = 60;
        int max_num_days         = 90;
        auto max_num_persons     = 400000;

        auto t0   = mio::abm::TimePoint(0); // Start time per simulation
        auto tmax = mio::abm::TimePoint(0) + mio::abm::days(max_num_days); // End time per simulation

        // Determine inital infection state distribution
        restart_timer(timer, "time for initial setup");
        auto initial_infection = determine_initial_infection_states_world(input_dir, start_date, dark_figure);
        restart_timer(timer, "time for determine_initial_infection_states_world");
        auto vacc_map = prepare_vaccination_state(mio::offset_date_by_days(start_date, (int)tmax.days()),
                                                  (input_dir / "pydata/Germany/vacc_county_ageinf_ma7.json").string());
        restart_timer(timer, "time for vaccinaiton state");
        for (int j = 0; j < num_runs; j++) {

            // Loop over a number of runs
            auto world = mio::abm::World(num_age_groupss);
            auto run_rng_counter =
                mio::rng_totalsequence_counter<uint64_t>(static_cast<uint32_t>(j), mio::Counter<uint32_t>(0));
            rng.set_counter(run_rng_counter);
            world.get_rng() = rng;

            create_sampled_world(world, input_dir, t0, max_num_persons, start_date, perc_easter_event,
                                 initial_infection, vacc_map, j);

            restart_timer(timer, "time taken for create sampled world");
            auto sim = mio::abm::Simulation(t0, std::move(world));

            //Logger
            mio::History<mio::abm::TimeSeriesWriter, LogInfectionPerLocationTypePerAgeGroup>
                historyInfectionPerLocationType{
                    Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups())};
            mio::History<mio::abm::TimeSeriesWriter, LogInfectionStatePerAgeGroup> historyInfectionStatePerAgeGroup{
                Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups())};
            mio::History<mio::abm::TimeSeriesWriter, LogCumulativeDetectedInfectionsPerAgeGroup>
                historyCumulativeDetectedInfectionsPerAgeGroup{
                    Eigen::Index(sim.get_world().parameters.get_num_groups())};

            // / NPIS//

            const auto location_it = sim.get_world().get_locations();

            // 1. Add testing strategies

            add_npi_testing_strategies_to_world(sim, tmax, testing_probability_sympt, ratio_asympt_to_sympt,
                                                lockdown_test_prob, after_lockdown_test_prob);

            // 2. Mask schemes for all locations
            // First set all locations to have mask usage, we need ffp2 masks
            for (auto& location : location_it) {
                location.set_required_mask(mio::abm::MaskType::FFP2);
                if (location.get_type() == mio::abm::LocationType::Home) {
                    location.set_npi_active(false);
                }
                else {
                    location.set_npi_active(true);
                }
            }

            // 3. Dampings everywhere except home
            for (auto& location : location_it) {
                if (location.get_type() == mio::abm::LocationType::School) {
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(0).seconds()), 0.5); // from 2021-03-01
                    location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= (0.5);
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(date_of_lockdown).seconds()),
                                         0.00); // from 2021-03-15
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(end_date_of_lockdown).seconds()),
                                         0.5); // from 2021-04-12 till 2021-05-30
                }
                if (location.get_type() == mio::abm::LocationType::BasicsShop) {
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(0).seconds()), 0.8); // from 2021-03-15
                    // location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= (0.5);
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(date_of_lockdown).seconds()),
                                         damping_community_lockdown); // from 2021-03-15
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(end_date_of_lockdown).seconds()),
                                         0.8); // from 2021-03-15
                }
                if (location.get_type() == mio::abm::LocationType::SocialEvent) {
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(0).seconds()), 0.8); // from 2021-03-15
                    // location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= (0.5);
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(date_of_lockdown).seconds()),
                                         damping_community_lockdown); // from 2021-03-15
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(end_date_of_lockdown).seconds()),
                                         0.8); // from 2021-03-15
                }
                if (location.get_type() == mio::abm::LocationType::Work) {
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(0).seconds()), 0.75); // from 2021-03-15
                    // location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= (0.75 * 0.75);
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(date_of_lockdown).seconds()),
                                         0.7); // from 2021-03-15
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(end_date_of_lockdown).seconds()),
                                         0.75); // from 2021-03-15
                }
            }

            sim.get_world().parameters.get<mio::abm::InfectionRateFromViralShed>() = viral_shedding_rate;
            sim.get_world().parameters.get<mio::abm::MaskProtection>()             = masks;
            sim.get_world().parameters.get<mio::abm::QuarantineEffectiveness>()    = quarantine_effectiveness;
            sim.get_world().parameters.get<mio::abm::QuarantineDuration>()         = quarantine_duration;

            sim.advance(mio::abm::TimePoint(mio::abm::days(date_of_lockdown).seconds()),
                        historyInfectionPerLocationType, historyInfectionStatePerAgeGroup,
                        historyCumulativeDetectedInfectionsPerAgeGroup);
            sim.get_world().parameters.get<mio::abm::InfectionRateFromViralShed>() =
                viral_shedding_rate * seasonality_april;

            for (auto& location : location_it) {
                if (location.get_type() != mio::abm::LocationType::Home) {
                    location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= contact_red_lockdown;
                }
            }

            sim.advance(mio::abm::TimePoint(mio::abm::days(end_date_of_lockdown).seconds()),
                        historyInfectionPerLocationType, historyInfectionStatePerAgeGroup,
                        historyCumulativeDetectedInfectionsPerAgeGroup);

            for (auto& location : location_it) {
                if (location.get_type() != mio::abm::LocationType::Home) {
                    location.get_infection_parameters().get<mio::abm::ContactRates>().array() *=
                        (1 / (contact_red_lockdown)) * after_lockdown_contact_reduction;
                }
            }
            sim.get_world().parameters.get<mio::abm::InfectionRateFromViralShed>() =
                viral_shedding_rate * seasonality_may;

            sim.advance(mio::abm::TimePoint(tmax.seconds()), historyInfectionPerLocationType,
                        historyInfectionStatePerAgeGroup, historyCumulativeDetectedInfectionsPerAgeGroup);

            ////Advance till here
            // Stop the clock after sim.advance and calculate the duration

            // TODO: update result of the simulation to be a vector of location result.
            auto temp_sim_infection_per_loc_tpye =
                std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionPerLocationType.get_log())};
            auto temp_sim_infection_state_per_age_group =
                std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionStatePerAgeGroup.get_log())};
            auto temp_sim_cumulative_detected_infections_per_age_group = std::vector<mio::TimeSeries<ScalarType>>{
                std::get<0>(historyCumulativeDetectedInfectionsPerAgeGroup.get_log())};

            auto rmse = calculate_rmse_from_results(input_dir, temp_sim_infection_state_per_age_group[0],
                                                    temp_sim_cumulative_detected_infections_per_age_group[0],
                                                    max_num_days, start_date);
            //write rmse as debug output
            std::cout << "RMSE: " << rmse << std::endl;
            rmse_results_per_grid_point.at(i) += rmse;
        }
        rmse_results_per_grid_point.at(i) /= num_runs;
#pragma omp critical
        {
            write_grid_search_prematurely_to_file(rank, result_dir,
                                                  std::make_pair(grid_search_rank[i], rmse_results_per_grid_point[i]));
        }
    }

    // make the gathered results available to all ranks
    std::vector<std::pair<std::vector<double>, double>> my_results;
    for (size_t i = 0; i < grid_search_rank.size(); i++) {
        my_results.push_back(std::make_pair(grid_search_rank[i], rmse_results_per_grid_point[i]));
    }
#ifdef MEMILIO_ENABLE_MPI
    get_grid_search_results_and_write_them_to_file(rank, num_procs, result_dir, my_results);
#endif

    printf("done.\n");
    return mio::success();
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
                        std::vector<std::vector<double>> parameter_values, mio::RandomNumberGenerator rng,
                        bool save_single_runs = false)
{
    int num_procs, rank;
#ifdef MEMILIO_ENABLE_MPI
    MPI_Comm_size(mio::mpi::get_world(), &num_procs);
    MPI_Comm_rank(mio::mpi::get_world(), &rank);
#else
    num_procs = 1;
    rank      = 0;
#endif
    std::vector<double> rmse_results_per_grid_point;
    rmse_results_per_grid_point.resize(parameter_values.size());
    for (size_t par_i = 0; par_i < parameter_values.size(); par_i++) {
        auto params = parameter_values[par_i];

        std::cout << "Parameter values: ";
        for (size_t j = 0; j < params.size(); j++) {
            std::cout << params.at(j) << " ";
        }
        std::cout << std::endl;

        auto run_distribution = distribute_runs(num_runs, num_procs);
        auto start_run_idx =
            std::accumulate(run_distribution.begin(), run_distribution.begin() + size_t(rank), size_t(0));
        auto end_run_idx = start_run_idx + run_distribution[size_t(rank)];

        const double viral_shedding_rate        = params[0];
        const double dark_figure                = params[1];
        const double contact_red_lockdown       = params[2];
        const double damping_community_lockdown = 0.5;
        // const double testing_probability_sympt  = 0.036;
        const double testing_probability_sympt = params[3];

        const double lockdown_test_prob       = 1.2;
        const double after_lockdown_test_prob = 1.0;

        const auto seasonality_april = 0.95;
        const auto seasonality_may   = 0.85;

        const double masks                            = params[7];
        const double after_lockdown_contact_reduction = 0.50;
        // const double ratio_asympt_to_sympt            = 20.0;
        const double ratio_asympt_to_sympt = params[4];
        const double perc_easter_event     = 0.2;
        // const auto quarantine_duration        = mio::abm::days(10);
        // const double quarantine_effectiveness = 0.5;
        const auto quarantine_duration        = mio::abm::days(params[5]);
        const double quarantine_effectiveness = params[6];

        mio::Date start_date{2021, 3, 1};
        int date_of_lockdown     = 29;
        int end_date_of_lockdown = 60;
        int max_num_days         = 90;
        auto max_num_persons     = 400000;
        bool npis_on             = true;

        auto t0   = mio::abm::TimePoint(0); // Start time per simulation
        auto tmax = mio::abm::TimePoint(0) + mio::abm::days(max_num_days); // End time per simulation

        auto ensemble_infection_per_loc_type_per_age_group = std::vector<
            std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of infection per location type per age group results
        ensemble_infection_per_loc_type_per_age_group.reserve(size_t(num_runs));

        auto ensemble_infection_state_per_age_group =
            std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of infection state per age group results
        ensemble_infection_state_per_age_group.reserve(size_t(num_runs));

        auto ensemble_test_per_loc_type_per_age_group = std::vector<
            std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of tests per location type per age group results
        ensemble_test_per_loc_type_per_age_group.reserve(size_t(num_runs));

        auto ensemble_positive_test_per_loc_type_per_age_group = std::vector<std::vector<
            mio::TimeSeries<ScalarType>>>{}; // Vector of positive tests per location type per age group results
        ensemble_positive_test_per_loc_type_per_age_group.reserve(size_t(num_runs));

        auto ensemble_cumulative_detected_infections =
            std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of cumulative detected infections
        ensemble_cumulative_detected_infections.reserve(size_t(num_runs));

        auto ensemble_estimated_reproduction_number =
            std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of estimated reproduction number
        ensemble_estimated_reproduction_number.reserve(size_t(num_runs));

        auto ensemble_estimated_new_detected_infections =
            std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of estimated reproduction number
        ensemble_estimated_new_detected_infections.reserve(size_t(num_runs));

        auto ensemble_params = std::vector<std::vector<mio::abm::World>>{}; // Vector of all worlds
        ensemble_params.reserve(size_t(num_runs));

        // Determine inital infection state distribution
        restart_timer(timer, "time for initial setup");
        auto initial_infection = determine_initial_infection_states_world(input_dir, start_date, dark_figure);
        restart_timer(timer, "time for determine_initial_infection_states_world");
        auto vacc_map = prepare_vaccination_state(mio::offset_date_by_days(start_date, (int)tmax.days()),
                                                  (input_dir / "pydata/Germany/vacc_county_ageinf_ma7.json").string());
        restart_timer(timer, "time for vaccinaiton state");
        // Loop over a number of runs
        for (size_t run_idx = start_run_idx; run_idx < end_run_idx; run_idx++) {

            auto world = mio::abm::World(num_age_groupss);
            auto run_rng_counter =
                mio::rng_totalsequence_counter<uint64_t>(static_cast<uint32_t>(run_idx), mio::Counter<uint32_t>(0));
            rng.set_counter(run_rng_counter);
            world.get_rng() = rng;

            create_sampled_world(world, input_dir, t0, max_num_persons, start_date, perc_easter_event,
                                 initial_infection, vacc_map, (uint32_t)run_idx);

            restart_timer(timer, "time taken for create sampled world");
            auto sim = mio::abm::Simulation(t0, std::move(world));

            //Logger
            mio::History<mio::abm::TimeSeriesWriter, LogInfectionPerLocationTypePerAgeGroup>
                historyInfectionPerLocationTypePerAgeGroup{
                    Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups())};
            mio::History<mio::abm::TimeSeriesWriter, LogInfectionStatePerAgeGroup> historyInfectionStatePerAgeGroup{
                Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups())};
            mio::History<mio::abm::TimeSeriesWriter, LogTestPerLocationTypePerAgeGroup>
                historyTestPerLocationTypePerAgeGroup{
                    Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups())};
            mio::History<mio::abm::TimeSeriesWriter, LogPositiveTestPerLocationTypePerAgeGroup>
                historyPositiveTestPerLocationTypePerAgeGroup{
                    Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups())};
            mio::History<mio::abm::TimeSeriesWriter, LogCumulativeDetectedInfectionsPerAgeGroup>
                historyCumulativeDetectedInfectionsPerAgeGroup{
                    Eigen::Index(sim.get_world().parameters.get_num_groups())};
            mio::History<mio::abm::TimeSeriesWriter, LogEstimatedReproductionNumber> historyEstimatedReproductionNumber{
                Eigen::Index(1)};

            // / NPIS//
            if (npis_on) {

                const auto location_it = sim.get_world().get_locations();
                //Testing strategies
                add_npi_testing_strategies_to_world(sim, tmax, testing_probability_sympt, ratio_asympt_to_sympt,
                                                    lockdown_test_prob, after_lockdown_test_prob);

                // 2. Mask schemes for all locations
                // First set all locations to have mask usage, we need ffp2 masks
                for (auto& location : location_it) {
                    location.set_required_mask(mio::abm::MaskType::FFP2);
                    if (location.get_type() == mio::abm::LocationType::Home) {
                        location.set_npi_active(false);
                    }
                    else {
                        location.set_npi_active(true);
                    }
                }

                // 3. Dampings everywhere except home
                for (auto& location : location_it) {
                    if (location.get_type() == mio::abm::LocationType::School) {
                        location.add_damping(mio::abm::TimePoint(mio::abm::days(0).seconds()), 0.5); // from 2021-03-01
                        location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= (0.5);
                        location.add_damping(mio::abm::TimePoint(mio::abm::days(date_of_lockdown).seconds()),
                                             0.00); // from 2021-03-15
                        location.add_damping(mio::abm::TimePoint(mio::abm::days(end_date_of_lockdown).seconds()),
                                             0.5); // from 2021-04-12 till 2021-05-30
                    }
                    if (location.get_type() == mio::abm::LocationType::BasicsShop) {
                        location.add_damping(mio::abm::TimePoint(mio::abm::days(0).seconds()), 0.8); // from 2021-03-15
                        // location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= (0.5);
                        location.add_damping(mio::abm::TimePoint(mio::abm::days(date_of_lockdown).seconds()),
                                             damping_community_lockdown); // from 2021-03-15
                        location.add_damping(mio::abm::TimePoint(mio::abm::days(end_date_of_lockdown).seconds()),
                                             0.8); // from 2021-03-15
                    }
                    if (location.get_type() == mio::abm::LocationType::SocialEvent) {
                        location.add_damping(mio::abm::TimePoint(mio::abm::days(0).seconds()), 0.8); // from 2021-03-15
                        // location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= (0.5);
                        location.add_damping(mio::abm::TimePoint(mio::abm::days(date_of_lockdown).seconds()),
                                             damping_community_lockdown); // from 2021-03-15
                        location.add_damping(mio::abm::TimePoint(mio::abm::days(end_date_of_lockdown).seconds()),
                                             0.8); // from 2021-03-15
                    }
                    if (location.get_type() == mio::abm::LocationType::Work) {
                        location.add_damping(mio::abm::TimePoint(mio::abm::days(0).seconds()), 0.75); // from 2021-03-15
                        // location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= (0.75 * 0.75);
                        location.add_damping(mio::abm::TimePoint(mio::abm::days(date_of_lockdown).seconds()),
                                             0.7); // from 2021-03-15
                        location.add_damping(mio::abm::TimePoint(mio::abm::days(end_date_of_lockdown).seconds()),
                                             0.75); // from 2021-03-15
                    }
                }

                sim.get_world().parameters.get<mio::abm::InfectionRateFromViralShed>() = viral_shedding_rate;
                sim.get_world().parameters.get<mio::abm::MaskProtection>()             = masks;
                sim.get_world().parameters.get<mio::abm::QuarantineEffectiveness>()    = quarantine_effectiveness;
                sim.get_world().parameters.get<mio::abm::QuarantineDuration>()         = quarantine_duration;

                sim.advance(mio::abm::TimePoint(mio::abm::days(date_of_lockdown).seconds()),
                            historyInfectionPerLocationTypePerAgeGroup, historyInfectionStatePerAgeGroup,
                            historyTestPerLocationTypePerAgeGroup, historyPositiveTestPerLocationTypePerAgeGroup,
                            historyCumulativeDetectedInfectionsPerAgeGroup, historyEstimatedReproductionNumber);
                sim.get_world().parameters.get<mio::abm::InfectionRateFromViralShed>() =
                    viral_shedding_rate * seasonality_april;

                for (auto& location : location_it) {
                    if (location.get_type() != mio::abm::LocationType::Home) {
                        location.get_infection_parameters().get<mio::abm::ContactRates>().array() *=
                            contact_red_lockdown;
                    }
                }

                sim.advance(mio::abm::TimePoint(mio::abm::days(end_date_of_lockdown).seconds()),
                            historyInfectionPerLocationTypePerAgeGroup, historyInfectionStatePerAgeGroup,
                            historyTestPerLocationTypePerAgeGroup, historyPositiveTestPerLocationTypePerAgeGroup,
                            historyCumulativeDetectedInfectionsPerAgeGroup, historyEstimatedReproductionNumber);

                for (auto& location : location_it) {
                    if (location.get_type() != mio::abm::LocationType::Home) {
                        location.get_infection_parameters().get<mio::abm::ContactRates>().array() *=
                            ((1 / contact_red_lockdown) * after_lockdown_contact_reduction);
                    }
                }
                sim.get_world().parameters.get<mio::abm::InfectionRateFromViralShed>() =
                    viral_shedding_rate * seasonality_may;

                sim.advance(mio::abm::TimePoint(tmax.seconds()), historyInfectionPerLocationTypePerAgeGroup,
                            historyInfectionStatePerAgeGroup, historyTestPerLocationTypePerAgeGroup,
                            historyPositiveTestPerLocationTypePerAgeGroup,
                            historyCumulativeDetectedInfectionsPerAgeGroup, historyEstimatedReproductionNumber);
            }
            else {
                sim.get_world().parameters.get<mio::abm::InfectionRateFromViralShed>() = viral_shedding_rate;
                sim.advance(mio::abm::TimePoint(mio::abm::days(20).seconds()),
                            historyInfectionPerLocationTypePerAgeGroup, historyInfectionStatePerAgeGroup);
            }
            ////Advance till here
            // Stop the clock after sim.advance and calculate the duration
            restart_timer(timer, "time taken for simulation end");
            // TODO: update result of the simulation to be a vector of location result.
            auto temp_sim_infection_per_loc_type_per_age_group = std::vector<mio::TimeSeries<ScalarType>>{
                std::get<0>(historyInfectionPerLocationTypePerAgeGroup.get_log())};
            auto temp_sim_infection_state_per_age_group =
                std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionStatePerAgeGroup.get_log())};
            auto temp_sim_test_per_loc_type_per_age_group =
                std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyTestPerLocationTypePerAgeGroup.get_log())};
            auto temp_sim_positive_test_per_loc_type_per_age_group = std::vector<mio::TimeSeries<ScalarType>>{
                std::get<0>(historyPositiveTestPerLocationTypePerAgeGroup.get_log())};
            auto temp_sim_cumulative_detected_infections_per_age_group = std::vector<mio::TimeSeries<ScalarType>>{
                std::get<0>(historyCumulativeDetectedInfectionsPerAgeGroup.get_log())};
            auto temp_sim_estimated_new_detected_infections = std::vector<mio::TimeSeries<ScalarType>>{
                get_new_detected_infections(temp_sim_cumulative_detected_infections_per_age_group[0])};
            auto temp_sim_estimated_reproduction_number =
                std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyEstimatedReproductionNumber.get_log())};

            // Push result of the simulation back to the result vector
            ensemble_infection_per_loc_type_per_age_group.emplace_back(temp_sim_infection_per_loc_type_per_age_group);
            ensemble_infection_state_per_age_group.emplace_back(temp_sim_infection_state_per_age_group);
            ensemble_test_per_loc_type_per_age_group.emplace_back(temp_sim_test_per_loc_type_per_age_group);
            ensemble_positive_test_per_loc_type_per_age_group.emplace_back(
                temp_sim_positive_test_per_loc_type_per_age_group);
            ensemble_cumulative_detected_infections.emplace_back(temp_sim_cumulative_detected_infections_per_age_group);
            ensemble_estimated_new_detected_infections.emplace_back(temp_sim_estimated_new_detected_infections);
            ensemble_estimated_reproduction_number.emplace_back(temp_sim_estimated_reproduction_number);

            std::cout << "Run " << run_idx + 1 << " of " << num_runs << " finished." << std::endl;

            // calculate RMSE
            int number_of_days_for_rmse = 0;
            if (npis_on) {
                number_of_days_for_rmse = 90;
            }
            else {
                number_of_days_for_rmse = 20;
            }
            auto rmse = calculate_rmse_from_results(input_dir, temp_sim_infection_state_per_age_group[0],
                                                    temp_sim_cumulative_detected_infections_per_age_group[0],
                                                    number_of_days_for_rmse, start_date);

            rmse_results_per_grid_point.at(par_i) += rmse;

            std::cout << "RMSE: " << rmse << std::endl;

            //HACK since // gather_results(rank, num_procs, num_runs, ensemble_params);
            //for now this doesnt work, but we can still save the results of the last world since the
            //parameters are the same for each run
            if (rank == 0 && run_idx == end_run_idx - 1) {
                for (size_t i = 0; i < num_runs; i++) {
                    ensemble_params.emplace_back(std::vector<mio::abm::World>{sim.get_world()});
                }
            }
        }

        rmse_results_per_grid_point.at(par_i) /= num_runs;
        printf("RMSE: %f\n, par_i: %ld\n", rmse_results_per_grid_point.at(par_i), par_i);

#ifdef MEMILIO_ENABLE_MPI

        //gather results
        auto final_ensemble_infection_per_loc_type_per_age_group =
            gather_results(rank, num_procs, num_runs, ensemble_infection_per_loc_type_per_age_group);
        auto final_ensemble_infection_state_per_age_group =
            gather_results(rank, num_procs, num_runs, ensemble_infection_state_per_age_group);
        auto final_ensemble_test_per_loc_type_per_age_group =
            gather_results(rank, num_procs, num_runs, ensemble_test_per_loc_type_per_age_group);
        auto final_ensemble_positive_test_per_loc_type_per_age_group =
            gather_results(rank, num_procs, num_runs, ensemble_positive_test_per_loc_type_per_age_group);
        auto final_ensemble_cumulative_detected_infections =
            gather_results(rank, num_procs, num_runs, ensemble_cumulative_detected_infections);
        auto final_ensemble_estimated_new_detected_infections =
            gather_results(rank, num_procs, num_runs, ensemble_estimated_new_detected_infections);
        auto final_ensemble_estimated_reproduction_number =
            gather_results(rank, num_procs, num_runs, ensemble_estimated_reproduction_number);

        if (rank == 0) {
            BOOST_OUTCOME_TRY(save_results(
                final_ensemble_infection_per_loc_type_per_age_group, ensemble_params, {0},
                result_dir / "infection_per_location_type_per_age_group" / std::to_string(par_i), save_single_runs));
            BOOST_OUTCOME_TRY(save_results(final_ensemble_infection_state_per_age_group, ensemble_params, {0},
                                           result_dir / "infection_state_per_age_group" / std::to_string(par_i),
                                           save_single_runs));
            BOOST_OUTCOME_TRY(save_results(final_ensemble_test_per_loc_type_per_age_group, ensemble_params, {0},
                                           result_dir / "test_per_location_type_per_age_group" / std::to_string(par_i),
                                           save_single_runs));
            BOOST_OUTCOME_TRY(
                save_results(final_ensemble_positive_test_per_loc_type_per_age_group, ensemble_params, {0},
                             result_dir / "positive_test_per_location_type_per_age_group" / std::to_string(par_i),
                             save_single_runs));
            BOOST_OUTCOME_TRY(save_results(final_ensemble_cumulative_detected_infections, ensemble_params, {0},
                                           result_dir / "cumulative_detected_infections" / std::to_string(par_i),
                                           save_single_runs));
            BOOST_OUTCOME_TRY(save_results(final_ensemble_estimated_new_detected_infections, ensemble_params, {0},
                                           result_dir / "new_detected_infections" / std::to_string(par_i),
                                           save_single_runs));
            BOOST_OUTCOME_TRY(save_results(final_ensemble_estimated_reproduction_number, ensemble_params, {0},
                                           result_dir / "estimated_reproduction_number" / std::to_string(par_i),
                                           save_single_runs, true, true));
        }
#else
        BOOST_OUTCOME_TRY(save_results(ensemble_infection_per_loc_type_per_age_group, ensemble_params, {0},
                                       result_dir / "infection_per_location_type_per_age_group" / std::to_string(par_i),
                                       save_single_runs));
        BOOST_OUTCOME_TRY(save_results(ensemble_infection_state_per_age_group, ensemble_params, {0},
                                       result_dir / "infection_state_per_age_group" / std::to_string(par_i),
                                       save_single_runs));
        BOOST_OUTCOME_TRY(save_results(ensemble_test_per_loc_type_per_age_group, ensemble_params, {0},
                                       result_dir / "test_per_location_type_per_age_group" / std::to_string(par_i),
                                       save_single_runs));
        BOOST_OUTCOME_TRY(save_results(
            ensemble_positive_test_per_loc_type_per_age_group, ensemble_params, {0},
            result_dir / "positive_test_per_location_type_per_age_group" / std::to_string(par_i), save_single_runs));
        BOOST_OUTCOME_TRY(save_results(ensemble_cumulative_detected_infections, ensemble_params, {0},
                                       result_dir / "cumulative_detected_infections" / std::to_string(par_i),
                                       save_single_runs));
        BOOST_OUTCOME_TRY(save_results(ensemble_estimated_reproduction_number, ensemble_params, {0},
                                       result_dir / "estimated_reproduction_number" / std::to_string(par_i),
                                       save_single_runs, 1));
        BOOST_OUTCOME_TRY(save_results(ensemble_estimated_new_detected_infections, ensemble_params, {0},
                                       result_dir / "new_detected_infections" / std::to_string(par_i), save_single_runs,
                                       1));

#endif
        restart_timer(timer, "time taken for data gathering and saving results");

        printf("done.\n");
    }
    return mio::success();
}

// From: https://en.cppreference.com/w/cpp/chrono/c/strftime
const std::string currentDateTime()
{
    // Example of the very popular RFC 3339 format UTC time
    std::time_t time = std::time({});
    char timeString[std::size("yyyy-mm-dddddddd")];
    std::strftime(std::data(timeString), std::size(timeString), "%F%H%M%S", std::gmtime(&time));
    return timeString;
}

mio::IOResult<bool> create_result_folders(std::string const& result_dir, int n_params = 0, bool grid_search = false)
{
    std::string inf_p_loc_t_p_ag      = result_dir + "/infection_per_location_type_per_age_group/";
    std::string inf_state_p_ag        = result_dir + "/infection_state_per_age_group/";
    std::string test_p_loc_t_p_ag     = result_dir + "/test_per_location_type_per_age_group/";
    std::string pos_test_p_loc_t_p_ag = result_dir + "/positive_test_per_location_type_per_age_group/";
    std::string cum_det_inf           = result_dir + "/cumulative_detected_infections/";
    std::string est_rep_num           = result_dir + "/estimated_reproduction_number/";
    std::string new_det_inf           = result_dir + "/new_detected_infections/";

    BOOST_OUTCOME_TRY(mio::create_directory(result_dir));
    BOOST_OUTCOME_TRY(mio::create_directory(inf_p_loc_t_p_ag));
    BOOST_OUTCOME_TRY(mio::create_directory(inf_state_p_ag));
    BOOST_OUTCOME_TRY(mio::create_directory(test_p_loc_t_p_ag));
    BOOST_OUTCOME_TRY(mio::create_directory(pos_test_p_loc_t_p_ag));
    BOOST_OUTCOME_TRY(mio::create_directory(cum_det_inf));
    BOOST_OUTCOME_TRY(mio::create_directory(new_det_inf));
    BOOST_OUTCOME_TRY(mio::create_directory(est_rep_num));
    if (n_params > 0) {
        // we create n_param folders in each of the subfolders
        for (int i = 0; i < n_params; i++) {
            BOOST_OUTCOME_TRY(mio::create_directory(inf_p_loc_t_p_ag + "/" + std::to_string(i)));
            BOOST_OUTCOME_TRY(mio::create_directory(inf_state_p_ag + "/" + std::to_string(i)));
            BOOST_OUTCOME_TRY(mio::create_directory(test_p_loc_t_p_ag + "/" + std::to_string(i)));
            BOOST_OUTCOME_TRY(mio::create_directory(pos_test_p_loc_t_p_ag + "/" + std::to_string(i)));
            BOOST_OUTCOME_TRY(mio::create_directory(cum_det_inf + "/" + std::to_string(i)));
            BOOST_OUTCOME_TRY(mio::create_directory(est_rep_num + "/" + std::to_string(i)));
            BOOST_OUTCOME_TRY(mio::create_directory(new_det_inf + "/" + std::to_string(i)));
        }
    }

    if (grid_search) {
        std::string grid_search_dir = result_dir + "/grid_search/";
        BOOST_OUTCOME_TRY(mio::create_directory(grid_search_dir));
    }
    return mio::success();
}

mio::IOResult<bool> copy_result_folder(std::string const& from_dir, std::string const& to_dir)
{
    BOOST_OUTCOME_TRY(mio::create_directory(to_dir));
    fs::copy(from_dir, to_dir, fs::copy_options::overwrite_existing | fs::copy_options::recursive);
    return mio::success();
}

int main(int argc, char** argv)
{

    mio::set_log_level(mio::LogLevel::err);
    auto start = std::chrono::system_clock::now();

    int num_procs, rank;
#ifdef MEMILIO_ENABLE_MPI
    mio::mpi::init();
    MPI_Comm_size(mio::mpi::get_world(), &num_procs);
    MPI_Comm_rank(mio::mpi::get_world(), &rank);
#else
    num_procs = 1;
    rank      = 0;
#endif
    // we need one seed
    std::initializer_list<uint32_t> seeds = {14159265u, 35897932u};
    auto rng                              = mio::RandomNumberGenerator();
    rng.seed(seeds);
    rng.synchronize();

    std::string input_dir = "/p/project1/loki/memilio/memilio/data";
    // std::string input_dir = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data";
    // std::string input_dir = "/Users/david/Documents/HZI/memilio/data";
    // std::string input_dir       = "C:/Users/korf_sa/Documents/rep/data";

#ifdef MEMILIO_ENABLE_MPI
    // we need to send every rank the same folder
    std::string result_dir;
    if (rank == 0)
        result_dir = input_dir + "/results_" + currentDateTime();
    int line_size = result_dir.size();
    MPI_Bcast(&line_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0)
        result_dir.resize(line_size);
    MPI_Bcast(const_cast<char*>(result_dir.data()), line_size, MPI_CHAR, 0, MPI_COMM_WORLD);
#else
    std::string result_dir = input_dir + "/results_" + currentDateTime();
#endif

    size_t num_runs;
    bool run_grid_search = false;

    if (argc == 2) {
        num_runs        = std::stoi(argv[1]);
        run_grid_search = false;
        printf("Running with number of runs = %d.\n", (int)num_runs);
        printf("Saving results to \"%s\".\n", result_dir.c_str());
    }
    else if (argc == 3) {
        num_runs        = 11;
        run_grid_search = true;
        printf("running with grid search\n");
        printf("Running with number of runs %d.\n", (int)num_runs);
        printf("Saving results to \"%s\".\n", result_dir.c_str());
    }
    else {
        num_runs = 1;
        printf("Running with number of runs = %d.\n", (int)num_runs);
        printf("Saving results to \"%s\".\n", result_dir.c_str());
    }

    timer = TIME_NOW;

    if (run_grid_search) {
        // grid search for parameters:
        // 1: Viral Shed
        // 2: Dark figure
        // 3: Contact reduction lockdown
        // 4: testing prob symptomatic
        // 5: odds ratio for asymptomatic to symptomatic

        // std::vector<std::pair<double, double>> grid_boundaries = {
        //     {1.4, 2.0}, {2.5, 5.5}, {0.35, 0.85}, {0.01, 0.045}, {2, 15}};
        // std::vector<int> points_per_dim = {6, 6, 6, 6, 6};

        std::vector<double> grid_boundaries = {1.52, 4.3, 0.85, 0.024, 4.6};
        // std::vector<double> grid_boundaries = {1.52, 4.3, 0.75, 0.024, 4.6};
        // std::vector<double> grid_boundaries = {1.64, 4.3, 0.55, 0.038, 12.4};
        std::vector<int> points_per_dim = {6, 6, 6, 6, 6};

        // std::vector<std::pair<double, double>> grid_boundaries = {
        // {1.79, 1.83}, {3.28, 3.29}, {0.52, 0.56}, {0.03, 0.04}, {1.0, 15.0}};
        // std::vector<int> points_per_dim = {4, 4, 4, 9, 9};

        auto grid = grid_points(grid_boundaries, points_per_dim);
        if (rank == 0) {
            auto created = create_result_folders(result_dir, 0, run_grid_search);
            if (!created) {
                std::cout << created.error().formatted_message();
                return created.error().code().value();
            }
        }
        auto result = run_with_grid_search(input_dir, result_dir, num_runs, grid, rng);
    }
    else {
        // std::vector<std::vector<double>> parameters = {{1.596}, {4.171}, {0.7275}, {0.02472}, {4.83}, {10.0}, {0.5}};
        // std::vector<std::vector<double>> parameters = {{1.52}, {4.3}, {0.85}, {0.024}, {4.6}, {10.0}, {0.5}};

        // std::vector<std::vector<double>> parameters = {{1.596}, {4.171}, {0.7125}, {0.02472}, {4.83}, {10.0}, {0.5}};

        // std::vector<std::vector<double>> parameters = {{1.596}, {4.171}, {1.0}, {0.02472, 0.07416, 0.1236},
        //                                                {4.83},  {10.0},  {0.5}};

        // std::vector<std::vector<double>> parameters = {
        //     {1.596}, {4.171}, {0.7125}, {0.012, 0.024, 0.036, 0.048, 0.06}, {2.0, 5.0, 8.0, 11.0, 14.0}, {10.0}, {0.5}};
        // std::vector<std::vector<double>> parameters = {
        //     {1.596}, {4.171}, {0.7125}, {0.012, 0.024, 0.036, 0.048, 0.06}, {4.83}, {2, 5, 8, 11, 14}, {0.5}};
        // std::vector<std::vector<double>> parameters = {
        //     {1.596}, {4.171}, {0.7125}, {0.02472}, {4.83}, {2, 5, 8, 11, 14}, {0.0, 0.25, 0.5, 0.75, 1.0}};
        std::vector<std::vector<double>> parameters = {
            {1.596}, {4.171}, {0.7125}, {0.012, 0.024, 0.036, 0.048, 0.06}, {4.83}, {10.0}, {0.5},{0.2, 0.225, 0.25, 0.275, 0.3}};

        auto every_combination = every_combination_of_parameters(parameters);
        if (rank == 0) {
            auto created = create_result_folders(result_dir, every_combination.size(), run_grid_search);
            if (!created) {
                std::cout << created.error().formatted_message();
                return created.error().code().value();
            }
        }
        auto result = run(input_dir, result_dir, num_runs, every_combination, rng);
    }

    // copy results into a fixed name folder to have easier access
    std::string last_run_dir = input_dir + "/results_last_run";
    if (rank == 0)
        auto copied = copy_result_folder(result_dir, last_run_dir);

    mio::mpi::finalize();

    auto end                                      = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    //print time taken
    std::cout << "Time taken: " << elapsed_seconds.count() << "s\n";

    return 0;
}