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

mio::CustomIndexArray<double, mio::AgeGroup, mio::osecir::InfectionState> initial_infection_distribution{
    {mio::AgeGroup(num_age_groupss), mio::osecir::InfectionState::Count}, 0.5};

std::map<mio::Date, std::vector<std::pair<uint32_t, uint32_t>>> vacc_map;

/**
 * Determine initial distribution of infection states.
*/
void determine_initial_infection_states_world(const fs::path& input_dir, const mio::Date date, double sclaling_infected)
{
    // estimate intial population by ODE compartiments
    auto initial_graph                     = get_graph(date, 1, input_dir, sclaling_infected);
    const size_t braunschweig_id           = 16; // Braunschweig has ID 16
    auto braunschweig_node                 = initial_graph.value()[braunschweig_id];
    initial_infection_distribution.array() = braunschweig_node.populations.array().cast<double>();
}

/**
 * Assign an infection state to each person according to real world data read in through the ODE secir model.
 * Infections are set with probabilities computed by the values in the rows in initial_infection_distribution.
 */
void assign_infection_state_prob(mio::abm::World& world, mio::abm::TimePoint t)
{
    // convert initial population to ABM initial infections
    for (auto& person : world.get_persons()) {
        if (person.get_should_be_logged()) {
            auto rng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);

            auto infection_state = mio::osecir::InfectionState(mio::DiscreteDistribution<size_t>::get_instance()(
                rng, initial_infection_distribution.slice(person.get_age()).as_array().array()));

            if (infection_state != mio::osecir::InfectionState::Susceptible) {
                person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Alpha, person.get_age(),
                                                             world.parameters, t,
                                                             infection_state_map.at(infection_state)));
            }
        }
    }
}

/**
 * Assign an infection state to each person according to real world data read in through the ODE secir model.
 * Infections are set with the rounded values in the rows in initial_infection_distribution.
 * Only works if enough persons in the all age groups exist.
 * The number of agents in the model should fit to the sum of the rows in initial_infection_distribution,
 * otherwise many agents will be susceptible.
 */
void assign_infection_state(mio::abm::World& world, mio::abm::TimePoint t)
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

void prepare_vaccination_state(mio::Date simulation_end, const std::string& filename)
{
    // for saving previous day of vaccination
    std::vector<std::pair<uint32_t, uint32_t>> vacc_vector_prev(num_age_groupss);
    //inizialize the vector with 0
    for (size_t i = 0; i < num_age_groupss; ++i) {
        vacc_vector_prev[i] = std::make_pair(0, 0);
    }

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
}

/**
 * @brief assign an vaccination state to each person according to real world data read in through the ODE secir model.
 * 
 * @param input 
 * @return int 
 */
void assign_vaccination_state(mio::abm::World& world, mio::Date simulation_beginning)
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
            for (uint32_t i = 0; i < vacc_entry.second[age].second; ++i) {
                if (vaccinated_persons[age].size() == 0) {
                    // mio::log_error("Not enough vaccinated people to vacc 2nd time! ");
                }
                else {
                    // select random already vaccinated person and assign Vaccination
                    uint32_t id_rnd = vaccinated_persons[age][mio::UniformIntDistribution<size_t>::get_instance()(
                        world.get_rng(), 0U, vaccinated_persons[age].size() - 1)];
                    mio::abm::Person& person = world.get_person(id_rnd);
                    auto timePoint           = mio::abm::TimePoint(
                        mio::get_offset_in_days(vacc_entry.first, simulation_beginning) * 24 * 60 * 60);
                    person.add_new_vaccination(
                        mio::abm::Vaccination(mio::abm::ExposureType::GenericVaccine, timePoint));
                    vaccinated_persons[age].erase(
                        std::remove(vaccinated_persons[age].begin(), vaccinated_persons[age].end(), id_rnd),
                        vaccinated_persons[age].end());
                }
            }
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

void create_world_from_data(mio::abm::World& world, const std::string& filename, const int max_number_persons)
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
            auto& person = world.add_person(home, determine_age_group(age));
            person.set_mask_preferences({0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.2, 0.2});
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

    //Set testing parameters
    auto pcr_test_values                                          = mio::abm::TestParameters{0.9, 0.99};
    auto antigen_test_values                                      = mio::abm::TestParameters{0.69, 0.95};
    auto generic_test_values                                      = mio::abm::TestParameters{0.7, 0.95};
    params.get<mio::abm::TestData>()[mio::abm::TestType::PCR]     = pcr_test_values;
    params.get<mio::abm::TestData>()[mio::abm::TestType::Antigen] = antigen_test_values;
    params.get<mio::abm::TestData>()[mio::abm::TestType::Generic] = generic_test_values;

    // Set percentage parameters
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_0_to_4}]   = 0.5;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_5_to_14}]  = 0.55;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_15_to_34}] = 0.6;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_35_to_59}] = 0.69;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_60_to_79}] = 0.825;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_80_plus}]  = 0.9;

    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_0_to_4}]   = 0.0180;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_5_to_14}]  = 0.0237;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_15_to_34}] = 0.0373;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_35_to_59}] = 0.07;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_60_to_79}] = 0.17;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_80_plus}]  = 0.2374;

    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_0_to_4}]   = 0.1;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_5_to_14}]  = 0.11;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_15_to_34}] = 0.12;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_35_to_59}] = 0.13;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_60_to_79}] = 0.33;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_80_plus}]  = 0.62;

    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_0_to_4}]   = 0.12;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_5_to_14}]  = 0.13;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_15_to_34}] = 0.15;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_35_to_59}] = 0.28;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_60_to_79}] = 0.43;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_80_plus}]  = 0.48;

    // Set infection parameters

    // params.get<mio::abm::InfectionRateFromViralShed>()[{mio::abm::VirusVariant::Alpha}] = 3.5;

    // Set protection level from high viral load. Information based on: https://doi.org/10.1093/cid/ciaa886
    params.get<mio::abm::HighViralLoadProtectionFactor>() = [](ScalarType days) -> ScalarType {
        return mio::linear_interpolation_of_data_set<ScalarType, ScalarType>(
            {{0, 0.00863}, {1, 0.00969}, {7, 0.0029}, {10, 0.002}, {14, 0.0014}, {21, 0}}, days);
    };

    // Set protection level against an severe infection. Information based on: https://doi.org/10.1093/cid/ciaa886
    params.get<mio::abm::SeverityProtectionFactor>() = [](ScalarType days) -> ScalarType {
        return mio::linear_interpolation_of_data_set<ScalarType, ScalarType>({{0, 0.01967},
                                                                              {30, 0.01975},
                                                                              {60, 0.01977},
                                                                              {90, 0.01974},
                                                                              {120, 0.01963},
                                                                              {150, 0.01947},
                                                                              {180, 0.0193},
                                                                              {210, 0.01929},
                                                                              {240, 0.01923},
                                                                              {270, 0.01908},
                                                                              {300, 0.01893},
                                                                              {330, 0.01887},
                                                                              {360, 0.01887},
                                                                              {360, 0.015}},
                                                                             days);
    };

    //Set other parameters
    params.get<mio::abm::MaskProtection>()           = 0.55; //all masks have a 0.66 protection factor for now
    params.get<mio::abm::AerosolTransmissionRates>() = 0.0;
}

// set location specific parameters
void set_local_parameters(mio::abm::World& world)
{
    const int n_age_groups = (int)world.parameters.get_num_groups();

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
        // # we assume that a 20:4:4:1:1 ratio of contacts is made at home:school:work:SocialEv:basicshopping
        // First line: Sclaing according to which Timespan the contatcs have to be made
        // Second Line: Normalizing to 24 hours
        // Third Line: Scaling according intensity of contacts
        switch (loc.get_type()) {
        case mio::abm::LocationType::Home:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_home;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1.5;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.66; // 2/5z
            break;
        case mio::abm::LocationType::School:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_school;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 7.5;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1;
            break;
        case mio::abm::LocationType::Work:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_work;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 7.5;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.3; // 2/5z
            break;
        case mio::abm::LocationType::SocialEvent:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 15;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 5;
            break;
        case mio::abm::LocationType::BasicsShop:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 15;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.2; // 1/7
            break;
        case mio::abm::LocationType::Event:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 4;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 4; // 1/7
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
void add_testing_strategies(mio::abm::World& world, bool symptomatic, bool social_event)
{
    if (symptomatic) {
        std::cout << "Adding symptomatic testing strategy" << std::endl;
        auto testing_min_time_symptomatic = mio::abm::days(7);
        auto probability_symptomatic      = 1.0;
        auto start_date_test_symptomatic  = mio::abm::TimePoint(mio::abm::days(0).seconds()); // 2021-04-12
        auto end_date_test_symptomatic    = mio::abm::TimePoint(mio::abm::days(90).seconds()); // 2021-05-30
        auto test_type_symptomatic        = mio::abm::TestType::Antigen; // Antigen test
        auto test_parameters = world.parameters.get<mio::abm::TestData>()[test_type_symptomatic]; // Test parameters
        auto testing_criteria_symptomatic = mio::abm::TestingCriteria();
        testing_criteria_symptomatic.add_infection_state(mio::abm::InfectionState::InfectedSymptoms);
        auto testing_scheme_symptomatic = mio::abm::TestingScheme(
            testing_criteria_symptomatic, testing_min_time_symptomatic, start_date_test_symptomatic,
            end_date_test_symptomatic, test_parameters, probability_symptomatic);
        world.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Home, testing_scheme_symptomatic);
    }
    if (social_event) {
        auto testing_min_time_socev = mio::abm::days(7);
        auto probability_socev      = 0.8;
        auto start_date_test_socev  = mio::abm::TimePoint(mio::abm::days(0).seconds());
        auto end_date_test_socev    = mio::abm::TimePoint(mio::abm::days(90).seconds());
        auto test_type_socev        = mio::abm::TestType::Antigen; // Antigen test
        auto test_parameters        = world.parameters.get<mio::abm::TestData>()[test_type_socev]; // Test parameters
        auto testing_criteria_socev = mio::abm::TestingCriteria();
        auto testing_scheme_socev =
            mio::abm::TestingScheme(testing_criteria_socev, testing_min_time_socev, start_date_test_socev,
                                    end_date_test_socev, test_parameters, probability_socev);
        world.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::SocialEvent, testing_scheme_socev);
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
    auto date_need = mio::offset_date_by_days(start_date, -19);
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

double calculate_rmse_from_results(const fs::path& data_dir, mio::TimeSeries<ScalarType> sim_inf_states,
                                   int max_num_days, mio::Date start_date)
{
    // We need to read in the results from the results directory
    auto real_data_dead_path =
        mio::path_join((data_dir / "pydata" / "Germany").string(), "cases_all_county_age_ma7.json");
    auto real_data_icu_path = mio::path_join((data_dir / "pydata" / "Germany").string(), "county_divi_ma7.json");
    auto divi_data          = mio::read_divi_data(real_data_icu_path);
    auto death_data         = mio::read_confirmed_cases_data(real_data_dead_path);
    // read in the real data
    auto real_data_dead_vec = read_in_deaths(death_data.value(), start_date, max_num_days);
    auto real_data_icu_vec  = read_in_icu(divi_data.value(), start_date, max_num_days);

    // Simulated data
    std::vector<int> sim_data_vec_icu;
    std::vector<int> sim_data_vec_dead;
    for (int i = 0; i < max_num_days; i++) {
        int number_of_persons_in_icu = 0;
        int number_of_persons_dead   = 0;
        for (size_t j = 0; j < (size_t)num_age_groupss; j++) {
            auto index_icu = (((size_t)(mio::abm::InfectionState::Count)) * (j)) +
                             ((uint32_t)(mio::abm::InfectionState::InfectedCritical));
            auto index_dead =
                (((size_t)(mio::abm::InfectionState::Count)) * (j)) + ((uint32_t)(mio::abm::InfectionState::Dead));
            number_of_persons_in_icu += sim_inf_states[i * 24][index_icu];
            number_of_persons_dead += sim_inf_states[i * 24][index_dead];
        }
        sim_data_vec_icu.push_back(number_of_persons_in_icu);
        sim_data_vec_dead.push_back(number_of_persons_dead);
    }

    // now we calculate the RMSE
    double rmse_dead = 0;
    double rmse_icu  = 0;
    for (size_t i = 0; i < real_data_dead_vec.size(); i++) {
        rmse_dead += pow(real_data_dead_vec[i] - sim_data_vec_dead[i], 2);
        rmse_icu += pow((int)(sim_data_vec_icu[i] * 0.5) - real_data_icu_vec[i], 2);
    }
    rmse_dead = sqrt(rmse_dead / real_data_dead_vec.size());
    rmse_icu  = sqrt(rmse_icu / real_data_icu_vec.size());

    return rmse_dead + rmse_icu;
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
        for (int j = 0; j < number_of_points.at(i) - 1; j++) {
            temp.push_back(parameter_boundaries[i].first + j * step);
        }
        grid.push_back(temp);
    }
    return grid;
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
    std::vector<int> counter_per_dimension(grid.size(), 0);
    for (int i = 0; i < number_of_points; i++) {
        std::vector<double> temp;
        for (size_t j = 0; j < grid.size(); j++) {
            temp.push_back(grid[j][counter_per_dimension[j]]);
        }
        grid_search.push_back(temp);
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

    //shuffle this vector
    std::random_device rd;
    std::mt19937 g(rd());
    for (size_t i = 0; i < person_per_age_group.size(); i++) {
        std::shuffle(person_per_age_group[i].begin(), person_per_age_group[i].end(), g);
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

            if (j == 2) {
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
void create_sampled_world(mio::abm::World& world, const fs::path& input_dir, const mio::abm::TimePoint& t0,
                          int max_num_persons, mio::Date start_date_sim, double perc_easter_event)
{
    //Set global infection parameters (similar to infection parameters in SECIR model) and initialize the world

    set_parameters(world.parameters);

    restart_timer(timer, "time taken for setting up parameters and local parameters");

    // Create the world object from statistical data.
    create_world_from_data(world, (input_dir / "mobility/braunschweig_result_ffa8_modified2.csv").generic_string(),
                           max_num_persons);
    world.use_migration_rules(false);
    restart_timer(timer, "time taken for braunschweig trip input");

    create_easter_social_event(world, perc_easter_event);
    restart_timer(timer, "time taken for setting up easter social ebent");

    // Assign an infection state to each person.
    assign_infection_state(world, t0);
    restart_timer(timer, "time taken for assigning infection state");

    // Assign vaccination status to each person.
    assign_vaccination_state(world, start_date_sim);
    restart_timer(timer, "time taken for assigning vaccination state");
    set_local_parameters(world);
    // add_testing_strategies(world, true, false);
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
                if ((p.get_time_of_last_test() > prev_time)) {
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
                if (p.get_time_of_last_test() > prev_time &&
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

struct LogCumulativeDetectedInfections : mio::LogAlways {
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
                if (p.get_infection().is_detected()) {
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
     * @param[in] sim The simulation of the abm.
     * @return A pair of the TimePoint and the TimeSeries of the estimated reproduction number.
     */
    static Type log(const mio::abm::Simulation& sim)
    {
        Eigen::VectorXd estimation = Eigen::VectorXd::Zero(Eigen::Index(1));

        // time period to take into account for estimating the reproduction number
        // longer periods lead to more averaged results
        mio::abm::TimeSpan time_frame = mio::abm::days(1);
        const auto t                  = sim.get_time();
        const auto persons            = sim.get_world().get_persons();
        const auto virus              = mio::abm::VirusVariant::Alpha;

        // PRAGMA_OMP(parallel for)
        int number_newly_infected = 0;
        int number_infectious     = 0;
        for (auto i = size_t(0); i < persons.size(); ++i) {
            auto& p = persons[i];
            if (p.get_should_be_logged()) {
                // PRAGMA_OMP(atomic)
                if (p.is_infected(t) && !p.is_infected(t - time_frame)) {
                    number_newly_infected += 1;
                }
                // count infectious people at the midpoint of the time period as an estimation
                if (p.is_infected(t - time_frame / 2)) {
                    if (p.get_infection().get_viral_shed(t - time_frame / 2) != 0) {
                        number_infectious += 1;
                    }
                }
            }
        }
        // if the infection dies out, the reproduction number has no meaningful value
        if (number_infectious == 0) {
            return std::make_pair(t, estimation);
        }

        // assume equal parameters for each age
        auto vl_params = sim.get_world().parameters.get<mio::abm::ViralLoadDistributions>()[{virus, age_group_0_to_4}];

        // Assume uniform distribution of parameters. Taking mean value of uniform distribution. Factor 1/2 cancels in the division.
        double average_infectious_period =
            (vl_params.viral_load_peak.params.a() + vl_params.viral_load_peak.params.b()) /
                (vl_params.viral_load_incline.params.a() + vl_params.viral_load_incline.params.b()) -
            (vl_params.viral_load_peak.params.a() + vl_params.viral_load_peak.params.b()) /
                (vl_params.viral_load_decline.params.a() + vl_params.viral_load_decline.params.b());

        estimation[0] += (number_newly_infected * average_infectious_period) / (number_infectious * time_frame.days());
        return std::make_pair(t, estimation);
    }
};

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
        file.close();
    }
}

#endif

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
                                         std::vector<std::vector<double>> grid_points)
{
    int num_procs, rank;
#ifdef MEMILIO_ENABLE_MPI
    MPI_Comm_size(mio::mpi::get_world(), &num_procs);
    MPI_Comm_rank(mio::mpi::get_world(), &rank);
#else
    num_procs = 1;
    rank      = 0;
#endif
    mio::unused(num_runs);

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

PRAGMA_OMP(parallel for)
for (size_t i = 0; i < grid_search_rank.size(); i++) {
    auto params = grid_search_rank[i];

    // grid search for parameters:
    // 1: Viral Shed
    // 2: Seasonality April
    // 3: Seasonality May
    // 4: Perc Easter Event
    // 6: Dark Figure
    // 7.: Contact rate forst ssocial ebents closure
    // 8.: Contact intensity for social events
    // 9.: Contact rate for home
    auto viral_shedding_rate = params[0];
    auto seasonality_april   = params[1];
    auto seasonality_may     = params[2];
    auto perc_easter_event   = params[3];
    auto dark_figure         = params[4];
    auto contact_rate_ssc    = params[5];
    auto masks               = params[6];

    mio::Date start_date{2021, 3, 1};
    int max_num_days     = 90;
    auto max_num_persons = std::numeric_limits<int>::max();

    auto t0   = mio::abm::TimePoint(0); // Start time per simulation
    auto tmax = mio::abm::TimePoint(0) + mio::abm::days(max_num_days); // End time per simulation

    // Determine inital infection state distribution
    restart_timer(timer, "time for initial setup");
    determine_initial_infection_states_world(input_dir, start_date, dark_figure);
    restart_timer(timer, "time for determine_initial_infection_states_world");
    prepare_vaccination_state(mio::offset_date_by_days(start_date, (int)tmax.days()),
                              (input_dir / "pydata/Germany/vacc_county_ageinf_ma7.json").string());
    restart_timer(timer, "time for vaccinaiton state");

    // Loop over a number of runs
    auto world = mio::abm::World(num_age_groupss);

    create_sampled_world(world, input_dir, t0, max_num_persons, start_date, perc_easter_event);

    restart_timer(timer, "time taken for create sampled world");
    auto sim = mio::abm::Simulation(t0, std::move(world));

    //Logger
    mio::History<mio::abm::TimeSeriesWriter, LogInfectionPerLocationTypePerAgeGroup> historyInfectionPerLocationType{
        Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups())};
    mio::History<mio::abm::TimeSeriesWriter, LogInfectionStatePerAgeGroup> historyInfectionStatePerAgeGroup{
        Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups())};
    mio::History<mio::abm::TimeSeriesWriter, LogCumulativeDetectedInfections> historyCumulativeDetectedInfections{
        Eigen::Index(sim.get_world().parameters.get_num_groups())};

    // / NPIS//

    const auto location_it = sim.get_world().get_locations();

    // 1. testing schemes in schools
    auto testing_min_time_school = mio::abm::days(7);
    auto probability_school      = 1.0;
    auto start_date_test_school  = mio::abm::TimePoint(mio::abm::days(42).seconds()); // 2021-04-12
    auto end_date_test_school    = mio::abm::TimePoint(tmax); // 2021-05-30
    auto test_type_school        = mio::abm::TestType::Antigen; // Antigen test
    auto test_parameters = sim.get_world().parameters.get<mio::abm::TestData>()[test_type_school]; // Test parameters
    auto testing_criteria_school = mio::abm::TestingCriteria();
    auto testing_scheme_school =
        mio::abm::TestingScheme(testing_criteria_school, testing_min_time_school, start_date_test_school,
                                end_date_test_school, test_parameters, probability_school);
    sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::School, testing_scheme_school);

    // 2. testing schemes in work places for 35% of random workplaces
    std::vector<uint32_t> work_location_ids;
    for (auto& location : location_it) {
        if (location.get_type() == mio::abm::LocationType::Work) {
            work_location_ids.push_back(location.get_index());
        }
    }
    //take 35% of work locations
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(work_location_ids.begin(), work_location_ids.end(), g);
    auto num_work_locations = (int)(0.35 * work_location_ids.size());
    std::vector<uint32_t> work_location_ids_35(work_location_ids.begin(),
                                               work_location_ids.begin() + num_work_locations);
    auto testing_min_time_work = mio::abm::days(1);
    auto probability_work      = 1.0;
    auto start_date_test_work  = mio::abm::TimePoint(mio::abm::days(72).seconds());
    auto end_date_test_work    = mio::abm::TimePoint(tmax);
    auto test_type_work        = mio::abm::TestType::Antigen; // Antigen test
    auto test_parameters_work = sim.get_world().parameters.get<mio::abm::TestData>()[test_type_work]; // Test parameters
    auto testing_criteria_work = mio::abm::TestingCriteria();
    auto testing_scheme_work =
        mio::abm::TestingScheme(testing_criteria_work, testing_min_time_work, start_date_test_work, end_date_test_work,
                                test_parameters_work, probability_work);
    for (auto& location_id : work_location_ids_35) {
        sim.get_world().get_testing_strategy().add_testing_scheme(
            mio::abm::LocationId{location_id, mio::abm::LocationType::Work}, testing_scheme_work);
    }

    // 2.5 plus testing schemes at 20 % of basics shops
    std::vector<uint32_t> basics_shop_location_ids;
    for (auto& location : location_it) {
        if (location.get_type() == mio::abm::LocationType::BasicsShop) {
            basics_shop_location_ids.push_back(location.get_index());
        }
    }
    //take 20% of basics shop locations
    std::shuffle(basics_shop_location_ids.begin(), basics_shop_location_ids.end(), g);
    auto num_basics_shop_locations = (int)(0.2 * basics_shop_location_ids.size());
    std::vector<uint32_t> basics_shop_location_ids_20(basics_shop_location_ids.begin(),
                                                      basics_shop_location_ids.begin() + num_basics_shop_locations);
    auto testing_min_time_basics_shop = mio::abm::days(2);
    auto probability_basics_shop      = 1.0;
    auto start_date_test_basics_shop  = mio::abm::TimePoint(mio::abm::days(14).seconds());
    auto end_date_test_basics_shop    = mio::abm::TimePoint(tmax);
    auto test_type_basics_shop        = mio::abm::TestType::Antigen; // Antigen test
    auto test_parameters_basics_shop =
        sim.get_world().parameters.get<mio::abm::TestData>()[test_type_basics_shop]; // Test parameters
    auto testing_criteria_basics_shop = mio::abm::TestingCriteria();
    auto testing_scheme_basics_shop =
        mio::abm::TestingScheme(testing_criteria_basics_shop, testing_min_time_basics_shop, start_date_test_basics_shop,
                                end_date_test_basics_shop, test_parameters_basics_shop, probability_basics_shop);
    for (auto& location_id : basics_shop_location_ids_20) {
        sim.get_world().get_testing_strategy().add_testing_scheme(
            mio::abm::LocationId{location_id, mio::abm::LocationType::BasicsShop}, testing_scheme_basics_shop);
    }

    // 3. Mask schemes for all locations
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

    // // 4. Dampings for schools and Basic shops
    for (auto& location : location_it) {
        if (location.get_type() == mio::abm::LocationType::School) {
            location.add_damping(mio::abm::TimePoint(mio::abm::days(0).seconds()), 0.5); // from 2021-03-01
            location.add_damping(mio::abm::TimePoint(mio::abm::days(14).seconds()), 0.00); // from 2021-03-15
            location.add_damping(mio::abm::TimePoint(mio::abm::days(42).seconds()),
                                 0.5); // from 2021-04-12 till 2021-05-30
        }
        if (location.get_type() == mio::abm::LocationType::BasicsShop) {
            location.add_damping(mio::abm::TimePoint(mio::abm::days(14).seconds()), 0.8); // from 2021-03-15
        }
        if (location.get_type() == mio::abm::LocationType::Work) {
            location.add_damping(mio::abm::TimePoint(mio::abm::days(0).seconds()), 0.75); // from 2021-03-01
            location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.75;
        }
    }

    // 5. add capacity limits to some locations
    // first we need two lists, one for 50% of random social event locations and the other list for the other 50%
    // 1. -> Restrict trips to a SocialEvent Location to maximum of 10 (,5,2) for the times ['2021-03-01 to 2021-03-14'], ['2021-03-15 to 2021-05-09'], ['2021-05-10 to 2021-05-30'], as this is the percentage outside this has to be randomly taken from an 50/50 split between inside and outside events, see https://de.statista.com/statistik/daten/studie/171168/umfrage/haeufig-betriebene-freizeitaktivitaeten/
    // 2. -> For the other 50%, we do : -> full closure for 70 days ['2021-03-01 to 2021-05-09'], partial closure (10%) for the remaining days ['2021-05-10 to 2021-05-31']
    // ----> Divide Social Event locations into a 50/50 split. First 50% gƒet the restrictive capacity
    std::vector<int> social_event_location_ids_small;
    std::vector<int> social_event_location_ids_big;
    for (auto& location : location_it) {
        if (location.get_type() == mio::abm::LocationType::SocialEvent) {
            social_event_location_ids_small.push_back(location.get_index());
        }
    }
    //take 50% of social event locations
    std::shuffle(social_event_location_ids_small.begin(), social_event_location_ids_small.end(), g);
    auto num_social_event_locations_big = (int)(0.25 * social_event_location_ids_small.size());
    social_event_location_ids_big.insert(social_event_location_ids_big.end(), social_event_location_ids_small.begin(),
                                         social_event_location_ids_small.begin() + num_social_event_locations_big);
    social_event_location_ids_small.erase(social_event_location_ids_small.begin(),
                                          social_event_location_ids_small.begin() + num_social_event_locations_big);

    // add capacity limits on day one
    for (auto& location : location_it) {
        if (std::find(social_event_location_ids_small.begin(), social_event_location_ids_small.end(),
                      location.get_index()) != social_event_location_ids_small.end()) {
            location.set_capacity(10, 0);
        }
        if (std::find(social_event_location_ids_big.begin(), social_event_location_ids_big.end(),
                      location.get_index()) != social_event_location_ids_big.end()) {
            location.set_capacity(0, 0);
        }
    }

    restart_timer(timer, "till advance 14");
    sim.get_world().parameters.get<mio::abm::MaskProtection>()             = masks;
    sim.get_world().parameters.get<mio::abm::InfectionRateFromViralShed>() = viral_shedding_rate;
    sim.advance(mio::abm::TimePoint(mio::abm::days(14).seconds()), historyInfectionStatePerAgeGroup,
                historyInfectionPerLocationType);
    mio::unused(seasonality_april, seasonality_may, contact_rate_ssc);
    std::cout << "day 14 finished" << std::endl;

    // small social events to capacity 5
    for (auto& location : location_it) {
        if (std::find(social_event_location_ids_small.begin(), social_event_location_ids_small.end(),
                      location.get_index()) != social_event_location_ids_small.end()) {
            location.set_capacity(5, 0);
            location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= contact_rate_ssc;
        }
        if (location.get_type() == mio::abm::LocationType::BasicsShop) {
            location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= contact_rate_ssc;
        }
    }

    restart_timer(timer, "till advance 31 (march ends)");
    sim.advance(mio::abm::TimePoint(mio::abm::days(31).seconds()), historyInfectionStatePerAgeGroup,
                historyInfectionPerLocationType);
    sim.get_world().parameters.get<mio::abm::InfectionRateFromViralShed>() = viral_shedding_rate * seasonality_april;

    restart_timer(timer, "till advance 42");
    sim.advance(mio::abm::TimePoint(mio::abm::days(42).seconds()), historyInfectionStatePerAgeGroup,
                historyInfectionPerLocationType);
    std::cout << "day 42 finished" << std::endl; // 2date 2021-04-12

    for (auto& location : location_it) {
        if (location.get_type() != mio::abm::LocationType::School) {
            location.set_npi_active(false);
        }
    }

    restart_timer(timer, "till advance 61");
    sim.advance(mio::abm::TimePoint(mio::abm::days(61).seconds()), historyInfectionStatePerAgeGroup,
                historyInfectionPerLocationType);
    std::cout << "day 61 finished (date 2021-05-01)" << std::endl;
    sim.get_world().parameters.get<mio::abm::InfectionRateFromViralShed>() = viral_shedding_rate * seasonality_may;

    restart_timer(timer, "till advance 72");
    sim.advance(mio::abm::TimePoint(mio::abm::days(70).seconds()), historyInfectionStatePerAgeGroup,
                historyInfectionPerLocationType);
    std::cout << "day 72 finished (date 2021-05-10)" << std::endl;

    for (auto& location : location_it) {
        if (std::find(social_event_location_ids_small.begin(), social_event_location_ids_small.end(),
                      location.get_index()) != social_event_location_ids_small.end()) {
            location.set_capacity(2, 0);
            location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.5;
        }
        //50% of big social events get reopened and capacity will be limited to xx
        int number_of_big_social_events = (int)(0.9 * social_event_location_ids_big.size());
        if (std::find(social_event_location_ids_big.begin(), social_event_location_ids_big.end(),
                      location.get_index()) != social_event_location_ids_big.end()) {
            number_of_big_social_events--;
            if (number_of_big_social_events >= 0) {
                location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.35;
                location.set_capacity(15, 0);
            }
        }
    }
    for (auto& location : location_it) {
        if (location.get_type() != mio::abm::LocationType::School &&
            location.get_type() != mio::abm::LocationType::Home) {
            location.set_npi_active(true);
        }
    }
    restart_timer(timer, "till advance tmax");
    sim.advance(tmax, historyInfectionStatePerAgeGroup, historyInfectionPerLocationType);
    std::cout << "day 90 finished" << std::endl;

    ////Advance till here
    // Stop the clock after sim.advance and calculate the duration
    restart_timer(timer, "time taken for simulation end");
    // TODO: update result of the simulation to be a vector of location result.
    auto temp_sim_infection_per_loc_tpye =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionPerLocationType.get_log())};
    auto temp_sim_infection_state_per_age_group =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionStatePerAgeGroup.get_log())};

    auto rmse =
        calculate_rmse_from_results(input_dir, temp_sim_infection_state_per_age_group[0], max_num_days, start_date);
    rmse_results_per_grid_point.at(i) = rmse;
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

    auto viral_shedding_rate = 5.0;
    auto seasonality_april   = 0.85;
    auto seasonality_may     = 0.5;
    auto perc_easter_event   = 0.4;
    auto dark_figure         = 2.0;
    auto contact_rate_ssc    = 0.3;
    auto masks               = 0.4;

    mio::Date start_date{2021, 3, 1};
    int max_num_days     = 90;
    auto max_num_persons = 400000;
    bool npis_on         = false;

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

    auto ensemble_positive_test_per_loc_type_per_age_group = std::vector<
        std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of positive tests per location type per age group results
    ensemble_positive_test_per_loc_type_per_age_group.reserve(size_t(num_runs));

    auto ensemble_cumulative_detected_infections =
        std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of cumulative detected infections
    ensemble_cumulative_detected_infections.reserve(size_t(num_runs));

    auto ensemble_estimated_reproduction_number =
        std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of estimated reproduction number
    ensemble_estimated_reproduction_number.reserve(size_t(num_runs));

    auto ensemble_params = std::vector<std::vector<mio::abm::World>>{}; // Vector of all worlds
    ensemble_params.reserve(size_t(num_runs));

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
    restart_timer(timer, "time for initial setup");
    determine_initial_infection_states_world(input_dir, start_date, dark_figure);
    restart_timer(timer, "time for determine_initial_infection_states_world");
    prepare_vaccination_state(mio::offset_date_by_days(start_date, (int)tmax.days()),
                              (input_dir / "pydata/Germany/vacc_county_ageinf_ma7.json").string());
    restart_timer(timer, "time for vaccinaiton state");
    // Loop over a number of runs
    for (size_t run_idx = start_run_idx; run_idx < end_run_idx; run_idx++) {
        auto world = mio::abm::World(num_age_groupss);

        create_sampled_world(world, input_dir, t0, max_num_persons, start_date, perc_easter_event);

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
        mio::History<mio::abm::TimeSeriesWriter, LogCumulativeDetectedInfections> historyCumulativeDetectedInfections{
            Eigen::Index(sim.get_world().parameters.get_num_groups())};
        mio::History<mio::abm::TimeSeriesWriter, LogEstimatedReproductionNumber> historyEstimatedReproductionNumber{
            Eigen::Index(1)};

        // / NPIS//
        if (npis_on) {

            const auto location_it = sim.get_world().get_locations();

            // 1. testing schemes in schools
            auto testing_min_time_school = mio::abm::days(7);
            auto probability_school      = 1.0;
            auto start_date_test_school  = mio::abm::TimePoint(mio::abm::days(42).seconds()); // 2021-04-12
            auto end_date_test_school    = mio::abm::TimePoint(tmax); // 2021-05-30
            auto test_type_school        = mio::abm::TestType::Antigen; // Antigen test
            auto test_parameters =
                sim.get_world().parameters.get<mio::abm::TestData>()[test_type_school]; // Test parameters
            auto testing_criteria_school = mio::abm::TestingCriteria();
            auto testing_scheme_school =
                mio::abm::TestingScheme(testing_criteria_school, testing_min_time_school, start_date_test_school,
                                        end_date_test_school, test_parameters, probability_school);
            sim.get_world().get_testing_strategy().add_testing_scheme(mio::abm::LocationType::School,
                                                                      testing_scheme_school);

            // 2. testing schemes in work places for 35% of random workplaces
            std::vector<uint32_t> work_location_ids;
            for (auto& location : location_it) {
                if (location.get_type() == mio::abm::LocationType::Work) {
                    work_location_ids.push_back(location.get_index());
                }
            }
            //take 35% of work locations
            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(work_location_ids.begin(), work_location_ids.end(), g);
            auto num_work_locations = (int)(0.35 * work_location_ids.size());
            std::vector<uint32_t> work_location_ids_35(work_location_ids.begin(),
                                                       work_location_ids.begin() + num_work_locations);
            auto testing_min_time_work = mio::abm::days(1);
            auto probability_work      = 1.0;
            auto start_date_test_work  = mio::abm::TimePoint(mio::abm::days(72).seconds());
            auto end_date_test_work    = mio::abm::TimePoint(tmax);
            auto test_type_work        = mio::abm::TestType::Antigen; // Antigen test
            auto test_parameters_work =
                sim.get_world().parameters.get<mio::abm::TestData>()[test_type_work]; // Test parameters
            auto testing_criteria_work = mio::abm::TestingCriteria();
            auto testing_scheme_work =
                mio::abm::TestingScheme(testing_criteria_work, testing_min_time_work, start_date_test_work,
                                        end_date_test_work, test_parameters_work, probability_work);
            for (auto& location_id : work_location_ids_35) {
                sim.get_world().get_testing_strategy().add_testing_scheme(
                    mio::abm::LocationId{location_id, mio::abm::LocationType::Work}, testing_scheme_work);
            }

            // 2.5 plus testing schemes at 20 % of basics shops
            std::vector<uint32_t> basics_shop_location_ids;
            for (auto& location : location_it) {
                if (location.get_type() == mio::abm::LocationType::BasicsShop) {
                    basics_shop_location_ids.push_back(location.get_index());
                }
            }
            //take 20% of basics shop locations
            std::shuffle(basics_shop_location_ids.begin(), basics_shop_location_ids.end(), g);
            auto num_basics_shop_locations = (int)(0.2 * basics_shop_location_ids.size());
            std::vector<uint32_t> basics_shop_location_ids_20(
                basics_shop_location_ids.begin(), basics_shop_location_ids.begin() + num_basics_shop_locations);
            auto testing_min_time_basics_shop = mio::abm::days(2);
            auto probability_basics_shop      = 1.0;
            auto start_date_test_basics_shop  = mio::abm::TimePoint(mio::abm::days(14).seconds());
            auto end_date_test_basics_shop    = mio::abm::TimePoint(tmax);
            auto test_type_basics_shop        = mio::abm::TestType::Antigen; // Antigen test
            auto test_parameters_basics_shop =
                sim.get_world().parameters.get<mio::abm::TestData>()[test_type_basics_shop]; // Test parameters
            auto testing_criteria_basics_shop = mio::abm::TestingCriteria();
            auto testing_scheme_basics_shop   = mio::abm::TestingScheme(
                testing_criteria_basics_shop, testing_min_time_basics_shop, start_date_test_basics_shop,
                end_date_test_basics_shop, test_parameters_basics_shop, probability_basics_shop);
            for (auto& location_id : basics_shop_location_ids_20) {
                sim.get_world().get_testing_strategy().add_testing_scheme(
                    mio::abm::LocationId{location_id, mio::abm::LocationType::BasicsShop}, testing_scheme_basics_shop);
            }

            // 3. Mask schemes for all locations
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

            // // 4. Dampings for schools and Basic shops
            for (auto& location : location_it) {
                if (location.get_type() == mio::abm::LocationType::School) {
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(0).seconds()), 0.5); // from 2021-03-01
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(14).seconds()),
                                         0.00); // from 2021-03-15
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(42).seconds()),
                                         0.5); // from 2021-04-12 till 2021-05-30
                }
                if (location.get_type() == mio::abm::LocationType::BasicsShop) {
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(14).seconds()), 0.8); // from 2021-03-15
                }
                if (location.get_type() == mio::abm::LocationType::Work) {
                    location.add_damping(mio::abm::TimePoint(mio::abm::days(0).seconds()), 0.75); // from 2021-03-01
                    location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.7;
                }
            }

            // 5. add capacity limits to some locations
            // first we need two lists, one for 50% of random social event locations and the other list for the other 50%
            // 1. -> Restrict trips to a SocialEvent Location to maximum of 10 (,5,2) for the times ['2021-03-01 to 2021-03-14'], ['2021-03-15 to 2021-05-09'], ['2021-05-10 to 2021-05-30'], as this is the percentage outside this has to be randomly taken from an 50/50 split between inside and outside events, see https://de.statista.com/statistik/daten/studie/171168/umfrage/haeufig-betriebene-freizeitaktivitaeten/
            // 2. -> For the other 50%, we do : -> full closure for 70 days ['2021-03-01 to 2021-05-09'], partial closure (10%) for the remaining days ['2021-05-10 to 2021-05-31']
            // ----> Divide Social Event locations into a 50/50 split. First 50% gƒet the restrictive capacity
            std::vector<int> social_event_location_ids_small;
            std::vector<int> social_event_location_ids_big;
            for (auto& location : location_it) {
                if (location.get_type() == mio::abm::LocationType::SocialEvent) {
                    social_event_location_ids_small.push_back(location.get_index());
                }
            }
            //take 50% of social event locations
            std::shuffle(social_event_location_ids_small.begin(), social_event_location_ids_small.end(), g);
            auto num_social_event_locations_big = (int)(0.25 * social_event_location_ids_small.size());
            social_event_location_ids_big.insert(
                social_event_location_ids_big.end(), social_event_location_ids_small.begin(),
                social_event_location_ids_small.begin() + num_social_event_locations_big);
            social_event_location_ids_small.erase(social_event_location_ids_small.begin(),
                                                  social_event_location_ids_small.begin() +
                                                      num_social_event_locations_big);

            // add capacity limits on day one
            for (auto& location : location_it) {
                if (std::find(social_event_location_ids_small.begin(), social_event_location_ids_small.end(),
                              location.get_index()) != social_event_location_ids_small.end()) {
                    location.set_capacity(10, 0);
                }
                if (std::find(social_event_location_ids_big.begin(), social_event_location_ids_big.end(),
                              location.get_index()) != social_event_location_ids_big.end()) {
                    location.set_capacity(0, 0);
                }
            }

            restart_timer(timer, "till advance 14");
            sim.get_world().parameters.get<mio::abm::InfectionRateFromViralShed>() = viral_shedding_rate;
            sim.get_world().parameters.get<mio::abm::MaskProtection>()             = masks;
            sim.advance(mio::abm::TimePoint(mio::abm::days(14).seconds()), historyInfectionPerLocationTypePerAgeGroup,
                        historyInfectionStatePerAgeGroup, historyTestPerLocationTypePerAgeGroup,
                        historyPositiveTestPerLocationTypePerAgeGroup, historyCumulativeDetectedInfections,
                        historyEstimatedReproductionNumber);
            std::cout << "day 14 finished" << std::endl;

            // small social events to capacity 5
            for (auto& location : location_it) {
                if (std::find(social_event_location_ids_small.begin(), social_event_location_ids_small.end(),
                              location.get_index()) != social_event_location_ids_small.end()) {
                    location.set_capacity(5, 0);
                    location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= contact_rate_ssc;
                }
                if (location.get_type() == mio::abm::LocationType::BasicsShop) {
                    location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= contact_rate_ssc;
                }
            }

            restart_timer(timer, "till advance 31 (march ends)");
            sim.advance(mio::abm::TimePoint(mio::abm::days(31).seconds()), historyInfectionPerLocationTypePerAgeGroup,
                        historyInfectionStatePerAgeGroup, historyTestPerLocationTypePerAgeGroup,
                        historyPositiveTestPerLocationTypePerAgeGroup, historyCumulativeDetectedInfections,
                        historyEstimatedReproductionNumber);
            sim.get_world().parameters.get<mio::abm::InfectionRateFromViralShed>() =
                viral_shedding_rate * seasonality_april;

            restart_timer(timer, "till advance 42");
            sim.advance(mio::abm::TimePoint(mio::abm::days(42).seconds()), historyInfectionPerLocationTypePerAgeGroup,
                        historyInfectionStatePerAgeGroup, historyTestPerLocationTypePerAgeGroup,
                        historyPositiveTestPerLocationTypePerAgeGroup, historyCumulativeDetectedInfections,
                        historyEstimatedReproductionNumber);
            std::cout << "day 42 finished" << std::endl; // 2date 2021-04-12

            for (auto& location : location_it) {
                if (location.get_type() != mio::abm::LocationType::School) {
                    location.set_npi_active(false);
                }
            }

            restart_timer(timer, "till advance 62");
            sim.advance(mio::abm::TimePoint(mio::abm::days(61).seconds()), historyInfectionPerLocationTypePerAgeGroup,
                        historyInfectionStatePerAgeGroup, historyTestPerLocationTypePerAgeGroup,
                        historyPositiveTestPerLocationTypePerAgeGroup, historyCumulativeDetectedInfections,
                        historyEstimatedReproductionNumber);
            std::cout << "day 62 finished (date 2021-05-01)" << std::endl;
            sim.get_world().parameters.get<mio::abm::InfectionRateFromViralShed>() =
                viral_shedding_rate * seasonality_may;

            restart_timer(timer, "till advance 72");
            sim.advance(mio::abm::TimePoint(mio::abm::days(71).seconds()), historyInfectionPerLocationTypePerAgeGroup,
                        historyInfectionStatePerAgeGroup, historyTestPerLocationTypePerAgeGroup,
                        historyPositiveTestPerLocationTypePerAgeGroup, historyCumulativeDetectedInfections,
                        historyEstimatedReproductionNumber);
            std::cout << "day 72 finished (date 2021-05-10)" << std::endl;

            for (auto& location : location_it) {
                if (std::find(social_event_location_ids_small.begin(), social_event_location_ids_small.end(),
                              location.get_index()) != social_event_location_ids_small.end()) {
                    location.set_capacity(2, 0);
                    location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.5;
                }
                //50% of big social events get reopened and capacity will be limited to xx
                int number_of_big_social_events = (int)(0.9 * social_event_location_ids_big.size());
                if (std::find(social_event_location_ids_big.begin(), social_event_location_ids_big.end(),
                              location.get_index()) != social_event_location_ids_big.end()) {
                    number_of_big_social_events--;
                    if (number_of_big_social_events >= 0) {
                        location.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.35;
                        location.set_capacity(15, 0);
                    }
                }
            }
            for (auto& location : location_it) {
                if (location.get_type() != mio::abm::LocationType::School &&
                    location.get_type() != mio::abm::LocationType::Home) {
                    location.set_npi_active(true);
                }
            }
            restart_timer(timer, "till advance tmax");
            sim.advance(tmax, historyInfectionPerLocationTypePerAgeGroup, historyInfectionStatePerAgeGroup,
                        historyTestPerLocationTypePerAgeGroup, historyPositiveTestPerLocationTypePerAgeGroup,
                        historyCumulativeDetectedInfections, historyEstimatedReproductionNumber);
            std::cout << "day 90 finished" << std::endl;
        }
        else {
            sim.advance(mio::abm::TimePoint(mio::abm::days(20).seconds()), historyInfectionPerLocationTypePerAgeGroup,
                        historyInfectionStatePerAgeGroup, historyTestPerLocationTypePerAgeGroup,
                        historyPositiveTestPerLocationTypePerAgeGroup, historyCumulativeDetectedInfections,
                        historyEstimatedReproductionNumber);
        }
        ////Advance till here
        // Stop the clock after sim.advance and calculate the duration
        restart_timer(timer, "time taken for simulation end");
        // TODO: update result of the simulation to be a vector of location result.
        auto temp_sim_infection_per_loc_type_per_age_group =
            std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionPerLocationTypePerAgeGroup.get_log())};
        auto temp_sim_infection_state_per_age_group =
            std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionStatePerAgeGroup.get_log())};
        auto temp_sim_test_per_loc_type_per_age_group =
            std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyTestPerLocationTypePerAgeGroup.get_log())};
        auto temp_sim_positive_test_per_loc_type_per_age_group = std::vector<mio::TimeSeries<ScalarType>>{
            std::get<0>(historyPositiveTestPerLocationTypePerAgeGroup.get_log())};
        auto temp_sim_cumulative_detected_infections_per_age_group =
            std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyCumulativeDetectedInfections.get_log())};
        auto temp_sim_estimated_reproduction_number =
            std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyEstimatedReproductionNumber.get_log())};
        // Push result of the simulation back to the result vector
        ensemble_infection_per_loc_type_per_age_group.emplace_back(temp_sim_infection_per_loc_type_per_age_group);
        ensemble_infection_state_per_age_group.emplace_back(temp_sim_infection_state_per_age_group);
        ensemble_test_per_loc_type_per_age_group.emplace_back(temp_sim_test_per_loc_type_per_age_group);
        ensemble_positive_test_per_loc_type_per_age_group.emplace_back(
            temp_sim_positive_test_per_loc_type_per_age_group);
        ensemble_cumulative_detected_infections.emplace_back(temp_sim_cumulative_detected_infections_per_age_group);
        ensemble_estimated_reproduction_number.emplace_back(temp_sim_estimated_reproduction_number);

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
    auto final_ensemble_estimated_reproduction_number =
        gather_results(rank, num_procs, num_runs, ensemble_estimated_reproduction_number);
    if (rank == 0) {
        BOOST_OUTCOME_TRY(save_results(final_ensemble_infection_per_loc_type_per_age_group, ensemble_params, {0},
                                       result_dir / "infection_per_location_type_per_age_group/", save_single_runs));
        BOOST_OUTCOME_TRY(save_results(final_ensemble_infection_state_per_age_group, ensemble_params, {0},
                                       result_dir / "infection_state_per_age_group/", save_single_runs));
        BOOST_OUTCOME_TRY(save_results(final_ensemble_test_per_loc_type_per_age_group, ensemble_params, {0},
                                       result_dir / "test_per_location_type_per_age_group/", save_single_runs));
        BOOST_OUTCOME_TRY(save_results(final_ensemble_positive_test_per_loc_type_per_age_group, ensemble_params, {0},
                                       result_dir / "positive_test_per_location_type_per_age_group/",
                                       save_single_runs));
        BOOST_OUTCOME_TRY(save_results(final_ensemble_cumulative_detected_infections, ensemble_params, {0},
                                       result_dir / "cumulative_detected_infections/", save_single_runs));
        BOOST_OUTCOME_TRY(save_results(final_ensemble_estimated_reproduction_number, ensemble_params, {0},
                                       result_dir / "estimated_reproduction_number/", save_single_runs));
    }
#else

    BOOST_OUTCOME_TRY(save_results(ensemble_infection_state_per_age_group, ensemble_params, {0},
                                   result_dir / "infection_state_per_age_group/", save_single_runs));
    BOOST_OUTCOME_TRY(save_results(ensemble_test_per_loc_type_per_age_group, ensemble_params, {0},
                                   result_dir / "test_per_location_type_per_age_group/", save_single_runs));
    BOOST_OUTCOME_TRY(save_results(ensemble_positive_test_per_loc_type_per_age_group, ensemble_params, {0},
                                   result_dir / "positive_test_per_location_type_per_age_group/", save_single_runs));
    BOOST_OUTCOME_TRY(save_results(ensemble_cumulative_detected_infections, ensemble_params, {0},
                                   result_dir / "cumulative_detected_infections/", save_single_runs));
    BOOST_OUTCOME_TRY(save_results(ensemble_estimated_reproduction_number, ensemble_params, {0},
                                   result_dir / "estimated_reproduction_number/", save_single_runs));

#endif
    restart_timer(timer, "time taken for data gathering and saving results");

    printf("done.\n");
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

mio::IOResult<bool> create_result_folders(std::string const& result_dir, bool grid_search = false)
{
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir + "/infection_per_location_type_per_age_group/"));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir + "/infection_state_per_age_group/"));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir + "/test_per_location_type_per_age_group/"));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir + "/positive_test_per_location_type_per_age_group/"));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir + "/cumulative_detected_infections/"));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir + "/estimated_reproduction_number/"));

    if (grid_search) {
        BOOST_OUTCOME_TRY(mio::create_directory(result_dir + "/grid_search/"));
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
#ifdef MEMILIO_ENABLE_MPI
    mio::mpi::init();
#endif

    std::string input_dir = "/p/project1/loki/memilio/memilio/data";
    // std::string input_dir = "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data";
    // std::string input_dir = "/Users/david/Documents/HZI/memilio/data";
    // std::string input_dir       = "C:/Users/korf_sa/Documents/rep/data";
    std::string precomputed_dir = input_dir + "/results";
    std::string result_dir      = input_dir + "/results_" + currentDateTime();

    size_t num_runs;
    bool run_grid_search = false;

    if (argc == 2) {
        num_runs        = std::stoi(argv[1]);
        run_grid_search = false;
        printf("Running with number of runs = %d.\n", (int)num_runs);
        printf("Saving results to \"%s\".\n", result_dir.c_str());
    }
    else if (argc == 3) {
        num_runs        = 1;
        run_grid_search = true;
        printf("running with grid search\n");
        printf("Running with number of runs = 1\n");
    }
    else {
        num_runs = 1;
        printf("Running with number of runs = %d.\n", (int)num_runs);
        printf("Saving results to \"%s\".\n", result_dir.c_str());
    }

    auto created = create_result_folders(result_dir, run_grid_search);
    if (!created) {
        std::cout << created.error().formatted_message();
        return created.error().code().value();
    }
    timer = TIME_NOW;

    if (run_grid_search) {
        // grid search for parameters:
        // 1: Viral Shed
        // 2: Seasonality April
        // 3: Seasonality May
        // 4: Perc Easter Event
        // 6: Dark Figure
        // 7.: Contact rate forst ssocial ebents closure
        // 8.: Masks
        std::vector<std::pair<double, double>> grid_boundaries = {{4.0, 7.0}, {0.8, 0.95}, {0.5, 0.7}, {0.3, 0.6},
                                                                  {1.0, 5.0}, {0.2, 0.6},  {0.3, 0.5}};

        std::vector<int> points_per_dim = {5, 3, 3, 3, 7, 5, 3};
        auto grid                       = grid_points(grid_boundaries, points_per_dim);
        std::cout << "Grid size: " << grid.size() << std::endl;
        auto result = run_with_grid_search(input_dir, result_dir, num_runs, grid);
    }
    else {
        auto result = run(input_dir, result_dir, num_runs);
    }

    // copy results into a fixed name folder to have easier access
    std::string last_run_dir = input_dir + "/results_last_run";
    auto copied              = copy_result_folder(result_dir, last_run_dir);

    mio::mpi::finalize();

    auto end                                      = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    //print time taken
    std::cout << "Time taken: " << elapsed_seconds.count() << "s\n";

    return 0;
}
