/*
* Copyright (C) 2020-2025 MEmilio
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
#include "abm/common_abm_loggers.h"
#include "abm/location_id.h"
#include "abm/lockdown_rules.h"
#include "abm/parameters.h"
#include "abm/person.h"
#include "abm/simulation.h"
#include "abm/model.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/io/io.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/uncertain_value.h"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"

#include <fstream>
#include <vector>
#include <iostream>

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
 * Set a value and distribution of an UncertainValue.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution.
 * @param p uncertain value to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void assign_uniform_distribution(mio::UncertainValue<>& p, ScalarType min, ScalarType max)
{
    p = mio::UncertainValue<>(0.5 * (max + min));
    p.set_distribution(mio::ParameterDistributionUniform(min, max));
}

/**
 * Determine the infection state of a person at the beginning of the simulation.
 * The infection states are chosen randomly. They are distributed according to the probabilites set in the example.
 * @return random infection state
 */
mio::abm::InfectionState determine_infection_state(mio::abm::PersonalRandomNumberGenerator& rng, ScalarType exposed,
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
void assign_infection_state(mio::abm::Model& model, mio::abm::TimePoint t, double exposed_prob,
                            double infected_no_symptoms_prob, double infected_symptoms_prob, double recovered_prob)
{
    auto persons = model.get_persons();
    for (auto& person : persons) {
        auto rng             = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        auto infection_state = determine_infection_state(rng, exposed_prob, infected_no_symptoms_prob,
                                                         infected_symptoms_prob, recovered_prob);
        if (infection_state != mio::abm::InfectionState::Susceptible)
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model.parameters, t, infection_state));
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

void create_model_from_data(mio::abm::Model& model, const std::string& filename, const mio::abm::TimePoint t0,
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

    std::map<uint32_t, mio::abm::LocationId> locations        = {};
    std::map<uint32_t, mio::abm::PersonId> pids_data_to_model = {};
    std::map<uint32_t, uint32_t> person_ids                   = {};
    std::map<uint32_t, std::pair<uint32_t, int>> locations_before;
    std::map<uint32_t, std::pair<uint32_t, int>> locations_after;

    // For the model we need: Hospitals, ICUs (for both we just create one for now), Homes for each unique householdID, One Person for each person_id with respective age and home_id.

    // We assume that no person goes to an hospital, altough e.g. "Sonstiges" could be a hospital
    auto hospital = model.add_location(mio::abm::LocationType::Hospital);
    model.get_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    model.get_location(hospital).set_capacity(std::numeric_limits<uint32_t>::max(),
                                              std::numeric_limits<uint32_t>::max());
    auto icu = model.add_location(mio::abm::LocationType::ICU);
    model.get_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    model.get_location(icu).set_capacity(std::numeric_limits<uint32_t>::max(), std::numeric_limits<uint32_t>::max());

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

    // Add all locations to the model
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
            home = model.add_location(mio::abm::LocationType::Home, 1);
            locations.insert({home_id, home});
            mio::abm::GeographicalLocation location_long_lat_home = {(double)row[index["lon_start"]] / 1e+5,
                                                                     (double)row[index["lat_start"]] / 1e+5};
            model.get_location(home).set_geographical_location(location_long_lat_home);
        }
        else {
            home = it_home->second;
        }

        mio::abm::LocationId location;
        auto it_location = locations.find(
            target_location_id); // Check if location already exists also for home which have the same id (home_id = target_location_id)
        if (it_location == locations.end()) {
            location = model.add_location(
                get_location_type(activity_end),
                1); // Assume one place has one activity, this may be untrue but not important for now(?)
            locations.insert({target_location_id, location});
            model.get_location(location).set_geographical_location(location_long_lat);
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

        uint32_t person_data_id = row[index["puid"]];
        if (person_ids.find(person_data_id) == person_ids.end())
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

        auto pid_itr = pids_data_to_model.find(person_data_id);

        if (pid_itr == pids_data_to_model.end()) { // person has not been added to model yet
            auto it_first_location_id = locations_before.find(person_data_id);
            if (it_first_location_id == locations_before.end()) {
                it_first_location_id = locations_after.find(person_data_id);
            }
            auto first_location_id = it_first_location_id->second.first;
            auto first_location    = locations.find(first_location_id)->second;
            auto person_model_id   = model.add_person(first_location, determine_age_group(age));
            auto home              = locations.find(home_id)->second;
            model.assign_location(person_model_id, home);
            model.assign_location(person_model_id, hospital);
            model.assign_location(person_model_id, icu);
            pid_itr = pids_data_to_model.insert_or_assign(person_data_id, person_model_id).first;
        }

        model.assign_location(
            pid_itr->second,
            target_location); //This assumes that we only have in each tripchain only one location type for each person
        if (locations.find(start_location_id) == locations.end()) {
            // For trips where the start location is not known use Home instead
            start_location = model.get_person(pid_itr->second).get_assigned_location(mio::abm::LocationType::Home);
        }
        model.get_trip_list().add_trip(mio::abm::Trip(
            pid_itr->second, mio::abm::TimePoint(0) + mio::abm::minutes(trip_start), target_location, start_location,
            mio::abm::TransportMode(transport_mode), mio::abm::ActivityType(acticity_end)));
    }
    model.get_trip_list().use_weekday_trips_on_weekend();
}

void set_parameters(mio::abm::Parameters params)
{
    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    params.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 and 35-59)
    params.get<mio::abm::AgeGroupGotoWork>().set_multiple({age_group_15_to_34, age_group_35_to_59}, true);

    params.set<mio::abm::IncubationPeriod>({{mio::abm::VirusVariant::Count, mio::AgeGroup(num_age_groups)}, 4.});

    //0-4
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]  = 0.276;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = 0.092;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.142;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]      = 0.001;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]             = 0.186;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]              = 0.015;
    params.get<mio::abm::SevereToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]                  = 0.001;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]           = 0.143;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]                = 0.001;

    //5-14
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = 0.276;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] =
        0.092;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = 0.142;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]    = 0.001;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]           = 0.186;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]            = 0.015;
    params.get<mio::abm::SevereToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]                = 0.001;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]         = 0.143;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]              = 0.001;
    params.get<mio::abm::RecoveredToSusceptible>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]      = 0.;

    //15-34
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        0.315;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        0.079;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.139;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}]    = 0.003;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}]           = 0.157;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}]            = 0.013;
    params.get<mio::abm::SevereToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}]                = 0.021;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}]         = 0.126;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}]              = 0.021;

    //35-59
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] =
        0.315;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] =
        0.079;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = 0.136;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]    = 0.009;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]           = 0.113;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]            = 0.02;
    params.get<mio::abm::SevereToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]                = 0.008;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]         = 0.05;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]              = 0.008;
    params.get<mio::abm::RecoveredToSusceptible>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]      = 0.;

    //60-79
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] =
        0.315;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] =
        0.079;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.123;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]    = 0.024;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]           = 0.083;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]            = 0.035;
    params.get<mio::abm::SevereToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]                = 0.023;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]         = 0.035;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]              = 0.023;
    params.get<mio::abm::RecoveredToSusceptible>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]      = 0.;

    //80+
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] = 0.315;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] =
        0.079;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] = 0.115;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]    = 0.033;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]           = 0.055;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]            = 0.036;
    params.get<mio::abm::SevereToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]                = 0.052;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]         = 0.035;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]              = 0.052;
    params.get<mio::abm::RecoveredToSusceptible>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]      = 0.;

    // Set each parameter for vaccinated people including personal infection and vaccine protection levels.
    // Summary: https://doi.org/10.1038/s41577-021-00550-x,

    //0-4
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]  = 0.161;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = 0.132;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.143;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]      = 0.001;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]             = 0.186;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]              = 0.015;
    params.get<mio::abm::SevereToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]                  = 0.001;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]           = 0.143;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]                = 0.001;
    params.get<mio::abm::RecoveredToSusceptible>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]        = 0.0;
    // Protection of reinfection is the same for all age-groups, based on:
    // https://doi.org/10.1016/S0140-6736(22)02465-5, https://doi.org/10.1038/s41591-021-01377-8
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::NaturalInfection, age_group_0_to_4,
                                                       mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.852},
         {180, 0.852},
         {210, 0.845},
         {240, 0.828},
         {270, 0.797},
         {300, 0.759},
         {330, 0.711},
         {360, 0.661},
         {390, 0.616},
         {420, 0.580},
         {450, 0.559},
         {450, 0.550}}};

    // Information is based on: https://doi.org/10.1016/S0140-6736(21)02183-8
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_0_to_4,
                                                       mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.5}, {30, 0.91}, {60, 0.92}, {90, 0.88}, {120, 0.84}, {150, 0.81}, {180, 0.88}, {450, 0.5}}};

    // Set up age-related severe protection levels, based on:
    // https://doi.org/10.1016/S0140-6736(22)02465-5
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::NaturalInfection, age_group_0_to_4,
                                                      mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.967},
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
         {450, 0.5}}};
    // Information is based on: https://doi.org/10.1016/S0140-6736(21)02183-8
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_0_to_4,
                                                      mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.5}, {30, 0.88}, {60, 0.91}, {90, 0.98}, {120, 0.94}, {150, 0.88}, {450, 0.5}}};

    //5-14
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = 0.161;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] =
        0.132;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = 0.143;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]    = 0.001;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]           = 0.186;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]            = 0.015;
    params.get<mio::abm::SevereToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]                = 0.001;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]         = 0.143;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]              = 0.001;
    params.get<mio::abm::RecoveredToSusceptible>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]      = 0.0;
    // Protection of reinfection is the same for all age-groups, based on:
    // https://doi.org/10.1016/S0140-6736(22)02465-5, https://doi.org/10.1038/s41591-021-01377-8
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::NaturalInfection, age_group_5_to_14,
                                                       mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.852},
         {180, 0.852},
         {210, 0.845},
         {240, 0.828},
         {270, 0.797},
         {300, 0.759},
         {330, 0.711},
         {360, 0.661},
         {390, 0.616},
         {420, 0.580},
         {450, 0.559},
         {450, 0.550}}};
    // Information is based on: https://doi.org/10.1016/S0140-6736(21)02183-8
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_5_to_14,
                                                       mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.5}, {30, 0.91}, {60, 0.92}, {90, 0.88}, {120, 0.84}, {150, 0.81}, {180, 0.88}, {450, 0.5}}};

    // Set up age-related severe protection levels, based on:
    // https://doi.org/10.1016/S0140-6736(22)02465-5
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::NaturalInfection, age_group_5_to_14,
                                                      mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.967},
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
         {450, 0.5}}};
    // Information is based on: https://doi.org/10.1016/S0140-6736(21)02183-8
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_5_to_14,
                                                      mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.5}, {30, 0.88}, {60, 0.91}, {90, 0.98}, {120, 0.94}, {150, 0.88}, {450, 0.5}}};

    //15-34
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        0.179;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        0.126;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.142;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}]    = 0.001;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}]           = 0.157;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}]            = 0.013;
    params.get<mio::abm::SevereToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}]                = 0.021;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}]         = 0.126;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}]              = 0.021;
    params.get<mio::abm::RecoveredToSusceptible>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}]      = 0.0;
    // Set up personal infection and vaccine protection levels, based on: https://doi.org/10.1038/s41577-021-00550-x, https://doi.org/10.1038/s41591-021-01377-8
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::NaturalInfection, age_group_15_to_34,
                                                       mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.852},
         {180, 0.852},
         {210, 0.845},
         {240, 0.828},
         {270, 0.797},
         {300, 0.759},
         {330, 0.711},
         {360, 0.661},
         {390, 0.616},
         {420, 0.580},
         {450, 0.559},
         {450, 0.550}}};
    // Information is based on: https://doi.org/10.1016/S0140-6736(21)02183-8
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_15_to_34,
                                                       mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.5}, {30, 0.89}, {60, 0.84}, {90, 0.78}, {120, 0.68}, {150, 0.57}, {180, 0.39}, {450, 0.1}}};
    // Set up age-related severe protection levels, based on:
    // https://doi.org/10.1016/S0140-6736(22)02465-5
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::NaturalInfection, age_group_15_to_34,
                                                      mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.967},
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
         {450, 0.5}}};
    // Information is from: https://doi.org/10.1016/S0140-6736(21)02183-8
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_15_to_34,
                                                      mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.5}, {30, 0.88}, {60, 0.91}, {90, 0.98}, {120, 0.94}, {150, 0.88}, {180, 0.90}, {450, 0.5}}};

    //35-59
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] =
        0.179;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] =
        0.126;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = 0.141;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]    = 0.003;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]           = 0.113;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]            = 0.02;
    params.get<mio::abm::SevereToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]                = 0.008;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]         = 0.05;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]              = 0.008;
    params.get<mio::abm::RecoveredToSusceptible>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}]      = 0.0;
    // Protection of reinfection is the same for all age-groups, based on:
    // https://doi.org/10.1016/S0140-6736(22)02465-5, https://doi.org/10.1038/s41591-021-01377-8
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::NaturalInfection, age_group_35_to_59,
                                                       mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.852},
         {180, 0.852},
         {210, 0.845},
         {240, 0.828},
         {270, 0.797},
         {300, 0.759},
         {330, 0.711},
         {360, 0.661},
         {390, 0.616},
         {420, 0.580},
         {450, 0.559},
         {450, 0.550}}};
    // Information is based on: https://doi.org/10.1016/S0140-6736(21)02183-8
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_35_to_59,
                                                       mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.5}, {30, 0.89}, {60, 0.84}, {90, 0.78}, {120, 0.68}, {150, 0.57}, {180, 0.39}, {450, 0.1}}};
    // Set up age-related severe protection levels, based on:
    // https://doi.org/10.1016/S0140-6736(22)02465-5
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::NaturalInfection, age_group_35_to_59,
                                                      mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.967},
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
         {450, 0.5}}};
    // Information is from: https://doi.org/10.1016/S0140-6736(21)02183-8
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_35_to_59,
                                                      mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.5}, {30, 0.88}, {60, 0.91}, {90, 0.98}, {120, 0.94}, {150, 0.88}, {180, 0.90}, {450, 0.5}}};
    //60-79
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] =
        0.179;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] =
        0.126;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.136;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]    = 0.009;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]           = 0.083;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]            = 0.035;
    params.get<mio::abm::SevereToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]                = 0.023;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]         = 0.035;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]              = 0.023;
    params.get<mio::abm::RecoveredToSusceptible>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}]      = 0.0;
    // Protection of reinfection is the same for all age-groups, based on:
    // https://doi.org/10.1016/S0140-6736(22)02465-5, https://doi.org/10.1038/s41591-021-01377-8
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::NaturalInfection, age_group_60_to_79,
                                                       mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.852},
         {180, 0.852},
         {210, 0.845},
         {240, 0.828},
         {270, 0.797},
         {300, 0.759},
         {330, 0.711},
         {360, 0.661},
         {390, 0.616},
         {420, 0.580},
         {450, 0.559},
         {450, 0.550}}};
    // Information is based on: https://doi.org/10.1016/S0140-6736(21)02183-8
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_60_to_79,
                                                       mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.5}, {30, 0.87}, {60, 0.85}, {90, 0.78}, {120, 0.67}, {150, 0.61}, {180, 0.50}, {450, 0.1}}};
    // Set up personal severe protection levels.
    // Protection of severe infection of age group 65 + is different from other age group, based on:
    // https://doi.org/10.1016/S0140-6736(22)02465-5
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::NaturalInfection, age_group_60_to_79,
                                                      mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.967},
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
         {360, 0.5}}};
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_60_to_79,
                                                      mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.5}, {30, 0.91}, {60, 0.86}, {90, 0.91}, {120, 0.94}, {150, 0.95}, {180, 0.90}, {450, 0.5}}};

    //80+
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] = 0.179;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] =
        0.126;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] = 0.133;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]    = 0.012;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]           = 0.055;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]            = 0.036;
    params.get<mio::abm::SevereToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]                = 0.052;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]         = 0.035;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]              = 0.052;
    params.get<mio::abm::RecoveredToSusceptible>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]      = 0.0;
    // Protection of reinfection is the same for all age-groups, based on:
    // https://doi.org/10.1016/S0140-6736(22)02465-5, https://doi.org/10.1038/s41591-021-01377-8
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::NaturalInfection, age_group_80_plus,
                                                       mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.852},
         {180, 0.852},
         {210, 0.845},
         {240, 0.828},
         {270, 0.797},
         {300, 0.759},
         {330, 0.711},
         {360, 0.661},
         {390, 0.616},
         {420, 0.580},
         {450, 0.559},
         {450, 0.550}}};
    // Information is from: https://doi.org/10.1016/S0140-6736(21)02183-8
    params.get<mio::abm::InfectionProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_80_plus,
                                                       mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.5}, {30, 0.80}, {60, 0.79}, {90, 0.75}, {120, 0.56}, {150, 0.49}, {180, 0.43}, {450, 0.1}}};
    // Set up personal severe protection levels.
    // Protection of severe infection of age group 65 + is different from other age group, based on:
    // https://doi.org/10.1016/S0140-6736(22)02465-5
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::NaturalInfection, age_group_0_to_4,
                                                      mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.967},
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
         {360, 0.5}}};
    // Information is based on: https://doi.org/10.1016/S0140-6736(21)02183-8
    params.get<mio::abm::SeverityProtectionFactor>()[{mio::abm::ProtectionType::GenericVaccine, age_group_80_plus,
                                                      mio::abm::VirusVariant::Wildtype}] = {
        mio::TimeSeriesFunctorType::LinearInterpolation,
        {{0, 0.5}, {30, 0.84}, {60, 0.88}, {90, 0.89}, {120, 0.86}, {150, 0.85}, {180, 0.83}, {450, 0.5}}};
}

/**
 * Create a sampled simulation with start time t0.
 * @param t0 The start time of the Simulation.
 */
mio::abm::Simulation create_sampled_simulation(const std::string& input_file, const mio::abm::TimePoint& t0,
                                               int max_num_persons)
{
    // Assumed percentage of infection state at the beginning of the simulation.
    ScalarType exposed_prob = 0.005, infected_no_symptoms_prob = 0.001, infected_symptoms_prob = 0.001,
               recovered_prob = 0.0;

    //Set global infection parameters (similar to infection parameters in SECIR model) and initialize the model
    auto model = mio::abm::Model(num_age_groups);

    set_parameters(model.parameters);

    // Create the model object from statistical data.
    create_model_from_data(model, input_file, t0, max_num_persons);
    model.use_mobility_rules(false);

    // Assign an infection state to each person.
    assign_infection_state(model, t0, exposed_prob, infected_no_symptoms_prob, infected_symptoms_prob, recovered_prob);

    auto sim = mio::abm::Simulation(t0, std::move(model));
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

    auto mobility_data = std::get<0>(history.get_log());
    std::ofstream myfile3("mobility_data.txt");
    myfile3 << "agent_id, trip_id, start_location, end_location, start_time, end_time, transport_mode, activity, "
               "infection_state \n";
    int trips_id = 0;
    for (uint32_t mobility_data_index = 2; mobility_data_index < mobility_data.size(); ++mobility_data_index) {
        myfile3 << "timestep Nr.: " << mobility_data_index - 1 << "\n";
        for (uint32_t trip_index = 0; trip_index < mobility_data[mobility_data_index].size(); trip_index++) {
            auto agent_id = std::get<0>(mobility_data[mobility_data_index][trip_index]);

            int start_index = mobility_data_index - 1;
            using Type      = std::tuple<mio::abm::PersonId, mio::abm::LocationId, mio::abm::TimePoint,
                                    mio::abm::TransportMode, mio::abm::ActivityType, mio::abm::InfectionState>;
            while (!std::binary_search(std::begin(mobility_data[start_index]), std::end(mobility_data[start_index]),
                                       mobility_data[mobility_data_index][trip_index],
                                       [](const Type& v1, const Type& v2) {
                                           return std::get<0>(v1) < std::get<0>(v2);
                                       })) {
                start_index--;
            }
            auto start_location_iterator =
                std::lower_bound(std::begin(mobility_data[start_index]), std::end(mobility_data[start_index]),
                                 mobility_data[mobility_data_index][trip_index], [](const Type& v1, const Type& v2) {
                                     return std::get<0>(v1) < std::get<0>(v2);
                                 });
            auto start_location = (int)std::get<1>(*start_location_iterator).get();

            auto end_location = (int)std::get<1>(mobility_data[mobility_data_index][trip_index]).get();

            auto start_time = (int)std::get<2>(mobility_data[mobility_data_index][trip_index]).seconds();
            auto end_time   = (int)std::get<2>(mobility_data[mobility_data_index][trip_index]).seconds();

            auto transport_mode  = (int)std::get<3>(mobility_data[mobility_data_index][trip_index]);
            auto activity        = (int)std::get<4>(mobility_data[mobility_data_index][trip_index]);
            auto infection_state = (int)std::get<5>(mobility_data[mobility_data_index][trip_index]);
            myfile3 << agent_id << ", " << trips_id << ", " << start_location << " , " << end_location << " , "
                    << start_time << " , " << end_time << " , " << transport_mode << " , " << activity << " , "
                    << infection_state << "\n";
            trips_id++;
        }
    }
    myfile3.close();
}

mio::IOResult<void> run(const std::string& input_file, const fs::path& result_dir, size_t num_runs,
                        bool save_single_runs = true)
{

    auto t0               = mio::abm::TimePoint(0); // Start time per simulation
    auto tmax             = mio::abm::TimePoint(0) + mio::abm::days(2); // End time per simulation
    auto ensemble_results = std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of collected results
    ensemble_results.reserve(size_t(num_runs));
    auto run_idx            = size_t(1); // The run index
    auto save_result_result = mio::IOResult<void>(mio::success()); // Variable informing over successful IO operations
    auto max_num_persons    = 1000;

    // Loop over a number of runs
    while (run_idx <= num_runs) {

        // Create the sampled simulation with start time t0.
        auto sim = create_sampled_simulation(input_file, t0, max_num_persons);
        //output object
        mio::History<mio::DataWriterToMemory, mio::abm::LogLocationInformation, mio::abm::LogPersonInformation,
                     mio::abm::LogDataForMobility>
            historyPersonInf;
        mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
            Eigen::Index(mio::abm::InfectionState::Count)};
        mio::History<mio::abm::DataWriterToMemoryDelta, mio::abm::LogDataForMobility> historyPersonInfDelta;
        // Collect the id of location in model.
        std::vector<int> loc_ids;
        for (auto& location : sim.get_model().get_locations()) {
            loc_ids.push_back(location.get_id().get());
        }
        // Advance the model to tmax
        sim.advance(tmax, historyPersonInf, historyTimeSeries, historyPersonInfDelta);
        // TODO: update result of the simulation to be a vector of location result.
        auto temp_sim_result = std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyTimeSeries.get_log())};
        // Push result of the simulation back to the result vector
        ensemble_results.push_back(temp_sim_result);
        // Option to save the current run result to file
        if (save_result_result && save_single_runs) {
            auto result_dir_run = result_dir / ("abm_result_run_" + std::to_string(run_idx) + ".h5");
            save_result_result  = save_result(ensemble_results.back(), loc_ids, 1, result_dir_run.string());
        }
        write_log_to_file_person_and_location_data(historyPersonInf);
        write_log_to_file_trip_data(historyPersonInfDelta);
        ++run_idx;
    }
    BOOST_OUTCOME_TRY(save_result_result);
    return mio::success();
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::warn);

    std::string result_dir = ".";
    std::string input_file = "";
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
