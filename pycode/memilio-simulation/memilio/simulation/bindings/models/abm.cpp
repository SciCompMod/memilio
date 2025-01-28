/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Martin Siggel, Daniel Abele, Martin J. Kuehn, Jan Kleinert, Khoa Nguyen
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
#include "abm/infection_state.h"
#include "abm/location_id.h"
#include "abm/location_type.h"
#include "abm/time.h"
#include "memilio/io/io.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "pybind_util.h"
#include "utils/custom_index_array.h"
#include "utils/parameter_set.h"
#include "utils/index.h"
#include "abm/simulation.h"
#include "abm/household.h"
#include "abm/personal_rng.h"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "memilio/io/hdf5_cpp.h"

#include "pybind11/attr.h"
#include "pybind11/cast.h"
#include "pybind11/pybind11.h"
#include "pybind11/operators.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <map>
#include <string>
#include <type_traits>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <filesystem>
#include <iostream>

namespace py = pybind11;

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
    using Type = std::vector<std::tuple<mio::abm::LocationId, mio::abm::LocationType>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        std::vector<std::tuple<mio::abm::LocationId, mio::abm::LocationType>> location_ids{};
        for (auto&& location : sim.get_model().get_locations()) {
            location_ids.push_back(std::make_tuple(location.get_id(), location.get_type()));
        }
        return location_ids;
    }
};

//AgentId logger
struct LogAgentIds : mio::LogOnce {
    using Type = std::vector<mio::abm::PersonId>;
    static Type log(const mio::abm::Simulation& sim)
    {
        std::vector<mio::abm::PersonId> agent_ids{};
        for (auto&& person : sim.get_model().get_persons()) {
            agent_ids.push_back(person.get_id());
        }
        return agent_ids;
    }
};

//agent logger
struct LogPersonsPerLocationAndInfectionTime : mio::LogAlways {
    using Type = std::vector<std::tuple<mio::abm::LocationId, mio::abm::LocationType, mio::abm::PersonId,
                                        mio::abm::TimeSpan, mio::abm::InfectionState>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        std::vector<std::tuple<mio::abm::LocationId, mio::abm::LocationType, mio::abm::PersonId, mio::abm::TimeSpan,
                               mio::abm::InfectionState>>
            location_ids_person{};
        for (auto&& person : sim.get_model().get_persons()) {
            location_ids_person.push_back(std::make_tuple(person.get_location(), person.get_location_type(),
                                                          person.get_id(), person.get_time_since_transmission(),
                                                          person.get_infection_state(sim.get_time())));
        }
        return location_ids_person;
    }
};

std::pair<double, double> get_my_and_sigma(double mean, double std)
{
    double my    = log(mean * mean / sqrt(mean * mean + std * std));
    double sigma = sqrt(log(1 + std * std / (mean * mean)));
    return {my, sigma};
}

mio::AgeGroup determine_age_group(uint32_t age)
{
    if (age <= 4) {
        return mio::AgeGroup(0);
    }
    else if (age <= 15) {
        return mio::AgeGroup(1);
    }
    else if (age <= 34) {
        return mio::AgeGroup(2);
    }
    else if (age <= 59) {
        return mio::AgeGroup(3);
    }
    else if (age <= 79) {
        return mio::AgeGroup(4);
    }
    else if (age > 79) {
        return mio::AgeGroup(5);
    }
    else {
        return mio::AgeGroup(0);
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

void write_infection_paths(std::string filename, mio::abm::Model& model, mio::abm::TimePoint tmax)
{
    auto file = fopen(filename.c_str(), "w");
    if (file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", filename);
    }
    else {
        fprintf(file, "Agent_id S E I_ns I_sy I_sev I_cri R D\n");
        for (auto& person : model.get_persons()) {
            fprintf(file, "%d ", person.get_id().get());
            if (person.get_infection_state(tmax) == mio::abm::InfectionState::Susceptible) {
                fprintf(file, "%.14f ", tmax.hours());
                for (auto i = 0; i < static_cast<int>(mio::abm::InfectionState::Count); ++i) {
                    fprintf(file, "0 ");
                }
            }
            else {
                auto time_S = std::max(
                    {person.get_infection().get_infection_start() - mio::abm::TimePoint(0), mio::abm::TimeSpan(0)});
                auto time_E    = person.get_infection().get_time_in_state(mio::abm::InfectionState::Exposed);
                auto time_INS  = person.get_infection().get_time_in_state(mio::abm::InfectionState::InfectedNoSymptoms);
                auto time_ISy  = person.get_infection().get_time_in_state(mio::abm::InfectionState::InfectedSymptoms);
                auto time_ISev = person.get_infection().get_time_in_state(mio::abm::InfectionState::InfectedSevere);
                auto time_ICri = person.get_infection().get_time_in_state(mio::abm::InfectionState::InfectedCritical);
                auto time_R    = mio::abm::TimePoint(0);
                auto time_D    = mio::abm::TimePoint(0);
                auto t_Infected = time_E + time_INS + time_ISy + time_ISev + time_ICri;
                if (person.get_infection_state(tmax) == mio::abm::InfectionState::Recovered) {
                    if (time_S.hours() == 0) {
                        time_R =
                            tmax - t_Infected + (person.get_infection().get_infection_start() - mio::abm::TimePoint(0));
                    }
                    else {
                        time_R = tmax - time_S - t_Infected;
                    }
                }
                else if (person.get_infection_state(tmax) == mio::abm::InfectionState::Dead) {
                    if (time_S.hours() == 0) {
                        time_D =
                            tmax - t_Infected + (person.get_infection().get_infection_start() - mio::abm::TimePoint(0));
                    }
                    else {
                        time_D = tmax - time_S - t_Infected;
                    }
                }
                fprintf(file, "%.14f ", time_S.hours());
                fprintf(file, "%.14f ", time_E.hours());
                fprintf(file, "%.14f ", time_INS.hours());
                fprintf(file, "%.14f ", time_ISy.hours());
                fprintf(file, "%.14f ", time_ISev.hours());
                fprintf(file, "%.14f ", time_ICri.hours());
                fprintf(file, "%.14f ", time_R.hours());
                fprintf(file, "%.14f ", time_D.hours());
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }
}

void write_compartments(std::string filename, mio::abm::Model& model,
                        mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                     LogPersonsPerLocationAndInfectionTime, LogAgentIds>& history)
{
    auto file = fopen(filename.c_str(), "w");
    if (file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", filename);
    }
    else {
        auto log = history.get_log();
        auto tps = std::get<0>(log);
        fprintf(file, "t S E Ins Isy Isev Icri R D\n");
        for (auto t = size_t(0); t < tps.size(); ++t) {
            auto tp = mio::abm::TimePoint(0) + mio::abm::hours(t);
            fprintf(file, "%.14f ", tps[t]);
            std::vector<int> comps(static_cast<size_t>(mio::abm::InfectionState::Count));
            for (auto& person : model.get_persons()) {
                auto state = person.get_infection_state(tp);
                comps[static_cast<size_t>(state)] += 1;
            }
            for (auto c : comps) {
                fprintf(file, "%d ", c);
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }
}

std::string loc_to_string(mio::abm::LocationId loc_id, mio::abm::LocationType type)
{
    return "0" + std::to_string(static_cast<int>(type)) + std::to_string(loc_id.get());
}

int time_since_transmission(mio::abm::TimeSpan t)
{
    return t.days() > 1000 ? -1 : t.hours();
}

#ifdef MEMILIO_HAS_HDF5
mio::IOResult<void> write_h5(std::string filename,
                             mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                          LogPersonsPerLocationAndInfectionTime, LogAgentIds>& history)
{
    // create hdf5 file
    mio::H5File file{H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(file.id, mio::StatusCode::FileNotFound, filename);
    // get agents ids
    auto log          = history.get_log();
    auto agent_ids    = std::get<3>(log)[0];
    auto logPerPerson = std::get<2>(log);
    for (auto& id : agent_ids) {
        auto group_name = std::to_string(id.get());
        mio::H5Group agent_h5group{H5Gcreate(file.id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
        MEMILIO_H5_CHECK(agent_h5group.id, mio::StatusCode::UnknownError,
                         "H5Group could not be created for agent (" + group_name + ")");
        const int num_timepoints = static_cast<int>(logPerPerson.size());
        // Dimension for both data sets
        hsize_t dims_t[] = {static_cast<hsize_t>(num_timepoints)};
        // DataSpace for both data sets
        mio::H5DataSpace dspace_t{H5Screate_simple(1, dims_t, NULL)};
        MEMILIO_H5_CHECK(dspace_t.id, mio::StatusCode::UnknownError, "DataSpace could not be created.");
        std::vector<char*> LocIds(num_timepoints);
        std::vector<int> tsm(num_timepoints);
        for (int t = 0; t < num_timepoints; ++t) {
            std::string s =
                loc_to_string(std::get<0>(logPerPerson[t][id.get()]), std::get<1>(logPerPerson[t][id.get()]));
            LocIds[t] = new char[s.size()];
            strcpy(LocIds[t], s.c_str());
            tsm[t] = time_since_transmission(std::get<3>(logPerPerson[t][id.get()]));
        }
        // string dim
        hsize_t str_dim[1] = {LocIds.size()};
        mio::H5DataSpace dspace_str{H5Screate_simple(1, str_dim, NULL)};

        hid_t datatype = H5Tcopy(H5T_C_S1);
        H5Tset_size(datatype, H5T_VARIABLE);

        mio::H5DataSet dset_LocIds{
            H5Dcreate(agent_h5group.id, "loc_ids", datatype, dspace_str.id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
        MEMILIO_H5_CHECK(dset_LocIds.id, mio::StatusCode::UnknownError, "LocId DataSet could not be created (LocIds).");
        MEMILIO_H5_CHECK(H5Dwrite(dset_LocIds.id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, LocIds.data()),
                         mio::StatusCode::UnknownError, "LocId data could not be written.");
        for (int t = 0; t < num_timepoints; ++t) {
            delete[] LocIds[t];
        }
        mio::H5DataSet dset_tsm{H5Dcreate(agent_h5group.id, "time_since_transm", H5T_NATIVE_INT, dspace_t.id,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
        MEMILIO_H5_CHECK(dset_tsm.id, mio::StatusCode::UnknownError,
                         "TST DataSet could not be created (timeSinceTransmission).");
        MEMILIO_H5_CHECK(H5Dwrite(dset_tsm.id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tsm.data()),
                         mio::StatusCode::UnknownError, "TST data could not be written.");
    }
    return mio::success();
}
#endif

void write_mapping_to_file(std::string filename, std::map<int, std::vector<std::string>>& mapping)
{
    auto file = fopen(filename.c_str(), "w");
    if (file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", filename);
    }
    else {
        for (auto it = mapping.begin(); it != mapping.end(); it++) {
            fprintf(file, "%d", it->first);
            for (auto s : it->second) {
                fprintf(file, " %s", s.c_str());
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }
}

void initialize_model(mio::abm::Model& model, std::string person_file, size_t num_hospitals, std::string outfile)
{
    // Mapping of ABM locations to traffic areas/cells
    // - each traffic area is mapped to a vector containing strings with LocationType and LocationId
    std::map<int, std::vector<std::string>> loc_area_mapping;
    // Mapping of traffic data location ids to ABM location ids
    std::map<int, mio::abm::LocationId> home_locations;
    std::map<int, mio::abm::LocationId> shop_locations;
    std::map<int, mio::abm::LocationId> event_locations;
    std::map<int, mio::abm::LocationId> school_locations;
    std::map<int, mio::abm::LocationId> work_locations;
    std::vector<mio::abm::LocationId> hospitals;
    std::vector<mio::abm::LocationId> icus;

    const boost::filesystem::path p = person_file;
    if (!boost::filesystem::exists(p)) {
        mio::log_error("Cannot read in data. File does not exist.");
    }
    // File pointer
    std::fstream fin;

    // Open an existing file
    fin.open(person_file, std::ios::in);
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

    // Create hospitals and ICU
    for (auto i = size_t(0); i < num_hospitals; ++i) {
        auto hosp = model.add_location(mio::abm::LocationType::Hospital);
        auto icu  = model.add_location(mio::abm::LocationType::ICU);
        hospitals.push_back(hosp);
        icus.push_back(icu);
    }

    while (std::getline(fin, line)) {
        row.clear();

        // read columns in this row
        split_line(line, &row);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        uint32_t age = row[index["age"]];

        int home_id   = row[index["home_id"]];
        int home_zone = row[index["home_zone"]];

        mio::abm::LocationId home;

        auto iter_home = home_locations.find(home_id);
        // check whether home location already exists in model
        if (iter_home == home_locations.end()) {
            // if home location does not exists yet, create new location and insert it to mapping
            home = model.add_location(mio::abm::LocationType::Home);
            home_locations.insert({home_id, home});
            std::string loc =
                "0" + std::to_string(static_cast<int>(mio::abm::LocationType::Home)) + std::to_string(home.get());
            auto zone_iter = loc_area_mapping.find(home_zone);
            if (zone_iter == loc_area_mapping.end()) {
                loc_area_mapping.insert({home_zone, {loc}});
            }
            else {
                loc_area_mapping[home_zone].push_back(loc);
            }
        }
        else {
            home = home_locations[home_id];
        }
        // Add person to model and assign home location to it
        auto pid     = model.add_person(home, determine_age_group(age));
        auto& person = model.get_person(pid);
        person.set_assigned_location(mio::abm::LocationType::Home, home);

        int shop_id   = row[index["shop_id"]];
        int shop_zone = row[index["shop_zone"]];

        mio::abm::LocationId shop;

        auto iter_shop = shop_locations.find(shop_id);
        // Check whether shop location already exists in model
        if (iter_shop == shop_locations.end()) {
            // Create shop location and add it to mapping
            shop = model.add_location(mio::abm::LocationType::BasicsShop);
            // Shops with ids -1 are individual locations each
            if (shop_id != -1) {
                shop_locations.insert({shop_id, shop});
            }
            std::string loc =
                "0" + std::to_string(static_cast<int>(mio::abm::LocationType::BasicsShop)) + std::to_string(shop.get());
            auto zone_iter = loc_area_mapping.find(shop_zone);
            if (zone_iter == loc_area_mapping.end()) {
                loc_area_mapping.insert({shop_zone, {loc}});
            }
            else {
                loc_area_mapping[shop_zone].push_back(loc);
            }
        }
        else {
            shop = shop_locations[shop_id];
        }
        // Assign shop to person
        person.set_assigned_location(mio::abm::LocationType::BasicsShop, shop);

        int event_id   = row[index["event_id"]];
        int event_zone = row[index["event_zone"]];

        mio::abm::LocationId event;

        auto iter_event = event_locations.find(event_id);
        // Check whether event location already exists in model
        if (iter_event == event_locations.end()) {
            //Create event location and add it to mapping
            event = model.add_location(mio::abm::LocationType::SocialEvent);
            // Events with id -1 are individual locations each
            if (event_id != -1) {
                event_locations.insert({event_id, event});
            }
            std::string loc = "0" + std::to_string(static_cast<int>(mio::abm::LocationType::SocialEvent)) +
                              std::to_string(event.get());
            auto zone_iter = loc_area_mapping.find(event_zone);
            if (zone_iter == loc_area_mapping.end()) {
                loc_area_mapping.insert({event_zone, {loc}});
            }
            else {
                loc_area_mapping[event_zone].push_back(loc);
            }
        }
        else {
            event = event_locations[event_id];
        }
        // Assign event location to person
        person.set_assigned_location(mio::abm::LocationType::SocialEvent, event);

        // Check if person is school-aged
        if (person.get_age() == mio::AgeGroup(1)) {
            int school_id   = row[index["school_id"]];
            int school_zone = row[index["school_zone"]];

            mio::abm::LocationId school;

            auto iter_school = school_locations.find(school_id);
            // Check whether school location is already in model
            if (iter_school == school_locations.end()) {
                // Add schools locations to model and insert it in mapping
                school = model.add_location(mio::abm::LocationType::School);
                // schools with id -1 are individual locations each
                if (school_id != -1) {
                    school_locations.insert({school_id, school});
                }
                std::string loc = "0" + std::to_string(static_cast<int>(mio::abm::LocationType::School)) +
                                  std::to_string(school.get());
                auto zone_iter = loc_area_mapping.find(school_zone);
                if (zone_iter == loc_area_mapping.end()) {
                    loc_area_mapping.insert({school_zone, {loc}});
                }
                else {
                    loc_area_mapping[school_zone].push_back(loc);
                }
            }
            else {
                school = school_locations[school_id];
            }
            // Assign school location to person
            person.set_assigned_location(mio::abm::LocationType::School, school);
        }
        // Check if person is work-aged
        if (person.get_age() == mio::AgeGroup(2) || person.get_age() == mio::AgeGroup(3)) {
            int work_id   = row[index["work_id"]];
            int work_zone = row[index["work_zone"]];

            mio::abm::LocationId work;

            auto iter_work = work_locations.find(work_id);
            // Check whether work location already exists in model
            if (iter_work == work_locations.end()) {
                // Add work location to model and insert it in mapping
                work = model.add_location(mio::abm::LocationType::Work);
                // Locations with id -1 are individual locations each
                if (work_id != -1) {
                    work_locations.insert({work_id, work});
                }
                std::string loc =
                    "0" + std::to_string(static_cast<int>(mio::abm::LocationType::Work)) + std::to_string(work.get());
                auto zone_iter = loc_area_mapping.find(work_zone);
                if (zone_iter == loc_area_mapping.end()) {
                    loc_area_mapping.insert({work_zone, {loc}});
                }
                else {
                    loc_area_mapping[work_zone].push_back(loc);
                }
            }
            else {
                work = work_locations[work_id];
            }
            // Assign work location to person
            person.set_assigned_location(mio::abm::LocationType::Work, work);
        }
        // Assign Hospital and ICU
        std::vector<mio::abm::LocationId>::iterator randItHosp = hospitals.begin();
        std::advance(randItHosp, std::rand() % hospitals.size());
        person.set_assigned_location(mio::abm::LocationType::Hospital, *randItHosp);
        std::vector<mio::abm::LocationId>::iterator randItIcu = icus.begin();
        std::advance(randItIcu, std::rand() % icus.size());
        person.set_assigned_location(mio::abm::LocationType::ICU, *randItIcu);
    }

    // Add hospitals to Mapping
    for (size_t i = size_t(0); i < hospitals.size(); ++i) {
        auto it = loc_area_mapping.begin();
        std::advance(it, rand() % loc_area_mapping.size());
        loc_area_mapping[it->first].push_back("0" + std::to_string(static_cast<int>(mio::abm::LocationType::Hospital)) +
                                              std::to_string(hospitals[i].get()));
        loc_area_mapping[it->first].push_back("0" + std::to_string(static_cast<int>(mio::abm::LocationType::ICU)) +
                                              std::to_string(icus[i].get()));
    }
    write_mapping_to_file(outfile, loc_area_mapping);
}

PYBIND11_MODULE(_simulation_abm, m)
{
    pymio::iterable_enum<mio::abm::InfectionState>(m, "InfectionState")
        .value("Susceptible", mio::abm::InfectionState::Susceptible)
        .value("Exposed", mio::abm::InfectionState::Exposed)
        .value("InfectedNoSymptoms", mio::abm::InfectionState::InfectedNoSymptoms)
        .value("InfectedSymptoms", mio::abm::InfectionState::InfectedSymptoms)
        .value("InfectedSevere", mio::abm::InfectionState::InfectedSevere)
        .value("InfectedCritical", mio::abm::InfectionState::InfectedCritical)
        .value("Recovered", mio::abm::InfectionState::Recovered)
        .value("Dead", mio::abm::InfectionState::Dead)
        .value("Count", mio::abm::InfectionState::Count);

    pymio::iterable_enum<mio::abm::ProtectionType>(m, "ProtectionType")
        .value("NoProtection", mio::abm::ProtectionType::NoProtection)
        .value("NaturalInfection", mio::abm::ProtectionType::NaturalInfection)
        .value("GenericVaccine", mio::abm::ProtectionType::GenericVaccine);

    pymio::iterable_enum<mio::abm::VirusVariant>(m, "VirusVariant").value("Wildtype", mio::abm::VirusVariant::Wildtype);

    pymio::iterable_enum<mio::abm::LocationType>(m, "LocationType")
        .value("Home", mio::abm::LocationType::Home)
        .value("School", mio::abm::LocationType::School)
        .value("Work", mio::abm::LocationType::Work)
        .value("SocialEvent", mio::abm::LocationType::SocialEvent)
        .value("BasicsShop", mio::abm::LocationType::BasicsShop)
        .value("Hospital", mio::abm::LocationType::Hospital)
        .value("ICU", mio::abm::LocationType::ICU)
        .value("Car", mio::abm::LocationType::Car)
        .value("PublicTransport", mio::abm::LocationType::PublicTransport)
        .value("TransportWithoutContact", mio::abm::LocationType::TransportWithoutContact)
        .value("Cemetery", mio::abm::LocationType::Cemetery);

    pymio::iterable_enum<mio::abm::TestType>(m, "TestType")
        .value("Generic", mio::abm::TestType::Generic)
        .value("Antigen", mio::abm::TestType::Antigen)
        .value("PCR", mio::abm::TestType::PCR);

    pymio::bind_class<mio::abm::TimeSpan, pymio::EnablePickling::Never>(m, "TimeSpan")
        .def(py::init<int>(), py::arg("seconds") = 0)
        .def_property_readonly("seconds", &mio::abm::TimeSpan::seconds)
        .def_property_readonly("hours", &mio::abm::TimeSpan::hours)
        .def_property_readonly("days", &mio::abm::TimeSpan::days)
        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self -= py::self)
        .def(py::self * int{})
        .def(py::self *= int{})
        .def(py::self / int{})
        .def(py::self /= int{})
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self < py::self)
        .def(py::self <= py::self)
        .def(py::self > py::self)
        .def(py::self <= py::self);

    m.def("seconds", &mio::abm::seconds);
    m.def("minutes", &mio::abm::minutes);
    m.def("hours", &mio::abm::hours);
    m.def("days", py::overload_cast<int>(&mio::abm::days));

    pymio::bind_class<mio::abm::TimePoint, pymio::EnablePickling::Never>(m, "TimePoint")
        .def(py::init<int>(), py::arg("seconds") = 0)
        .def_property_readonly("seconds", &mio::abm::TimePoint::seconds)
        .def_property_readonly("days", &mio::abm::TimePoint::days)
        .def_property_readonly("hours", &mio::abm::TimePoint::hours)
        .def_property_readonly("day_of_week", &mio::abm::TimePoint::day_of_week)
        .def_property_readonly("hour_of_day", &mio::abm::TimePoint::hour_of_day)
        .def_property_readonly("time_since_midnight", &mio::abm::TimePoint::time_since_midnight)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self < py::self)
        .def(py::self <= py::self)
        .def(py::self > py::self)
        .def(py::self >= py::self)
        .def(py::self - py::self)
        .def(py::self + mio::abm::TimeSpan{})
        .def(py::self += mio::abm::TimeSpan{})
        .def(py::self - mio::abm::TimeSpan{})
        .def(py::self -= mio::abm::TimeSpan{});

    pymio::bind_class<mio::abm::TestParameters, pymio::EnablePickling::Never>(m, "TestParameters")
        .def(py::init<double, double, mio::abm::TimeSpan, mio::abm::TestType>())
        .def_readwrite("sensitivity", &mio::abm::TestParameters::sensitivity)
        .def_readwrite("specificity", &mio::abm::TestParameters::specificity)
        .def_readwrite("required_time", &mio::abm::TestParameters::required_time)
        .def_readwrite("type", &mio::abm::TestParameters::type);

    pymio::bind_CustomIndexArray<mio::UncertainValue<double>, mio::abm::VirusVariant, mio::AgeGroup>(
        m, "_AgeParameterArray");
    pymio::bind_CustomIndexArray<mio::abm::TestParameters, mio::abm::TestType>(m, "_TestData");
    pymio::bind_Index<mio::abm::ProtectionType>(m, "ProtectionTypeIndex");
    pymio::bind_ParameterSet<mio::abm::ParametersBase, pymio::EnablePickling::Never>(m, "ParametersBase");
    pymio::bind_class<mio::abm::Parameters, pymio::EnablePickling::Never, mio::abm::ParametersBase>(m, "Parameters")
        .def(py::init<int>())
        .def("check_constraints", &mio::abm::Parameters::check_constraints);

    pymio::bind_ParameterSet<mio::abm::LocalInfectionParameters, pymio::EnablePickling::Never>(
        m, "LocalInfectionParameters")
        .def(py::init<size_t>());

    pymio::bind_class<mio::abm::LocationId, pymio::EnablePickling::Never>(m, "LocationId")
        .def(py::init<uint32_t>(), py::arg("id"))
        .def("index", &mio::abm::LocationId::get);

    pymio::bind_class<mio::abm::PersonId, pymio::EnablePickling::Never>(m, "PersonId")
        .def(py::init<uint32_t>(), py::arg("id"))
        .def("index", &mio::abm::PersonId::get);

    pymio::bind_class<mio::abm::Person, pymio::EnablePickling::Never>(m, "Person")
        .def("set_assigned_location",
             py::overload_cast<mio::abm::LocationType, mio::abm::LocationId>(&mio::abm::Person::set_assigned_location))
        .def("add_new_infection",
             [](mio::abm::Person& self, mio::abm::Infection& infection, mio::abm::TimePoint t) {
                 self.add_new_infection(std::move(infection), t);
             })
        .def("assigned_location",
             [](mio::abm::Person& self, mio::abm::LocationType type) {
                 return self.get_assigned_location(type);
             })
        .def("infection_state",
             [](mio::abm::Person& self, mio::abm::TimePoint t) {
                 return self.get_infection_state(t);
             })
        .def_property_readonly("infection", py::overload_cast<>(&mio::abm::Person::get_infection, py::const_))
        .def_property_readonly("location", py::overload_cast<>(&mio::abm::Person::get_location, py::const_))
        .def_property_readonly("age", &mio::abm::Person::get_age)
        .def_property_readonly("id", &mio::abm::Person::get_id)
        .def_property_readonly("is_in_quarantine", &mio::abm::Person::is_in_quarantine);

    pymio::bind_class<mio::abm::HouseholdMember, pymio::EnablePickling::Never>(m, "HouseholdMember")
        .def(py::init<int>())
        .def("set_age_weight", &mio::abm::HouseholdMember::set_age_weight);

    pymio::bind_class<mio::abm::Household, pymio::EnablePickling::Never>(m, "Household")
        .def(py::init<>())
        .def("add_members", &mio::abm::Household::add_members);

    m.def("add_household_group_to_model", &mio::abm::add_household_group_to_model);

    pymio::bind_class<mio::abm::HouseholdGroup, pymio::EnablePickling::Never>(m, "HouseholdGroup")
        .def(py::init<>())
        .def("add_households", &mio::abm::HouseholdGroup::add_households);

    pymio::bind_class<mio::abm::TestingCriteria, pymio::EnablePickling::Never>(m, "TestingCriteria")
        .def(py::init<const std::vector<mio::AgeGroup>&, const std::vector<mio::abm::InfectionState>&>(),
             py::arg("age_groups"), py::arg("infection_states"));

    pymio::bind_class<mio::abm::TestingScheme, pymio::EnablePickling::Never>(m, "TestingScheme")
        .def(py::init<const mio::abm::TestingCriteria&, mio::abm::TimeSpan, mio::abm::TimePoint, mio::abm::TimePoint,
                      const mio::abm::TestParameters&, double>(),
             py::arg("testing_criteria"), py::arg("testing_validity_period"), py::arg("start_date"),
             py::arg("end_date"), py::arg("test_parameters"), py::arg("probability"))
        .def_property_readonly("active", &mio::abm::TestingScheme::is_active);

    pymio::bind_class<mio::abm::ProtectionEvent, pymio::EnablePickling::Never>(m, "ProtectionEvent")
        .def(py::init<mio::abm::ProtectionType, mio::abm::TimePoint>(), py::arg("type"), py::arg("time"))
        .def_readwrite("type", &mio::abm::ProtectionEvent::type)
        .def_readwrite("time", &mio::abm::ProtectionEvent::time);

    pymio::bind_class<mio::abm::TestingStrategy, pymio::EnablePickling::Never>(m, "TestingStrategy")
        .def(py::init<const std::vector<mio::abm::TestingStrategy::LocalStrategy>&>());

    pymio::bind_class<mio::abm::Infection, pymio::EnablePickling::Never>(m, "Infection")
        .def(py::init([](mio::abm::Model& model, mio::abm::Person& person, mio::abm::VirusVariant variant,
                         mio::abm::TimePoint start_date, mio::abm::InfectionState start_state, bool detected) {
            auto rng = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
            return mio::abm::Infection(rng, variant, person.get_age(), model.parameters, start_date, start_state,
                                       person.get_latest_protection(), detected);
        }))
        .def("get_infection_start", &mio::abm::Infection::get_infection_start)
        .def("get_time_in_state", [](mio::abm::Infection& self, mio::abm::InfectionState state) {
            return self.get_time_in_state(state);
        });

    pymio::bind_class<mio::abm::Location, pymio::EnablePickling::Never>(m, "Location")
        .def("set_capacity", &mio::abm::Location::set_capacity)
        .def_property_readonly("type", &mio::abm::Location::get_type)
        .def_property_readonly("id", &mio::abm::Location::get_id)
        .def_property("infection_parameters",
                      py::overload_cast<>(&mio::abm::Location::get_infection_parameters, py::const_),
                      [](mio::abm::Location& self, mio::abm::LocalInfectionParameters params) {
                          self.get_infection_parameters() = params;
                      });

    //copying and moving of ranges enabled below, see PYMIO_IGNORE_VALUE_TYPE
    pymio::bind_Range<decltype(std::declval<const mio::abm::Model>().get_locations())>(m, "_ModelLocationsRange");
    pymio::bind_Range<decltype(std::declval<const mio::abm::Model>().get_persons())>(m, "_ModelPersonsRange");

    pymio::bind_class<mio::abm::Trip, pymio::EnablePickling::Never>(m, "Trip")
        .def(py::init<uint32_t, mio::abm::TimePoint, mio::abm::LocationId, mio::abm::LocationId,
                      std::vector<uint32_t>>(),
             py::arg("person_id"), py::arg("time"), py::arg("destination"), py::arg("origin"),
             py::arg("cells") = std::vector<uint32_t>())
        .def_readwrite("person_id", &mio::abm::Trip::person_id)
        .def_readwrite("time", &mio::abm::Trip::time)
        .def_readwrite("destination", &mio::abm::Trip::destination)
        .def_readwrite("origin", &mio::abm::Trip::origin)
        .def_readwrite("cells", &mio::abm::Trip::cells);

    pymio::bind_class<mio::abm::TripList, pymio::EnablePickling::Never>(m, "TripList")
        .def(py::init<>())
        .def("add_trip", &mio::abm::TripList::add_trip, py::arg("trip"), py::arg("weekend") = false)
        .def("next_trip", &mio::abm::TripList::get_next_trip, py::arg("weekend") = false)
        .def("num_trips", &mio::abm::TripList::num_trips, py::arg("weekend") = false);

    pymio::bind_class<mio::abm::Model, pymio::EnablePickling::Never>(m, "Model")
        .def(py::init<int32_t>())
        .def("add_location", &mio::abm::Model::add_location, py::arg("location_type"), py::arg("num_cells") = 1)
        .def("add_person", py::overload_cast<mio::abm::LocationId, mio::AgeGroup>(&mio::abm::Model::add_person),
             py::arg("location_id"), py::arg("age_group"))
        .def("assign_location", &mio::abm::Model::assign_location, py::arg("person_id"), py::arg("location_id"))
        .def("get_location", py::overload_cast<mio::abm::LocationId>(&mio::abm::Model::get_location, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_rng", &mio::abm::Model::get_rng, py::return_value_policy::reference_internal)
        .def_property_readonly("locations", py::overload_cast<>(&mio::abm::Model::get_locations, py::const_),
                               py::keep_alive<1, 0>{}) //keep this model alive while contents are referenced in ranges
        .def_property_readonly("persons", py::overload_cast<>(&mio::abm::Model::get_persons, py::const_),
                               py::keep_alive<1, 0>{})
        .def_property(
            "trip_list", py::overload_cast<>(&mio::abm::Model::get_trip_list),
            [](mio::abm::Model& self, const mio::abm::TripList& list) {
                self.get_trip_list() = list;
            },
            py::return_value_policy::reference_internal)
        .def_property("use_mobility_rules", py::overload_cast<>(&mio::abm::Model::use_mobility_rules, py::const_),
                      py::overload_cast<bool>(&mio::abm::Model::use_mobility_rules))
        .def_readwrite("parameters", &mio::abm::Model::parameters)
        .def_property(
            "testing_strategy", py::overload_cast<>(&mio::abm::Model::get_testing_strategy, py::const_),
            [](mio::abm::Model& self, mio::abm::TestingStrategy strategy) {
                self.get_testing_strategy() = strategy;
            },
            py::return_value_policy::reference_internal);

    pymio::bind_class<mio::abm::Simulation, pymio::EnablePickling::Never>(m, "Simulation")
        .def(py::init<mio::abm::TimePoint, size_t>())
        .def("advance",
             &mio::abm::Simulation::advance<mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                                         LogPersonsPerLocationAndInfectionTime, LogAgentIds>>)
        .def("advance",
             static_cast<void (mio::abm::Simulation::*)(mio::abm::TimePoint)>(&mio::abm::Simulation::advance),
             py::arg("tmax"))
        .def_property_readonly("model", py::overload_cast<>(&mio::abm::Simulation::get_model));

    pymio::bind_class<mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                   LogPersonsPerLocationAndInfectionTime, LogAgentIds>,
                      pymio::EnablePickling::Never>(m, "History")
        .def(py::init<>())
        .def_property_readonly("log", [](mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                                      LogPersonsPerLocationAndInfectionTime, LogAgentIds>& self) {
            return self.get_log();
        });

    m.def(
        "set_viral_load_parameters",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double min_peak,
           double max_peak, double min_incline, double max_incline, double min_decline, double max_decline) {
            infection_params.get<mio::abm::ViralLoadDistributions>()[{variant, age}] =
                mio::abm::ViralLoadDistributionsParameters{
                    {min_peak, max_peak}, {min_incline, max_incline}, {min_decline, max_decline}};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_incubationPeriod",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto incubation_period_params                                      = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::IncubationPeriod>()[{variant, age}] = {incubation_period_params.first,
                                                                                  incubation_period_params.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedNoSymptomsToSymptoms",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedNoSymptomsToSymptoms = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{variant, age}] = {
                TimeInfectedNoSymptomsToSymptoms.first, TimeInfectedNoSymptomsToSymptoms.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedNoSymptomsToRecovered",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedNoSymptomsToRecovered = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{variant, age}] = {
                TimeInfectedNoSymptomsToRecovered.first, TimeInfectedNoSymptomsToRecovered.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedSymptomsToSevere",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedSymptomsToSevere = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedSymptomsToSevere>()[{variant, age}] = {
                TimeInfectedSymptomsToSevere.first, TimeInfectedSymptomsToSevere.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedSymptomsToRecovered",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedSymptomsToRecovered = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedSymptomsToRecovered>()[{variant, age}] = {
                TimeInfectedSymptomsToRecovered.first, TimeInfectedSymptomsToRecovered.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedSevereToRecovered",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedSevereToRecovered = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedSevereToRecovered>()[{variant, age}] = {
                TimeInfectedSevereToRecovered.first, TimeInfectedSevereToRecovered.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedSevereToCritical",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedSevereToCritical = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedSevereToCritical>()[{variant, age}] = {
                TimeInfectedSevereToCritical.first, TimeInfectedSevereToCritical.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedCriticalToRecovered",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedCriticalToRecovered = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedCriticalToRecovered>()[{variant, age}] = {
                TimeInfectedCriticalToRecovered.first, TimeInfectedCriticalToRecovered.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedCriticalToDead",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedCriticalToDead                                              = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedCriticalToDead>()[{variant, age}] = {
                TimeInfectedCriticalToDead.first, TimeInfectedCriticalToDead.second};
        },
        py::return_value_policy::reference_internal);

    m.def("set_AgeGroupGoToSchool", [](mio::abm::Parameters& infection_params, mio::AgeGroup age) {
        infection_params.get<mio::abm::AgeGroupGotoSchool>()[age] = true;
    });

    m.def("set_AgeGroupGoToWork", [](mio::abm::Parameters& infection_params, mio::AgeGroup age) {
        infection_params.get<mio::abm::AgeGroupGotoWork>()[age] = true;
    });

    m.def(
        "set_infectivity_parameters",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double min_alpha,
           double max_alpha, double min_beta, double max_beta) {
            infection_params.get<mio::abm::InfectivityDistributions>()[{variant, age}] =
                mio::abm::InfectivityDistributionsParameters{{min_alpha, max_alpha}, {min_beta, max_beta}};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_seeds",
        [](mio::abm::Model& model, int seed) {
            auto rng = mio::RandomNumberGenerator();
            rng.seed({static_cast<uint32_t>(seed)});
            model.get_rng() = rng;
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_log_level_warn",
        []() {
            mio::set_log_level(mio::LogLevel::warn);
        },
        py::return_value_policy::reference_internal);

    m.def("initialize_model", &initialize_model, py::return_value_policy::reference_internal);
    m.def("save_infection_paths", &write_infection_paths, py::return_value_policy::reference_internal);
    m.def("save_comp_output", &write_compartments, py::return_value_policy::reference_internal);

#ifdef MEMILIO_HAS_HDF5
    m.def(
        "write_h5",
        [](std::string filename, mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                              LogPersonsPerLocationAndInfectionTime, LogAgentIds>& history) {
            auto res = write_h5(filename, history);
            if (res == mio::success()) {
                return 0;
            }
            else {
                return 1;
            }
        },
        py::return_value_policy::reference_internal);
#endif

    m.attr("__version__") = "dev";
}

PYMIO_IGNORE_VALUE_TYPE(decltype(std::declval<mio::abm::Model>().get_locations()))
PYMIO_IGNORE_VALUE_TYPE(decltype(std::declval<mio::abm::Model>().get_persons()))
