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
#include "H5Tpublic.h"
#include "H5public.h"
#include "abm/infection_state.h"
#include "abm/location_id.h"
#include "abm/location_type.h"
#include "abm/time.h"
#include "memilio/io/io.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"
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
#include "munich_postprocessing/output_processing.h"

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
#include <tuple>
#include <type_traits>
#include <utility>
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
                                        mio::abm::TimeSpan, mio::abm::InfectionState, int>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        std::vector<std::tuple<mio::abm::LocationId, mio::abm::LocationType, mio::abm::PersonId, mio::abm::TimeSpan,
                               mio::abm::InfectionState, int>>
            location_ids_person{};
        for (auto&& person : sim.get_model().get_persons()) {
            int ww_id = sim.get_model().get_location(person.get_location()).get_wastewater_id();
            location_ids_person.push_back(std::make_tuple(person.get_location(), person.get_location_type(),
                                                          person.get_id(), person.get_time_since_transmission(),
                                                          person.get_infection_state(sim.get_time()), ww_id));
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
    return t.days() > 500 ? -1 : t.hours();
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

mio::IOResult<void> write_h5_v2(std::string filename,
                                mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                             LogPersonsPerLocationAndInfectionTime, LogAgentIds>& history)
{
    //Option 1:
    // One hdf5 group per agent. For every agent there three datasets.
    // 1st dataset: is a vector where even indices a the time point and odd indices the corresponding time since transmission.
    // 2nd dataset: is a vector with the location change tps
    // 3rd dataset: is a vector with (new) location ids

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
        // get time since transmission
        std::vector<int> tsm(num_timepoints);
        std::transform(logPerPerson.begin(), logPerPerson.end(), tsm.begin(), [&id](const auto& entry) {
            return time_since_transmission(std::get<3>(entry[id.get()]));
        });
        std::vector<int> tsm_result_vec;
        bool is_infected = false;
        // Add time since transmission at the beginning of the simulation
        tsm_result_vec.push_back(0);
        tsm_result_vec.push_back(tsm[0]);
        if (tsm[0] > -1) {
            is_infected = true;
        }
        // Add time points of transmission and recovery
        for (size_t t = 1; t < tsm.size(); ++t) {
            if (!is_infected) {
                //time point of transmission
                if (tsm[t] > -1) {
                    tsm_result_vec.push_back(t);
                    tsm_result_vec.push_back(tsm[t]);
                    is_infected = true;
                }
            }
            else {
                //time point of recovery
                if (tsm[t] == -1) {
                    tsm_result_vec.push_back(t);
                    tsm_result_vec.push_back(tsm[t]);
                    is_infected = false;
                }
            }
        }
        // Dimension for 1st data set
        hsize_t dims_d1[] = {static_cast<hsize_t>(tsm_result_vec.size())};
        // Get locations
        std::vector<int> loc_tps_result_vec;
        std::vector<char*> loc_ids_result_vec;
        // Add location at the beginning of the simulation
        loc_tps_result_vec.push_back(0);
        std::string loc = loc_to_string(std::get<0>(logPerPerson[0][id.get()]), std::get<1>(logPerPerson[0][id.get()]));
        loc_ids_result_vec.push_back(new char[loc.size()]);
        strcpy(loc_ids_result_vec[loc_ids_result_vec.size() - 1], loc.c_str());
        // Iterate over all locations and only add change points
        for (size_t t = 1; t < num_timepoints; ++t) {
            std::string loc_new =
                loc_to_string(std::get<0>(logPerPerson[t][id.get()]), std::get<1>(logPerPerson[t][id.get()]));
            if (loc_new != loc) {
                loc_tps_result_vec.push_back(t);
                loc_ids_result_vec.push_back(new char[loc_new.size()]);
                strcpy(loc_ids_result_vec[loc_ids_result_vec.size() - 1], loc_new.c_str());
            }
            loc = loc_new;
        }
        // Dimension for 2nd and 3rd data set
        hsize_t dims_d2d3[1] = {loc_ids_result_vec.size()};
        // DataSpace for data set 1
        mio::H5DataSpace dspace_d1{H5Screate_simple(1, dims_d1, NULL)};
        MEMILIO_H5_CHECK(dspace_d1.id, mio::StatusCode::UnknownError, "V2: DataSpace 1 could not be created.");
        // DataSpace for data set 2
        mio::H5DataSpace dspace_d2{H5Screate_simple(1, dims_d2d3, NULL)};
        MEMILIO_H5_CHECK(dspace_d2.id, mio::StatusCode::UnknownError, "V2: DataSpace 2 could not be created.");
        // DataSpace for data set 3
        mio::H5DataSpace dspace_d3{H5Screate_simple(1, dims_d2d3, NULL)};
        MEMILIO_H5_CHECK(dspace_d3.id, mio::StatusCode::UnknownError, "V2: DataSpace 3 could not be created.");
        // Data type for loc ids
        hid_t datatype = H5Tcopy(H5T_C_S1);
        H5Tset_size(datatype, H5T_VARIABLE);
        // Add data set 1 to HDF5 file
        mio::H5DataSet dset_tsm{H5Dcreate(agent_h5group.id, "time_since_transm", H5T_NATIVE_INT, dspace_d1.id,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
        MEMILIO_H5_CHECK(dset_tsm.id, mio::StatusCode::UnknownError,
                         "V2: Time since transmission DataSet could not be created.");
        MEMILIO_H5_CHECK(H5Dwrite(dset_tsm.id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tsm_result_vec.data()),
                         mio::StatusCode::UnknownError, "V2: Time since transmission data could not be written.");
        // Add data set 2 to HDF5 file
        mio::H5DataSet dset_loc_tps{H5Dcreate(agent_h5group.id, "loc_tps", H5T_NATIVE_INT, dspace_d2.id, H5P_DEFAULT,
                                              H5P_DEFAULT, H5P_DEFAULT)};
        MEMILIO_H5_CHECK(dset_loc_tps.id, mio::StatusCode::UnknownError, "V2: Loc tps DataSet could not be created.");
        MEMILIO_H5_CHECK(
            H5Dwrite(dset_loc_tps.id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, loc_tps_result_vec.data()),
            mio::StatusCode::UnknownError, "V2: Loc tps data could not be written.");
        // Add data set 3 to HDF5 file
        mio::H5DataSet dset_loc_ids{
            H5Dcreate(agent_h5group.id, "loc_ids", datatype, dspace_d3.id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
        MEMILIO_H5_CHECK(dset_loc_ids.id, mio::StatusCode::UnknownError, "V2: Loc ids DataSet could not be created.");
        MEMILIO_H5_CHECK(H5Dwrite(dset_loc_ids.id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, loc_ids_result_vec.data()),
                         mio::StatusCode::UnknownError, "V2: Loc ids data could not be written.");
    }
    return mio::success();
}

mio::IOResult<void> write_h5_v3(std::string filename,
                                mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                             LogPersonsPerLocationAndInfectionTime, LogAgentIds>& history)
{
    //Option 2:
    // hdf5 file contains two datasets (no groups)
    // 1st dataset: #agents x 2 matrix with entry (i,0) tp of transmission and entry (i,1) tp of recovery
    // 2nd dataset: #agents x #tps with entry (i,t) loc at t

    // create hdf5 file
    mio::H5File file{H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(file.id, mio::StatusCode::FileNotFound, filename);

    // get agents ids
    auto log          = history.get_log();
    auto agent_ids    = std::get<3>(log)[0];
    auto logPerPerson = std::get<2>(log);

    //Create group for all data sets
    mio::H5Group h5group{H5Gcreate(file.id, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(h5group.id, mio::StatusCode::UnknownError, "H5Group \"data\" could not be created ");

    const int num_agents = static_cast<int>(agent_ids.size());
    // Dimension and matrix for 1st data set
    hsize_t dims_d1[] = {static_cast<hsize_t>(num_agents), static_cast<hsize_t>(2)};
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tsm_matrix(num_agents, 2);

    const int num_timepoints = static_cast<int>(logPerPerson.size());
    // Dimension and matrix for 2nd data set
    hsize_t dims_d2[] = {static_cast<hsize_t>(num_agents), static_cast<hsize_t>(num_timepoints)};
    Eigen::Matrix<char*, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> loc_matrix(num_agents, num_timepoints);
    // Fill matrices
    for (auto& id : agent_ids) {
        // initialize time point transmission and recovery with Null
        tsm_matrix(id.get(), 0)  = 100000;
        tsm_matrix(id.get(), 1)  = 100000;
        bool tp_transmission_set = false;
        bool tp_recovery_set     = false;
        for (int t = 0; t < num_timepoints; ++t) {
            // set time point of transmission
            if (!tp_transmission_set && time_since_transmission(std::get<3>(logPerPerson[t][id.get()])) > -1) {
                if (t == 0) {
                    tsm_matrix(id.get(), 0) = -time_since_transmission(std::get<3>(logPerPerson[t][id.get()]));
                }
                else {
                    tsm_matrix(id.get(), 0) = t + time_since_transmission(std::get<3>(logPerPerson[t][id.get()]));
                }
                tp_transmission_set = true;
            }
            // set time point of recovery
            if (tp_transmission_set && !tp_recovery_set &&
                time_since_transmission(std::get<3>(logPerPerson[t][id.get()])) == -1) {
                tsm_matrix(id.get(), 1) = t;
                tp_recovery_set         = true;
            }
            std::string s =
                loc_to_string(std::get<0>(logPerPerson[t][id.get()]), std::get<1>(logPerPerson[t][id.get()]));
            loc_matrix(id.get(), t) = new char[s.size()];
            strcpy(loc_matrix(id.get(), t), s.c_str());
        }
    }
    // DataSpace for data set 1
    mio::H5DataSpace dspace_d1{H5Screate_simple(2, dims_d1, NULL)};
    MEMILIO_H5_CHECK(dspace_d1.id, mio::StatusCode::UnknownError, "V3: DataSpace 1 could not be created.");
    // DataSpace for data set 2
    mio::H5DataSpace dspace_d2{H5Screate_simple(2, dims_d2, NULL)};
    MEMILIO_H5_CHECK(dspace_d2.id, mio::StatusCode::UnknownError, "V3: DataSpace 2 could not be created.");
    // Add data set 1 to HDF5 file
    mio::H5DataSet dset_t_r_tp{H5Dcreate(h5group.id, "transm_recovery_tp", H5T_NATIVE_DOUBLE, dspace_d1.id, H5P_DEFAULT,
                                         H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(dset_t_r_tp.id, mio::StatusCode::UnknownError,
                     "V3: Transmission and recovery tp dataSet could not be created.");
    MEMILIO_H5_CHECK(H5Dwrite(dset_t_r_tp.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tsm_matrix.data()),
                     mio::StatusCode::UnknownError, "V3: Transmission and recovery tp data could not be written.");
    // Add data set 2 to HDF5 file
    hid_t datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, H5T_VARIABLE);
    mio::H5DataSet dset_loc_ids{
        H5Dcreate(h5group.id, "loc_ids", datatype, dspace_d2.id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(dset_loc_ids.id, mio::StatusCode::UnknownError, "V3: Loc ids DataSet could not be created.");
    MEMILIO_H5_CHECK(H5Dwrite(dset_loc_ids.id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, loc_matrix.data()),
                     mio::StatusCode::UnknownError, "V3: Loc ids data could not be written.");
    return mio::success();
}

mio::IOResult<void> write_h5_v4(std::string filename,
                                mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                             LogPersonsPerLocationAndInfectionTime, LogAgentIds>& history)
{
    //Option 2:
    // hdf5 file contains two datasets (no groups)
    // 1st dataset: #agents x 2 matrix with entry (i,0) tp of transmission and entry (i,1) tp of recovery
    // 2nd dataset: #agents x #tps with entry (i,t) TAN area ID at t. For locations that are not in any TAN area, the ID is 0.

    // create hdf5 file
    mio::H5File file{H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(file.id, mio::StatusCode::FileNotFound, filename);

    // get agents ids
    auto log          = history.get_log();
    auto agent_ids    = std::get<3>(log)[0];
    auto logPerPerson = std::get<2>(log);

    //Create group for all data sets
    mio::H5Group h5group{H5Gcreate(file.id, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(h5group.id, mio::StatusCode::UnknownError, "H5Group \"data\" could not be created ");

    const int num_agents = static_cast<int>(agent_ids.size());
    // Dimension and matrix for 1st data set
    hsize_t dims_d1[] = {static_cast<hsize_t>(num_agents), static_cast<hsize_t>(2)};
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tsm_matrix(num_agents, 2);

    const int num_timepoints = static_cast<int>(logPerPerson.size());
    // Dimension and matrix for 2nd data set
    hsize_t dims_d2[] = {static_cast<hsize_t>(num_agents), static_cast<hsize_t>(num_timepoints)};
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> loc_matrix(num_agents, num_timepoints);
    // Fill matrices
    for (auto& id : agent_ids) {
        // initialize time point transmission and recovery with Null
        tsm_matrix(id.get(), 0)  = 100000;
        tsm_matrix(id.get(), 1)  = 100000;
        bool tp_transmission_set = false;
        bool tp_recovery_set     = false;
        for (int t = 0; t < num_timepoints; ++t) {
            // set time point of transmission
            if (!tp_transmission_set && time_since_transmission(std::get<3>(logPerPerson[t][id.get()])) > -1) {
                if (t == 0) {
                    tsm_matrix(id.get(), 0) = -time_since_transmission(std::get<3>(logPerPerson[t][id.get()]));
                }
                else {
                    tsm_matrix(id.get(), 0) = t + time_since_transmission(std::get<3>(logPerPerson[t][id.get()]));
                }
                tp_transmission_set = true;
            }
            // set time point of recovery
            if (tp_transmission_set && !tp_recovery_set &&
                time_since_transmission(std::get<3>(logPerPerson[t][id.get()])) == -1) {
                tsm_matrix(id.get(), 1) = t;
                tp_recovery_set         = true;
            }
            loc_matrix(id.get(), t) = std::get<5>(logPerPerson[t][id.get()]);
        }
    }
    // DataSpace for data set 1
    mio::H5DataSpace dspace_d1{H5Screate_simple(2, dims_d1, NULL)};
    MEMILIO_H5_CHECK(dspace_d1.id, mio::StatusCode::UnknownError, "V3: DataSpace 1 could not be created.");
    // DataSpace for data set 2
    mio::H5DataSpace dspace_d2{H5Screate_simple(2, dims_d2, NULL)};
    MEMILIO_H5_CHECK(dspace_d2.id, mio::StatusCode::UnknownError, "V3: DataSpace 2 could not be created.");
    // Add data set 1 to HDF5 file
    mio::H5DataSet dset_t_r_tp{H5Dcreate(h5group.id, "transm_recovery_tp", H5T_NATIVE_DOUBLE, dspace_d1.id, H5P_DEFAULT,
                                         H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(dset_t_r_tp.id, mio::StatusCode::UnknownError,
                     "V3: Transmission and recovery tp dataSet could not be created.");
    MEMILIO_H5_CHECK(H5Dwrite(dset_t_r_tp.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tsm_matrix.data()),
                     mio::StatusCode::UnknownError, "V3: Transmission and recovery tp data could not be written.");
    // Add data set 2 to HDF5 file
    mio::H5DataSet dset_loc_ids{
        H5Dcreate(h5group.id, "loc_ids", H5T_NATIVE_INT, dspace_d2.id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(dset_loc_ids.id, mio::StatusCode::UnknownError, "V3: Loc ids DataSet could not be created.");
    MEMILIO_H5_CHECK(H5Dwrite(dset_loc_ids.id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, loc_matrix.data()),
                     mio::StatusCode::UnknownError, "V3: Loc ids data could not be written.");
    return mio::success();
}

mio::IOResult<void> write_h5_v5(std::string filename,
                                mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                             LogPersonsPerLocationAndInfectionTime, LogAgentIds>& history)
{
    // hdf5 file contains two datasets (no groups)
    // 1st dataset: #agents x 2 matrix with entry (i,0) tp of transmission and entry (i,1) tp of recovery
    // 2nd dataset: #agents x #tps with entry (i,t) LocationType at t.

    // create hdf5 file
    mio::H5File file{H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(file.id, mio::StatusCode::FileNotFound, filename);

    // get agents ids
    auto log          = history.get_log();
    auto agent_ids    = std::get<3>(log)[0];
    auto logPerPerson = std::get<2>(log);

    //Create group for all data sets
    mio::H5Group h5group{H5Gcreate(file.id, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(h5group.id, mio::StatusCode::UnknownError, "H5Group \"data\" could not be created ");

    const int num_agents = static_cast<int>(agent_ids.size());
    // Dimension and matrix for 1st data set
    hsize_t dims_d1[] = {static_cast<hsize_t>(num_agents), static_cast<hsize_t>(2)};
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tsm_matrix(num_agents, 2);

    const int num_timepoints = static_cast<int>(logPerPerson.size());
    // Dimension and matrix for 2nd data set
    hsize_t dims_d2[] = {static_cast<hsize_t>(num_agents), static_cast<hsize_t>(num_timepoints)};
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> loc_matrix(num_agents, num_timepoints);
    // Fill matrices
    for (auto& id : agent_ids) {
        // initialize time point transmission and recovery with Null
        tsm_matrix(id.get(), 0)  = 100000;
        tsm_matrix(id.get(), 1)  = 100000;
        bool tp_transmission_set = false;
        bool tp_recovery_set     = false;
        for (int t = 0; t < num_timepoints; ++t) {
            // set time point of transmission
            if (!tp_transmission_set && time_since_transmission(std::get<3>(logPerPerson[t][id.get()])) > -1) {
                if (t == 0) {
                    tsm_matrix(id.get(), 0) = -time_since_transmission(std::get<3>(logPerPerson[t][id.get()]));
                }
                else {
                    tsm_matrix(id.get(), 0) = t + time_since_transmission(std::get<3>(logPerPerson[t][id.get()]));
                }
                tp_transmission_set = true;
            }
            // set time point of recovery
            if (tp_transmission_set && !tp_recovery_set &&
                time_since_transmission(std::get<3>(logPerPerson[t][id.get()])) == -1) {
                tsm_matrix(id.get(), 1) = t;
                tp_recovery_set         = true;
            }
            loc_matrix(id.get(), t) = static_cast<int>(std::get<1>(logPerPerson[t][id.get()]));
        }
    }
    // DataSpace for data set 1
    mio::H5DataSpace dspace_d1{H5Screate_simple(2, dims_d1, NULL)};
    MEMILIO_H5_CHECK(dspace_d1.id, mio::StatusCode::UnknownError, "V3: DataSpace 1 could not be created.");
    // DataSpace for data set 2
    mio::H5DataSpace dspace_d2{H5Screate_simple(2, dims_d2, NULL)};
    MEMILIO_H5_CHECK(dspace_d2.id, mio::StatusCode::UnknownError, "V3: DataSpace 2 could not be created.");
    // Add data set 1 to HDF5 file
    mio::H5DataSet dset_t_r_tp{H5Dcreate(h5group.id, "transm_recovery_tp", H5T_NATIVE_DOUBLE, dspace_d1.id, H5P_DEFAULT,
                                         H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(dset_t_r_tp.id, mio::StatusCode::UnknownError,
                     "V3: Transmission and recovery tp dataSet could not be created.");
    MEMILIO_H5_CHECK(H5Dwrite(dset_t_r_tp.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tsm_matrix.data()),
                     mio::StatusCode::UnknownError, "V3: Transmission and recovery tp data could not be written.");
    // Add data set 2 to HDF5 file
    mio::H5DataSet dset_loc_ids{
        H5Dcreate(h5group.id, "loc_ids", H5T_NATIVE_INT, dspace_d2.id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(dset_loc_ids.id, mio::StatusCode::UnknownError, "V3: Loc ids DataSet could not be created.");
    MEMILIO_H5_CHECK(H5Dwrite(dset_loc_ids.id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, loc_matrix.data()),
                     mio::StatusCode::UnknownError, "V3: Loc ids data could not be written.");
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

void set_wastewater_ids(std::string mapping_file, mio::abm::Model& model)
{
    // read in mapping from loc to TAN area id and save it in std::map
    std::map<std::tuple<mio::abm::LocationType, mio::abm::LocationId>, int> loc_to_tan_id;

    const boost::filesystem::path p = mapping_file;
    if (!boost::filesystem::exists(p)) {
        mio::log_error("Cannot read in mapping data. File does not exist.");
    }
    // File pointer
    std::fstream fin;

    // Open an existing file
    fin.open(mapping_file, std::ios::in);
    std::string line;

    while (std::getline(fin, line)) {
        std::stringstream ss(line);
        std::vector<std::string> v;
        while (getline(ss, line, ' ')) {
            // store token string in the vector
            v.push_back(line);
        }
        int tan_id   = std::stoi(v[1]);
        int loc_type = std::stoi(v[0].substr(0, 2));
        int loc_id   = std::stoi(v[0].substr(2));
        loc_to_tan_id.insert({std::make_tuple(mio::abm::LocationType(loc_type), mio::abm::LocationId(loc_id)), tan_id});
    }

    // add cemetary
    loc_to_tan_id.insert({std::make_tuple(mio::abm::LocationType::Cemetery, mio::abm::LocationId(0)), 0});

    model.set_wastewater_ids(loc_to_tan_id);
}

void initialize_model(mio::abm::Model& model, std::string person_file, std::string hosp_file, std::string outfile,
                      size_t max_work_size, size_t max_school_size)
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
    std::map<std::pair<mio::abm::LocationType, mio::abm::LocationId>, mio::abm::LocationId> hosp_to_icu;
    std::vector<double> hospital_weights;
    std::vector<double> icu_weights;
    // Mapping of assigned agents to school and work locations
    std::map<int, std::map<mio::abm::LocationId, size_t>> school_sizes;
    std::map<int, std::map<mio::abm::LocationId, size_t>> work_sizes;

    // Read in hospitals
    const boost::filesystem::path h = hosp_file;
    if (!boost::filesystem::exists(h)) {
        mio::log_error("Cannot read in data. Hospital file does not exist.");
    }
    // File pointer
    std::fstream fin_hosp;
    // Open an existing file
    fin_hosp.open(hosp_file, std::ios::in);
    std::vector<int32_t> row_hosp;
    std::vector<std::string> row_string_hosp;
    std::string line_hosp;
    // Read the Titles from the Data file
    std::getline(fin_hosp, line_hosp);
    line_hosp.erase(std::remove(line_hosp.begin(), line_hosp.end(), '\r'), line_hosp.end());
    std::vector<std::string> titles_hosp;
    boost::split(titles_hosp, line_hosp, boost::is_any_of(","));
    uint32_t count_of_titles_hosp              = 0;
    std::map<std::string, uint32_t> index_hosp = {};
    for (auto const& title : titles_hosp) {
        index_hosp.insert({title, count_of_titles_hosp});
        row_string_hosp.push_back(title);
        count_of_titles_hosp++;
    }
    while (std::getline(fin_hosp, line_hosp)) {
        row_hosp.clear();

        // read columns in this row
        split_line(line_hosp, &row_hosp);
        line_hosp.erase(std::remove(line_hosp.begin(), line_hosp.end(), '\r'), line_hosp.end());

        int beds          = row_hosp[index_hosp["beds"]];
        int icu_beds      = row_hosp[index_hosp["icu_beds"]];
        int hospital_zone = row_hosp[index_hosp["hospital_zone"]];
        auto hosp         = model.add_location(mio::abm::LocationType::Hospital);
        hospitals.push_back(hosp);
        hospital_weights.push_back(beds);
        std::string loc =
            "0" + std::to_string(static_cast<int>(mio::abm::LocationType::Hospital)) + std::to_string(hosp.get());
        auto zone_iter = loc_area_mapping.find(hospital_zone);
        if (zone_iter == loc_area_mapping.end()) {
            loc_area_mapping.insert({hospital_zone, {loc}});
        }
        else {
            loc_area_mapping[hospital_zone].push_back(loc);
        }
        // Add icu if there is one
        if (icu_beds > 0) {
            auto icu = model.add_location(mio::abm::LocationType::ICU);
            icus.push_back(icu);
            icu_weights.push_back(icu_beds);
            hosp_to_icu.insert({std::make_pair(mio::abm::LocationType::Hospital, hosp), icu});
            std::string loc_icu =
                "0" + std::to_string(static_cast<int>(mio::abm::LocationType::ICU)) + std::to_string(icu.get());
            zone_iter = loc_area_mapping.find(hospital_zone);
            if (zone_iter == loc_area_mapping.end()) {
                loc_area_mapping.insert({hospital_zone, {loc_icu}});
            }
            else {
                loc_area_mapping[hospital_zone].push_back(loc_icu);
            }
        }
    }

    // Read in persons
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
        model.get_location(home).increase_size();

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
        model.get_location(shop).increase_size();

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
        model.get_location(event).increase_size();

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
                    // Add school to size map
                    school_sizes[school_id].insert({school, 1});
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
                if (school_sizes[school_id][school] == max_school_size) {
                    // Check if a new school has to be open or if there still is a school that has capacity
                    bool found = false;
                    for (auto const& [key, val] : school_sizes[school_id]) {
                        if (val < max_school_size) {
                            found = true;
                            school_sizes[school_id][key] += 1;
                            school = key;
                            break;
                        }
                    }
                    if (!found) {
                        // Create new school
                        school = model.add_location(mio::abm::LocationType::School);
                        school_locations.insert({school_id, school});
                        // Add school to size map
                        school_sizes[school_id].insert({school, 1});
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
                }
                else {
                    school_sizes[school_id][school] += 1;
                }
            }
            // Assign school location to person
            person.set_assigned_location(mio::abm::LocationType::School, school);
            model.get_location(school).increase_size();
        }
        // Check if person is work-aged
        if (person.get_age() == mio::AgeGroup(2) || person.get_age() == mio::AgeGroup(3)) {
            int work_id   = row[index["work_id"]];
            int work_zone = row[index["work_zone"]];

            if (work_zone == -2) {
                mio::log_error("Person with id {} has work age but no work zone", row[index["puid"]]);
            }

            mio::abm::LocationId work;

            auto iter_work = work_locations.find(work_id);
            // Check whether work location already exists in model
            if (iter_work == work_locations.end()) {
                // Add work location to model and insert it in mapping
                work = model.add_location(mio::abm::LocationType::Work);
                // Locations with id -1 are individual locations each
                if (work_id != -1) {
                    work_locations.insert({work_id, work});
                    // Add work to size map
                    work_sizes[work_id].insert({work, 1});
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
                if (work_sizes[work_id][work] == max_work_size) {
                    // Check if a new work has to be opened or if there still is a school that has capacity
                    bool found = false;
                    for (auto const& [key, val] : work_sizes[work_id]) {
                        if (val < max_work_size) {
                            found = true;
                            work_sizes[work_id][key] += 1;
                            work = key;
                            break;
                        }
                    }
                    if (!found) {
                        // Create new work
                        work = model.add_location(mio::abm::LocationType::Work);
                        work_locations.insert({work_id, work});
                        // Add work to size map
                        work_sizes[work_id].insert({work, 1});
                        std::string loc = "0" + std::to_string(static_cast<int>(mio::abm::LocationType::Work)) +
                                          std::to_string(work.get());
                        auto zone_iter = loc_area_mapping.find(work_zone);
                        if (zone_iter == loc_area_mapping.end()) {
                            loc_area_mapping.insert({work_zone, {loc}});
                        }
                        else {
                            loc_area_mapping[work_zone].push_back(loc);
                        }
                    }
                }
                else {
                    work_sizes[work_id][work] += 1;
                }
            }
            // Assign work location to person
            person.set_assigned_location(mio::abm::LocationType::Work, work);
            model.get_location(work).increase_size();
        }
        // Assign Hospital and ICU
        size_t hosp = mio::DiscreteDistribution<size_t>::get_instance()(model.get_rng(), hospital_weights);
        person.set_assigned_location(mio::abm::LocationType::Hospital, hospitals[hosp]);
        model.get_location(hospitals[hosp]).increase_size();
        if (hosp_to_icu.count(std::make_pair(mio::abm::LocationType::Hospital, hospitals[hosp])) > 0) {
            person.set_assigned_location(
                mio::abm::LocationType::ICU,
                hosp_to_icu[std::make_pair(mio::abm::LocationType::Hospital, hospitals[hosp])]);
            model.get_location(hosp_to_icu[std::make_pair(mio::abm::LocationType::Hospital, hospitals[hosp])])
                .increase_size();
        }
        else {
            size_t icu = mio::DiscreteDistribution<size_t>::get_instance()(model.get_rng(), icu_weights);
            person.set_assigned_location(mio::abm::LocationType::ICU, icus[icu]);
            model.get_location(icus[icu]).increase_size();
        }
    }

    write_mapping_to_file(outfile, loc_area_mapping);
}

void write_size_per_location(std::string out_file, mio::abm::Model& model)
{
    std::map<std::string, size_t> size_per_loc;
    // Count number of assigned persons for each location
    for (auto& a : model.get_persons()) {
        for (auto& loc_id : a.get_assigned_locations()) {
            auto& loc = model.get_location(loc_id);
            if (loc_id != mio::abm::LocationId::invalid_id() && loc.get_type() != mio::abm::LocationType::Cemetery) {
                std::string loc_string =
                    "0" + std::to_string(static_cast<int>(loc.get_type())) + std::to_string(loc.get_id().get());
                auto string_iter = size_per_loc.find(loc_string);
                if (string_iter == size_per_loc.end()) {
                    size_per_loc.insert({loc_string, 1});
                }
                else {
                    size_per_loc[loc_string] += 1;
                }
            }
        }
    }
    //write map to file
    auto file = fopen(out_file.c_str(), "w");
    if (file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", out_file);
    }
    else {
        for (auto it = size_per_loc.begin(); it != size_per_loc.end(); it++) {
            fprintf(file, "%s", (it->first).c_str());
            fprintf(file, " %d", int(it->second));
            fprintf(file, "\n");
        }
        fclose(file);
    }
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
    pymio::bind_CustomIndexArray<double, mio::abm::VirusVariant>(m, "_InfectionRateArray");
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
                         mio::abm::TimePoint start_date, mio::abm::InfectionState start_state, bool detected,
                         bool shift_init, double shift_rate) {
            auto rng = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
            return mio::abm::Infection(rng, variant, person.get_age(), model.parameters, start_date, start_state,
                                       person.get_latest_protection(), detected, shift_init, shift_rate);
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
        .def("add_infection_rate_damping", &mio::abm::Model::add_infection_rate_damping, py::arg("t"),
             py::arg("factor"))
        .def("add_location_closure", &mio::abm::Model::add_location_closure, py::arg("t"), py::arg("loc_type"),
             py::arg("percentage"))
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
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double my,
           double sigma) {
            infection_params.get<mio::abm::IncubationPeriod>()[{variant, age}] = {my, sigma};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedNoSymptomsToSymptoms",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double my,
           double sigma) {
            infection_params.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{variant, age}] = {my, sigma};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedNoSymptomsToRecovered",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double my,
           double sigma) {
            infection_params.get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{variant, age}] = {my, sigma};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedSymptomsToSevere",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double my,
           double sigma) {
            infection_params.get<mio::abm::TimeInfectedSymptomsToSevere>()[{variant, age}] = {my, sigma};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedSymptomsToRecovered",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double my,
           double sigma) {
            infection_params.get<mio::abm::TimeInfectedSymptomsToRecovered>()[{variant, age}] = {my, sigma};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedSevereToRecovered",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double my,
           double sigma) {
            infection_params.get<mio::abm::TimeInfectedSevereToRecovered>()[{variant, age}] = {my, sigma};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedSevereToCritical",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double my,
           double sigma) {
            infection_params.get<mio::abm::TimeInfectedSevereToCritical>()[{variant, age}] = {my, sigma};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedCriticalToRecovered",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double my,
           double sigma) {
            infection_params.get<mio::abm::TimeInfectedCriticalToRecovered>()[{variant, age}] = {my, sigma};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedCriticalToDead",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double my,
           double sigma) {
            infection_params.get<mio::abm::TimeInfectedCriticalToDead>()[{variant, age}] = {my, sigma};
        },
        py::return_value_policy::reference_internal);

    m.def("set_AgeGroupGoToSchool", [](mio::abm::Parameters& infection_params, mio::AgeGroup age) {
        infection_params.get<mio::abm::AgeGroupGotoSchool>()[age] = true;
    });

    m.def("set_AgeGroupGoToWork", [](mio::abm::Parameters& infection_params, mio::AgeGroup age) {
        infection_params.get<mio::abm::AgeGroupGotoWork>()[age] = true;
    });

    m.def("set_AgeGroupGoToShop", [](mio::abm::Parameters& infection_params, mio::AgeGroup age) {
        infection_params.get<mio::abm::AgeGroupGotoShop>()[age] = true;
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
    m.def("set_wastewater_ids", &set_wastewater_ids, py::return_value_policy::reference_internal);
    m.def("write_size_per_location", &write_size_per_location, py::return_value_policy::reference_internal);

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

    m.def(
        "write_h5_v2",
        [](std::string filename, mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                              LogPersonsPerLocationAndInfectionTime, LogAgentIds>& history) {
            auto res = write_h5_v2(filename, history);
            if (res == mio::success()) {
                return 0;
            }
            else {
                return 1;
            }
        },
        py::return_value_policy::reference_internal);

    m.def(
        "write_h5_v3",
        [](std::string filename, mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                              LogPersonsPerLocationAndInfectionTime, LogAgentIds>& history) {
            auto res = write_h5_v3(filename, history);
            if (res == mio::success()) {
                return 0;
            }
            else {
                return 1;
            }
        },
        py::return_value_policy::reference_internal);

    m.def(
        "write_h5_v4",
        [](std::string filename, mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                              LogPersonsPerLocationAndInfectionTime, LogAgentIds>& history) {
            auto res = write_h5_v4(filename, history);
            if (res == mio::success()) {
                return 0;
            }
            else {
                return 1;
            }
        },
        py::return_value_policy::reference_internal);

    m.def(
        "write_h5_v5",
        [](std::string filename, mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                              LogPersonsPerLocationAndInfectionTime, LogAgentIds>& history) {
            auto res = write_h5_v5(filename, history);
            if (res == mio::success()) {
                return 0;
            }
            else {
                return 1;
            }
        },
        py::return_value_policy::reference_internal);

    m.def(
        "calculate_infections_per_quantity",
        [](std::string sim_result_folder, std::string save_folder, int num_sims, int start_sim,
           std::string person_file) {
            std::string save_file_loctype_infections    = save_folder + "num_agents_infections_loctype.txt";
            std::string save_file_ww_area_infections    = save_folder + "num_agents_infections_area.txt";
            std::string save_file_hh_size_ag_infections = save_folder + "num_agents_infections_hh_size_ag.txt";
            calculate_infections_per_quantity(sim_result_folder, save_file_loctype_infections, "locType", "v5",
                                              num_sims, start_sim);
            calculate_infections_per_quantity(sim_result_folder, save_file_ww_area_infections, "area", "v4", num_sims,
                                              start_sim);
            calculate_infections_per_hh_size(sim_result_folder, person_file, save_file_hh_size_ag_infections, num_sims,
                                             start_sim);
        },
        py::return_value_policy::reference_internal);

    m.def(
        "calculate_agents_per_quantity_age_groups",
        [](std::string sim_result_folder, std::string save_folder, int num_sims, int start_sim,
           std::string person_file) {
            std::string save_file_loctype = save_folder + "num_agents_loctype_ag.txt";
            std::string save_file_ww_area = save_folder + "num_agents_area_ag.txt";
            calculate_agents_per_quantity_age_groups(sim_result_folder, save_file_loctype, "locType", person_file,
                                                     num_sims, start_sim, "v5");
            calculate_agents_per_quantity_age_groups(sim_result_folder, save_file_ww_area, "area", person_file,
                                                     num_sims, start_sim, "v4");
        },
        py::return_value_policy::reference_internal);
#endif

    m.attr("__version__") = "dev";
}

PYMIO_IGNORE_VALUE_TYPE(decltype(std::declval<mio::abm::Model>().get_locations()))
PYMIO_IGNORE_VALUE_TYPE(decltype(std::declval<mio::abm::Model>().get_persons()))
