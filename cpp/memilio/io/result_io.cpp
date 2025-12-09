/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kuehn
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
#include "memilio/io/result_io.h"

#ifdef MEMILIO_HAS_HDF5

#include "memilio/io/hdf5_cpp.h"
#include "memilio/math/eigen_util.h"
#include "memilio/epidemiology/damping.h"

#include <vector>
#include <iostream>
#include <string>

namespace mio
{
IOResult<void> save_result(const std::vector<TimeSeries<ScalarType>>& results, const std::vector<int>& ids,
                           int num_groups, const std::string& filename)
{
    int region_idx = 0;
    H5File file{H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(file.id, StatusCode::FileNotFound, filename);
    for (auto& result : results) {
        auto h5group_name = "/" + std::to_string(ids[region_idx]);
        H5Group region_h5group{H5Gcreate(file.id, h5group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
        MEMILIO_H5_CHECK(region_h5group.id, StatusCode::UnknownError,
                         "Group could not be created (" + h5group_name + ")");

        const int num_timepoints      = static_cast<int>(result.get_num_time_points());
        const int num_infectionstates = (int)result.get_num_elements() / num_groups;

        hsize_t dims_t[] = {static_cast<hsize_t>(num_timepoints)};
        H5DataSpace dspace_t{H5Screate_simple(1, dims_t, NULL)};
        MEMILIO_H5_CHECK(dspace_t.id, StatusCode::UnknownError, "Time DataSpace could not be created.");
        H5DataSet dset_t{H5Dcreate(region_h5group.id, "Time", H5T_NATIVE_DOUBLE, dspace_t.id, H5P_DEFAULT, H5P_DEFAULT,
                                   H5P_DEFAULT)};
        MEMILIO_H5_CHECK(dset_t.id, StatusCode::UnknownError, "Time DataSet could not be created (Time).");
        auto values_t = std::vector<ScalarType>(result.get_times().begin(), result.get_times().end());
        MEMILIO_H5_CHECK(H5Dwrite(dset_t.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values_t.data()),
                         StatusCode::UnknownError, "Time data could not be written.");

        auto total = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(
                         num_timepoints, num_infectionstates)
                         .eval();

        for (int group_idx = 0; group_idx <= num_groups; ++group_idx) {
            auto group = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(
                             num_timepoints, num_infectionstates)
                             .eval();
            if (group_idx < num_groups) {
                for (Eigen::Index t_idx = 0; t_idx < result.get_num_time_points(); ++t_idx) {
                    auto v           = result[t_idx].transpose().eval();
                    auto group_slice = mio::slice(v, {group_idx * num_infectionstates, num_infectionstates});
                    mio::slice(group, {t_idx, 1}, {0, num_infectionstates}) = group_slice;
                    mio::slice(total, {t_idx, 1}, {0, num_infectionstates}) += group_slice;
                }
            }

            hsize_t dims_values[] = {static_cast<hsize_t>(num_timepoints), static_cast<hsize_t>(num_infectionstates)};
            H5DataSpace dspace_values{H5Screate_simple(2, dims_values, NULL)};
            MEMILIO_H5_CHECK(dspace_values.id, StatusCode::UnknownError, "Values DataSpace could not be created.");
            auto dset_name = group_idx == num_groups ? std::string("Total") : "Group" + std::to_string(group_idx + 1);
            H5DataSet dset_values{H5Dcreate(region_h5group.id, dset_name.c_str(), H5T_NATIVE_DOUBLE, dspace_values.id,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
            MEMILIO_H5_CHECK(dset_values.id, StatusCode::UnknownError, "Values DataSet could not be created.");

            MEMILIO_H5_CHECK(H5Dwrite(dset_values.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                      group_idx == num_groups ? total.data() : group.data()),
                             StatusCode::UnknownError, "Values data could not be written.");
        }
        region_idx++;
    }
    return success();
}

herr_t store_group_name(hid_t /*id*/, const char* name, const H5L_info_t* /*linfo*/, void* opdata)
{
    auto h5group_names = reinterpret_cast<std::vector<std::string>*>(opdata);
    h5group_names->push_back(name);
    return 0;
}

IOResult<std::vector<SimulationResult>> read_result(const std::string& filename)
{
    std::vector<SimulationResult> results;

    H5File file{H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)};
    MEMILIO_H5_CHECK(file.id, StatusCode::FileNotFound, filename);

    std::vector<std::string> h5group_names;
    MEMILIO_H5_CHECK(H5Literate(file.id, H5_INDEX_NAME, H5_ITER_INC, NULL, &store_group_name, &h5group_names),
                     StatusCode::UnknownError, "Group names could not be read.");

    std::sort(h5group_names.begin(), h5group_names.end(),
            [](const std::string &a, const std::string &b) {
                return std::stof(a) < std::stof(b);
            });

    for (auto& h5group_name : h5group_names) {
        H5Group region_h5group{H5Gopen(file.id, h5group_name.c_str(), H5P_DEFAULT)};
        MEMILIO_H5_CHECK(region_h5group.id, StatusCode::UnknownError,
                         "Group could not be opened (" + h5group_name + ")");

        std::vector<std::string> h5dset_names;
        MEMILIO_H5_CHECK(
            H5Literate(region_h5group.id, H5_INDEX_NAME, H5_ITER_INC, NULL, &store_group_name, &h5dset_names),
            StatusCode::UnknownError, "Dataset names could not be read.");

        std::string h5_key = std::any_of(h5dset_names.begin(), h5dset_names.end(),
                                         [](const std::string& str) {
                                             return str.find("Group") == 0;
                                         })
                                 ? "Group"
                                 : "End";

        auto num_groups = (Eigen::Index)std::count_if(h5dset_names.begin(), h5dset_names.end(), [&h5_key](auto&& str) {
            return str.find(h5_key) != std::string::npos;
        });

        H5DataSet dataset_t{H5Dopen(region_h5group.id, "Time", H5P_DEFAULT)};
        MEMILIO_H5_CHECK(dataset_t.id, StatusCode::UnknownError, "Time DataSet could not be read.");

        // dataset dimensions
        H5DataSpace dataspace_t{H5Dget_space(dataset_t.id)};
        MEMILIO_H5_CHECK(dataspace_t.id, StatusCode::UnknownError, "Time DataSpace could not be read.");
        const auto n_dims_t = 1;
        hsize_t dims_t[n_dims_t];
        H5Sget_simple_extent_dims(dataspace_t.id, dims_t, NULL);

        auto num_timepoints = Eigen::Index(dims_t[0]);
        auto time           = std::vector<ScalarType>(num_timepoints);
        MEMILIO_H5_CHECK(H5Dread(dataset_t.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time.data()),
                         StatusCode::UnknownError, "Time data could not be read.");

        auto dataset_name_total("/" + h5group_name + "/Total");
        H5DataSet dataset_total{H5Dopen(file.id, dataset_name_total.c_str(), H5P_DEFAULT)};
        MEMILIO_H5_CHECK(dataset_total.id, StatusCode::UnknownError, "Totals DataSet could not be read.");

        //read data space dimensions
        H5DataSpace dataspace_total{H5Dget_space(dataset_total.id)};
        MEMILIO_H5_CHECK(dataspace_total.id, StatusCode::UnknownError, "Totals DataSpace could not be read.");
        const auto n_dims_total = 2;
        hsize_t dims_total[n_dims_total];
        H5Sget_simple_extent_dims(dataspace_total.id, dims_total, NULL);
        if (num_timepoints != Eigen::Index(dims_total[0])) {
            return failure(StatusCode::InvalidFileFormat, "Number of time points does not match.");
        }
        auto num_infectionstates = Eigen::Index(dims_total[1]);

        auto total_values = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(
            num_timepoints, num_infectionstates);
        MEMILIO_H5_CHECK(
            H5Dread(dataset_total.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, total_values.data()),
            StatusCode::UnknownError, "Totals data could not be read");

        auto totals = TimeSeries<ScalarType>(num_infectionstates);
        totals.reserve(num_timepoints);
        for (auto t_idx = 0; t_idx < num_timepoints; ++t_idx) {
            totals.add_time_point(time[t_idx], slice(total_values, {t_idx, 1}, {0, num_infectionstates}).transpose());
        }

        auto groups = TimeSeries<ScalarType>(num_infectionstates * num_groups);
        groups.reserve(num_timepoints);
        for (Eigen::Index t_idx = 0; t_idx < num_timepoints; ++t_idx) {
            groups.add_time_point(time[t_idx]);
        }

        std::vector<int> h5_key_indices;
        // Extract group indices from h5dset_names
        for (const auto& name : h5dset_names) {
            if (name.find(h5_key) == 0) {
                h5_key_indices.push_back(std::stoi(name.substr(h5_key.size())));
            }
        }

        for (auto h5_key_indx = size_t(0); h5_key_indx < h5_key_indices.size(); h5_key_indx++) {
            auto group_name = "/" + h5group_name + "/" + h5_key + std::to_string(h5_key_indices[h5_key_indx]);
            H5DataSet dataset_values{H5Dopen(file.id, group_name.c_str(), H5P_DEFAULT)};
            MEMILIO_H5_CHECK(dataset_values.id, StatusCode::UnknownError, "Values DataSet could not be read.");

            //read data space dimensions
            H5DataSpace dataspace_values{H5Dget_space(dataset_values.id)};
            MEMILIO_H5_CHECK(dataspace_values.id, StatusCode::UnknownError, "Values DataSpace could not be read.");
            const auto n_dims_values = 2;
            hsize_t dims_values[n_dims_values];
            H5Sget_simple_extent_dims(dataspace_values.id, dims_values, NULL);
            if (num_timepoints != Eigen::Index(dims_values[0])) {
                return failure(StatusCode::InvalidFileFormat, "Number of time points does not match.");
            }
            if (num_infectionstates != Eigen::Index(dims_values[1])) {
                return failure(StatusCode::InvalidFileFormat, "Number of infection states does not match.");
            }

            auto group_values = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(
                num_timepoints, num_infectionstates);
            MEMILIO_H5_CHECK(
                H5Dread(dataset_values.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_values.data()),
                StatusCode::UnknownError, "Values data could not be read");

            for (Eigen::Index idx_t = 0; idx_t < num_timepoints; idx_t++) {
                for (Eigen::Index idx_c = 0; idx_c < num_infectionstates; idx_c++) {
                    groups[idx_t][num_infectionstates * h5_key_indx + idx_c] = group_values(idx_t, idx_c);
                }
            }
        }

        results.push_back(SimulationResult(groups, totals));
    }
    return success(results);
}

} // namespace mio

#endif //MEMILIO_HAS_HDF5
