/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Wadim Koslow, Daniel Abele
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
#include "secir_vaccine/secir_result_io.h"

#ifdef MEMILIO_HAS_HDF5

#include "memilio/io/hdf5_cpp.h"
#include "memilio/math/eigen_util.h"
#include "memilio/epidemiology/damping.h"

#include <vector>
#include <iostream>
#include <string>

namespace mio
{
namespace vaccinated
{

    IOResult<void> save_result(const std::vector<TimeSeries<double>>& results, const std::vector<int>& ids,
                               const std::string& filename)
    {
        int county = 0;
        H5File file{H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)};
        MEMILIO_H5_CHECK(file.id, StatusCode::FileNotFound, filename);
        for (auto& result : results) {
            auto group_name = "/" + std::to_string(ids[county]);
            H5Group county_group{H5Gcreate(file.id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
            MEMILIO_H5_CHECK(county_group.id, StatusCode::UnknownError,
                             "Group could not be created (" + group_name + ")");

            const int n_data    = static_cast<int>(result.get_num_time_points());
            const int n_compart = (int)InfectionState::Count;
            const int nb_groups = static_cast<int>(result[0].size()) / n_compart;

            hsize_t dims_t[] = {static_cast<hsize_t>(n_data)};
            H5DataSpace dspace_t{H5Screate_simple(1, dims_t, NULL)};
            MEMILIO_H5_CHECK(dspace_t.id, StatusCode::UnknownError, "Time DataSpace could not be created.");
            H5DataSet dset_t{H5Dcreate(county_group.id, "Time", H5T_NATIVE_DOUBLE, dspace_t.id, H5P_DEFAULT,
                                       H5P_DEFAULT, H5P_DEFAULT)};
            MEMILIO_H5_CHECK(dset_t.id, StatusCode::UnknownError, "Time DataSet could not be created (Time).");
            auto values_t = std::vector<double>(result.get_times().begin(), result.get_times().end());
            MEMILIO_H5_CHECK(H5Dwrite(dset_t.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values_t.data()),
                             StatusCode::UnknownError, "Time data could not be written.");

            auto total = std::vector<Eigen::Matrix<double, n_compart, 1>>(
                n_data, Eigen::Matrix<double, n_compart, 1>::Constant(0));

            for (int group = 0; group < nb_groups + 1; ++group) {
                auto dset = std::vector<Eigen::Matrix<double, n_compart, 1>>(n_data);
                if (group < nb_groups) {
                    for (size_t irow = 0; irow < static_cast<size_t>(result.get_num_time_points()); ++irow) {
                        auto v     = result[irow].eval();
                        auto slice = mio::slice(v, {group * n_compart, n_compart});
                        dset[irow] = slice;
                        total[irow] += slice;
                    }
                }

                hsize_t dims_values[] = {static_cast<hsize_t>(n_data), static_cast<hsize_t>(n_compart)};
                H5DataSpace dspace_values{H5Screate_simple(2, dims_values, NULL)};
                MEMILIO_H5_CHECK(dspace_values.id, StatusCode::UnknownError, "Values DataSpace could not be created.");
                auto dset_name = group == nb_groups ? std::string("Total") : "Group" + std::to_string(group + 1);
                H5DataSet dset_values{H5Dcreate(county_group.id, dset_name.c_str(), H5T_NATIVE_DOUBLE, dspace_values.id,
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
                MEMILIO_H5_CHECK(dset_values.id, StatusCode::UnknownError, "Values DataSet could not be created.");

                MEMILIO_H5_CHECK(H5Dwrite(dset_values.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                          group == nb_groups ? total.data() : dset.data()),
                                 StatusCode::UnknownError, "Values data could not be written.");
            }
            county++;
        }
        return success();
    }

    herr_t store_group_name(hid_t loc_id, const char* name, const H5L_info_t* linfo, void* opdata)
    {
        unused(linfo);
        hid_t group;
        auto group_names = reinterpret_cast<std::vector<std::string>*>(opdata);
        group            = H5Gopen2(loc_id, name, H5P_DEFAULT);
        group_names->push_back(name);
        H5Gclose(group);
        return 0;
    }

    IOResult<std::vector<SecirSimulationResult>> read_result(const std::string& filename, int nb_groups)
    {
        const int nb_compart = (int)InfectionState::Count;
        std::vector<SecirSimulationResult> results;

        H5File file{H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)};
        MEMILIO_H5_CHECK(file.id, StatusCode::FileNotFound, filename);

        std::vector<std::string> group_names;
        MEMILIO_H5_CHECK(H5Literate(file.id, H5_INDEX_NAME, H5_ITER_INC, NULL, &store_group_name, &group_names),
                         StatusCode::UnknownError, "Group names could not be read.");

        for (auto& name : group_names) {
            auto dataset_name_t = "/" + name + "/Time";
            H5DataSet dataset_t{H5Dopen(file.id, dataset_name_t.c_str(), H5P_DEFAULT)};
            MEMILIO_H5_CHECK(dataset_t.id, StatusCode::UnknownError, "Time DataSet could not be read.");

            // dataset dimensions
            H5DataSpace dataspace_t{H5Dget_space(dataset_t.id)};
            MEMILIO_H5_CHECK(dataspace_t.id, StatusCode::UnknownError, "Time DataSpace could not be read.");
            const auto n_dims_t = 1;
            hsize_t dims_t[n_dims_t];
            H5Sget_simple_extent_dims(dataspace_t.id, dims_t, NULL);

            auto time = std::vector<double>(dims_t[0]);
            MEMILIO_H5_CHECK(H5Dread(dataset_t.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time.data()),
                             StatusCode::UnknownError, "Time data could not be read.");

            TimeSeries<double> groups(nb_compart * nb_groups), totals(nb_compart);
            groups.reserve(dims_t[0]);
            totals.reserve(dims_t[0]);
            for (size_t i = 0; i < dims_t[0]; i++) {
                groups.add_time_point(time[i]);
                totals.add_time_point(time[i]);
            }

            for (int i = 0; i < nb_groups; ++i) {
                auto group_name = "/" + name + "/Group" + std::to_string(i + 1);
                H5DataSet dataset_values{H5Dopen(file.id, group_name.c_str(), H5P_DEFAULT)};
                MEMILIO_H5_CHECK(dataset_values.id, StatusCode::UnknownError, "Values DataSet could not be read.");

                //read data space dimensions
                H5DataSpace dataspace_values{H5Dget_space(dataset_values.id)};
                MEMILIO_H5_CHECK(dataspace_values.id, StatusCode::UnknownError, "Values DataSpace could not be read.");
                const auto n_dims_values = 2;
                hsize_t dims_values[n_dims_values];
                H5Sget_simple_extent_dims(dataspace_values.id, dims_values, NULL);

                auto group_values = std::vector<double[nb_compart]>(dims_values[0]);
                MEMILIO_H5_CHECK(
                    H5Dread(dataset_values.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_values.data()),
                    StatusCode::UnknownError, "Values data could not be read");

                for (size_t j = 0; j < dims_values[0]; j++) {
                    for (size_t k = 0; k < nb_compart; k++) {
                        groups[j][nb_compart * i + k] = group_values[j][k];
                    }
                }
            }

            auto dataset_name_total("/" + name + "/Total");
            H5DataSet dataset_total{H5Dopen(file.id, dataset_name_total.c_str(), H5P_DEFAULT)};
            MEMILIO_H5_CHECK(dataset_total.id, StatusCode::UnknownError, "Totals DataSet could not be read.");

            //read data space dimensions
            H5DataSpace dataspace_total{H5Dget_space(dataset_total.id)};
            MEMILIO_H5_CHECK(dataspace_total.id, StatusCode::UnknownError, "Totals DataSpace could not be read.");
            const auto n_dims_total = 2;
            hsize_t dims_total[n_dims_total];
            H5Sget_simple_extent_dims(dataspace_total.id, dims_total, NULL);

            auto total_values = std::vector<Eigen::Matrix<double, nb_compart, 1>>(dims_total[0]);
            MEMILIO_H5_CHECK(
                H5Dread(dataset_total.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, total_values.data()),
                StatusCode::UnknownError, "Totals data could not be read");

            std::copy(begin(total_values), end(total_values), totals.begin());
            results.push_back(SecirSimulationResult(groups, totals));
        }
        return success(results);
    }

} // namespace vaccinated
} // namespace mio

#endif //MEMILIO_HAS_HDF5
