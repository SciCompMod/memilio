#include <epidemiology_io/secir_result_io.h>
#include "epidemiology_io/hdf5_cpp.h"
#include <epidemiology/utils/eigen_util.h>
#include <epidemiology/secir/secir.h>
#include <epidemiology/secir/damping.h>

#include <vector>
#include <iostream>
#include <string>

using namespace H5;

namespace epi
{

IOResult<void> save_result(const std::vector<TimeSeries<double>>& results, const std::vector<int>& ids,
                           const std::string& filename)
{
    try {
        int county = 0;
        H5File file(filename, H5F_ACC_TRUNC);
        for (auto& result : results) {

            H5::Group county_group = file.createGroup("/" + std::to_string(ids[county]));
            const int n_dims       = 2;

            const int n_data    = static_cast<int>(result.get_num_time_points());
            const int n_compart = (int)InfectionState::Count;
            const int nb_groups = static_cast<int>(result[0].size()) / n_compart;

            hsize_t dims_t[] = {static_cast<hsize_t>(n_data)};
            auto dspace_t    = DataSpace(1, dims_t);
            auto dset_t      = county_group.createDataSet("Time", PredType::NATIVE_DOUBLE, dspace_t);
            auto times       = std::vector<double>(result.get_times().begin(), result.get_times().end());
            dset_t.write(times.data(), PredType::NATIVE_DOUBLE);

            auto total = std::vector<Eigen::Matrix<double, n_compart, 1>>(
                n_data, Eigen::Matrix<double, n_compart, 1>::Constant(0));

            for (int group = 0; group < nb_groups + 1; ++group) {
                auto dset = std::vector<Eigen::Matrix<double, n_compart, 1>>(n_data);
                if (group < nb_groups) {
                    for (size_t irow = 0; irow < static_cast<size_t>(result.get_num_time_points()); ++irow) {
                        auto v     = result[irow].eval();
                        auto slice = epi::slice(v, {group * n_compart, n_compart});
                        dset[irow] = slice;
                        total[irow] += slice;
                    }
                }

                hsize_t dims_values[] = {static_cast<hsize_t>(n_data), static_cast<hsize_t>(n_compart)};
                auto dspace_values    = DataSpace(n_dims, dims_values);
                auto dset_name        = group == nb_groups ? std::string("Total") : "Group" + std::to_string(group + 1);
                auto dset_values      = county_group.createDataSet(dset_name, PredType::NATIVE_DOUBLE, dspace_values);

                if (group == nb_groups)
                    dset_values.write(total.data(), PredType::NATIVE_DOUBLE);
                else
                    dset_values.write(dset.data(), PredType::NATIVE_DOUBLE);
            }
            county++;
        }
        return success();
    }
    catch (H5::FileIException& ex) {
        return failure(StatusCode::FileNotFound, filename + ", " + ex.getDetailMsg());
    }
    catch (H5::Exception& ex) {
        return failure(StatusCode::UnknownError, filename + ", " + ex.getDetailMsg());
    }
}

herr_t file_info(hid_t loc_id, const char* name, const H5L_info_t* linfo, void* opdata)
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
    try {
        const H5std_string FILE_NAME(filename);
        const int nb_compart = (int)InfectionState::Count;
        std::vector<SecirSimulationResult> results;

        H5File file(FILE_NAME, H5F_ACC_RDONLY);
        hid_t file_it = file.getLocId();
        std::vector<std::string> group_names;
        herr_t idx = H5Literate(file_it, H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, &group_names);
        unused(idx);

        for (auto& name : group_names) {
            H5std_string DATASET_NAME_TIME("/" + name + "/Time");
            DataSet dataset_time = file.openDataSet(DATASET_NAME_TIME);

            DataSpace filespace_time = dataset_time.getSpace();

            const auto n_dims_time = 1;
            hsize_t dims_time[n_dims_time]; // dataset dimensions
            filespace_time.getSimpleExtentDims(dims_time);

            DataSpace mspace1(n_dims_time, dims_time);

            auto time = std::vector<double>(dims_time[0]);
            dataset_time.read(time.data(), PredType::NATIVE_DOUBLE, mspace1, filespace_time);

            TimeSeries<double> groups(nb_compart * nb_groups), totals(nb_compart);
            groups.reserve(dims_time[0]);
            totals.reserve(dims_time[0]);
            for (size_t i = 0; i < dims_time[0]; i++) {
                groups.add_time_point(time[i]);
                totals.add_time_point(time[i]);
            }

            for (int i = 0; i < nb_groups; ++i) {
                auto DATASET_NAME_GROUP = "/" + name + "/Group" + std::to_string(i + 1);
                auto dataset_group      = file.openDataSet(DATASET_NAME_GROUP);

                auto filespace_group = dataset_group.getSpace();

                const int n_dims_group = 2;
                hsize_t dims_group[n_dims_group]; // dataset dimensions
                filespace_group.getSimpleExtentDims(dims_group);

                DataSpace mspace2(n_dims_group, dims_group);

                auto group = std::vector<double[nb_compart]>(dims_group[0]);
                dataset_group.read(group.data(), PredType::NATIVE_DOUBLE, mspace2, filespace_group);

                for (size_t j = 0; j < dims_group[0]; j++) {
                    for (size_t k = 0; k < nb_compart; k++) {
                        groups[j][nb_compart * i + k] = group[j][k];
                    }
                }
            }

            H5std_string DATASET_NAME_TOTAL("/" + name + "/Total");
            DataSet dataset_total = file.openDataSet(DATASET_NAME_TOTAL);

            DataSpace filespace = dataset_total.getSpace();

            const int n_dims_total = 2;
            hsize_t dims_total[n_dims_total]; // dataset dimensions
            filespace.getSimpleExtentDims(dims_total);

            DataSpace mspace3(n_dims_total, dims_total);

            auto total_buf = std::vector<Eigen::Matrix<double, nb_compart, 1>>(dims_total[0]);
            dataset_total.read(total_buf.data(), PredType::NATIVE_DOUBLE, mspace3, filespace);

            std::copy(begin(total_buf), end(total_buf), totals.begin());
            results.push_back(SecirSimulationResult(groups, totals));
        }
        return results;
    }
    catch (H5::FileIException& ex) {
        return failure(StatusCode::FileNotFound, filename + ", " + ex.getDetailMsg());
    }
    catch (H5::Exception& ex) {
        return failure(StatusCode::InvalidFileFormat, filename + ", " + ex.getDetailMsg());
    }
}

} // namespace epi
