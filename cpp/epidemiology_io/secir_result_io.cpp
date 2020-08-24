#include <epidemiology_io/secir_result_io.h>
#include <epidemiology/secir.h>
#include <epidemiology/damping.h>
#include <epidemiology/eigen_util.h>
#include <vector>

#include <iostream>
#include <string>
#include <H5Cpp.h>
using namespace H5;

namespace epi
{

SecirSimulationResult::SecirSimulationResult()
{
}

SecirSimulationResult::SecirSimulationResult(std::vector<double> time, std::vector<Eigen::VectorXd> groups,
                                             std::vector<std::vector<double>> total)
    : m_time{time}
    , m_groups{groups}
    , m_total{total}
{
}

std::vector<double> SecirSimulationResult::get_time_vector() const&
{
    return m_time;
}

std::vector<Eigen::VectorXd> SecirSimulationResult::get_groups_vectors() const&
{
    return m_groups;
}

std::vector<std::vector<double>> SecirSimulationResult::get_totals_vector() const&
{
    return m_total;
}

void save_result(const std::vector<double>& times, const std::vector<Eigen::VectorXd>& secir,
                 const std::string& filename)
{
    const int n_dims = 2;

    const int n_data    = secir.size();
    const int n_compart = SecirCompartments::SecirCount;
    const int nb_groups = secir[0].size() / n_compart;

    H5File file(filename, H5F_ACC_TRUNC);

    hsize_t dims[] = {static_cast<hsize_t>(n_data)};
    auto dataspace = DataSpace(1, dims);

    auto DATASET_NAME = "Time";
    auto dataset      = file.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(times.data(), PredType::NATIVE_DOUBLE);

    auto total =
        std::vector<Eigen::Matrix<double, n_compart, 1>>(n_data, Eigen::Matrix<double, n_compart, 1>::Constant(0));

    for (int group = 0; group < nb_groups + 1; ++group) {
        auto dset = std::vector<Eigen::Matrix<double, n_compart, 1>>(n_data);
        if (group < nb_groups) {
            for (size_t irow = 0; irow < secir.size(); ++irow) {
                auto slice = epi::slice(secir[irow], {group * n_compart, n_compart});
                dset[irow] = slice;
                total[irow] += slice;
            }
        }

        hsize_t dims[] = {static_cast<hsize_t>(n_data), static_cast<hsize_t>(n_compart)};
        auto dataspace = DataSpace(n_dims, dims);

        auto DATASET_NAME = group == nb_groups ? std::string("Total") : "Group" + std::to_string(group + 1);
        auto dataset      = file.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace);

        if (group == nb_groups)
            dataset.write(total.data(), PredType::NATIVE_DOUBLE);
        else
            dataset.write(dset.data(), PredType::NATIVE_DOUBLE);
    }
}

SecirSimulationResult read_result(const std::string& filename, int nb_groups)
{
    const H5std_string FILE_NAME(filename);
    const int nb_compart = SecirCompartments::SecirCount;

    H5File file(FILE_NAME, H5F_ACC_RDONLY);
    H5std_string DATASET_NAME_TIME("Time");
    DataSet dataset_time = file.openDataSet(DATASET_NAME_TIME);

    DataSpace filespace_time = dataset_time.getSpace();

    const auto n_dims_time = 1;
    hsize_t dims_time[n_dims_time]; // dataset dimensions
    filespace_time.getSimpleExtentDims(dims_time);

    DataSpace mspace1(n_dims_time, dims_time);

    auto time = std::vector<double>(dims_time[0]);
    dataset_time.read(time.data(), PredType::NATIVE_DOUBLE, mspace1, filespace_time);

    std::vector<Eigen::VectorXd> groups;
    for (int i = 0; i < dims_time[0]; i++) {
        Eigen::VectorXd group_vector(nb_compart * nb_groups);
        groups.push_back(group_vector);
    }

    for (int i = 0; i < nb_groups; ++i) {
        auto DATASET_NAME_GROUP = "Group" + std::to_string(i + 1);
        auto dataset_group      = file.openDataSet(DATASET_NAME_GROUP);

        auto filespace_group = dataset_group.getSpace();

        const int n_dims_group = 2;
        hsize_t dims_group[n_dims_group]; // dataset dimensions
        filespace_group.getSimpleExtentDims(dims_group);

        DataSpace mspace2(n_dims_group, dims_group);

        auto group = std::vector<double[nb_compart]>(dims_group[0]);
        dataset_group.read(group.data(), PredType::NATIVE_DOUBLE, mspace2, filespace_group);

        for (int j = 0; j < dims_group[0]; j++) {
            for (int k = 0; k < nb_compart; k++) {
                groups[j][nb_compart * i + k] = group[j][k];
            }
        }
    }

    H5std_string DATASET_NAME_TOTAL("Total");
    DataSet dataset_total = file.openDataSet(DATASET_NAME_TOTAL);

    DataSpace filespace = dataset_total.getSpace();

    const int n_dims_total = 2;
    hsize_t dims_total[n_dims_total]; // dataset dimensions
    filespace.getSimpleExtentDims(dims_total);

    DataSpace mspace3(n_dims_total, dims_total);

    auto total_buf = std::vector<double[nb_compart]>(dims_total[0]);
    dataset_total.read(total_buf.data(), PredType::NATIVE_DOUBLE, mspace3, filespace);

    std::vector<std::vector<double>> total(total_buf.size());
    std::transform(begin(total_buf), end(total_buf), begin(total), [nb_compart](auto v) {
        return std::vector<double>(v, v + nb_compart);
    });

    return SecirSimulationResult(time, groups, total);
}

} // namespace epi
