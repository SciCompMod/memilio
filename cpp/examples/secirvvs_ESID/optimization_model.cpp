#include "optimization_model.h"

OptimizationModel::OptimizationModel(boost::filesystem::path data_directory, double t0, double tmax, double num_age_groups)
    : m_data_directory(std::move(data_directory))
    , m_t0(t0)
    , m_tmax(tmax)
    , m_num_age_groups(num_age_groups)
{
    std::cout << "OptimizationModel initialized with data directory: " << m_data_directory << std::endl;
}

const boost::filesystem::path& OptimizationModel::data_directory() const
{
    return m_data_directory;
}

double OptimizationModel::t0() const
{
    return m_t0;
}
double OptimizationModel::tmax() const
{
    return m_tmax;
}
double OptimizationModel::num_age_groups() const
{
    return m_num_age_groups;
}