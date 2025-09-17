#include "optimization_model.h"

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
int OptimizationModel::num_age_groups() const
{
    return m_num_age_groups;
}