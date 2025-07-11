#include "optimization_model.h"

OptimizationModel::OptimizationModel(const std::filesystem::path& data_directory, double t0, double tmax)
    : m_data_directory(data_directory)
    , m_t0(t0)
    , m_tmax(tmax)
{
    std::cout << "OptimizationModel initialized with data directory: " << m_data_directory << std::endl;
}

std::filesystem::path OptimizationModel::data_directory() const
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