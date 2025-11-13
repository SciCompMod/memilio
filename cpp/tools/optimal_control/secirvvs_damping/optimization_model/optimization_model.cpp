#include "optimization_model.h"

OptimizationModel::OptimizationModel(const boost::filesystem::path& data_directory, size_t simulation_days,
                                     size_t num_age_groups)
    : m_data_directory(data_directory)
    , m_simulation_days(simulation_days)
    , m_num_age_groups(num_age_groups)
    , m_cache(std::make_shared<Cache>())
{
    std::cout << "OptimizationModel initialized with data directory: " << m_data_directory << std::endl;
}

boost::filesystem::path OptimizationModel::data_directory() const
{
    return this->m_data_directory;
}

size_t OptimizationModel::simulation_days() const
{
    return this->m_simulation_days;
}
size_t OptimizationModel::num_age_groups() const
{
    return this->m_num_age_groups;
}
