#include "optimization_model.h"

OptimizationModel::OptimizationModel(const boost::filesystem::path& data_directory,
                                     const boost::filesystem::path& optimal_control_settings_directory,

                                     size_t simulation_days, size_t num_age_groups, double npis_duration,
                                     double npis_interval, double npis_base_value)
    : m_data_directory(data_directory)
    , m_optimal_control_settings_directory(optimal_control_settings_directory)
    , m_simulation_days(simulation_days)
    , m_num_age_groups(num_age_groups)
    , m_npis_duration(npis_duration)
    , m_npis_interval(npis_interval)
    , m_npis_base_value(npis_base_value)
    , m_cache(std::make_shared<Cache>())
{
    std::cout << "OptimizationModel initialized with data directory: " << m_data_directory << std::endl;
}

boost::filesystem::path OptimizationModel::data_directory() const
{
    return this->m_data_directory;
}

boost::filesystem::path OptimizationModel::optimal_control_settings_directory() const
{
    return this->m_optimal_control_settings_directory;
}

size_t OptimizationModel::simulation_days() const
{
    return this->m_simulation_days;
}
size_t OptimizationModel::num_age_groups() const
{
    return this->m_num_age_groups;
}

double OptimizationModel::npis_duration() const
{
    return this->m_npis_duration;
}

double OptimizationModel::npis_interval() const
{
    return this->m_npis_interval;
}

double OptimizationModel::npis_base_value() const
{
    return this->m_npis_base_value;
}
