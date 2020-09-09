#include "epidemiology/utils/uncertain_value.h"

namespace epi
{

void UncertainValue::set_distribution(const ParameterDistribution& dist)
{
    m_dist.reset(dist.clone());
}

observer_ptr<ParameterDistribution> UncertainValue::get_distribution()
{
    return m_dist.get();
}

observer_ptr<const ParameterDistribution> UncertainValue::get_distribution() const
{
    return m_dist.get();
}

double UncertainValue::draw_sample()
{
    if (m_dist) {
        m_value = m_dist->get_sample();
    }

    return m_value;
}

} // namespace epi
