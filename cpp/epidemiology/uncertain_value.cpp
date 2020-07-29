#include "uncertain_value.h"

namespace epi
{

void UncertainValue::set_distribution(const ParameterDistribution& dist)
{
    m_dist.reset(dist.clone());
}

observer_ptr<ParameterDistribution> UncertainValue::get_distribution() const
{
    return m_dist.get();
}

void UncertainValue::draw_sample()
{
    if (m_dist) {
        m_value = m_dist->get_rand_sample();
    }
}

UncertainValue::~UncertainValue()
{
}

} // namespace epi
