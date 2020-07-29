#ifndef UNCERTAINVALUE_H
#define UNCERTAINVALUE_H

#include <memory>
#include <epidemiology/memory.h>
#include <epidemiology/parameter_studies/parameter_distributions.h>

namespace epi
{

class UncertainValue
{
public:
    UncertainValue(double v = 0.)
        : m_value(v)
    {
    }

    UncertainValue(UncertainValue&& other) = default;

    UncertainValue(const UncertainValue& other)
        : m_value(other.m_value)
    {
        if (other.m_dist) {
            m_dist.reset(other.m_dist->clone());
        }
    }

    UncertainValue& operator=(const UncertainValue& other)
    {
        UncertainValue tmp(other);
        std::swap(*this, tmp);
        return *this;
    }

    operator double() const
    {
        return m_value;
    }

    operator double&()
    {
        return m_value;
    }

    UncertainValue& operator=(double v)
    {
        m_value = v;
        return *this;
    }

    /**
     * @brief Sets the distribution of the value.
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void set_distribution(const ParameterDistribution& dist);

    /**
     * @brief Returns the parameter distribution.
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_distribution() const;

    /**
     * @brief Sets the value by sampling from the distribution
     *
     * If no distribution is set, the value is not changed.
     */
    void draw_sample();

    ~UncertainValue();

private:
    double m_value;
    std::unique_ptr<ParameterDistribution> m_dist;
};

} // namespace epi

#endif // UNCERTAINVALUE_H
