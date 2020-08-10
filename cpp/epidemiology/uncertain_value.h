#ifndef UNCERTAINVALUE_H
#define UNCERTAINVALUE_H

#include <memory>
#include <epidemiology/memory.h>
#include <epidemiology/parameter_studies/parameter_distributions.h>

namespace epi
{

/**
 * @brief The UncertainValue class consists of a 
 *        scalar double value and a Distribution object
 * 
 * The UncertainValue class represents a model parameter that
 * can take a scalar value but that is subjected to a uncertainty.
 * The uncertainty is represented by a distribution object of kind
 * ParameterDistribution and the current scalar value can be 
 * replaced by drawing a new sample from the the distribution
 */
class UncertainValue
{
public:
    UncertainValue(double v = 0.)
        : m_value(v)
    {
    }

    UncertainValue(UncertainValue&& other) = default;

    /**
    * @brief Create an UncertainValue by cloning scalar value 
    *        and distribution of another UncertainValue
    */
    UncertainValue(const UncertainValue& other)
        : m_value(other.m_value)
    {
        if (other.m_dist) {
            m_dist.reset(other.m_dist->clone());
        }
    }

    /**
    * @brief Set an UncertainValue from another UncertainValue
    *        containing a double and a distribution
    */
    UncertainValue& operator=(const UncertainValue& other)
    {
        UncertainValue tmp(other);
        std::swap(*this, tmp);
        return *this;
    }

    /**
     * @brief Conversion to double by returning the double contained in UncertainValue
     */
    operator double() const
    {
        return m_value;
    }

    /**
     * @brief Conversion to double reference by returning the double contained in UncertainValue
     */
    operator double&()
    {
        return m_value;
    }

    /**
     * @brief Set an UncertainValue from a double, distribution remains unchanged.
     */
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
     *        and returns the new value
     *
     * If no distribution is set, the value is not changed.
     */
    double draw_sample();

    ~UncertainValue();

private:
    double m_value;
    std::unique_ptr<ParameterDistribution> m_dist;
};

} // namespace epi

#endif // UNCERTAINVALUE_H
