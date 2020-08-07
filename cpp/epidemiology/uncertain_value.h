#ifndef UNCERTAINVALUE_H
#define UNCERTAINVALUE_H

#include <memory>
#include <epidemiology/memory.h>
#include <epidemiology/secir.h>
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

/**
 * @brief The UncertainContactMatrix class consists of a
 *        ContactFrequencyMatrix and certain distributions that describe
 *        the relative and uncertain changes of diagonal and 
 *        offdiagonal values over time 
 * 
 * The UncertainContactMatrix class represents a matrix-style model parameter 
 * that can take a ContactFrequencyMatrix value but that is subjected to a uncertainty,
 * based on contact pattern changes realized by an uncertain number and shape of epi::Damping(s).
 * The uncertainty is represented by a several distributions.
 * The number of and days where contact pattern changes is characterized by two distributions.
 * The changes in diagonal entries are supposed to be at the same magnitude, its base 
 * value defined by one distribution. Their relative deviation from the base value is 
 * given by another distribution. Another distribution yields the relative deviation
 * of the offdiagonal entries with respect to their corresponding diagonals.
 */
class UncertainContactMatrix
{
public:
    UncertainContactMatrix(SecirParams::ContactFrequencyMatrix cont_freq = ContactFrequencyMatrix{1})
        : m_cont_freq(cont_freq)
    {
    }

    UncertainContactMatrix(SecirParams::ContactFrequencyMatrix&& cont_freq)
        : m_cont_freq(std::move(cont_freq))
    {
    }

    UncertainContactMatrix(UncertainContactMatrix&& other) = default;

    /**
    * @brief Create an UncertainContactMatrix by cloning ContactFrequencyMatrix 
    *        and distributions of another UncertainContactMatrix
    */
    UncertainContactMatrix(const UncertainContactMatrix& other)
        : m_cont_freq(other.m_cont_freq)
    {
        if (other.m_damp_diag_base) {
            m_damp_diag_base.reset(other.m_damp_diag_base->clone());
        }
        if (other.m_damp_diag_rel) {
            m_damp_diag_rel.reset(other.m_damp_diag_rel->clone());
        }
        if (other.m_damp_offdiag_rel) {
            m_damp_offdiag_rel.reset(other.m_damp_offdiag_rel->clone());
        }
        if (other.m_damp_nb) {
            m_damp_nb.reset(other.m_damp_nb->clone());
        }
        if (other.m_damp_days) {
            m_damp_days.reset(other.m_damp_days->clone());
        }
    }

    /**
    * @brief Set an UncertainContactMatrix from another UncertainContactMatrix, 
     *        containing a ContactFrequencyMatrix and related distributions
     *        for changes in contact patterns over time.
    */
    UncertainContactMatrix& operator=(const UncertainContactMatrix& other)
    {
        UncertainContactMatrix tmp(other);
        std::swap(*this, tmp);
        return *this;
    }

    /**
     * @brief Conversion to const ContactFrequencyMatrix reference by returning the 
     *        ContactFrequencyMatrix contained in UncertainContactMatrix
     */
    operator SecirParams::ContactFrequencyMatrix const &() const
    {
        return m_cont_freq;
    }

    /**
     * @brief Set an UncertainContactMatrix from a ContactFrequencyMatrix, 
     *        all distributions remain unchanged.
     */
    UncertainContactMatrix& operator=(SecirParams::ContactFrequencyMatrix cont_freq)
    {
        m_cont_freq = cont_freq;
        return *this;
    }

    /**
     * @brief Returns the ContactFrequencyMatrix reference 
     *        of the UncertainContactMatrix object
     */
    SecirParams::ContactFrequencyMatrix& get_cont_freq_mat()
    {
        return m_cont_freq;
    }    

    /**
     * @brief Returns the const ContactFrequencyMatrix reference 
     *        of the UncertainContactMatrix object
     */
    SecirParams::ContactFrequencyMatrix const& get_cont_freq_mat() const
    {
        return m_cont_freq;
    }

    /**
     * @brief Sets the random distribution for the number of 
     *        contact pattern changes over time
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void UncertainContactMatrix::set_dist_damp_nb(const ParameterDistribution& dist)
    {
        m_damp_nb.reset(dist.clone());
    }

    /**
     * @brief Sets the random distribution of the actual 
     *        points in time where the contact pattern changes
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void UncertainContactMatrix::set_dist_damp_days(const ParameterDistribution& dist)
    {
        m_damp_days.reset(dist.clone());
    }

    /**
     * @brief Sets the random distribution for a multiplicative base factor
     *        (changing the diagonal entries of the ContactFrequencyMatrix)
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void UncertainContactMatrix::set_dist_damp_diag_base(const ParameterDistribution& dist)
    {
        m_damp_diag_base.reset(dist.clone());
    }

    /**
     * @brief Sets the random distribution for the multiplicative factors
     *        changing each diagonal entry of the ContactFrequencyMatrix
     *        relative to the base value sampled by the distribution from
     *        get_dist_damp_diag_base()
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void UncertainContactMatrix::set_dist_damp_diag_rel(const ParameterDistribution& dist)
    {
        m_damp_diag_rel.reset(dist.clone());
    }

    /**
     * @brief Sets the random distribution for the multiplicative factors
     *        changing each offdiagonal entry of the ContactFrequencyMatrix
     *        relative to the corresponding diagonal entries by the distribution 
     *        from get_dist_damp_diag_rel()
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void UncertainContactMatrix::set_dist_damp_offdiag_rel(const ParameterDistribution& dist)
    {
        m_damp_offdiag_rel.reset(dist.clone());
    }

    /**
     * @brief Returns the random distribution for the number of 
     *        contact pattern changes over time
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_dist_damp_nb() const
    {
        return m_damp_nb.get();
    }

    /**
     * @brief Returns the random distribution of the actual 
     *        points in time where the contact pattern changes
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_dist_damp_days() const
    {
        return m_damp_days.get();
    }

    /**
     * @brief Returns the random distribution for a multiplicative base factor
     *        (changing the diagonal entries of the ContactFrequencyMatrix)
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_dist_damp_diag_base() const
    {
        return m_damp_diag_base.get();
    }

    /**
     * @brief Returns the random distribution for the multiplicative factors
     *        changing each diagonal entry of the ContactFrequencyMatrix
     *        relative to the base value sampled by the distribution from
     *        get_dist_damp_diag_base()
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_dist_damp_diag_rel() const
    {
        return m_damp_diag_rel.get();
    }

    /**
     * @brief Returns the random distribution for the multiplicative factors
     *        changing each offdiagonal entry of the ContactFrequencyMatrix
     *        relative to the corresponding diagonal entries by the distribution 
     *        from get_dist_damp_diag_rel()
     *        
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_dist_damp_offdiag_rel() const
    {
        return m_damp_offdiag_rel.get();
    }

    /**
     * @brief Sets the value by sampling from the distributions
     *
     * If no distribution is set, the value is not changed.
     * @param accum accumulating current and newly sampled dampings if true;
     *              default: false; removing all previously set dampings
     */
    SecirParams::ContactFrequencyMatrix draw_sample(bool accum = false)
    {
        if (accum) {
            m_cont_freq.clear_dampings();
        }

        if (m_damp_nb && m_damp_days && m_damp_diag_base && m_damp_diag_rel && m_damp_offdiag_rel) {

            int nb_dampings = (int)(m_damp_nb->get_sample() + 0.5);
            for (int i = 0; i < nb_dampings; i++) {

                double day            = m_damp_days->get_sample();
                double damp_diag_base = m_damp_diag_base->get_sample();

                // diagonal entries
                std::vector<double> damp_diag_val(m_cont_freq.get_size(), 0);
                for (int j = 0; j < m_cont_freq.get_size(); j++) {
                    damp_diag_val[j] = damp_diag_base * m_damp_diag_rel->get_sample();
                    m_cont_freq.add_damping(Damping(day, damp_diag_val[j]), j, j);
                }

                // offdiagonal entries
                for (int j = 0; j < m_cont_freq.get_size(); j++) {

                    for (int k = j + 1; k < m_cont_freq.get_size(); k++) {
                        double damp_offdiag_val = 0.5 * damp_diag_val[j] * m_damp_offdiag_rel->get_sample() +
                                                  0.5 * damp_diag_val[k] * m_damp_offdiag_rel->get_sample();
                        m_cont_freq.add_damping(Damping(day, damp_offdiag_val), j, k);
                    }
                }
            }
        }
        else {
            epi::log_warning("UncertainContactMatrix distributions not set, no sampling conducted.")
        }

        return m_cont_freq;
    }

    ~UncertainContactMatrix();

private:
    SecirParams::ContactFrequencyMatrix m_cont_freq;
    std::unique_ptr<ParameterDistribution>
        m_damp_nb; // random number of dampings (one damping is understood as nb_groups^2 many dampings at the same day)
    std::unique_ptr<ParameterDistribution> m_damp_days; // random number of day where to implement damping
    std::unique_ptr<ParameterDistribution>
        m_damp_diag_base; // random number of base value for the diagonal of the damping matrix
    std::unique_ptr<ParameterDistribution> m_damp_diag_rel; // random number of variation from base value for diagonal
    std::unique_ptr<ParameterDistribution>
        m_damp_offdiag_rel; // random number of variation from diagonal value for offdiagonal
};

} // namespace epi

#endif // UNCERTAINVALUE_H
