#ifndef UNCERTAINMATRIX_H
#define UNCERTAINMATRIX_H

#include "epidemiology/utils/memory.h"
#include "epidemiology/utils/parameter_distributions.h"
#include "epidemiology/secir/contact_matrix.h"
#include "epidemiology/secir/damping_sampling.h"

#include <memory>

namespace epi
{

/**
 * @brief The UncertainContactMatrix class consists of a
 *        ContactMatrix and certain distributions that describe
 *        the relative and uncertain changes of diagonal and 
 *        offdiagonal values over time 
 * 
 * The UncertainContactMatrix class represents a matrix-style model parameter 
 * that can take a ContactMatrix value but that is subjected to a uncertainty,
 * based on contact pattern changes realized by an uncertain number and shape of dampings.
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
    UncertainContactMatrix(size_t num_matrices = 1, Eigen::Index num_groups = 1);

    UncertainContactMatrix(const ContactMatrixGroup& cont_freq);

    /**
     * @brief Conversion to const ContactMatrix reference by returning the 
     *        ContactMatrix contained in UncertainContactMatrix
     */
    operator ContactMatrixGroup const &() const;

    /**
     * @brief Conversion to ContactMatrix reference by returning the 
     *        ContactMatrix contained in UncertainContactMatrix
     */
    operator ContactMatrixGroup&();

    /**
     * @brief Set an UncertainContactMatrix from a ContactMatrix, 
     *        all distributions remain unchanged.
     */
    UncertainContactMatrix& operator=(const ContactMatrixGroup& cont_freq);

    /**
     * @brief Returns the ContactMatrix reference 
     *        of the UncertainContactMatrix object
     */
    ContactMatrixGroup& get_cont_freq_mat();

    /**
     * @brief Returns the const ContactMatrix reference 
     *        of the UncertainContactMatrix object
     */
    ContactMatrixGroup const& get_cont_freq_mat() const;

    /**
     * @return list of damping samplings.
     */
    const std::vector<DampingSampling>& get_dampings() const
    {
        return m_dampings;
    }
    std::vector<DampingSampling>& get_dampings()
    {
        return m_dampings;
    }

    /**
     * @brief Sets the value by sampling from the distributions
     *
     * If no distribution is set, the value is not changed.
     * @param accum accumulating current and newly sampled dampings if true;
     *              default: false; removing all previously set dampings
     */
    ContactMatrixGroup draw_sample(bool accum = false);

private:
    ContactMatrixGroup m_cont_freq;
    std::vector<DampingSampling> m_dampings;
};

} // namespace epi

#endif // UNCERTAINMATRIX_H
