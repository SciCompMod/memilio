#ifndef UNCERTAINMATRIX_H
#define UNCERTAINMATRIX_H

#include "epidemiology/utils/memory.h"
#include "epidemiology/utils/parameter_distributions.h"
#include "epidemiology/secir/contact_matrix.h"

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
    * @brief Create an UncertainContactMatrix by cloning ContactMatrix 
    *        and distributions of another UncertainContactMatrix
    */
    UncertainContactMatrix(const UncertainContactMatrix& other);

    /**
    * @brief Set an UncertainContactMatrix from another UncertainContactMatrix, 
     *        containing a ContactMatrix and related distributions
     *        for changes in contact patterns over time.
    */
    UncertainContactMatrix& operator=(const UncertainContactMatrix& other);

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
     * @brief Sets the random distribution for the number of 
     *        contact pattern changes over time
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void set_distribution_damp_nb(const ParameterDistribution& dist);

    /**
     * @brief Sets the random distribution of the actual 
     *        points in time where the contact pattern changes
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void set_distribution_damp_days(const ParameterDistribution& dist);

    /**
     * @brief Sets the random distribution for a multiplicative base factor
     *        (changing the diagonal entries of the ContactMatrix)
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void set_distribution_damp_diag_base(const ParameterDistribution& dist);

    /**
     * @brief Sets the random distribution for the multiplicative factors
     *        changing each diagonal entry of the ContactMatrix
     *        relative to the base value sampled by the distribution from
     *        get_distribution_damp_diag_base()
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void set_distribution_damp_diag_rel(const ParameterDistribution& dist);

    /**
     * @brief Sets the random distribution for the multiplicative factors
     *        changing each offdiagonal entry of the ContactMatrix
     *        relative to the corresponding diagonal entries by the distribution 
     *        from get_distribution_damp_diag_rel()
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void set_distribution_damp_offdiag_rel(const ParameterDistribution& dist);

    /**
     * @brief Returns the random distribution for the number of 
     *        contact pattern changes over time
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_distribution_damp_nb();

    /**
     * @brief Returns the random distribution for the number of 
     *        contact pattern changes over time
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<const ParameterDistribution> get_distribution_damp_nb() const;

    /**
     * @brief Returns the random distribution of the actual 
     *        points in time where the contact pattern changes
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_distribution_damp_days();

    /**
     * @brief Returns the random distribution of the actual 
     *        points in time where the contact pattern changes
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<const ParameterDistribution> get_distribution_damp_days() const;

    /**
     * @brief Returns the random distribution for a multiplicative base factor
     *        (changing the diagonal entries of the ContactMatrix)
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_distribution_damp_diag_base();

    /**
     * @brief Returns the random distribution for a multiplicative base factor
     *        (changing the diagonal entries of the ContactMatrix)
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<const ParameterDistribution> get_distribution_damp_diag_base() const;

    /**
     * @brief Returns the random distribution for the multiplicative factors
     *        changing each diagonal entry of the ContactMatrix
     *        relative to the base value sampled by the distribution from
     *        get_distribution_damp_diag_base()
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_distribution_damp_diag_rel();

    /**
     * @brief Returns the random distribution for the multiplicative factors
     *        changing each diagonal entry of the ContactMatrix
     *        relative to the base value sampled by the distribution from
     *        get_distribution_damp_diag_base()
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<const ParameterDistribution> get_distribution_damp_diag_rel() const;

    /**
     * @brief Returns the random distribution for the multiplicative factors
     *        changing each offdiagonal entry of the ContactMatrix
     *        relative to the corresponding diagonal entries by the distribution 
     *        from get_distribution_damp_diag_rel()
     *        
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_distribution_damp_offdiag_rel();

    /**
     * @brief Returns the random distribution for the multiplicative factors
     *        changing each offdiagonal entry of the ContactMatrix
     *        relative to the corresponding diagonal entries by the distribution 
     *        from get_distribution_damp_diag_rel()
     *        
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<const ParameterDistribution> get_distribution_damp_offdiag_rel() const;

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

#endif // UNCERTAINMATRIX_H
