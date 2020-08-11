#ifndef UNCERTAINMATRIX_H
#define UNCERTAINMATRIX_H

#include <memory>
#include <epidemiology/memory.h>
#include <epidemiology/damping.h>
#include <epidemiology/parameter_studies/parameter_distributions.h>

namespace epi
{

/**
 * @brief Initializes a Contact Frequency matrix for 
 **/
class ContactFrequencyMatrix
{
public:
    /**
     * @brief Standard constructor of 1x1-matrix of contact frequencies (i.e., one contact frequency)
     */
    ContactFrequencyMatrix();

    /**
     * @brief Constructor of a nb_groups x nb_groups-contact frequencies matrix
     * @param[in] nb_groups number of groups in the model
     */
    ContactFrequencyMatrix(size_t nb_groups);

    /**
     * @brief returns the size of the contact frequency matrix
     */
    int get_size() const;

    /**
     * @brief sets the contact frequency in the SECIR model; in case of multiple groups, set the contact rate cr_ij=cr_ji=cont_freq
     * @param cont_freq contact rate/frequency in 1/day unit
     * @param self_group own group
     * @param contact_group group which gets in contact with own group
     */
    void set_cont_freq(double cont_freq, int self_group, int contact_group);

    /**
     * @brief returns the contact frequency set for the SECIR model in 1/day unit; in case of multiple groups, returns the contact rate cr_ij=cr_ji
     */
    double get_cont_freq(int self_group, int contact_group) const;

    /**
     * @brief sets the damping in the SECIR model; in case of multiple groups, set the contact rate d_ij=d_ji=cont_freq
     * @param damping dampings over the whole time line in day unit
     * @param self_group own group
     * @param contact_group group which gets in contact with own group
     */
    void set_dampings(Dampings const& damping, int self_group, int contact_group);

    /**
     * @brief returns the dampings set for the SECIR model in 1/day unit; in case of multiple groups, returns the damping d_ij=d_ji
     */
    const Dampings& get_dampings(int self_group, int contact_group) const;

    /**
     * @brief add damping to the dampings object specified by self_ and contact_group
     * @param damping one damping in day unit
     * @param self_group own group
     * @param contact_group group which gets in contact with own group
     */
    void add_damping(Damping const& damping, int self_group, int contact_group);

    /**
     * @brief removes all previously set dampings
     */
    void clear_dampings();

private:
    std::vector<std::vector<double>> m_cont_freq;
    // This defines a damping factor for a mitigation strategy for different points in time.
    std::vector<std::vector<Dampings>> m_dampings;
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
    UncertainContactMatrix();

    UncertainContactMatrix(ContactFrequencyMatrix cont_freq);

    /**
    * @brief Create an UncertainContactMatrix by cloning ContactFrequencyMatrix 
    *        and distributions of another UncertainContactMatrix
    */
    UncertainContactMatrix(const UncertainContactMatrix& other);

    /**
    * @brief Set an UncertainContactMatrix from another UncertainContactMatrix, 
     *        containing a ContactFrequencyMatrix and related distributions
     *        for changes in contact patterns over time.
    */
    UncertainContactMatrix& operator=(const UncertainContactMatrix& other);

    /**
     * @brief Conversion to const ContactFrequencyMatrix reference by returning the 
     *        ContactFrequencyMatrix contained in UncertainContactMatrix
     */
    operator ContactFrequencyMatrix const &() const;

    /**
     * @brief Conversion to ContactFrequencyMatrix reference by returning the 
     *        ContactFrequencyMatrix contained in UncertainContactMatrix
     */
    operator ContactFrequencyMatrix&();

    /**
     * @brief Set an UncertainContactMatrix from a ContactFrequencyMatrix, 
     *        all distributions remain unchanged.
     */
    UncertainContactMatrix& operator=(ContactFrequencyMatrix cont_freq);

    /**
     * @brief Returns the ContactFrequencyMatrix reference 
     *        of the UncertainContactMatrix object
     */
    ContactFrequencyMatrix& get_cont_freq_mat();

    /**
     * @brief Returns the const ContactFrequencyMatrix reference 
     *        of the UncertainContactMatrix object
     */
    ContactFrequencyMatrix const& get_cont_freq_mat() const;

    /**
     * @brief Sets the random distribution for the number of 
     *        contact pattern changes over time
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void set_dist_damp_nb(const ParameterDistribution& dist);

    /**
     * @brief Sets the random distribution of the actual 
     *        points in time where the contact pattern changes
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void set_dist_damp_days(const ParameterDistribution& dist);

    /**
     * @brief Sets the random distribution for a multiplicative base factor
     *        (changing the diagonal entries of the ContactFrequencyMatrix)
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void set_dist_damp_diag_base(const ParameterDistribution& dist);

    /**
     * @brief Sets the random distribution for the multiplicative factors
     *        changing each diagonal entry of the ContactFrequencyMatrix
     *        relative to the base value sampled by the distribution from
     *        get_dist_damp_diag_base()
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void set_dist_damp_diag_rel(const ParameterDistribution& dist);

    /**
     * @brief Sets the random distribution for the multiplicative factors
     *        changing each offdiagonal entry of the ContactFrequencyMatrix
     *        relative to the corresponding diagonal entries by the distribution 
     *        from get_dist_damp_diag_rel()
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void set_dist_damp_offdiag_rel(const ParameterDistribution& dist);

    /**
     * @brief Returns the random distribution for the number of 
     *        contact pattern changes over time
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_dist_damp_nb() const;

    /**
     * @brief Returns the random distribution of the actual 
     *        points in time where the contact pattern changes
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_dist_damp_days() const;

    /**
     * @brief Returns the random distribution for a multiplicative base factor
     *        (changing the diagonal entries of the ContactFrequencyMatrix)
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_dist_damp_diag_base() const;

    /**
     * @brief Returns the random distribution for the multiplicative factors
     *        changing each diagonal entry of the ContactFrequencyMatrix
     *        relative to the base value sampled by the distribution from
     *        get_dist_damp_diag_base()
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_dist_damp_diag_rel() const;

    /**
     * @brief Returns the random distribution for the multiplicative factors
     *        changing each offdiagonal entry of the ContactFrequencyMatrix
     *        relative to the corresponding diagonal entries by the distribution 
     *        from get_dist_damp_diag_rel()
     *        
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_dist_damp_offdiag_rel() const;

    /**
     * @brief Sets the value by sampling from the distributions
     *
     * If no distribution is set, the value is not changed.
     * @param accum accumulating current and newly sampled dampings if true;
     *              default: false; removing all previously set dampings
     */
    ContactFrequencyMatrix draw_sample(bool accum = false);

private:
    ContactFrequencyMatrix m_cont_freq;
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
