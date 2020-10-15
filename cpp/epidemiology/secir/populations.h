#ifndef POPULATIONS_H
#define POPULATIONS_H

#include <epidemiology/utils/uncertain_value.h>

#include <Eigen/Core>
#include <vector>

namespace epi
{

/**
 * @brief A class for compartment populations
 *
 * Populations can be split up into different categories, e.g. by
 * age group, yearly income group, gender etc. Compartmental models
 * introduce the additional category of infection type. For the SEIR
 * model these are Susceptible, Exposed, Infected and Removed.
 *
 * This class is a wrapper around a tensor that stores the individual populations
 * in each category and it provides some helper functions to get and set
 * the population of subgroups.
 *
 * It is implemented as a "flat tensor" of compartment populations
 *
 */
class Populations
{
public:
    /**
     * @brief Standard constructor of population parameters' class in the SECIR model
     */
    Populations();
    Populations(const Populations&) = default;
    Populations(Populations&&)      = default;
    Populations& operator=(const Populations&) = default;

    Populations(std::vector<size_t> const& category_sizes);

    /**
     * @brief get_num_compartments returns the number of compartments
     *
     * This corresponds to the product of the category sizes
     *
     * @return number of compartments
     */
    size_t get_num_compartments() const;

    /**
     * @brief get_category_sizes returns the vector of category sizes
     * @return vector of category sizes
     */
    std::vector<size_t> const& get_category_sizes() const;

    /**
     * @brief get_compartments returns a reference to the vector of populations
     * @return vector of populations
     */
    Eigen::VectorXd get_compartments() const;

    /**
     * @brief get returns the population of one compartment
     * @param indices a vector containing the indices for each category
     * @return the population of compartment
     */
    UncertainValue const& get(std::vector<size_t> const& indices) const;

    /**
     * @brief get returns the population of one compartment
     * @param indices a vector containing the indices for each category
     * @return the population of compartment
     */
    UncertainValue& get(std::vector<size_t> const& indices);

    /**
     * @brief get_group_total returns the total population of a group in one category
     * @param category_idx index of the category
     * @param group_idx index of the group
     * @return total population of the group
     */
    double get_group_total(size_t category_idx, size_t group_idx) const;

    /**
     * @brief get_total returns the total population of all compartments
     * @return total population
     */
    double get_total() const;

    /**
     * @brief set sets the scalar value of the population of one compartment
     * @param indices a vector containing the indices for each category
     * @param value the new value for the compartment's population
     */
    void set(std::vector<size_t> const& indices, UncertainValue const& value);

    /**
     * @brief set sets the scalar value of the population of one compartment
     * @param indices a vector containing the indices for each category
     * @param value the new value for the compartment's population
     */
    void set(std::vector<size_t> const& indices, double value);

    /**
     * @brief set sets the random distribution of the population of one compartment
     * @param indices a vector containing the indices for each category
     * @param value the new random distribution for the compartment's population
     */
    void set(std::vector<size_t> const& indices, ParameterDistribution const& dist);

    /**
     * @brief set_group_total sets the total population for a given group
     *
     * This function rescales all the compartments populations proportionally. If all compartments
     * have zero population, the total population gets distributed equally over all
     * compartments
     *
     * @param category_idx The index of the category of the group
     * @param group_idx The index of the group within the category
     * @param value the new value for the total population
     */
    void set_group_total(size_t category_idx, size_t group_idx, double value);

    /**
     * @brief set_difference_from_group_total sets the total population for a given group from a difference
     *
     * This function sets the population size 
     *
     * @param indices The indices of the compartment
     * @param category_idx The index of the category of the group
     * @param group_idx The index of the group within the category
     * @param total_population the new value for the total population
     */
    void set_difference_from_group_total(std::vector<size_t> const& indices, size_t category_idx, size_t group_idx,
                                         double total_population);

    /**
     * @brief set_total sets the total population
     *
     * This function rescales all the compartments populations proportionally. If all compartments
     * have zero population, the total population gets distributed equally over all
     * compartments. 
     *
     * @param value the new value for the total population
     */
    void set_total(double value);

    /**
     * @brief set_difference_from_total takes a total population as input and sets the compartment
     * of a given index to have the difference between the input total population and the current
     * population in all other compartments
     * @param indices the index of the compartment
     * @param total_population the new value for the total population
     */
    void set_difference_from_total(std::vector<size_t> const& indices, double total_population);

    /**
     * @brief get_flat_index returns the flat index into the stored population, given the
     * indices of each category
     * @param indices a vector of indices
     * @return a flat index into the data structure storing the compartment populations
     */
    size_t get_flat_index(std::vector<size_t> const& indices) const;

    /**
     * @brief checks whether the population Parameters satisfy their corresponding constraints and applys them
     */
    void apply_constraints();

    /**
     * @brief checks whether the population Parameters satisfy their corresponding constraints
     */
    void check_constraints() const;

private:
    // A vector storying the size of each category
    std::vector<size_t> m_category_sizes;

    // A vector containing the population of all compartments
    std::vector<UncertainValue> m_y;
};

} // namespace epi

#endif // POPULATIONS_H
