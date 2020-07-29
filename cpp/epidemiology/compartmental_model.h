#pragma once

#include <vector>
#include <numeric>
#include "Eigen/Core"

namespace epi
{

/**
 * @brief Base class for compartmental models
 *
 * A compartmental model consists of compartments of a population
 * and migrations between two compartments.
 *
 * A compartmental model can have several categories, such as
 * infection status, age group or region. Each category can have
 * a different number of groups. The total number of comparments
 * is the product of all the numbers of groups of each compartment
 *
 * This is implemented as a "flat tensor" of compartment populations
 *
 */
class CompartmentalModel
{
public:
    CompartmentalModel(std::vector<size_t> const& category_sizes);

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
    Eigen::VectorXd const& get_compartments() const;

    /**
     * @brief get returns the population of one compartment
     * @param indices an initializer list containing the indices for each category
     * @return the population of compartment
     */
    double get(std::vector<size_t> const& indices) const;

    /**
     * @brief get_total_population returns the total population of all compartments
     * @return total population
     */
    double get_total_population() const;

    /**
     * @brief set sets the population of one compartment
     * @param indices an initializer list containing the indices for each category
     * @param value the new value for the compartment's population
     */
    void set(std::initializer_list<size_t> const& indices, double value);

    /**
     * @brief set_total_populaton sets the total population
     *
     * It rescales all the compartments populations proportionally. If all compartments
     * have zero population, the total population gets distributed equally over all
     * compartments
     *
     * @param value the new value for the total population
     */
    void set_total_populaton(double value);

    // TODO: getters and setters for slices, ranges subsets etc., e.g. for all infected and
    // hospitalized in in the age group 18-25 in Koeln, Bonn and Rhein-Sieg-Kreis

private:
    // A vector storying the size of each category
    std::vector<size_t> m_category_sizes;

    // A vector containing the population of all compartments
    Eigen::VectorXd m_y;
};

} // namespace epi
