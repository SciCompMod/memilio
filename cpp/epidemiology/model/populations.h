#ifndef POPULATIONS_H
#define POPULATIONS_H

#include <epidemiology/utils/uncertain_value.h>
#include "epidemiology/utils/tensor_helpers.h"

#include <Eigen/Core>
#include <vector>
#include <array>
#include <numeric>

namespace
{

// Some metaprogramming to get the index of a given type in a parameter pack.
// Taken from https://stackoverflow.com/questions/26169198/how-to-get-the-index-of-a-type-in-a-variadic-type-pack
template <typename T, typename... Ts>
struct Index;

template <typename T, typename... Ts>
struct Index<T, T, Ts...> : std::integral_constant<std::size_t, 0> {
};

template <typename T, typename U, typename... Ts>
struct Index<T, U, Ts...> : std::integral_constant<std::size_t, 1 + Index<T, Ts...>::value> {
};

template <typename T, typename... Ts>
constexpr std::size_t Index_v = Index<T, Ts...>::value;

//calculate the product of a parameter pack
template <size_t... Ns>
constexpr size_t product()
{
    size_t p = 1;
    for (auto n : {Ns...})
        p *= n;
    return p;
}

} // namespace

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

template <class... CATEGORIES>
class Populations
{
public:
    using Index = std::tuple<CATEGORIES...>;

    Populations()
    {
    }

    /**
     * @brief get_num_compartments returns the number of compartments
     *
     * This corresponds to the product of the category sizes
     *
     * @return number of compartments
     */
    size_t get_num_compartments() const
    {
        return product<static_cast<size_t>(CATEGORIES::Count)...>();
    }

    /**
     * @brief get_compartments returns a copy of the vector of populations. This can be used
     * as initial conditions for the ODE solver
     * @return Eigen::VectorXd  of populations
     */
    Eigen::VectorXd get_compartments() const
    {
        Eigen::VectorXd m_y_eigen(m_y.size());
        for (auto i = 0; i < m_y.size(); i++) {
            m_y_eigen[i] = m_y[i];
        }

        return m_y_eigen;
    }

    /**
     * @brief get returns the population of one compartment
     * @param indices a vector containing the indices for each category
     * @return the population of compartment
     */
    UncertainValue& get(CATEGORIES... Cats)
    {
        return m_y[get_flat_index(Cats...)];
    }

    /**
     * @brief get returns the population of one compartment
     * @param indices a vector containing the indices for each category
     * @return the population of compartment
     */
    UncertainValue const& get(CATEGORIES... Cats) const
    {
        return m_y[get_flat_index(Cats...)];
    }

    template <typename Arr>
    static auto get_from(Arr& y, CATEGORIES... Cats)
    {
        return y[get_flat_index(Cats...)];
    }

    template <typename Arr>
    static auto const get_from(Arr const& y, CATEGORIES... Cats)
    {
        return y[get_flat_index(Cats...)];
    }

    /**
     * @brief get_group_total returns the total population of a group in one category
     * @param T enum class for the category
     * @param group_idx enum value of the group
     * @return total population of the group
     */
    template <class T>
    double get_group_total(T group_idx) const
    {
        //TODO maybe implement an iterator/visitor pattern rather than calculating indices?
        size_t idx          = static_cast<size_t>(group_idx);
        size_t category_idx = Index_v<T, CATEGORIES...>;

        double sum   = 0;
        auto indices = get_slice_indices(category_idx, idx, dimensions);
        for (auto i : indices) {
            sum += m_y[i];
        }
        return sum;
    }

    /**
     * @brief get_total returns the total population of all compartments
     * @return total population
     */
    double get_total() const
    {
        double sum = 0;
        for (auto i = 0; i < m_y.size(); i++) {
            sum += m_y[i];
        }
        return sum;
    }

    /**
     * @brief set sets the scalar value of the population of one compartment
     * @param indices a vector containing the indices for each category
     * @param value the new value for the compartment's population
     */
    void set(UncertainValue const& value, CATEGORIES... Cats)
    {
        m_y[get_flat_index(Cats...)] = value;
    }

    /**
     * @brief set sets the scalar value of the population of one compartment
     * @param indices a vector containing the indices for each category
     * @param value the new value for the compartment's population
     */
    void set(double value, CATEGORIES... Cats)
    {
        m_y[get_flat_index(Cats...)] = value;
    }

    /**
     * @brief set sets the random distribution of the population of one compartment
     * @param indices a vector containing the indices for each category
     * @param value the new random distribution for the compartment's population
     */
    void set(ParameterDistribution const& dist, CATEGORIES... Cats)
    {
        m_y[get_flat_index(Cats...)].set_distribution(dist);
    }

    /**
     * @brief set_group_total sets the total population for a given group
     *
     * This function rescales all the compartments populations proportionally. If all compartments
     * have zero population, the total population gets distributed equally over all
     * compartments
     *
     * @param T The enum class of the category of the group
     * @param group_idx The enum value of the group within the category
     * @param value the new value for the total population
     */
    template <class T>
    void set_group_total(double value, T group_idx)
    {
        double current_population = get_group_total<T>(group_idx);

        size_t idx          = static_cast<size_t>(group_idx);
        size_t category_idx = Index_v<T, CATEGORIES...>;

        //TODO slice indices are calcualated twice...
        auto indices = get_slice_indices(category_idx, idx, dimensions);

        if (fabs(current_population) < 1e-12) {
            for (auto i : indices) {
                m_y[i] = value / indices.size();
            }
        }
        else {
            for (auto i : indices) {
                m_y[i] *= value / current_population;
            }
        }
    }

    /**
     * @brief set_difference_from_group_total sets the total population for a given group from a difference
     *
     * This function sets the population size 
     *
     * @param total_population the new value for the total population
     * @param indices The indices of the compartment
     * @param T The enum class of the category of the group
     * @param group_idx The enum of the group within the category
     */
    template <class T>
    void set_difference_from_group_total(double total_group_population, T group_idx, CATEGORIES... Cats)

    {
        // TODO: https://stackoverflow.com/questions/20162903/template-parameter-packs-access-nth-type-and-nth-element
        // // is the given index part of the group?
        // assert(indices[category_idx] == group_idx);

        double current_population = get_group_total<T>(group_idx);
        size_t idx                = get_flat_index(Cats...);
        current_population -= m_y[idx];

        assert(current_population <= total_group_population);

        m_y[idx] = total_group_population - current_population;
    }

    /**
     * @brief set_total sets the total population
     *
     * This function rescales all the compartments populations proportionally. If all compartments
     * have zero population, the total population gets distributed equally over all
     * compartments. 
     *
     * @param value the new value for the total population
     */
    void set_total(double value)
    {
        double current_population = get_total();
        if (fabs(current_population) < 1e-12) {
            double ysize = double(m_y.size());
            for (auto i = 0; i < m_y.size(); i++) {
                m_y[i] = value / ysize;
            }
        }
        else {
            for (auto i = 0; i < m_y.size(); i++) {
                m_y[i] *= value / current_population;
            }
        }
    }

    /**
     * @brief set_difference_from_total takes a total population as input and sets the compartment
     * of a given index to have the difference between the input total population and the current
     * population in all other compartments
     * @param indices the index of the compartment
     * @param total_population the new value for the total population
     */
    void set_difference_from_total(double total_population, CATEGORIES... Cats)
    {
        double current_population = get_total();
        size_t idx                = get_flat_index(Cats...);
        current_population -= m_y[idx];

        assert(current_population <= total_population);

        m_y[idx] = total_population - current_population;
    }

    /**
     * @brief get_flat_index returns the flat index into the stored population, given the
     * indices of each category
     * @param indices a vector of indices
     * @return a flat index into the data structure storing the compartment populations
     */
    static size_t get_flat_index(CATEGORIES... Cats)
    {
        return flatten_index({static_cast<size_t>(Cats)...}, dimensions);
    }

    /**
     * @brief checks whether the population Parameters satisfy their corresponding constraints and applys them
     */
    void apply_constraints()
    {
        for (auto i = 0; i < m_y.size(); i++) {
            if (m_y[i] < 0) {
                log_warning("Constraint check: Compartment size {:d} changed from {:.4f} to {:d}", i, m_y[i], 0);
                m_y[i] = 0;
            }
        }
    }

    /**
     * @brief checks whether the population Parameters satisfy their corresponding constraints
     */
    void check_constraints() const
    {
        for (auto i = 0; i < m_y.size(); i++) {
            if (m_y[i] < 0) {
                log_error("Constraint check: Compartment size {:d} is {:.4f} and smaller {:d}", i, m_y[i], 0);
            }
        }
    }

    // An array storying the size of each category
    static std::array<size_t, sizeof...(CATEGORIES)> dimensions;

private:
    // A vector containing the population of all compartments
    std::array<UncertainValue, product<static_cast<size_t>(CATEGORIES::Count)...>()> m_y{};
};

// initialize array storying the size of each category
template <class... CATEGORIES>
std::array<size_t, sizeof...(CATEGORIES)> Populations<CATEGORIES...>::dimensions = {
    static_cast<size_t>(CATEGORIES::Count)...};

} // namespace epi

#endif // POPULATIONS_H
