#ifndef POPULATIONS_H
#define POPULATIONS_H

#include <epidemiology/utils/uncertain_value.h>
#include "epidemiology/utils/tensor_helpers.h"
#include "epidemiology/utils/ScalarType.h"

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

//calculate the product of a size_t parameter pack
template <size_t...>
struct product;

template <>
struct product<> {
    static constexpr size_t value = 1;
};

template <size_t i, size_t... tail>
struct product<i, tail...> {
    static constexpr size_t value = i * product<tail...>::value;
};

} // namespace

namespace epi
{

/**
 * @brief A class template for compartment populations
 *
 * Populations can be split up into different categories, e.g. by
 * age group, yearly income group, gender etc. Compartmental models
 * introduce the additional category of infection type. For the SEIR
 * model these are Susceptible, Exposed, Infected and Removed. Every category
 * is assumed to contain a finite number of groups.
 *
 * This template can be instantiated given an arbitrary set of categories.
 * The categories are assumed to be an enum class with a member Count refering to
 * the number of elements in the enum.
 *
 * The class created from this template contains a "flat array" of compartment
 * populations and some functions for retrieving or setting the populations.
 *
 */

template <class... CATEGORIES>
class Populations
{
public:
    using Type = UncertainValue;

    // This type can be used by other classes to refer to a concrete compartment
    using Index = std::tuple<CATEGORIES...>;

    /**
     * @brief Populations default constructor
     */
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
    static size_t constexpr get_num_compartments()
    {
        return size;
    }

    /**
     * @brief get_compartments returns an Eigen copy of the vector of populations. This can be used
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
     * @param Cats enum values for each category
     * @return the population of compartment
     */
    Type& get(CATEGORIES... Cats)
    {
        return m_y[get_flat_index(Cats...)];
    }

    /**
     * @brief get returns the population of one compartment
     * @param Cats enum values for each category
     * @return the population of compartment
     */
    Type const& get(CATEGORIES... Cats) const
    {
        return m_y[get_flat_index(Cats...)];
    }

    /**
     * @brief get_grom returns the value of a flat container
     * at the flat index corresponding to set of enum values.
     * It is the same as get, except that it takes the values
     * from an outside reference flat container, rather than the
     * initial values stored within this class
     * @param y a reference to a flat container
     * @param Cats emi, va;ies fpr eacj category
     * @return the population of compartment
     */
    template <typename Arr>
    static auto get_from(Arr& y, CATEGORIES... Cats)
    {
        return y[get_flat_index(Cats...)];
    }

    /**
     * @brief get_grom returns the value of a flat container
     * at the flat index corresponding to set of enum values.
     * It is the same as get, except that it takes the values
     * from an outside reference flat container, rather than the
     * initial values stored within this class
     * @param y a reference to a flat container
     * @param Cats emi, va;ies fpr eacj category
     * @return the population of compartment
     */
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
    ScalarType get_group_total(T group_idx) const
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
    ScalarType get_total() const
    {
        return std::accumulate(m_y.begin(), m_y.end(), ScalarType(0.));
    }

    /**
     * @brief set sets the scalar value of the population of one compartment
     * @param indices a vector containing the indices for each category
     * @param value the new value for the compartment's population
     */
    void set(Type const& value, CATEGORIES... Cats)
    {
        m_y[get_flat_index(Cats...)] = value;
    }

    /**
     * @brief set sets the scalar value of the population of one compartment
     * @param indices a vector containing the indices for each category
     * @param value the new value for the compartment's population
     */
    void set(ScalarType value, CATEGORIES... Cats)
    {
        m_y[get_flat_index(Cats...)] = value;
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
    void set_group_total(ScalarType value, T group_idx)
    {
        ScalarType current_population = get_group_total(group_idx);

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
    void set_difference_from_group_total(ScalarType total_group_population, T group_idx, CATEGORIES... Cats)

    {
        ScalarType current_population = get_group_total(group_idx);
        size_t idx                    = get_flat_index(Cats...);
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
    void set_total(ScalarType value)
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
    static size_t constexpr size = product<static_cast<size_t>(CATEGORIES::Count)...>::value;

private:
    // A vector containing the population of all compartments
    std::array<Type, size> m_y{};
};

// initialize array storying the size of each category
template <class... CATEGORIES>
std::array<size_t, sizeof...(CATEGORIES)> Populations<CATEGORIES...>::dimensions = {
    static_cast<size_t>(CATEGORIES::Count)...};

} // namespace epi

#endif // POPULATIONS_H
