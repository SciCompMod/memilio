#ifndef CUSTOMINDEXARRAY_H
#define CUSTOMINDEXARRAY_H

#include <epidemiology/utils/uncertain_value.h>
#include "epidemiology/utils/tensor_helpers.h"
#include "epidemiology/utils/ScalarType.h"
#include "epidemiology/utils/eigen.h"

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
struct Index<T, T, Ts...> : std::integral_constant<int, 0> {
};

template <typename T, typename U, typename... Ts>
struct Index<T, U, Ts...> : std::integral_constant<int, 1 + Index<T, Ts...>::value> {
};

template <typename T, typename... Ts>
constexpr int Index_v = Index<T, Ts...>::value;

//calculate the product of a integer parameter pack
template <int...>
struct product;

template <>
struct product<> {
    static constexpr int value = 1;
};

template <int i, int... tail>
struct product<i, tail...> {
    static constexpr int value = i * product<tail...>::value;
};


//some metaprogramming to transform a tuple into a parameter pack and use it as
//an argument in a function.
//Taken from https://stackoverflow.com/questions/7858817/unpacking-a-tuple-to-call-a-matching-function-pointer/9288547#9288547

template <typename Function, typename Tuple, size_t... I>
decltype(auto) call(Function f, Tuple t, std::index_sequence<I...>)
{
    return f(std::get<I>(t)...);
}

template <typename Function, typename Tuple>
decltype(auto) call(Function f, Tuple t)
{
    static constexpr auto size = std::tuple_size<Tuple>::value;
    return call(f, t, std::make_index_sequence<size>{});
}

} // namespace

namespace epi
{

/**
 * @brief A class template for an array with custom indices
 *
 * This class stores an array of elements that can be queried using
 * a variadic number of custom index types. Each index type is associated
 * with a category, or dimension into a multidimensional array.
 *
 * Each custom index type must have a nested Count element that is castable
 * to an integer. Typically it is assumed, that custom indices are enums.
 *
 * Example:
 *
 * enum class AgeGroup
 * {
 *    Young,
 *    Old,
 *    Count = 2
 * };
 *
 * enum class Gender
 * {
 *    Female,
 *    Male,
 *    Diverse,
 *    Count = 3
 * };
 *
 * CustomIndexArray<size_t, AgeGroup, Gender> populations;
 *
 * Here, populations represents a 2x3 size_t array (though the data is stored contigously).
 * An element can be accessed using a flat row-major index or by using categories:
 *
 * auto x = populations.get(4);
 * auto y = populations.get(AgeGroup::Old, Gender::Male);
 * assert(x == y);
 *
 * @tparam Typ the type stored in the array
 * @tparam Categories The custom Index types (enums)
 *
 */

template <class Typ, class... Categories>
class CustomIndexArray
{
public:

    using Type  = Typ;
    using Index = std::tuple<Categories...>;
    using N     = product<static_cast<size_t>(Categories::Count)...>;
    using InternalArrayType = Eigen::Array<Type, N::value, 1>;

    // An array storying the size of each category
    static std::array<size_t, sizeof...(Categories)> dimensions;

    /**
     * @brief CustomIndexArray default constructor
     *
     * leaves entries uninitialized
     */
    CustomIndexArray() = default;

    /**
     * @brief CustomIndexArray constructor, that initializes the array
     * to constant instances of `CustsomIndexArray::Type`.
     *
     * It forwards the arguments for initializer_list construction of the
     * contained objects.
     *
     * @tparam Ts The argument types of the constructor arguments of Type
     * @param args The constructor arguments of Type
     */
    template <class... Ts,
              typename std::enable_if_t<std::is_constructible<Type, Ts...>::value>* = nullptr>
    CustomIndexArray(Ts&&... args)
        : m_y(InternalArrayType::Constant(N::value, 1, {std::forward<Ts>(args)...}))
    {
    }

    /**
     * @brief get_num_compartments returns the number of compartments
     *
     * This corresponds to the product of the category sizes
     *
     * @return number of compartments
     */
    static size_t constexpr size()
    {
        return N::value;
    }

    /**
     * @brief get_compartments returns an Eigen copy of the vector of populations. This can be used
     * as initial conditions for the ODE solver
     * @return Eigen::VectorXd  of populations
     */
    auto const& array() const
    {
        return m_y;
    }
    auto& array()
    {
        return m_y;
    }

    /**
     * @brief returns the entry of the array given a flat index index
     * @param index a flat index
     * @return the value at the index
     */
    Type& operator[](Index const& index) {
        return m_y[(Eigen::Index)call(get_flat_index, index)];
    }

    /**
     * @brief returns the entry of the array given a flat index index
     * @param index a flat index
     * @return the value at the index
     */
    Type const& operator[](Index const& index) const {
        return m_y[(Eigen::Index)call(get_flat_index, index)];
    }


    /**
     * @brief get_flat_index returns the flat index into the stored array, given the
     * indices of each category
     * @param indices the custom indices for each category
     * @return a flat index into the data structure storing the compartment populations
     */
    static size_t get_flat_index(Categories... cats)
    {
        return flatten_index({static_cast<size_t>(cats)...}, dimensions);
    }


protected:
    // An array containing the elements
    InternalArrayType m_y{};
};

// initialize array storying the size of each category
template <class Type, class... Categories>
std::array<size_t, sizeof...(Categories)> CustomIndexArray<Type, Categories...>::dimensions = {
    static_cast<size_t>(Categories::Count)...};

} // namespace epi

#endif // CUSTOMINDEXARRAY_H
