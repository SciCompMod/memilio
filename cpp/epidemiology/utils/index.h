#ifndef INDEX_H
#define INDEX_H

#include "epidemiology/utils/type_safe.h"

namespace epi
{

template  <typename... CategoryTags> class Index;

/**
 * @brief An Index with a single template parameter is a typesafe wrapper for size_t
 * that is associated with a Tag. It is used to index into a CustomIndexArray
 *
 * CustomIndexArray<Tag1, Tag2> a;
 * a[{Index<Tag1>(0), Index<Tag2>(14)}]
 *
 * Will retrieve the element associated with the indices 0 and 14 for the Tags Tag1
 * and Tag2 respectively.
 *
 * Optionally, the tag can be derived from Index to shorten the notation:
 *
 * struct Tag1 : public Index<Tag1>{
 *    Tag1(size_t val) : Index<Tag1>(val) {}
 * };
 * struct Tag2 : public Index<Tag2>{
 *    Tag2(size_t val) : Index<Tag2>(val) {}
 * };
 * CustomIndexArray<Tag1, Tag2> a;
 * a[{Tag1(0), Tag2(14)}]
 *
 * @tparam CategoryTag A tag for the typesafe index
 *
 */
template <typename CategoryTag>
class Index<CategoryTag> : public TypeSafe<size_t, Index<CategoryTag>>
                         , public OperatorComparison<Index<CategoryTag>>
                         , public OperatorAdditionSubtraction<Index<CategoryTag>>
{
public:
    using TypeSafe<size_t, Index<CategoryTag>>::TypeSafe;

    static size_t constexpr size = 1;

    /**
     * @brief Constructor from enum, if CategoryTag is an enum
     */
    template <typename Dummy = CategoryTag,
                  std::enable_if_t<std::is_enum<Dummy>::value, void>* = nullptr>
    Index(Dummy val) : TypeSafe<size_t, Index<CategoryTag>>((size_t)val) {}

    /**
     * @brief Constructor from size_t
     * @param val
     */
    explicit Index(size_t val) : TypeSafe<size_t, Index<CategoryTag>>(val) {}

};


/**
 * @brief An Index with more than one template parameter combines several Index objects.
 * It is used to index into a multidimensional CustomIndexArray
 *
 * @tparam CategoryTag Variadic template parameter for the Tags used in the MultiIndex
 */
template <typename... CategoryTag>
class Index
{
public:
    static size_t constexpr size = sizeof...(CategoryTag);

    // constructor from Indices
    Index(Index<CategoryTag> const&..._indices) : indices{_indices...} {}

    // comparison operators
    bool operator==(Index const& other) const
    {
        return indices == other.indices;
    }

    bool operator!=(Index const& other) const
    {
        return !(this == other);
    }

    std::tuple<Index<CategoryTag>...> indices;
};


// retrieves the Index at the Ith position for a Index with more than one Tag
template <size_t I, typename... CategoryTags,
          std::enable_if_t<(sizeof...(CategoryTags) > 1), void>* = nullptr>
constexpr typename std::tuple_element<I, std::tuple<Index<CategoryTags>...> >::type& get(Index<CategoryTags...>& i) noexcept
{
    return std::get<I>(i.indices);
}

// retrieves the Index at the Ith position for a Index with one Tag, equals identity function
template <size_t I, typename... CategoryTags,
          std::enable_if_t<(sizeof...(CategoryTags) == 1), void>* = nullptr>
constexpr typename std::tuple_element<I, std::tuple<Index<CategoryTags>...> >::type& get(Index<CategoryTags...>& i) noexcept
{
    static_assert(I==0, "I must be equal to zero for an Index with just one template parameter");
    return i;
}

// retrieves the Index at the Ith position for a Index with more than one Tag const version
template <size_t I, typename... CategoryTags,
          std::enable_if_t<(sizeof...(CategoryTags) > 1), void>* = nullptr>
constexpr typename std::tuple_element<I, std::tuple<Index<CategoryTags>...> >::type const& get(Index<CategoryTags...> const& i) noexcept
{
    return std::get<I>(i.indices);
}

// retrieves the Index at the Ith position for a Index with one Tag, equals identity function const version
template <size_t I, typename... CategoryTags,
          std::enable_if_t<(sizeof...(CategoryTags) == 1), void>* = nullptr>
constexpr typename std::tuple_element<I, std::tuple<Index<CategoryTags>...> >::type const& get(Index<CategoryTags...> const& i) noexcept
{
    static_assert(I==0, "I must be equal to zero for an Index with just one template parameter");
    return i;
}

// retrieves the Index for the tag Tag of a Index with more than one Tag
template <typename Tag, typename... CategoryTags,
          std::enable_if_t<(sizeof...(CategoryTags) > 1), void>* = nullptr>
constexpr Index<Tag>& get(Index<CategoryTags...>& i) noexcept
{
    return std::get<Index<Tag>>(i.indices);
}

// retrieves the Index for the tag Tag of a Index with one Tag, equals identity function
template <typename Tag, typename... CategoryTags,
          std::enable_if_t<(sizeof...(CategoryTags) == 1), void>* = nullptr>
constexpr Index<Tag>& get(Index<CategoryTags...>& i) noexcept
{
    using IndexTag = std::tuple_element_t<0, std::tuple<CategoryTags...>>;
    static_assert(std::is_same<Tag, IndexTag>::value, "Tags must match for an Index with just one template parameter");
    return i;
}

// retrieves the Index for the tag Tag for a Index with more than one Tag const version
template <typename Tag, typename... CategoryTags,
          std::enable_if_t<(sizeof...(CategoryTags) > 1), void>* = nullptr>
constexpr Index<Tag> const& get(Index<CategoryTags...> const& i) noexcept
{
    return std::get<Index<Tag>>(i.indices);
}

// retrieves the Index for the tag Tag for a Index with one Tag, equals identity function const version
template <typename Tag, typename... CategoryTags,
          std::enable_if_t<(sizeof...(CategoryTags) == 1), void>* = nullptr>
constexpr Index<Tag> const& get(Index<CategoryTags...> const& i) noexcept
{
    using IndexTag = std::tuple_element_t<0, std::tuple<CategoryTags...>>;
    static_assert(std::is_same<Tag, IndexTag>::value, "Tags must match for an Index with just one template parameter");
    return i;
}


} // namespace epi

#endif
