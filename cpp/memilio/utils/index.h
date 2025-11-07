/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef INDEX_H
#define INDEX_H

#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/type_safe.h"

namespace mio
{

template <typename... CategoryTags>
class Index;

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
class MEMILIO_ENABLE_EBO Index<CategoryTag> : public TypeSafe<size_t, Index<CategoryTag>>,
                                              public OperatorComparison<Index<CategoryTag>>,
                                              public OperatorAdditionSubtraction<Index<CategoryTag>>,
                                              public OperatorScalarMultiplicationDivision<Index<CategoryTag>, size_t>
{
public:
    using TypeSafe<size_t, Index<CategoryTag>>::TypeSafe;

    static size_t constexpr size = 1;

    static Index constexpr Zero()
    {
        return Index((size_t)0);
    }

    /**
     * @brief Constructor from enum, if CategoryTag is an enum
     */
    template <typename Dummy = CategoryTag, std::enable_if_t<std::is_enum<Dummy>::value, void>* = nullptr>
    Index(Dummy val)
        : TypeSafe<size_t, Index<CategoryTag>>((size_t)val)
    {
    }

    /**
     * @brief Constructor from size_t
     * @param val
     */
    explicit Index(size_t val)
        : TypeSafe<size_t, Index<CategoryTag>>(val)
    {
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        mio::serialize(io, size_t(*this));
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Index> deserialize(IOContext& io)
    {
        BOOST_OUTCOME_TRY(auto&& i, mio::deserialize(io, Tag<size_t>{}));
        return success(Index(i));
    }
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

    static Index constexpr Zero()
    {
        return Index(Index<CategoryTag>::Zero()...);
    }

    // constructor from Indices
    Index(Index<CategoryTag> const&... _indices)
        : indices{_indices...}
    {
    }

private:
    Index(const std::tuple<Index<CategoryTag>...>& _indices)
        : indices(_indices)
    {
    }

public:
    // comparison operators
    bool operator==(Index const& other) const
    {
        return indices == other.indices;
    }

    bool operator!=(Index const& other) const
    {
        return !(this == other);
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("MultiIndex");
        obj.add_element("Indices", indices);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Index> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("MultiIndex");
        auto tup = obj.expect_element("Indices", Tag<decltype(indices)>{});
        return mio::apply(
            io,
            [](auto&& tup_) {
                return Index(tup_);
            },
            tup);
    }

    std::tuple<Index<CategoryTag>...> indices;
};

/// Specialization of type_at_index for Index. @see type_at_index.
template <size_t Tag, class... CategoryTags>
struct type_at_index<Tag, ::mio::Index<CategoryTags...>> : public type_at_index<Tag, CategoryTags...> {
};

/// Specialization of index_of_type for Index. @see index_of_type.
template <class Tag, class... CategoryTags>
struct index_of_type<Tag, ::mio::Index<CategoryTags...>> : public index_of_type<Tag, CategoryTags...> {
};

/// Specialization of index_of_type for Index. Resolves ambiguity when using Index%s as items. @see index_of_type.
template <class... CategoryTags>
struct index_of_type<Index<CategoryTags...>, Index<CategoryTags...>> {
    static constexpr std::size_t value = 0;
};

// retrieves the Index at the Ith position for a Index with more than one Tag
template <size_t I, typename... CategoryTags, std::enable_if_t<(sizeof...(CategoryTags) > 1), void>* = nullptr>
constexpr typename std::tuple_element<I, std::tuple<Index<CategoryTags>...>>::type&
get(Index<CategoryTags...>& i) noexcept
{
    return std::get<I>(i.indices);
}

// retrieves the Index at the Ith position for a Index with one Tag, equals identity function
template <size_t I, typename... CategoryTags, std::enable_if_t<(sizeof...(CategoryTags) == 1), void>* = nullptr>
constexpr typename std::tuple_element<I, std::tuple<Index<CategoryTags>...>>::type&
get(Index<CategoryTags...>& i) noexcept
{
    static_assert(I == 0, "I must be equal to zero for an Index with just one template parameter");
    return i;
}

// retrieves the Index at the Ith position for a Index with more than one Tag const version
template <size_t I, typename... CategoryTags, std::enable_if_t<(sizeof...(CategoryTags) > 1), void>* = nullptr>
constexpr typename std::tuple_element<I, std::tuple<Index<CategoryTags>...>>::type const&
get(Index<CategoryTags...> const& i) noexcept
{
    return std::get<I>(i.indices);
}

// retrieves the Index at the Ith position for a Index with one Tag, equals identity function const version
template <size_t I, typename... CategoryTags, std::enable_if_t<(sizeof...(CategoryTags) == 1), void>* = nullptr>
constexpr typename std::tuple_element<I, std::tuple<Index<CategoryTags>...>>::type const&
get(Index<CategoryTags...> const& i) noexcept
{
    static_assert(I == 0, "I must be equal to zero for an Index with just one template parameter");
    return i;
}

// retrieves the Index for the tag Tag of a Index with more than one Tag
template <typename Tag, typename... CategoryTags, std::enable_if_t<(sizeof...(CategoryTags) > 1), void>* = nullptr>
constexpr Index<Tag>& get(Index<CategoryTags...>& i) noexcept
{
    return std::get<Index<Tag>>(i.indices);
}

// retrieves the Index for the tag Tag of a Index with one Tag, equals identity function
template <typename Tag, typename... CategoryTags, std::enable_if_t<(sizeof...(CategoryTags) == 1), void>* = nullptr>
constexpr Index<Tag>& get(Index<CategoryTags...>& i) noexcept
{
    using IndexTag = std::tuple_element_t<0, std::tuple<CategoryTags...>>;
    static_assert(std::is_same<Tag, IndexTag>::value, "Tags must match for an Index with just one template parameter");
    return i;
}

// retrieves the Index for the tag Tag for a Index with more than one Tag const version
template <typename Tag, typename... CategoryTags, std::enable_if_t<(sizeof...(CategoryTags) > 1), void>* = nullptr>
constexpr Index<Tag> const& get(Index<CategoryTags...> const& i) noexcept
{
    return std::get<Index<Tag>>(i.indices);
}

// retrieves the Index for the tag Tag for a Index with one Tag, equals identity function const version
template <typename Tag, typename... CategoryTags, std::enable_if_t<(sizeof...(CategoryTags) == 1), void>* = nullptr>
constexpr Index<Tag> const& get(Index<CategoryTags...> const& i) noexcept
{
    using IndexTag = std::tuple_element_t<0, std::tuple<CategoryTags...>>;
    static_assert(std::is_same<Tag, IndexTag>::value, "Tags must match for an Index with just one template parameter");
    return i;
}

namespace details
{
/// @brief Extracts CategoryTags from the tagged Index and returns a subindex of SuperIndex with the given categories.
template <class... CategoryTags, class SuperIndex>
inline Index<CategoryTags...> reduce_index_impl(const SuperIndex& i, mio::Tag<Index<CategoryTags...>>)
{
    // the subindex may not be trivially constructible, so we pass its type using mio::Tag
    // the type has to be passed as an argument to determine its CategoryTags

    // below, we use get<index_of_type<>> instead of get<> directly to handle categories that are not unique
    // (that is, `get<CategoryTags>(i)...` fails to compile for SuperIndex=Index<T, T>)
    return Index<CategoryTags...>{get<index_of_type_v<CategoryTags, SuperIndex>>(i)...};
}

/// @brief Creates and returns a SuperIndex from SubIndex, using entries from the given SubIndex or fill_value.
template <class... CategoryTags, class... Subset>
inline Index<CategoryTags...> extend_index_impl(const Index<Subset...>& i, const size_t fill_value,
                                                mio::Tag<Index<CategoryTags...>>)
{
    using SuperIndex = Index<CategoryTags...>;
    using SubIndex   = Index<Subset...>;
    // The superindex may not be trivially constructible, so we pass its type using mio::Tag.
    // The type has to be passed as an argument to determine its CategoryTags.

    return SuperIndex{[&]() {
        // This is an IIFE, which is invoked for each category (note the '...' after the function call).
        // So CategoryTags without a '...' is seen by each IIFE as exactly one category from this variadic template.
        if constexpr (is_type_in_list_v<CategoryTags, Subset...>) {
            // We use get<index_of_type<>> instead of get<> directly to handle categories that are not unique
            // (that is, `get<CategoryTags>(i)...` fails to compile for SuperIndex=Index<T, T>)
            return get<index_of_type_v<CategoryTags, SubIndex>>(i);
        }
        else {
            return Index<CategoryTags>(fill_value);
        }
    }()...};
}
} // namespace details

/**
 * @brief Create a SubIndex by copying values from SuperIndex.
 * If a type T is contained multiple times in SuperIndex, only the first occurance of T is used.
 * For example, `reduce_index<Index<T, T>>(Index<T, T, T>{1,2,3})` returns `{1,1}`.
 * @param[in] index An instance of SuperIndex
 * @tparam SubIndex An Index that contains a subset of the categories from SuperIndex.
 * @tparam SuperIndex Any Index.
 * @return A (sub)index with the given categories and values from index.
 */
template <class SubIndex, class SuperIndex>
SubIndex reduce_index(const SuperIndex& index)
{
    return details::reduce_index_impl(index, mio::Tag<SubIndex>{});
}

/**
 * @brief Create a SuperIndex by copying values from SubIndex, filling new categories with fill_value.
 * If a type T is contained multiple times in SubIndex, only the first occurance of T is used.
 * For example, `extend_index<Index<T, T, T>>(Index<T, T>{1,2})` returns `{1,1,1}`.
 * @param[in] index An instance of SubIndex
 * @param[in] fill_value The value to use for categories not in SubIndex.
 * @tparam SuperIndex Any Index.
 * @tparam SubIndex An Index that contains a subset of the categories from SuperIndex.
 * @return A (super)index with the given categories and values from index.
 */
template <class SuperIndex, class SubIndex>
SuperIndex extend_index(const SubIndex& index, size_t fill_value = 0)
{
    return details::extend_index_impl(index, fill_value, mio::Tag<SuperIndex>{});
}

} // namespace mio

#endif
