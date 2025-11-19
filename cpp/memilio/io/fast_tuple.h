#ifndef MIO_FAST_TUPLE_H
#define MIO_FAST_TUPLE_H

#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/metaprogramming.h"

#include <cstddef>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

namespace mio
{

namespace details
{

template <size_t Index, class Type, bool = std::is_empty_v<Type>>
struct FastTupleContainer;

template <size_t Index, class Type>
struct MEMILIO_ENABLE_EBO FastTupleContainer<Index, Type, true> : Type {
    constexpr FastTupleContainer()
        : Type()
    {
    }

    constexpr FastTupleContainer(const Type& t)
        : Type(t)
    {
    }

    template <class T>
    constexpr FastTupleContainer(T&& t)
        : Type(std::forward<T>(t))
    {
    }
};
template <size_t Index, class Type>
struct FastTupleContainer<Index, Type, false> {
    constexpr FastTupleContainer()
        : data()
    {
    }

    constexpr FastTupleContainer(const Type& t)
        : data(t)
    {
    }

    template <class T>
    constexpr FastTupleContainer(T&& t)
        : data(std::forward<T>(t))
    {
    }

    Type data;
};

/**
 * @brief Get the underlying type from a FastTupleContainer.
 * Using that the FastTupleImpl inherits from FastTupleContainer%s with unique indices, this getter only needs an index
 * to directly obtain the correct container from a FastTupleImpl.
 */
template <size_t Index, class Type>
constexpr Type& get_ftc(FastTupleContainer<Index, Type, true>& c)
{
    return c;
}
template <size_t Index, class Type>
constexpr const Type& get_ftc(const FastTupleContainer<Index, Type, true>& c)
{
    return c;
}
template <size_t Index, class Type>
constexpr Type&& get_ftc(FastTupleContainer<Index, Type, true>&& c)
{
    return c;
}

template <size_t Index, class Type>
constexpr Type& get_ftc(FastTupleContainer<Index, Type, false>& c)
{
    return c.data;
}
template <size_t Index, class Type>
constexpr const Type& get_ftc(const FastTupleContainer<Index, Type, false>& c)
{
    return c.data;
}
template <size_t Index, class Type>
constexpr Type&& get_ftc(FastTupleContainer<Index, Type, false>&& c)
{
    return c.data;
}

template <class Indices, class... Types>
struct MEMILIO_ENABLE_EBO FastTupleImpl;

template <size_t... Indices, class... Types>
struct MEMILIO_ENABLE_EBO
    FastTupleImpl<std::index_sequence<Indices...>, Types...> : FastTupleContainer<Indices, Types>... {

    constexpr FastTupleImpl() = default;

    template <typename... Args>
        requires(sizeof...(Args) == sizeof...(Types) && (std::is_same_v<Types, std::decay_t<Args>> && ...))
    explicit constexpr FastTupleImpl(Args&&... args)
        : FastTupleContainer<Indices, Types>(std::forward<Args>(args))...
    {
    }
};

} // namespace details

template <class... Types>
struct MEMILIO_ENABLE_EBO FastTuple : details::FastTupleImpl<std::make_index_sequence<sizeof...(Types)>, Types...> {
    using Base = details::FastTupleImpl<std::make_index_sequence<sizeof...(Types)>, Types...>;

    static constexpr size_t size()
    {
        return sizeof...(Types);
    }

    constexpr FastTuple() = default;

    explicit constexpr FastTuple(const Types&... types)
        requires(std::copy_constructible<Types> && ...)
        : Base(types...)
    {
    }

    explicit constexpr FastTuple(Types&&... types)
        requires(std::move_constructible<Types> && ...)
        : Base(std::move(types)...)
    {
    }
};

template <>
struct FastTuple<> {
    static constexpr size_t size()
    {
        return 0;
    }

    constexpr FastTuple() = default;

    template <size_t I>
    constexpr decltype(auto) get();
};

template <size_t I, class... Types>
constexpr decltype(auto) get(FastTuple<Types...>& t)
{
    return details::get_ftc<I>(t);
}
template <size_t I, class... Types>
constexpr decltype(auto) get(const FastTuple<Types...>& t)
{
    return details::get_ftc<I>(t);
}
template <size_t I, class... Types>
constexpr decltype(auto) get(FastTuple<Types...>&& t)
{
    return details::get_ftc<I>(t);
}

template <class Type, class... Types>
constexpr decltype(auto) get(FastTuple<Types...>& t)
{
    return details::get_ftc<index_of_type_v<Type, Types...>>(t);
}
template <class Type, class... Types>
constexpr decltype(auto) get(const FastTuple<Types...>& t)
{
    return details::get_ftc<index_of_type_v<Type, Types...>>(t);
}
template <class Type, class... Types>
constexpr decltype(auto) get(FastTuple<Types...>&& t)
{
    return details::get_ftc<index_of_type_v<Type, Types...>>(t);
}

namespace details
{

template <class T, class S, size_t... It, size_t... Is>
decltype(auto) cat_impl(T t, std::index_sequence<It...>, S s, std::index_sequence<Is...>)
{
    return FastTuple(get<It>(t)..., get<Is>(s)...);
}

template <class T, class S, size_t... It>
decltype(auto) cat_impl(T t, std::index_sequence<It...>, S, std::index_sequence<>)
{
    return t;
}

template <class T, class S, size_t... Is>
decltype(auto) cat_impl(T, std::index_sequence<>, S s, std::index_sequence<Is...>)
{
    return s;
}

template <class T, class S>
decltype(auto) cat_impl(T t, std::index_sequence<>, S, std::index_sequence<>)
{
    return t;
}

} // namespace details

template <class... T, class... S>
constexpr FastTuple<T..., S...> cat(const FastTuple<T...>& t, const FastTuple<S...>& s)
{
    return details::cat_impl(t, std::make_index_sequence<FastTuple<T...>::size()>{}, s,
                             std::make_index_sequence<FastTuple<S...>::size()>{});
}

template <class... T, class... S>
constexpr FastTuple<T..., S...> cat(FastTuple<T...>&& t, FastTuple<S...>&& s)
{
    return details::cat_impl(t, std::make_index_sequence<FastTuple<T...>::size()>{}, s,
                             std::make_index_sequence<FastTuple<S...>::size()>{});
}

namespace details
{

template <class F, class Tuple, std::size_t... I>
constexpr decltype(auto) apply_impl(F&& f, Tuple&& t, std::index_sequence<I...>) // exposition only
{
    return std::invoke(std::forward<F>(f), get<I>(std::forward<Tuple>(t))...);
}

} // namespace details

template <class F, class... Types>
constexpr decltype(auto) apply(F&& f, FastTuple<Types...>& t)
{
    return details::apply_impl(std::forward<F>(f), t, std::make_index_sequence<FastTuple<Types...>::size()>{});
}
template <class F, class... Types>
constexpr decltype(auto) apply(F&& f, const FastTuple<Types...>& t)
{
    return details::apply_impl(std::forward<F>(f), t, std::make_index_sequence<FastTuple<Types...>::size()>{});
}
template <class F, class... Types>
constexpr decltype(auto) apply(F&& f, FastTuple<Types...>&& t)
{
    return details::apply_impl(std::forward<F>(f), std::move(t),
                               std::make_index_sequence<FastTuple<Types...>::size()>{});
}

namespace details
{

template <class Tuple, size_t... I>
constexpr bool equal_impl(Tuple&& t, Tuple&& s, std::index_sequence<I...>)
{
    return ((get_ftc<I>(t) == get_ftc<I>(s)) && ...);
}

} // namespace details

template <class... Types>
constexpr bool operator==(const FastTuple<Types...>& t, const FastTuple<Types...>& s)
{
    return details::equal_impl(t, s, std::make_index_sequence<FastTuple<Types...>::size()>{});
}

/// Specialization of type_at_index for FastTuple. @see type_at_index.
template <size_t Index, class... Types>
struct type_at_index<Index, FastTuple<Types...>> : public type_at_index<Index, Types...> {
};

/// Specialization of index_of_type for FastTuple. @see index_of_type.
template <class Type, class... Types>
struct index_of_type<Type, FastTuple<Types...>> : public index_of_type<Type, Types...> {
};

/// Specialization of index_of_type for FastTuple. Resolves ambiguity when using FastTuple as items. @see index_of_type.
template <class... Types>
struct index_of_type<FastTuple<Types...>, FastTuple<Types...>> {
    static constexpr std::size_t value = 0;
};

} // namespace mio

#endif // MIO_FAST_TUPLE_H
