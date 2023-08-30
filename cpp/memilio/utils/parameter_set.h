/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Jan Kleinert, Daniel Abele
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
#ifndef EPI_UTILS_PARAMETER_SET_H
#define EPI_UTILS_PARAMETER_SET_H

#include "memilio/io/io.h"
#include "memilio/utils/stl_util.h"
#include <tuple>
#include <utility>

namespace mio
{

namespace details
{

//helpers for get_default
template <class X, class = void, class... Args>
struct has_get_default_member_function : std::false_type {
};

template <class T, class... Args>
struct has_get_default_member_function<T, void_t<decltype(T::get_default(std::declval<Args>()...))>, Args...>
    : std::true_type {
};

} // namespace details

/**
 * @brief check whether a get_default function exists
 * @tparam The type to check for the existence of the member function
 */
template <class T, class... Args>
using has_get_default_member_function = details::has_get_default_member_function<T, void, Args...>;

/**
 * @brief the properties of a parameter
 * @tparam Tag the parameter type
 */
template <class Tag>
struct ParameterTagTraits {
    using Type = typename Tag::Type;

    //get_default either with get_default member function or fallback on default constructor
    template <class Dummy = Tag, class... Ts>
    static std::enable_if_t<has_get_default_member_function<Dummy, Ts...>::value, Type> get_default(Ts&&... args)
    {
        return Tag::get_default(std::forward<Ts>(args)...);
    }

    template <class Dummy = Tag, class... Ts>
    static std::enable_if_t<
        !has_get_default_member_function<Dummy, Ts...>::value && std::is_default_constructible<Type>::value, Type>
    get_default(Ts&&...)
    {
        return Type{};
    }
};

namespace details
{
//stores a parameter and tags it so it can be unambiguously found in a tuple
template <class TagT>
class TaggedParameter
{
public:
    using Tag    = TagT;
    using Traits = ParameterTagTraits<Tag>;
    using Type   = typename Traits::Type;

    template <class... Ts, class Dummy1 = void,
              class = std::enable_if_t<std::is_constructible<Type, Ts...>::value, Dummy1>>
    TaggedParameter(Ts&&... args)
        : m_value(std::forward<Ts>(args)...)
    {
    }

    operator Type&()
    {
        return get();
    }

    operator const Type&() const
    {
        return get();
    }

    const Type& get() const
    {
        return m_value;
    }

    Type& get()
    {
        return m_value;
    }

    template <class T>
    bool operator==(const TaggedParameter<T>& other) const
    {
        return m_value == other.m_value;
    }

    template <class T>
    bool operator!=(const TaggedParameter<T>& other) const
    {
        return m_value != other.m_value;
    }

private:
    Type m_value;
};

//defines value = true if predicate defines value=true for all types in parameter pack
template <template <class...> class Pred, class... Tail>
struct AllOf;

template <template <class...> class Pred>
struct AllOf<Pred> : public std::true_type {
};

template <template <class...> class Pred, class Head, class... Tail>
struct AllOf<Pred, Head, Tail...> {
    static const constexpr bool value = Pred<Head>::value && AllOf<Pred, Tail...>::value;
};

//defines value = true if predicate defines value=true for any type in parameter pack
template <template <class...> class Pred, class... Tail>
struct AnyOf;

template <template <class...> class Pred>
struct AnyOf<Pred> : public std::false_type {
};

template <template <class...> class Pred, class Head, class... Tail>
struct AnyOf<Pred, Head, Tail...> {
    static const constexpr bool value = Pred<Head>::value || AnyOf<Pred, Tail...>::value;
};

//for X = template<T1, T2> X => BindTail<X, A>::type<B> = X<B, A>
template <template <class...> class F, class... Tail>
struct BindTail {
    template <class... Head>
    struct type : F<Head..., Tail...> {
    };
    //according to the standard, this must be a real type, can't be an alias of F.
    //An alias is immediately replaced and discarded when the compiler sees it, but this could leave the template F with some
    //parameters bound and some free, which is not allowed. A struct that derives from F is persistent during compilation.
};

//for template<T1, T2> X => BindHead<X, A>::type<B> = X<A, B>
template <template <class...> class F, class... Head>
struct BindHead {
    template <class... Tail>
    struct type : F<Head..., Tail...> {
    };
    //can't be an alias, see BindTail
};

} // namespace details

/**
 * @brief A tag used for tag-dispatching the Constructor of ParameterSet,
 * triggering default initialization of all parameters using the get_default
 * member function.
 */
struct NoDefaultInit {
};

template <typename T>
using is_no_default_init_tag = std::is_same<NoDefaultInit, T>;

/**
 * @brief a set of parameters defined at compile time
 *
 * parameters added as template parameters (tags)
 *
 * example:
 *
 *     struct FooParamTag { using type = X; ... };
 *     ParameterSet<FooParamTag>
 *
 * @tparam Tags All parameter types contained in this set. The types should be unique.
 */
template <class... Tags>
class ParameterSet
{
public:
    /**
     * @brief Non-initializing default constructor.
     * This constructor exists if all parameters are default constructible
     * and it can be used by calling the constructor with an empty object NoDefaultInit.
     * It serves in cases where the get_default() functions of the parameters are very
     * costly and should not be called as parameters will set not non-default values
     * anyway.
     */
    template <
        class Dummy = void,
        class = std::enable_if_t<details::AllOf<std::is_default_constructible, typename Tags::Type...>::value, Dummy>>
    explicit ParameterSet(NoDefaultInit)
    {
    }

    /**
     * @brief default initializing constructor
     * Initializes each parameter using either the get_default function defined in the parameter tag or the default constructor.
     * this constructor exists if all parameters have get_default() without arguments or a default constructor.
     */
    template <class Dummy = void,
              class       = std::enable_if_t<
                  details::AllOf<has_get_default_member_function, ParameterTagTraits<Tags>...>::value, Dummy>>
    ParameterSet()
        : m_tup(ParameterTagTraits<Tags>::get_default()...)
    {
    }

    /**
     * @brief default initializing constructor.
     * Initializes each parameter using either the get_default function defined in the parameter tag or the default constructor.
     * this constructor exists if all parameters have get_default(args...) with the same number of arguments or a default constructor.
     * Arguments get forwarded to get_default of parameters.
     @tparam T1 First argument of get_default(...) function called on the parameters, e.g., T1=AgeGroup&.
     @tparam TN Further arguments passed to get_default(...) functions. 
     */
    template <class T1, class... TN,
              class = std::enable_if_t<
              // Avoid erroneous template deduction for T1=ParameterSet as this constructor could falsely be considered
              // as a copy constructor for non-const lvalue references.
                  conjunction_v<negation<std::is_same<std::decay_t<T1>, ParameterSet>>, 
                        details::AllOf<details::BindTail<has_get_default_member_function, T1, TN...>::template type,
                                 ParameterTagTraits<Tags>...>>>>
    explicit ParameterSet(T1&& arg1, TN&&... argn)
        : m_tup(ParameterTagTraits<Tags>::get_default(arg1, argn...)...)
    {
    }

    /**
     * @brief get value of a parameter
     * @tparam Tag the queried parameter
     * @return The value of the parameter
     */
    template <class Tag>
    const typename ParameterTagTraits<Tag>::Type& get() const
    {
        return std::get<details::TaggedParameter<Tag>>(m_tup).get();
    }

    /**
     * @brief get value of a parameter
     * @tparam Tag the queried parameter
     * @return The value of the parameter
     */
    template <class Tag>
    typename ParameterTagTraits<Tag>::Type& get()
    {
        return std::get<details::TaggedParameter<Tag>>(m_tup).get();
    }

    /**
     * @brief set value of a parameter
     * @tparam Tag the parameter
     */
    template <class Tag>
    void set(const typename ParameterTagTraits<Tag>::Type& value)
    {
        get<Tag>() = value;
    }

    /**
     * @brief set value of a parameter
     * @tparam Tag the parameter
     */
    template <class Tag, class T>
    void set(T&& arg)
    {
        get<Tag>() = std::forward<T>(arg);
    }

    /**
     * @brief (re)set parameter to its default value
     *
     * only exists if parameter defines get_default
     *
     * @tparam Tag the parameter
     */
    template <class Tag, class... T>
    std::enable_if_t<has_get_default_member_function<ParameterTagTraits<Tag>, T...>::value, void> set_default(T&&... ts)
    {
        get<Tag>() = ParameterTagTraits<Tag>::get_default(std::forward<T>(ts)...);
    }

    /**
     * @brief returns the number of parameters
     * @return //number of parameters
     */
    static constexpr size_t size()
    {
        return sizeof...(Tags);
    }

    bool operator==(const ParameterSet& b) const
    {
        return m_tup == b.m_tup;
    }

    bool operator!=(const ParameterSet& b) const
    {
        return m_tup != b.m_tup;
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("ParameterSet");
        foreach (*this, [&obj](auto& p, auto t) mutable {
            using Tag = decltype(t);
            obj.add_element(Tag::name(), p);
        });
    }

private:
    ParameterSet(const typename Tags::Type&... t)
        : m_tup(t...)
    {
    }

    //entry to recursively deserialize all parameters in the ParameterSet
    //IOContext: serializer
    //IOObject: object that stores the serialized ParameterSet
    //Rs: IOResult<T> for each Parameter Tag that has already been deserialized
    template<class IOContext, class IOObject, class... Rs, std::enable_if_t<(sizeof...(Rs) < sizeof...(Tags)), void*> = nullptr>
    static IOResult<ParameterSet> deserialize_recursive(IOContext& io, IOObject& obj, Rs&&... rs)
    {
        //read current parameter, append result to results of previous parameters, recurse to next parameter
        const size_t I = sizeof...(Rs);
        using TaggedParameter = std::tuple_element_t<I, decltype(ParameterSet::m_tup)>;
        auto r = obj.expect_element(TaggedParameter::Tag::name(), mio::Tag<typename TaggedParameter::Type>{});
        return deserialize_recursive(io, obj, std::forward<Rs>(rs)..., std::move(r));
    }

    //end of recursion to deserialize parameters in the ParameterSet
    template<class IOContext, class IOObject, class... Rs, std::enable_if_t<(sizeof...(Rs) == sizeof...(Tags)), void*> = nullptr>
    static IOResult<ParameterSet> deserialize_recursive(IOContext& io, IOObject& /*obj*/, Rs&&... rs)
    {
        //one result for each parameters, so no more parameters to read
        //expand results, build finished ParameterSet, stop recursion
        return mio::apply(io, [](const typename Tags::Type&... t) {
            return ParameterSet(t...);
        }, std::forward<Rs>(rs)...);
    }

public:
    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<ParameterSet> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("ParameterSet");
        return deserialize_recursive(io, obj);
    }

private:
    std::tuple<details::TaggedParameter<Tags>...> m_tup;
};

namespace details
{
//helpers for ParameterTag
template <size_t I, class ParamSet>
struct ParameterTag;

template <size_t I, class Head, class... Tail>
struct ParameterTag<I, ParameterSet<Head, Tail...>> : public ParameterTag<I - 1, ParameterSet<Tail...>> {
};

template <class Head, class... Tail>
struct ParameterTag<0, ParameterSet<Head, Tail...>> {
    using Type = Head;
};
} // namespace details

/**
 * @brief get the the tag of the I-th parameter in a set
 */
template <size_t I, class ParamSet>
using ParameterTag = details::ParameterTag<I, ParamSet>;

template <size_t I, class ParamSet>
using ParameterTagT = typename ParameterTag<I, ParamSet>::Type;

namespace details
{
//helpers for foreach
template <class... Tail, class Params, class F>
std::enable_if_t<sizeof...(Tail) == 0, void> foreach_impl(Params&&, F)
{
}

template <class Head, class... Tail, class Params, class F>
void foreach_impl(Params&& p, F f)
{
    f(p.template get<Head>(), Head{});
    foreach_impl<Tail...>(p, f);
}

template <class Params, size_t... Tail, class F>
std::enable_if_t<sizeof...(Tail) == 0, void> foreach_tag_impl(F, std::index_sequence<Tail...>)
{
}

template <class Params, size_t Head, size_t... Tail, class F>
void foreach_tag_impl(F f, std::index_sequence<Head, Tail...>)
{
    f(ParameterTagT<Head, Params>{});
    foreach_tag_impl<Params, Tail...>(f, std::index_sequence<Tail...>{});
}
} // namespace details

/**
 * @brief call f(t) for all parameters in a ParameterSet with
 * t a default constructed parameter tag
 *
 * @tparam Params a ParameterSet
 * @tparam F The function type of f
 * @param f The function to be called
 */
template <class Params, class F>
void foreach_tag(F f)
{
    details::foreach_tag_impl<Params>(f, std::make_index_sequence<Params::size()>{});
}

/**
 * @brief call f(p, t) for all parameters in a ParameterSet with
 * p the value of the parameter
 * t a default constructed parameter tag
 *
 * @tparam F The function type of f
 * @tparam Tags the parameters
 * @param p the ParameterSet
 * @param f The function to be called
 */
template <class F, class... Tags>
void foreach (const ParameterSet<Tags...>& p, F f)
{
    details::foreach_impl<Tags...>(p, f);
}

template <class F, class... Tags>
void foreach (ParameterSet<Tags...>& p, F f)
{
    details::foreach_impl<Tags...>(p, f);
}

} // namespace mio

#endif //EPI_UTILS_PARAMETER_SET_H
