#ifndef EPI_UTILS_PARAMETER_SET_H
#define EPI_UTILS_PARAMETER_SET_H

#include <utility>
#include <tuple>

namespace epi
{

namespace details
{
    //see std::void_t (c++ 17)
    template <typename... Ts>
    struct make_void {
        typedef void type;
    };
    template <typename... Ts>
    using void_t = typename make_void<Ts...>::type;

    //helpers for get_default
    template <class X, class = void, class... Args>
    struct has_get_default_member_function
        : std::false_type
    {};

    template <class T, class... Args>
    struct has_get_default_member_function<T, void_t<decltype(T::get_default(std::declval<Args>()...))>, Args...>
        : std::true_type
    {};

    //helpers for check_constraint
    template <class T, class X = void>
    struct has_check_constraints_member_function
        : std::false_type
    {};

    template <class T>
    struct has_check_constraints_member_function<T ,details::void_t<decltype(T::check_constraints(std::declval<typename T::Type const&>()))>>
        : std::true_type
    {};

    //helpers for apply_constraints
    template <class T, class X = void>
    struct has_apply_constraints_member_function
        : std::false_type
    {};

    template <class T>
    struct has_apply_constraints_member_function<T, details::void_t<decltype(std::declval<T>().apply_constraints())>>
        : std::true_type
    {};
} // namespace details

/**
 * @brief check whether a get_default function exists
 * @tparam The type to check for the existence of the member function
 */
template <class T, class... Args>
using has_get_default_member_function = details::has_get_default_member_function<T, void, Args...>;

/**
 * @brief check whether a check_constraints function exists
 * @tparam The type to check for the existence of the member function
 */
template <class T>
using has_check_constraints_member_function = details::has_check_constraints_member_function<T>;

/**
 * @brief check whether a apply_constraints function exists
 * @tparam The type to check for the existence of the member function
 */
template <class T>
using has_apply_constraints_member_function = details::has_apply_constraints_member_function<T>;

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

        //check_constraints either with check_constraints member function or fallback on (void)()
        template <class Dummy = Tag>
        std::enable_if_t<has_check_constraints_member_function<Dummy>::value, void> check_constraints() const
        {
            Tag::check_constraints(m_value);
        }
        template <class Dummy = Tag>
        std::enable_if_t<
            !has_check_constraints_member_function<Dummy>::value, void>
        check_constraints() const
        {
        }

        template<class T>
        bool operator==(const TaggedParameter<T>& other)
        {
            return m_value == other.m_value;
        }

        template<class T>
        bool operator!=(const TaggedParameter<T>& other)
        {
            return m_value != other.m_value;
        }

    private:
        Type m_value;
    };

    //defines value = true if predicate defines value=true for all types in parameter pack
    template <template <class> class Pred, class... Tail>
    struct AllOf;

    template <template <class> class Pred>
    struct AllOf<Pred> : public std::true_type {
    };

    template <template <class> class Pred, class Head, class... Tail>
    struct AllOf<Pred, Head, Tail...> {
        static const constexpr bool value = Pred<Head>::value && AllOf<Pred, Tail...>::value;
    };

    //defines value = true if predicate defines value=true for any type in parameter pack
    template <template <class> class Pred, class... Tail>
    struct AnyOf;

    template <template <class> class Pred>
    struct AnyOf<Pred> : public std::false_type {
    };

    template <template <class> class Pred, class Head, class... Tail>
    struct AnyOf<Pred, Head, Tail...> {
        static const constexpr bool value = Pred<Head>::value || AnyOf<Pred, Tail...>::value;
    };

    // call std::get<i>(tup).check_constraints() for all i.
    //
    // In C++17 this is just a fold expression
    template<std::size_t I = 0, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type
    check_constraints_for_each_parameter(std::tuple<Tp...> const&)
    {
    }

    template<std::size_t I = 0, typename... Tp>
    inline typename std::enable_if<I < sizeof...(Tp), void>::type
    check_constraints_for_each_parameter(std::tuple<Tp...> const& t)
    {
        std::get<I>(t).check_constraints();
        check_constraints_for_each_parameter<I + 1, Tp...>(t);
    }

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
     * @brief default constructor
     *
     * exists if all parameters are default constructible
     */
    template <
        class Dummy = void,
        class = std::enable_if_t<details::AllOf<std::is_default_constructible, typename Tags::Type...>::value, Dummy>>
    explicit ParameterSet(NoDefaultInit)
    {
    }

    /**
     * @brief default initializing constructor
     * exists if all parameters have get_default.
     * Arguments get forwarded to get_default of parameters
     */
    template <class... T,
              class = std::enable_if_t<(sizeof...(T)==0) || !details::AnyOf<is_no_default_init_tag, T...>::value, void>>
    ParameterSet(T&&... args)
        : m_tup(ParameterTagTraits<Tags>::get_default(std::forward<T>(args)...)...)
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
    template <class Tag,
              class... T>
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

    void check_constraints() const
    {
        details::check_constraints_for_each_parameter(m_tup);
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

} // namespace epi

#endif //EPI_UTILS_PARAMETER_SET_H
