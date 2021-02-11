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
    template <class T, class X = void>
    struct has_get_default_member_function {
        static const bool value = false;
    };
    template <class T>
    struct has_get_default_member_function<T, details::void_t<decltype(T::get_default())>> {
        static const bool value = true;
    };
} // namespace details

//check whether a get_default function exists
template <class T>
using has_get_default_member_function = details::has_get_default_member_function<T>;

//the properties of a parameter
//can be specialized for further customization
template <class Tag>
struct ParameterTagTraits {
    using Type = typename Tag::Type;

    //get_default either with get_default member function or fallback on default constructor
    template <class Dummy = Tag>
    static std::enable_if_t<has_get_default_member_function<Dummy>::value, decltype(Dummy::get_default())> get_default()
    {
        return Tag::get_default();
    }
    template <class Dummy = Tag>
    static std::enable_if_t<
        !has_get_default_member_function<Dummy>::value && std::is_default_constructible<Type>::value, Type>
    get_default()
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
} // namespace details

struct DefaultInit {
};

//a set of parameters defined at compile time
//parameters added as template parameters (tags)
//example: struct FooParamTag { using type = X; ... }; ParameterSet<FooParamTag>
template <class... Tags>
class ParameterSet
{
public:
    //default constructor
    //exists if all parameters are default constructible
    template <
        class Dummy = void,
        class = std::enable_if_t<details::AllOf<std::is_default_constructible, typename Tags::Type...>::value, Dummy>>
    ParameterSet()
    {
    }

    //default initializing constructor
    //exists if all parameters have get_default
    template <class Dummy = void,
              class       = std::enable_if_t<
                  details::AllOf<has_get_default_member_function, ParameterTagTraits<Tags>...>::value, Dummy>>
    explicit ParameterSet(DefaultInit)
        : m_tup(ParameterTagTraits<Tags>::get_default()...)
    {
    }

    //explicit initializing constructor
    //initializes the n-th parameter using the n-th argument
    template <class... T, class = std::enable_if_t<
                              (sizeof...(T) >= 1 &&
                               std::is_constructible<std::tuple<details::TaggedParameter<Tags>...>, T...>::value),
                              void>>
    explicit ParameterSet(T&&... args)
        : m_tup(std::forward<T>(args)...)
    {
    }

    //get value of parameter
    template <class Tag>
    const typename ParameterTagTraits<Tag>::Type& get() const
    {
        return std::get<details::TaggedParameter<Tag>>(m_tup).get();
    }

    template <class Tag>
    typename ParameterTagTraits<Tag>::Type& get()
    {
        return std::get<details::TaggedParameter<Tag>>(m_tup).get();
    }

    //set value of parameter
    template <class Tag>
    void set(const typename ParameterTagTraits<Tag>::Type& value)
    {
        get<Tag>() = value;
    }

    template <class Tag, class T>
    void set(T&& arg)
    {
        get<Tag>() = std::forward<T>(arg);
    }

    //(re)set parameter to its default value
    //only exists if parameter defines get_default
    template <class Tag, class = std::enable_if_t<has_get_default_member_function<ParameterTagTraits<Tag>>::value, Tag>>
    void set_default()
    {
        get<Tag>() = ParameterTagTraits<Tag>::get_default();
    }

    //number of parameters
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

//get the the tag of the I-th parameter in a set
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

//call f(t) for all parameters with
//t a default constructed parameter tag
template <class Params, class F>
void foreach_tag(F f)
{
    details::foreach_tag_impl<Params>(f, std::make_index_sequence<Params::size()>{});
}

//call f(p, t) for all parameters with
//p the value of the parameter
//t a default constructed parameter tag
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