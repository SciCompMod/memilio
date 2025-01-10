/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding
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
#ifndef MIO_IO_DEFAULT_SERIALIZE_H_
#define MIO_IO_DEFAULT_SERIALIZE_H_

#include "memilio/io/io.h"
#include "memilio/utils/metaprogramming.h"

#include <tuple>
#include <type_traits>
#include <utility>

namespace mio
{

/**
 * @brief A pair of name and reference.
 * 
 * Used for default (de)serialization.
 * This object holds a char pointer to a name and reference to value. Mind their lifetime!
 * @tparam ValueType The (non-cv, non-reference) type of the value. 
 */
template <class ValueType>
struct NamedRef {
    using Reference = ValueType&;

    const char* name;
    Reference value;

    /**
     * @brief Create a named reference.
     *
     * @param n A string literal. 
     * @param v A non-const lvalue reference to the value. 
     */
    explicit NamedRef(const char* n, Reference v)
        : name(n)
        , value(v)
    {
    }
};

namespace details
{

/**
 * @brief Helper type to detect whether T has a default_serialize member function.
 * Use has_default_serialize.
 * @tparam T Any type.
 */
template <class T>
using default_serialize_expr_t = decltype(std::declval<T>().default_serialize());

/// Add a name-value pair to an io object.
template <class IOObject, class Member>
void add_named_ref(IOObject& obj, const NamedRef<Member> named_ref)
{
    obj.add_element(named_ref.name, named_ref.value);
}

/// Unpack all name-value pairs from the tuple and add them to a new io object with the given name.
template <class IOContext, class... Members>
void default_serialize_impl(IOContext& io, const char* name, const NamedRef<Members>... named_refs)
{
    auto obj = io.create_object(name);
    (add_named_ref(obj, named_refs), ...);
}

/// Retrieve a name-value pair from an io object.
template <class IOObject, class Member>
IOResult<Member> expect_named_ref(IOObject& obj, const NamedRef<Member> named_ref)
{
    return obj.expect_element(named_ref.name, Tag<Member>{});
}

/// Read an io object and its members from the io context using the given names and assign the values to a.
template <class IOContext, class DefaultSerializable, class... Members>
IOResult<DefaultSerializable> default_deserialize_impl(IOContext& io, DefaultSerializable& a, const char* name,
                                                       NamedRef<Members>... named_refs)
{
    auto obj = io.expect_object(name);

    // we cannot use expect_named_ref directly in apply, as function arguments have no guarantueed order of evaluation
    std::tuple<IOResult<Members>...> results{expect_named_ref(obj, named_refs)...};

    return apply(
        io,
        [&a, &named_refs...](const Members&... result_values) {
            // if all results are successfully deserialized, they are unpacked into result_values
            // then all class variables are overwritten (via the named_refs) with these values
            ((named_refs.value = result_values), ...);
            return a;
        },
        results);
}

} // namespace details

/**
 * @brief List of a class's members.
 * 
 * Used for default (de)serialization.
 * Holds a char pointer to the class name as well as a tuple of NamedRefs with all added class members.
 * Initially, the template parameter pack should be left empty. It will be filled by calling Members::add.
 * @tparam ValueTypes The (non-cv, non-reference) types of member variables.
 */
template <class... ValueTypes>
struct Members {
    // allow other Members access to the private constructor
    template <class...>
    friend struct Members;

    /**
     * @brief Initialize Members with a class name. Use the member function `add` to specify the class's variables.
     * @param[in] class_name Name of a class.
     */
    Members(const char* class_name)
        : name(class_name)
        , named_refs()
    {
    }

    /**
     * @brief Add a class member.
     *
     * Use this function consecutively for all members, e.g. `Members("class").add("a", a).add("b", b).add...`.
     *
     * @param[in] member_name The name used for serialization. Should be the same as or similar to the class member.
     * For example, a good option a private class member `m_time` is simply `"time"`.
     * @param[in] member A class member. Always pass this variable directly, do not use getters or accessors.
     * @return A Members object with all previous class members and the newly added one.  
     */
    template <class T>
    [[nodiscard]] Members<ValueTypes..., T> add(const char* member_name, T& member)
    {
        return Members<ValueTypes..., T>{name, std::tuple_cat(named_refs, std::tuple(NamedRef{member_name, member}))};
    }

    const char* name; ///< Name of the class.
    std::tuple<NamedRef<ValueTypes>...> named_refs; ///<  Names and references to members of the class.

private:
    /**
     * @brief Initialize Members directly. Used by the add function.
     * @param[in] class_name Name of a class.
     * @param[in] named_references Tuple of added class Members.
     */
    Members(const char* class_name, std::tuple<NamedRef<ValueTypes>...> named_references)
        : name(class_name)
        , named_refs(named_references)
    {
    }
};

/**
 * @brief Creates an instance of T for later initialization.
 *
 * The default implementation uses the default constructor of T, if available. If there is no default constructor, this
 * class can be spezialized to provide the method `static T create()`. If there is a default constructor, but it is
 * private, DefaultFactory<T> can be marked as friend instead.
 *
 * The state of the object retured by `create()` is completely arbitrary, and may be invalid. Make sure to set it to a
 * valid state before using it further.
 *
 * @tparam T The type to create.
 */
template <class T>
struct DefaultFactory {
    /// @brief Creates a new instance of T.
    static T create()
    {
        return T{};
    }
};

/**
 * @brief Detect whether T has a default_serialize member function.
 * @tparam T Any type.
 */
template <class T>
using has_default_serialize = is_expression_valid<details::default_serialize_expr_t, T>;

/**
 * @brief Serialization implementation for the default serialization feature.
 * Disables itself (SFINAE) if there is no default_serialize member or if a serialize member is present.
 * Generates the serialize method depending on the NamedRefs given by default_serialize.
 * @tparam IOContext A type that models the IOContext concept.
 * @tparam DefaultSerializable A type that can be default serialized.
 * @param io An IO context.
 * @param a An instance of DefaultSerializable to be serialized.
 */
template <class IOContext, class DefaultSerializable,
          std::enable_if_t<has_default_serialize<DefaultSerializable>::value &&
                               !has_serialize<IOContext, DefaultSerializable>::value,
                           DefaultSerializable*> = nullptr>
void serialize_internal(IOContext& io, const DefaultSerializable& a)
{
    // Note that the following cons_cast is only safe if we do not modify members.
    const auto members = const_cast<DefaultSerializable&>(a).default_serialize();
    // unpack members and serialize
    std::apply(
        [&io, &members](auto... named_refs) {
            details::default_serialize_impl(io, members.name, named_refs...);
        },
        members.named_refs);
}

/**
 * @brief Deserialization implementation for the default serialization feature.
 * Disables itself (SFINAE) if there is no default_serialize member or if a deserialize meember is present.
 * Generates the deserialize method depending on the NamedRefs given by default_serialize.
 * @tparam IOContext A type that models the IOContext concept.
 * @tparam DefaultSerializable A type that can be default serialized.
 * @param io An IO context.
 * @param tag Defines the type of the object that is to be deserialized (i.e. DefaultSerializble).
 * @return The restored object if successful, an error otherwise.
 */
template <class IOContext, class DefaultSerializable,
          std::enable_if_t<has_default_serialize<DefaultSerializable>::value &&
                               !has_deserialize<IOContext, DefaultSerializable>::value,
                           DefaultSerializable*> = nullptr>
IOResult<DefaultSerializable> deserialize_internal(IOContext& io, Tag<DefaultSerializable> tag)
{
    mio::unused(tag);
    DefaultSerializable a = DefaultFactory<DefaultSerializable>::create();
    auto members          = a.default_serialize();
    // unpack members and deserialize
    return std::apply(
        [&io, &members, &a](auto... named_refs) {
            return details::default_deserialize_impl(io, a, members.name, named_refs...);
        },
        members.named_refs);
}

} // namespace mio

#endif // MIO_IO_DEFAULT_SERIALIZE_H_
