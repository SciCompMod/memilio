/*
* Copyright (C) 2020-2024 MEmilio
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
#ifndef MIO_IO_AUTO_SERIALIZE_H_
#define MIO_IO_AUTO_SERIALIZE_H_

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
 * Used for auto-(de)serialization.
 * This object holds a pointer to a name and reference to value. Mind their lifetime!
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

template <class... Ts>
struct Members {

    Members(const char* class_name)
        : name(class_name)
        , name_value_pairs()
    {
    }

    Members(const char* class_name, std::tuple<NamedRef<Ts>...> nvps)
        : name(class_name)
        , name_value_pairs(nvps)
    {
    }

    template <class T>
    [[nodiscard]] Members<Ts..., T> add(const char* member_name, T& member)
    {
        return Members<Ts..., T>{name, std::tuple_cat(name_value_pairs, std::tuple(NamedRef{member_name, member}))};
    }

    const char* name;
    std::tuple<NamedRef<Ts>...> name_value_pairs;
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

namespace details
{

/**
 * @brief Helper type to detect whether T has a auto_serialize member function.
 * @tparam T Any type.
 */
template <class T>
using auto_serialize_expr_t = decltype(std::declval<T>().auto_serialize());

/// Add a name-value pair to an io object.
template <class IOObject, class Member>
void add_nvp(IOObject& obj, const NamedRef<Member> nvp)
{
    obj.add_element(nvp.name, nvp.value);
}

/// Unpack all name-value pairs from the tuple and add them to a new io object with the given name.
template <class IOContext, class... Members>
void auto_serialize_impl(IOContext& io, const char* name, const NamedRef<Members>... nvps)
{
    auto obj = io.create_object(name);
    (add_nvp(obj, nvps), ...);
}

/// Retrieve a name-value pair from an io object.
template <class IOObject, class Member>
IOResult<Member> expect_nvp(IOObject& obj, const NamedRef<Member> nvp)
{
    return obj.expect_element(nvp.name, Tag<Member>{});
}

/// Read an io object and its members from the io context using the given names and assign the values to a.
template <class IOContext, class AutoSerializable, class... Members>
IOResult<AutoSerializable> auto_deserialize_impl(IOContext& io, AutoSerializable& a, const char* name,
                                                 NamedRef<Members>... nvps)
{
    auto obj = io.expect_object(name);

    return apply(
        io,
        [&a, &nvps...](const Members&... values) {
            ((nvps.value = values), ...);
            return a;
        },
        expect_nvp(obj, nvps)...);
}

} // namespace details

/**
 * @brief Detect whether T has a auto_serialize member function.
 * @tparam T Any type.
 */
template <class T>
using has_auto_serialize = is_expression_valid<details::auto_serialize_expr_t, T>;

/**
 * @brief Serialization implementation for the auto-serialization feature.
 * Disables itself (SFINAE) if there is no auto_serialize member or if a serialize member is present.
 * Generates the serialize method depending on the NVPs given by auto_serialize.
 * @tparam IOContext A type that models the IOContext concept.
 * @tparam AutoSerializable A type that can be auto-serialized.
 * @param io An IO context.
 * @param a An instance of AutoSerializable to be serialized.
 */
template <
    class IOContext, class AutoSerializable,
    std::enable_if_t<has_auto_serialize<AutoSerializable>::value && !has_serialize<IOContext, AutoSerializable>::value,
                     AutoSerializable*> = nullptr>
void serialize_internal(IOContext& io, const AutoSerializable& a)
{
    // Note that the following cons_cast is only safe if we do not modify members.
    const auto members = const_cast<AutoSerializable&>(a).auto_serialize();
    // unpack members and serialize
    std::apply(
        [&io, &members](auto... nvps) {
            details::auto_serialize_impl(io, members.name, nvps...);
        },
        members.name_value_pairs);
}

/**
 * @brief Deserialization implementation for the auto-serialization feature.
 * Disables itself (SFINAE) if there is no auto_serialize member or if a deserialize meember is present.
 * Generates the deserialize method depending on the NVPs given by auto_serialize.
 * @tparam IOContext A type that models the IOContext concept.
 * @tparam AutoSerializable A type that can be auto-serialized.
 * @param io An IO context.
 * @param tag Defines the type of the object that is to be deserialized (i.e. AutoSerializble).
 * @return The restored object if successful, an error otherwise.
 */
template <class IOContext, class AutoSerializable,
          std::enable_if_t<has_auto_serialize<AutoSerializable>::value &&
                               !has_deserialize<IOContext, AutoSerializable>::value,
                           AutoSerializable*> = nullptr>
IOResult<AutoSerializable> deserialize_internal(IOContext& io, Tag<AutoSerializable> tag)
{
    mio::unused(tag);
    AutoSerializable a = DefaultFactory<AutoSerializable>::create();
    auto members       = a.auto_serialize();
    // unpack members and deserialize
    return std::apply(
        [&io, &members, &a](auto... nvps) {
            return details::auto_deserialize_impl(io, a, members.name, nvps...);
        },
        members.name_value_pairs);
}

} // namespace mio

#endif // MIO_IO_AUTO_SERIALIZE_H_
