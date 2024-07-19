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

#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

namespace mio
{

/**
 * @brief Pair of name and value used for auto-(de)serialization.
 *
 * This object holds a view of a name and reference of a value. Mind their lifetime!
 * @tparam ValueType The (non-cv, non-reference) type of the value. 
 */
template <class ValueType>
struct NVP {
    using Type = ValueType&;
    /**
     * @brief Create a (name, value) pair.
     *
     * @param n A view of the name. 
     * @param v A non-const lvalue reference to the value. 
     */
    explicit NVP(const std::string_view n, Type v)
        : name(n)
        , value(v)
    {
    }
    const std::string_view name;
    Type value;

    NVP()                      = delete;
    NVP(const NVP&)            = delete;
    NVP(NVP&&)                 = default;
    NVP& operator=(const NVP&) = delete;
    NVP& operator=(NVP&&)      = delete;
};

/**
 * @brief Provide names and values for auto-(de)serialization.
 *
 * This function packages the class name and a name-value pair for each class member together to define both a
 * serialize and deserialize function (with limited features).
 *
 * Note that auto-serialization requires that all class members participate in serialization, and that
 * each class member is (auto-)serializable and assignable.
 * 
 * @tparam Targets List of each class member's type.
 * @param class_name The name of the class to auto-serialize.
 * @param class_members A name-value pair (NVP) for each class member.
 * @return Collection of all name views and value references used for auto-(de)serialization. 
 */
template <class... Targets>
[[nodiscard]] inline auto make_auto_serialization(const std::string_view class_name, NVP<Targets>&&... class_members)
{
    return std::make_pair(class_name, std::make_tuple(std::move(class_members)...));
}

/**
 * @brief Creates an instance of AutoSerializable for auto-deserialization.
 *
 * The default implementation uses the default constructor of AutoSerializable, if available. If there is no default
 * constructor, this class must be spezialized to provide the method `static AutoSerializable create()`. If there is
 * a default constructor, but it is private, AutoSerializableFactory<AutoSerializable> can be marked as friend instead.
 *
 * The state of the object retured by `create()` is completely arbitrary, as it is expected that auto-deserialization
 * will overwrite the value of each class member.
 *
 * @tparam AutoSerializable A type with an auto_serialize member.
 */
template <class AutoSerializable>
struct AutoSerializableFactory {
    static AutoSerializable create()
    {
        return AutoSerializable{};
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

template <class IOObject, class Target>
void add_nvp(IOObject& obj, NVP<Target> const&& nvp)
{
    obj.add_element(std::string{nvp.name}, nvp.value);
}

template <class IOContext, class... Targets>
void auto_serialize_impl(IOContext& io, const std::string_view name, std::tuple<NVP<Targets>...> const&& targets)
{
    auto obj = io.create_object(std::string{name});

    std::apply(
        [&obj](NVP<Targets> const&&... nvps) {
            (add_nvp(obj, std::move(nvps)), ...);
        },
        std::move(targets));
}

template <class IOObject, class Target>
IOResult<Target> expect_nvp(IOObject& obj, NVP<Target>&& nvp)
{
    return obj.expect_element(std::string{nvp.name}, Tag<Target>{});
}

template <class IOContext, class AutoSerializable, class... Targets>
IOResult<AutoSerializable> auto_deserialize_impl(IOContext& io, AutoSerializable& a, std::string_view name,
                                                 std::tuple<NVP<Targets>...>&& targets)
{
    auto obj = io.expect_object(std::string{name});

    auto unpacked_apply = [&io, &a, &obj](NVP<Targets>... nvps) {
        return apply(
            io,
            [&a, &nvps...](const Targets&... values) {
                ((nvps.value = values), ...);
                return a;
            },
            expect_nvp(obj, std::move(nvps))...);
    };

    return std::apply(unpacked_apply, std::move(targets));
}

} // namespace details

/**
 * @brief Detect whether T has a auto_serialize member function.
 * @tparam T Any type.
 */
template <class T>
using has_auto_serialize = is_expression_valid<details::auto_serialize_expr_t, T>;

// disables itself if a deserialize member is present or if there is no auto_serialize member
// generates serialize method depending on NVPs given by auto_serialize
template <class IOContext, class AutoSerializable,
          std::enable_if_t<has_auto_serialize<AutoSerializable>::value &&
                               not has_serialize<IOContext, AutoSerializable>::value,
                           AutoSerializable*> = nullptr>
void serialize_internal(IOContext& io, const AutoSerializable& t)
{
    // Note that this cast is only safe if we do not modify targets.
    const auto targets = const_cast<AutoSerializable*>(&t)->auto_serialize();
    details::auto_serialize_impl(io, targets.first, std::move(targets.second));
}

// disables itself if a deserialize member is present or if there is no auto_serialize member
// generates deserialize method depending on NVPs given by auto_serialize
template <class IOContext, class AutoSerializable,
          std::enable_if_t<has_auto_serialize<AutoSerializable>::value &&
                               not has_deserialize<IOContext, AutoSerializable>::value,
                           AutoSerializable*> = nullptr>
IOResult<AutoSerializable> deserialize_internal(IOContext& io, Tag<AutoSerializable>)
{
    AutoSerializable a = AutoSerializableFactory<AutoSerializable>::create();
    auto targets       = a.auto_serialize();
    return details::auto_deserialize_impl(io, a, targets.first, std::move(targets.second));
}

} // namespace mio

#endif // MIO_IO_AUTO_SERIALIZE_H_
