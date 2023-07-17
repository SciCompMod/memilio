/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#ifndef EPI_IO_JSON_SERIALIZER_H
#define EPI_IO_JSON_SERIALIZER_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "memilio/io/io.h"
#include "memilio/io/serializer_base.h"
#include "memilio/utils/metaprogramming.h"
#include "json/json.h"
#include <fstream>
#include <utility>
#include <limits>

namespace mio
{

/**
 * JsonType allows the conversion of basic types for serialization.
 * All types that are directly handled by the json serializer inherit 
 * from std::true_type. They also define the transform functions for 
 * conversion from and to Json::Value. All others inherit from 
 * std::false_type, external serialization functions need to be available.
 * @tparam T the type to be serialized.
 * @{
 */
template <class T, class = void>
struct JsonType : std::false_type {
};
//bool
template <>
struct JsonType<bool> : std::true_type {
    static IOResult<bool> transform(const Json::Value& js)
    {
        if (js.isBool()) {
            return success(js.asBool());
        }
        return failure(StatusCode::InvalidType, "Json value is not a bool.");
    }
    static Json::Value transform(bool b)
    {
        return Json::Value(b);
    }
};
//all integers less than 32 bit must be stored as int
template <class T>
using is_small_integral = std::integral_constant<bool, (std::is_integral<T>::value && sizeof(T) <= 4)>;
//signed small ints
template <class T>
struct JsonType<T, std::enable_if_t<conjunction_v<is_small_integral<T>, std::is_signed<T>>>> : std::true_type {
    static IOResult<T> transform(const Json::Value& js)
    {
        if (js.isInt()) {
            int i = js.asInt();
            if (i >= int(std::numeric_limits<T>::min()) && i <= int(std::numeric_limits<T>::max())) {
                return success(T(i));
            }
            return failure(StatusCode::OutOfRange, "Json value is not in range of the requested integer type.");
        }
        return failure(StatusCode::InvalidType, "Json value is not an integer.");
    }
    static Json::Value transform(T i)
    {
        return Json::Value(Json::Int(i));
    }
};
//unsigned small ints
template <class T>
struct JsonType<T, std::enable_if_t<conjunction_v<is_small_integral<T>, std::is_unsigned<T>>>> : std::true_type {
    static IOResult<T> transform(const Json::Value& js)
    {
        if (js.isUInt()) {
            unsigned int i = js.asUInt();
            if (i <= (unsigned int)std::numeric_limits<T>::max()) {
                return success(T(i));
            }
            return failure(StatusCode::OutOfRange,
                           "Json value is not in range of the requested unsigned integer type.");
        }
        return failure(StatusCode::InvalidType, "Json value is not an unsigned integer.");
    }
    static Json::Value transform(T i)
    {
        return Json::Value(Json::UInt(i));
    }
};
template <class T>
using is_64bit_integral = std::integral_constant<bool, (std::is_integral<T>::value && sizeof(T) == 8)>;
//signed big ints
template <class T>
struct JsonType<T, std::enable_if_t<conjunction_v<is_64bit_integral<T>, std::is_signed<T>>>> : std::true_type {
    static IOResult<T> transform(const Json::Value& js)
    {
        if (js.isInt64()) {
            return success(js.asInt64());
        }
        return failure(StatusCode::InvalidType, "Json value is not a 64 bit integer.");
    }
    static Json::Value transform(T i)
    {
        return Json::Value(Json::Int64(i));
    }
};
//unsigned big ints
template <class T>
struct JsonType<T, std::enable_if_t<conjunction_v<is_64bit_integral<T>, std::is_unsigned<T>>>> : std::true_type {
    static IOResult<T> transform(const Json::Value& js)
    {
        if (js.isUInt64()) {
            return success(js.asUInt64());
        }
        return failure(StatusCode::InvalidType, "Json value is not an unsigned 64bit integer.");
    }
    static Json::Value transform(T i)
    {
        return Json::Value(Json::UInt64(i));
    }
};
//double
template <>
struct JsonType<double> : std::true_type {

    static IOResult<double> transform(const Json::Value& js)
    {
        if (js.isDouble()) {
            return success(js.asDouble());
        }
        return failure(StatusCode::InvalidType, "Json value is not a double.");
    }
    static Json::Value transform(double d)
    {
        return Json::Value(d);
    }
};
//float
template <>
struct JsonType<float> : std::true_type {
    static IOResult<float> transform(const Json::Value& js)
    {
        if (js.isDouble()) {
            auto d = js.asDouble();
            if (d > double(std::numeric_limits<float>::max()) || d < -double(std::numeric_limits<float>::max())) {
                return failure(StatusCode::OutOfRange, "Json value is not in single precision floating point range.");
            }
            return success(float(d));
        }
        return failure(StatusCode::InvalidType, "Json value is not a float.");
    }
    static Json::Value transform(float f)
    {
        return Json::Value(double(f));
    }
};
//string
template <>
struct JsonType<std::string> : std::true_type {
    static IOResult<std::string> transform(const Json::Value& js)
    {
        if (js.isString()) {
            return success(js.asString());
        }
        return failure(StatusCode::InvalidType, "Json value is not a string.");
    }
    static Json::Value transform(const std::string& s)
    {
        return Json::Value(s);
    }
};
//string literals
template <>
struct JsonType<const char*> : std::true_type {
    //cannot be read, but may be written (e.g. string literal), so only one way transform
    static Json::Value transform(const char* s)
    {
        return Json::Value(s);
    }
};
/**@}*/

/**
 * Implementation of the IOObject concept for JSON format.
 */
class JsonObject : public SerializerBase
{
public:
    /**
     * Constructor to set the status, flags, and json value to store the data.
     * @param status status, shared with the parent IO context and objects.
     * @param value reference to the json value that will store the data.
     * @param flags flags to determine the behavior of serialization.
     */
    JsonObject(Json::Value& value, const std::shared_ptr<IOStatus>& status, int flags)
        : SerializerBase{status, flags}
        , m_value{value}
    {
    }

    /**
     * add element to json value.
     * @tparam T the type of the value to be serialized.
     * @param name name of the element.
     * @param value value of the element.
     * @{
     */
    template <class T, std::enable_if_t<JsonType<T>::value, void*> = nullptr>
    void add_element(const std::string& name, const T& value);
    template <class T, std::enable_if_t<!JsonType<T>::value, void*> = nullptr>
    void add_element(const std::string& name, const T& value);
    /**@}*/

    /**
     * add optional element to json value.
     * @tparam T the type of the value to be serialized.
     * @param name name of the element.
     * @param value pointer to value of the element, may be null.
     */
    template <class T>
    void add_optional(const std::string& name, const T* value);

    /**
     * add list of elementss to json value.
     * @tparam Iter type of the iterators that represent the list.
     * @param name name of the list.
     * @param b iterator to first element in the list.
     * @param e iterator to end of the list.
     */
    template <class Iter>
    void add_list(const std::string& name, Iter b, Iter e);

    /**
     * retrieve element from the json value.
     * @tparam T the type of value to be deserialized.
     * @param name name of the element.
     * @param tag define type of the element for overload resolution.
     * @return retrieved element if succesful, error otherwise.
     * @{
     */
    template <class T, std::enable_if_t<JsonType<T>::value, void*> = nullptr>
    IOResult<T> expect_element(const std::string& name, Tag<T> tag) const;
    template <class T, std::enable_if_t<!JsonType<T>::value, void*> = nullptr>
    IOResult<T> expect_element(const std::string& name, Tag<T> tag) const;
    /**@}*/

    /**
     * retrieve optional element from the json value.
     * Json does not distinguish between empty optionals and missing elements.
     * @tparam T the type of value to be deserialized.
     * @param name name of the element.
     * @param tag define type of the element for overload resolution.
     * @return retrieved element if name is found and can be deserialized, empty optional if not found, error otherwise.
     */
    template <class T>
    IOResult<boost::optional<T>> expect_optional(const std::string& name, mio::Tag<T> tag);

    /**
     * retrieve list of elements from the json value.
     * @tparam T the type of the elements in the list to be deserialized.
     * @param name name of the list.
     * @param tag define type of the list elements for overload resolution.
     * @param return vector of deserialized elements if succesful, error otherwise.
     */
    template <class T>
    IOResult<std::vector<T>> expect_list(const std::string& name, Tag<T> tag);

    /**
     * The json value that data is stored in.
     */
    const auto& value() const&
    {
        return m_value;
    }

private:
    Json::Value& m_value;
};

/**
 * Implemenetation of IOContext concept for JSON format.
 */
class JsonContext : public SerializerBase
{
public:
    /**
     * Create context for serialization, set status and flags.
     * Creates an empty json value that data can be added to.
     * @param status status of serialization, shared with parent IO contexts and objects.
     * @param flags flags to determine behavior of serialization.
     */
    JsonContext(const std::shared_ptr<IOStatus>& status, int flags)
        : SerializerBase{status, flags}
        , m_value{}
    {
    }

    /**
     * Create context for deserialization, set status, flags and json value that contains data.
     * Data for deserialization can be read from the given json value.
     * @param value json value that contains the data.
     * @param status status of serialization, shared with parent IO contexts and objects.
     * @param flags flags to determine behavior of serialization.
     */
    JsonContext(const Json::Value& value, const std::shared_ptr<IOStatus>& status, int flags)
        : SerializerBase(status, flags)
        , m_value(value)
    {
    }

    /**
     * Create a JsonObject that accepts serialization data.
     * The type of the object is currently ignored but might
     * be used for verification later.
     * @param type name of the type of the object.
     * @return new JsonObject for serialization.
     */
    JsonObject create_object(const std::string& type)
    {
        mio::unused(type);
        m_value = Json::Value(Json::objectValue);
        return {m_value, m_status, m_flags};
    }

    /**
     * Create a JsonObject that contains serialized data.
     * The type of the object is currently ignored but might
     * be used for verification later.
     * @param type name of the type of the object.
     * @return new JsonObject for deserialization.
     */
    JsonObject expect_object(const std::string& type)
    {
        mio::unused(type);
        if (!m_value.isObject()) {
            *m_status = IOStatus(StatusCode::InvalidType, "Json value must be an object.");
        }
        return JsonObject(m_value, m_status, m_flags);
    }

    /**
     * The json value that contains the serialized data.
     * @{
     */
    const Json::Value& value() const&
    {
        return m_value;
    }
    /**
     * overload for rvalue references to allow moving the value out of this context.
     */
    Json::Value&& value() &&
    {
        return std::move(m_value);
    }
    /**@}*/

    /**
     * Serialize objects of basic types on their own.
     * Used to turn e.g. a single integer into json if it is not a member of an object or container.
     * @tparam T the type of value to be serialized.
     * @param io reference JsonContext.
     * @param t value to be serialized.
     */
    template <class T, std::enable_if_t<JsonType<T>::value, void*> = nullptr>
    friend void serialize_internal(JsonContext& io, const T& t)
    {
        io.m_value = JsonType<T>::transform(t);
    }

    /**
     * Deserialize objects of basic types on their own.
     * Used to restore e.g. a single integer from json if it is not a member of an object or container.
     * @tparam T the type of value to be deserialized.
     * @param io reference JsonContext.
     * @param t value to be serialized.
     */
    template <class T, std::enable_if_t<JsonType<T>::value, void*> = nullptr>
    friend IOResult<T> deserialize_internal(JsonContext& io, Tag<T>)
    {
        return JsonType<T>::transform(io.m_value);
    }

    /**
     * json specialization of serialization for containers.
     * Serialize containers as pure json array without a wrapping object.
     */
    template <class Container, std::enable_if_t<conjunction_v<is_container<Container>, negation<JsonType<Container>>,
                                                              negation<has_serialize<JsonContext, Container>>>,
                                                void*> = nullptr>
    friend void serialize_internal(JsonContext& io, const Container& v)
    {
        if (io.m_status->is_ok()) {
            auto array = Json::Value(Json::arrayValue);
            using std::begin;
            using std::end;
            for (auto it = begin(v); it != end(v); ++it) {
                auto ctxt = JsonContext(io.m_status, io.m_flags);
                mio::serialize(ctxt, *it);
                if (io.m_status) {
                    array.append(std::move(ctxt).value());
                }
            }
            io.m_value = std::move(array);
        }
    }

    /**
     * json specialization of deserialization for containers.
     * Deserialize containers from pure json arrays without a wrapping object.
     */
    template <class Container, std::enable_if_t<conjunction_v<is_container<Container>, negation<JsonType<Container>>,
                                                              negation<has_deserialize<JsonContext, Container>>>,
                                                void*> = nullptr>
    friend IOResult<Container> deserialize_internal(JsonContext& io, Tag<Container>)
    {
        const auto& array = io.m_value;
        if (array.isArray()) {
            Container v;
            for (auto&& el : array) {
                auto ctxt = JsonContext(el, io.m_status, io.m_flags);
                BOOST_OUTCOME_TRY(val, mio::deserialize(ctxt, Tag<typename Container::value_type>{}));
                v.insert(v.end(), val);
            }
            return success(std::move(v));
        }
        return failure(StatusCode::InvalidType, "Json value must be an array.");
    }

private:
    Json::Value m_value;
};

/**
 * Write the json value into a file.
 * @param path path of the file.
 * @param js_value json value to be written.
 * @return nothing if succesful, error code otherwise.
 */
inline IOResult<void> write_json(const std::string& path, const Json::Value& js_value)
{
    std::ofstream ofs(path);
    if (ofs.is_open()) {
        Json::StreamWriterBuilder swb;
        auto js_writer = std::unique_ptr<Json::StreamWriter>(swb.newStreamWriter());
        js_writer->write(js_value, &ofs);
        if (ofs) {
            return success();
        }
        else {
            return failure(StatusCode::UnknownError, "Unknown error writing json.");
        }
    }
    else {
        return failure(StatusCode::FileNotFound, path);
    }
}

/**
 * Serialize an object into json.
 * @tparam T the type of value to be serialized.
 * @param t the object to be serialized.
 * @param flags flags that determine the behavior of serialized; see mio::IOFlags.
 * @return json value if serialization is succesful, error code otherwise.
 */
template <class T>
IOResult<Json::Value> serialize_json(const T& t, int flags = IOF_None)
{
    JsonContext ctxt{std::make_shared<IOStatus>(), flags};
    serialize(ctxt, t);
    if (!ctxt.status()) {
        return failure(ctxt.status());
    }
    return success(ctxt.value());
}

/**
 * Serialize an object into json and write it into a file.
 * @tparam T the type of value to be serialized.
 * @param path the path of the file.
 * @param t the object to be serialized.
 * @param flags flags that determine the behavior of the serialization; see mio::IOFlags.
 * @return nothing if succesful, error code otherwise.
 */
template <class T>
IOResult<void> write_json(const std::string& path, const T& t, int flags = IOF_None)
{
    BOOST_OUTCOME_TRY(js, serialize_json(t, flags));
    return write_json(path, js);
}

/**
 * Read a json value from a file.
 * @param path path of the file.
 * @return a json value if succesful, error code otherwise.
 */
inline IOResult<Json::Value> read_json(const std::string& path)
{
    std::ifstream ifs(path);
    if (ifs.is_open()) {
        Json::CharReaderBuilder crb;
        std::string err_msg;
        Json::Value js_value;
        if (Json::parseFromStream(crb, ifs, &js_value, &err_msg)) {
            return success(js_value);
        }
        return failure(StatusCode::UnknownError, path + ", " + err_msg);
    }
    return failure(StatusCode::FileNotFound, path);
}

/**
 * Deserialize an object from json.
 * @tparam T the type of value to be deserialized.
 * @param js the json value.
 * @param tag defines the type of the object for overload resolution.
 * @param flags define behavior of serialization; see mio::IOFlags.
 * @return the deserialized object if succesful, error code otherwise.
 */
template <class T>
IOResult<T> deserialize_json(const Json::Value& js, Tag<T> tag, int flags = IOF_None)
{
    JsonContext ctxt{js, std::make_shared<IOStatus>(), flags};
    return deserialize(static_cast<JsonContext&>(ctxt), tag);
}

/**
 * Read a json value from a file and deserialize it into an object.
 * @tparam T the type of value to be deserialized.
 * @param path the path of the file.
 * @param tag defines the type of the object for overload resolution.
 * @param flags define behavior of serialization; see mio::IOFlags.
 * @return the deserialized object if succesful, error code otherwise.
 */
template <class T>
IOResult<T> read_json(const std::string& path, Tag<T> tag, int flags = IOF_None)
{
    BOOST_OUTCOME_TRY(js, read_json(path));
    return deserialize_json(js, tag, flags);
}

/////////////////////////////////////////////////////////////////
//Implementations for JsonContext/Object member functions below//
/////////////////////////////////////////////////////////////////

template <class T, std::enable_if_t<JsonType<T>::value, void*>>
void JsonObject::add_element(const std::string& name, const T& value)
{
    if (m_status->is_ok()) {
        m_value[name] = JsonType<T>::transform(value);
    }
}

template <class T, std::enable_if_t<!JsonType<T>::value, void*>>
void JsonObject::add_element(const std::string& name, const T& value)
{
    if (m_status->is_ok()) {
        auto ctxt = JsonContext(m_status, m_flags);
        mio::serialize(ctxt, value);
        if (m_status) {
            m_value[name] = std::move(ctxt).value();
        }
    }
}

template <class T>
void JsonObject::add_optional(const std::string& name, const T* value)
{
    if (value) {
        add_element(name, *value);
    }
}

template <class Iter>
void JsonObject::add_list(const std::string& name, Iter b, Iter e)
{
    if (m_status->is_ok()) {
        auto array = Json::Value(Json::arrayValue);
        for (auto it = b; it != e; ++it) {
            auto ctxt = JsonContext(m_status, m_flags);
            mio::serialize(ctxt, *it);
            if (m_status) {
                array.append(std::move(ctxt).value());
            }
        }
        m_value[name] = std::move(array);
    }
}

template <class T, std::enable_if_t<JsonType<T>::value, void*>>
IOResult<T> JsonObject::expect_element(const std::string& name, Tag<T> /*tag*/) const
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }

    const auto& element = m_value[name];
    if (element.isNull()) {
        return failure(StatusCode::KeyNotFound, name);
    }
    auto r = JsonType<T>::transform(element);
    if (r) {
        return r;
    }

    return failure(r.error().code(),
                   r.error().message() + " (" + name + ")"); //annotate type error message with element name
}

template <class T, std::enable_if_t<!JsonType<T>::value, void*>>
IOResult<T> JsonObject::expect_element(const std::string& name, Tag<T> tag) const
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }
    const auto& element = m_value[name];
    if (!element.isNull()) {
        auto ctxt = JsonContext(element, m_status, m_flags);
        auto r    = mio::deserialize(ctxt, tag);
        return r;
    }
    return failure(IOStatus{StatusCode::KeyNotFound, name});
}

template <class T>
IOResult<boost::optional<T>> JsonObject::expect_optional(const std::string& name, mio::Tag<T> tag)
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }
    const auto& element = m_value[name];
    if (element.isNull()) {
        return success(boost::optional<T>{});
    }
    auto r = expect_element(name, tag);
    if (r) {
        return success(r.value());
    }
    return failure(r.error());
}

template <class T>
IOResult<std::vector<T>> JsonObject::expect_list(const std::string& name, Tag<T> tag)
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }
    const auto& array = m_value[name];
    if (array.isArray()) {
        std::vector<T> v;
        v.reserve(array.size());
        for (auto&& el : array) {
            auto ctxt = JsonContext(el, m_status, m_flags);
            auto r    = mio::deserialize(ctxt, tag);
            if (r) {
                v.emplace_back(r.value());
            }
            else {
                return failure(r.error());
            }
        }
        return success(std::move(v));
    }
    return failure(StatusCode::KeyNotFound, name);
}

} // namespace mio

#endif //MEMILIO_HAS_JSONCPP

#endif //EPI_IO_JSON_SERIALIZER_H
