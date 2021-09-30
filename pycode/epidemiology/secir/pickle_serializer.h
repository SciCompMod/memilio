/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Martin Siggel, Daniel Abele, Martin J. Kuehn, Jan Kleinert
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
#ifndef PICKLE_SERIALIZER_H
#define PICKLE_SERIALIZER_H
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

#include "epidemiology/utils/io.h"
#include "epidemiology/utils/metaprogramming.h"
#include "boost/optional.hpp"
#include <fstream>
#include <utility>
#include <typeinfo>

namespace py = pybind11;

namespace epi
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
struct TupleType : std::false_type {
};
template <>
struct TupleType<bool> : std::true_type {
};


template <class T>
using is_small_integral = std::integral_constant<bool, (std::is_integral<T>::value && sizeof(T) <= 4)>;
//signed small ints
template <class T>
struct TupleType<T, std::enable_if_t<std::is_arithmetic<T>::value>> : std::true_type {
};

//unsigned small ints
template <class T>
struct TupleType<T, std::enable_if_t<conjunction_v<is_small_integral<T>, std::is_unsigned<T>>>> : std::true_type {
};

//signed big ints
template <>
struct TupleType<int64_t> : std::true_type {
};

//unsigned big ints
template <>
struct TupleType<uint64_t> : std::true_type {
};

//double
template <>
struct TupleType<double> : std::true_type {
};

//float
template <>
struct TupleType<float> : std::true_type {
};

//string
template <>
struct TupleType<std::string> : std::true_type {
};

//string literals
template <>
struct TupleType<const char*> : std::true_type {
};

template<typename T, class V>
static T transform_deserialize(const V& v)
{
    bool is_tuple = py::isinstance<py::tuple>(v);
    if (is_tuple && py::len(v) == 1)
    {
        return v[0].template cast<T>();
    }
    else
    {
        return v.template cast<T>();
    }
}
template<typename T>
static py::tuple transform_serialize(const T& v)
{
    return py::make_tuple(v);
}
/**@}*/

/**
 * Base class for implementations of serialization framework concepts.
 * Stores status and flags.
 */
class TupleBase
{
public:
    /**
     * Constructor that sets status and flags.
     */
    TupleBase(std::shared_ptr<IOStatus> status, int flags)
        : m_status(status)
        , m_flags(flags)
    {
        assert(status && "Status must not be null.");
    }

    /**
     * Flags that determine the behavior of serialization.
     * @see epi::IOFlags
     */
    int flags() const
    {
        return m_flags;
    }

    /**
     * Set flags that determine the behavior of serialization.
     * @see epi::IOFlags
     */
    void set_flags(int f)
    {
        m_flags = f;
    }

    /**
     * The current status of serialization.
     * Contains errors that occurred.
     */
    const IOStatus& status() const
    {
        return *m_status;
    }

    /**
     * Set the current status of serialization.
     */
    void set_error(const IOStatus& status)
    {
        if (*m_status) {
            *m_status = status;
        }
    }

protected:
    std::shared_ptr<IOStatus> m_status;
    int m_flags;
};

/**
 * Implementation of the IOObject concept for JSON format.
 */
class TupleObject : public TupleBase
{
public:
    /**
     * Constructor to set the status, flags, and json value to store the data.
     * @param status status, shared with the parent IO context and objects.
     * @param value reference to the json value that will store the data.
     * @param flags flags to determine the behavior of serialization.
     */
    TupleObject(const std::shared_ptr<IOStatus>& status, py::tuple& value, int flags)
        : TupleBase{status, flags}
        , m_value{value}
        , m_index(0)
    {
    }

    /**
     * add element to json value.
     * @tparam T the type of the value to be serialized.
     * @param name name of the element.
     * @param value value of the element.
     * @{
     */
    template <class T, std::enable_if_t<TupleType<T>::value, void*> = nullptr>
    void add_element( const std::string& name, const T& value);
    template <class T, std::enable_if_t<!TupleType<T>::value, void*> = nullptr>
    void add_element( const std::string& name, const T& value);
    /**@}*/

    /**
     * add optional element to json value.
     * @tparam T the type of the value to be serialized.
     * @param name name of the element.
     * @param value pointer to value of the element, may be null.
     */
    template <class T>
    void add_optional(const std::string& name,  const T* value);

    /**
     * add list of elementss to json value.
     * @tparam Iter type of the iterators that represent the list.
     * @param name name of the list.
     * @param b iterator to first element in the list.
     * @param e iterator to end of the list.
     * @{
     */
    template <class Iter, std::enable_if_t<TupleType<typename Iter::value_type>::value, void*> = nullptr>
    void add_list( const std::string& name, Iter b, Iter e);
    template <class Iter, std::enable_if_t<!TupleType<typename Iter::value_type>::value, void*> = nullptr>
    void add_list( const std::string& name, Iter b, Iter e);
    /**@}*/

    /**
     * retrieve element from the json value.
     * @tparam T the type of value to be deserialized.
     * @param name name of the element.
     * @param tag define type of the element for overload resolution.
     * @return retrieved element if succesful, error otherwise.
     * @{
     */
    template <class T, std::enable_if_t<TupleType<T>::value, void*> = nullptr>
    IOResult<T> expect_element(const std::string& name, Tag<T> tag);
    template <class T, std::enable_if_t<!TupleType<T>::value, void*> = nullptr>
    IOResult<T> expect_element(const std::string& name, Tag<T> tag);
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
    IOResult<boost::optional<T>> expect_optional(const std::string& name, epi::Tag<T> tag);

    /**
     * retrieve list of elements from the json value.
     * @tparam T the type of the elements in the list to be deserialized.
     * @param name name of the list.
     * @param tag define type of the list elements for overload resolution.
     * @param return vector of deserialized elements if succesful, error otherwise.
     * @{
     */
    template <class T, std::enable_if_t<TupleType<T>::value, void*> = nullptr>
    IOResult<std::vector<T>> expect_list(const std::string& name, Tag<T> tag);
    template <class T, std::enable_if_t<!TupleType<T>::value, void*> = nullptr>
    IOResult<std::vector<T>> expect_list(const std::string& name, Tag<T> tag);
    /**@}*/

    /**
     * The json value that data is stored in.
     */
    const auto& value() const&
    {
        return m_value;
    }

    size_t m_index;
private:
    py::tuple& m_value;
};

/**
 * Implemenetation of IOContext concept for JSON format.
 */
class TupleContext : public TupleBase
{
public:
    /**
     * Create context for serialization, set status and flags.
     * Creates an empty json value that data can be added to.
     * @param status status of serialization, shared with parent IO contexts and objects.
     * @param flags flags to determine behavior of serialization.
     */
    TupleContext(const std::shared_ptr<IOStatus>& status, int flags)
        : TupleBase{status, flags}
        , m_value{}
        , m_index(0) 
    {
    }

    /**
     * Create context for deserialization, set status, flags and json value that contains data.
     * Data for deserialization can be read from the given json value.
     * @param value json value that contains the data.
     * @param status status of serialization, shared with parent IO contexts and objects.
     * @param flags flags to determine behavior of serialization.
     */
    TupleContext(const py::tuple& value, const std::shared_ptr<IOStatus>& status, int flags)
        : TupleBase(status, flags)
        , m_value(value)
        , m_index(0)
    {
    }

    /**
     * Create a JsonObject that accepts serialization data.
     * The type of the object is currently ignored but might
     * be used for verification later.
     * @param type name of the type of the object.
     * @return new JsonObject for serialization.
     */
    TupleObject create_object(const std::string& type)
    {
        epi::unused(type);
        m_value = py::tuple();
        return {m_status, m_value, m_flags};
    }

    /**
     * Create a JsonObject that contains serialized data.
     * The type of the object is currently ignored but might
     * be used for verification later.
     * @param type name of the type of the object.
     * @return new JsonObject for deserialization.
     */
    TupleObject expect_object(const std::string& type)
    {
        epi::unused(type);
        return TupleObject(m_status, m_value, m_flags);
    }

    /**
     * The json value that contains the serialized data.
     * @{
     */
    const py::tuple& value() const&
    {
        return m_value;
    }
    /**
     * overload for rvalue references to allow moving the value out of this context.
     */
    py::tuple&& value() &&
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
    template <class T, std::enable_if_t<TupleType<T>::value, void*> = nullptr>
    friend void serialize_internal(TupleContext& io, const T& t)
    {
        io.m_value = transform_serialize(t);
    }

    /**
     * Deserialize objects of basic types on their own.
     * Used to restore e.g. a single integer from json if it is not a member of an object or container.
     * @tparam T the type of value to be deserialized.
     * @param io reference JsonContext.
     * @param t value to be serialized.
     */ 
    template <class T, std::enable_if_t<TupleType<T>::value, void*> = nullptr>
    friend IOResult<T> deserialize_internal(TupleContext& io, Tag<T>)
    {
        return epi::success(transform_deserialize<T>(io.m_value));
    }

private:
    size_t m_index;
    py::tuple m_value;
};

/**
 * Main class for (de-)serialization from/into json.
 * Root IOContext for json serialization, so always
 * starts with a status without any errors.
 */
class TupleSerializer : public TupleContext
{
public:
    /**
     * Constructor for deserialization, sets the flags and value that contains serialized data.
     * @param value value that contains serialized data.
     * @param flags flags that determine the behavior of serialization; see epi::IOFlags.
     */
    TupleSerializer(const py::tuple& value, int flags = IOF_None)
        : TupleContext(value, std::make_shared<IOStatus>(), flags)
    {
    }

    /**
     * Constructor for serialization, sets the flags.
     * Starts with an empty json value that will store the serialized data.
     * @param flags flags that determine the behavior of serialization; see epi::IOFlags.
     */
    TupleSerializer(int flags = IOF_None)
        : TupleSerializer(py::tuple{}, flags)
    {
    }
};

/**
 * Serialize an object into json.
 * @tparam T the type of value to be serialized.
 * @param t the object to be serialized.
 * @param flags flags that determine the behavior of serialized; see epi::IOFlags.
 * @return json value if serialization is succesful, error code otherwise.
 */
template <class T>
IOResult<py::tuple> serialize_tuple(const T& v, int flags = IOF_None)
{
    TupleSerializer t{flags};
    serialize(t, v);
    if (!t.status()) {
        return failure(t.status());
    }
    return success(t.value());
}


/**
 * Deserialize an object from json.
 * @tparam T the type of value to be deserialized.
 * @param js the json value.
 * @param tag defines the type of the object for overload resolution.
 * @param flags define behavior of serialization; see epi::IOFlags.
 * @return the deserialized object if succesful, error code otherwise.
 */
template <class T>
IOResult<T> deserialize_tuple(const py::tuple& t, Tag<T> tag, int flags = IOF_None)
{
    TupleSerializer ctxt{t, flags};
    return deserialize(ctxt, tag);
}


/////////////////////////////////////////////////////////////////
//Implementations for JsonContext/Object member functions below//
/////////////////////////////////////////////////////////////////

template <class T, std::enable_if_t<TupleType<T>::value, void*> = nullptr>
void TupleObject::add_element(const std::string& name, const T& value)
{
    if (m_status->is_ok()) {
        m_value += transform_serialize(value);
    }
}

template <class T, std::enable_if_t<!TupleType<T>::value, void*> = nullptr>
void TupleObject::add_element(const std::string& name, const T& value)
{
    if (m_status->is_ok()) {
        auto ctxt = TupleContext(m_status, m_flags);
        epi::serialize(ctxt, value);
        if (m_status) {
            m_value += py::make_tuple(std::move(ctxt).value());
        }
    }
}

template <class T>
void TupleObject::add_optional(const std::string& name, const T* value)
{
    if (value) {
        add_element(name, *value);
    }
    else {
        m_value += py::make_tuple(py::tuple(0));
    }
}

template <class Iter, std::enable_if_t<TupleType<typename Iter::value_type>::value, void*> = nullptr>
void TupleObject::add_list(const std::string& name, Iter b, Iter e)
{
    if (m_status->is_ok()) {
        auto new_tuple = py::tuple(0);
        for (auto it = b; it < e; ++it) {
            new_tuple += (transform_serialize(*it));
        }
        m_value += py::make_tuple(new_tuple);
    }
}

template <class Iter, std::enable_if_t<!TupleType<typename Iter::value_type>::value, void*> = nullptr>
void TupleObject::add_list(const std::string &name, Iter b, Iter e)
{
    if (m_status->is_ok()) {
        for (auto it = b; it < e; ++it) {
            auto ctxt = TupleContext(m_status, m_flags);
            epi::serialize(ctxt, *it);
            if (m_status) {
                m_value = m_value + py::make_tuple(std::move(ctxt).value());
            }
        }
    }
}

template <class T, std::enable_if_t<TupleType<T>::value, void*> = nullptr>
IOResult<T> TupleObject::expect_element(const std::string& name,Tag<T> /*tag*/)
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }
    if (py::len(m_value) <= m_index) {
        return failure(StatusCode::OutOfRange, "Index of " + name + " ot of range.");
    }

    return success(transform_deserialize<T>(m_value[m_index++]));
}

template <class T, std::enable_if_t<!TupleType<T>::value, void*> = nullptr>
IOResult<T> TupleObject::expect_element(const std::string& name,Tag<T> tag)
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }

    if (py::len(m_value) <= m_index) {
        return failure(StatusCode::OutOfRange, "Index of " + name + " ot of range.");
    }

    auto ctxt = TupleContext(m_value[m_index++], m_status, m_flags);
    return epi::deserialize(ctxt, tag);
}

template <class T>
IOResult<boost::optional<T>> TupleObject::expect_optional(const std::string &name, epi::Tag<T> tag)
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }
    if (py::len(m_value) <= m_index) {
        return failure(StatusCode::OutOfRange, "Index of " + name + " ot of range.");
    }
    const auto& element = m_value[m_index];
    bool is_tuple = py::isinstance<py::tuple>(element);
    if (is_tuple) {
        if (transform_deserialize<py::tuple>(element).is(py::tuple(0))) {
            m_index++;
            return success(boost::optional<T>{});
        }
    }
    auto r = expect_element(name, tag);

    if (r) {
        return success(r.value());
    }
    return failure(r.error());
}

template <class T, std::enable_if_t<TupleType<T>::value, void*> = nullptr>
IOResult<std::vector<T>> TupleObject::expect_list(const std::string &name, Tag<T> /*tag*/)
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }
    if (py::len(m_value) <= m_index) {
        return failure(StatusCode::OutOfRange, "Index of " + name + " ot of range.");
    }
    const auto& tuple = transform_deserialize<py::tuple>(m_value[m_index++]);
    std::vector<T> v;
    v.reserve(py::len(tuple));
    for (size_t i = 0; i < py::len(tuple); ++i) {
        v.emplace_back(transform_deserialize<T>(tuple[i]));
    }
    return success(std::move(v));
}

template <class T, std::enable_if_t<!TupleType<T>::value, void*> = nullptr>
IOResult<std::vector<T>> TupleObject::expect_list(const std::string &name, Tag<T> tag)
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }
    if (py::len(m_value) <= m_index) {
        return failure(StatusCode::OutOfRange, "Index of " + name + " ot of range.");
    }
    const auto& tuple = transform_deserialize<py::tuple>(m_value[m_index++]);
    std::vector<T> v;
    v.reserve(py::len(tuple));
    for (size_t i = 0; i < py::len(tuple); ++i) {
        auto ctxt = TupleContext(tuple[i], m_status, m_flags);
        auto r    = epi::deserialize(ctxt, tag);
        if (r) {
            v.emplace_back(std::move(r).value());
        }
        else {
            return failure(r.error());
        }
    }
    return success(std::move(v));
}

} // namespace epi

#endif
