/* 
* Copyright (C) 2020-2025 MEmilio
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

#include "memilio/io/io.h"
#include "memilio/utils/metaprogramming.h"
#include "boost/optional.hpp"
#include <fstream>
#include <utility>
#include <typeinfo>

namespace mio
{

/**
 * PickleType allows the conversion of basic types for serialization.
 * All types that are directly handled by the Pickle serializer inherit
 * from std::true_type.
 * std::false_type, external serialization functions need to be available.
 * @tparam T the type to be serialized.
 * @{
 */
template <class T, class = void>
struct PickleType : std::false_type {
};
template <>
struct PickleType<bool> : std::true_type {
};

template <class T>
using is_small_integral = std::integral_constant<bool, (std::is_integral<T>::value && sizeof(T) <= 4)>;

//small ints
template <class T>
struct PickleType<T, std::enable_if_t<is_small_integral<T>::value>> : std::true_type {
};

//signed big ints
template <>
struct PickleType<int64_t> : std::true_type {
};

//unsigned big ints
template <>
struct PickleType<uint64_t> : std::true_type {
};

//double
template <>
struct PickleType<double> : std::true_type {
};

//float
template <>
struct PickleType<float> : std::true_type {
};

//string
template <>
struct PickleType<std::string> : std::true_type {
};

//string literals
template <>
struct PickleType<const char*> : std::true_type {
};

/**
 * Transform functions for conversion from a pybind11::tuple element.
 * @tparam T the type to be deserialized.
 * @param v Parameter that needs to be deserialized
 */
template <typename T, class V>
static T from_tuple_element(const V& v)
{
    return v.template cast<T>();
}

/**
 * Transform functions for conversion to a pybind11::tuple.
 * @tparam T the type to be serialized.
 * @param v Parameter that needs to be serialized
 */
template <typename T>
static pybind11::tuple to_tuple(const T& v)
{
    static_assert(PickleType<T>::value, "v must be one of the types for which PickleType is true");
    return pybind11::make_tuple(v);
}
/**@}*/

/**
 * Base class for implementations of serialization framework concepts.
 * Stores status and flags.
 */
class PickleBase
{
public:
    /**
     * Constructor that sets status and flags.
     */
    PickleBase(std::shared_ptr<IOStatus> status, int flags)
        : m_status(status)
        , m_flags(flags)
    {
        assert(status && "Status must not be null.");
    }

    /**
     * Flags that determine the behavior of serialization.
     * @see mio::IOFlags
     */
    int flags() const
    {
        return m_flags;
    }

    /**
     * Set flags that determine the behavior of serialization.
     * @see mio::IOFlags
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
 * Implementation of the IOObject concept for Pickle format.
 */
class PickleObject : public PickleBase
{
public:
    /**
     * Constructor to set the status, flags, and tuple to store the data.
     * @param status status, shared with the parent IO context and objects.
     * @param value reference to the tuple that will store the data.
     * @param flags flags to determine the behavior of serialization.
     */
    PickleObject(const std::shared_ptr<IOStatus>& status, pybind11::tuple& value, int flags)
        : PickleBase{status, flags}
        , m_value{value}
        , m_index(0)
    {
    }

    /**
     * add element to the tuple.
     * @tparam T the type of the value to be serialized.
     * @param name name of the element.
     * @param value value of the element.
     * @{
     */
    template <class T, std::enable_if_t<PickleType<T>::value, void*> = nullptr>
    void add_element(const std::string& name, const T& value);
    template <class T, std::enable_if_t<!PickleType<T>::value, void*> = nullptr>
    void add_element(const std::string& name, const T& value);
    /**@}*/

    /**
     * add optional element to the tuple.
     * @tparam T the type of the value to be serialized.
     * @param name name of the element.
     * @param value pointer to value of the element, may be null.
     */
    template <class T>
    void add_optional(const std::string& name, const T* value);

    /**
     * add list of elementss to the tuple.
     * @tparam Iter type of the iterators that represent the list.
     * @param name name of the list.
     * @param b iterator to first element in the list.
     * @param e iterator to end of the list.
     * @{
     */
    template <class Iter,
              std::enable_if_t<PickleType<typename std::iterator_traits<Iter>::value_type>::value, void*> = nullptr>
    void add_list(const std::string& name, Iter b, Iter e);
    template <class Iter,
              std::enable_if_t<!PickleType<typename std::iterator_traits<Iter>::value_type>::value, void*> = nullptr>
    void add_list(const std::string& name, Iter b, Iter e);
    /**@}*/

    /**
     * retrieve element from the tuple.
     * @tparam T the type of value to be deserialized.
     * @param name name of the element.
     * @param tag define type of the element for overload resolution.
     * @return retrieved element if succesful, error otherwise.
     * @{
     */
    template <class T, std::enable_if_t<PickleType<T>::value, void*> = nullptr>
    IOResult<T> expect_element(const std::string& name, Tag<T> tag);
    template <class T, std::enable_if_t<!PickleType<T>::value, void*> = nullptr>
    IOResult<T> expect_element(const std::string& name, Tag<T> tag);
    /**@}*/

    /**
     * retrieve optional element from the tuple.
     * @tparam T the type of value to be deserialized.
     * @param name name of the element.
     * @param tag define type of the element for overload resolution.
     * @return retrieved element if name is found and can be deserialized, empty optional if not found, error otherwise.
     */
    template <class T>
    IOResult<boost::optional<T>> expect_optional(const std::string& name, mio::Tag<T> tag);

    /**
     * retrieve list of elements from the tuple.
     * @tparam T the type of the elements in the list to be deserialized.
     * @param name name of the list.
     * @param tag define type of the list elements for overload resolution.
     * @param return vector of deserialized elements if succesful, error otherwise.
     * @{
     */
    template <class T, std::enable_if_t<PickleType<T>::value, void*> = nullptr>
    IOResult<std::vector<T>> expect_list(const std::string& name, Tag<T> tag);
    template <class T, std::enable_if_t<!PickleType<T>::value, void*> = nullptr>
    IOResult<std::vector<T>> expect_list(const std::string& name, Tag<T> tag);
    /**@}*/

    /**
     * The tuple that data is stored in.
     */
    const auto& value() const&
    {
        return m_value;
    }

    size_t m_index;

private:
    pybind11::tuple& m_value;
};

/**
 * Implemenetation of IOContext concept for Pickle format.
 */
class PickleContext : public PickleBase
{
public:
    /**
     * Create context for serialization, set status and flags.
     * Creates an empty Tuple that data can be added to.
     * @param status status of serialization, shared with parent IO contexts and objects.
     * @param flags flags to determine behavior of serialization.
     */
    PickleContext(const std::shared_ptr<IOStatus>& status, int flags)
        : PickleBase{status, flags}
        , m_value{}
        , m_index(0)
    {
    }

    /**
     * Create context for deserialization, set status, flags and the tuple that contains data.
     * Data for deserialization can be read from the given tuple.
     * @param value tuple that contains the data.
     * @param status status of serialization, shared with parent IO contexts and objects.
     * @param flags flags to determine behavior of serialization.
     */
    PickleContext(const pybind11::tuple& value, const std::shared_ptr<IOStatus>& status, int flags)
        : PickleBase(status, flags)
        , m_value(value)
        , m_index(0)
    {
    }

    /**
     * Create a PickleObject that accepts serialization data.
     * The type of the object is currently ignored but might
     * be used for verification later.
     * @param type name of the type of the object.
     * @return new PickleObject for serialization.
     */
    PickleObject create_object(const std::string& type)
    {
        mio::unused(type);
        m_value = pybind11::tuple();
        return {m_status, m_value, m_flags};
    }

    /**
     * Create a PickleObject that contains serialized data.
     * The type of the object is currently ignored but might
     * be used for verification later.
     * @param type name of the type of the object.
     * @return new PickleObject for deserialization.
     */
    PickleObject expect_object(const std::string& type)
    {
        mio::unused(type);
        return PickleObject(m_status, m_value, m_flags);
    }

    /**
     * The tuple that contains the serialized data.
     * @{
     */
    const pybind11::tuple& value() const&
    {
        return m_value;
    }
    /**
     * overload for rvalue references to allow moving the value out of this context.
     */
    pybind11::tuple&& value() &&
    {
        return std::move(m_value);
    }

    /**@}*/

    /**
     * Serialize objects of basic types on their own.
     * Used to turn e.g. a single integer into tuple if it is not a member of an object or container.
     * @tparam T the type of value to be serialized.
     * @param io reference PickleContext.
     * @param t value to be serialized.
     */
    template <class T, std::enable_if_t<PickleType<T>::value, void*> = nullptr>
    friend void serialize_internal(PickleContext& io, const T& t)
    {
        io.m_value = to_tuple(t);
    }

    /**
     * Deserialize objects of basic types on their own.
     * Used to restore e.g. a single integer from a tuple if it is not a member of an object or container.
     * @tparam T the type of value to be deserialized.
     * @param io reference TupleContext.
     * @param t value to be serialized.
     */
    template <class T, std::enable_if_t<PickleType<T>::value, void*> = nullptr>
    friend IOResult<T> deserialize_internal(PickleContext& io, Tag<T>)
    {
        return mio::success(from_tuple_element<T>(io.m_value[0]));
    }

private:
    size_t m_index;
    pybind11::tuple m_value;
};

/**
 * Main class for (de-)serialization from/into pickle.
 * Root IOContext for pickle serialization, so always
 * starts with a status without any errors.
 */
class PickleSerializer : public PickleContext
{
public:
    /**
     * Constructor for deserialization, sets the flags and value that contains serialized data.
     * @param value value that contains serialized data.
     * @param flags flags that determine the behavior of serialization; see mio::IOFlags.
     */
    PickleSerializer(const pybind11::tuple& value, int flags = IOF_None)
        : PickleContext(value, std::make_shared<IOStatus>(), flags)
    {
    }

    /**
     * Constructor for serialization, sets the flags.
     * Starts with an empty tuple that will store the serialized data.
     * @param flags flags that determine the behavior of serialization; see mio::IOFlags.
     */
    PickleSerializer(int flags = IOF_None)
        : PickleSerializer(pybind11::tuple{}, flags)
    {
    }
};

/**
 * Serialize an object into a tuple.
 * @tparam T the type of value to be serialized.
 * @param t the object to be serialized.
 * @param flags flags that determine the behavior of serialized; see mio::IOFlags.
 * @return pybind11::tuple if serialization is succesful, error code otherwise.
 */
template <class T>
IOResult<pybind11::tuple> serialize_pickle(const T& v, int flags = IOF_None)
{
    PickleSerializer t{flags};
    serialize(t, v);
    if (!t.status()) {
        return failure(t.status());
    }
    return success(t.value());
}

/**
 * Deserialize an object from a tuple.
 * @tparam T the type of value to be deserialized.
 * @param t the tuple.
 * @param tag defines the type of the object for overload resolution.
 * @param flags define behavior of serialization; see mio::IOFlags.
 * @return the deserialized object if succesful, error code otherwise.
 */
template <class T>
IOResult<T> deserialize_pickle(const pybind11::tuple& t, Tag<T> tag, int flags = IOF_None)
{
    PickleSerializer ctxt{t, flags};
    return deserialize(ctxt, tag);
}

////////////////////////////////////////////////////////////////////
//Implementations for PickleContext/Object member functions below //
///////////////////////////////////////////////////////////////////

template <class T, std::enable_if_t<PickleType<T>::value, void*>>
void PickleObject::add_element(const std::string& name, const T& value)
{
    if (m_status->is_ok()) {
        m_value = m_value + to_tuple(value);
    }
}

template <class T, std::enable_if_t<!PickleType<T>::value, void*>>
void PickleObject::add_element(const std::string& name, const T& value)
{
    if (m_status->is_ok()) {
        auto ctxt = PickleContext(m_status, m_flags);
        mio::serialize(ctxt, value);
        if (m_status) {
            m_value = m_value + pybind11::make_tuple(std::move(ctxt).value());
        }
    }
}

template <class T>
void PickleObject::add_optional(const std::string& name, const T* value)
{
    auto end = value;
    if (value) {
        end += 1;
    }
    add_list(name, value, end);
}

template <class Iter, std::enable_if_t<PickleType<typename std::iterator_traits<Iter>::value_type>::value, void*>>
void PickleObject::add_list(const std::string& name, Iter b, Iter e)
{

    if (m_status->is_ok()) {
        auto new_tuple = pybind11::tuple(0);
        for (auto it = b; it < e; ++it) {
            new_tuple = new_tuple + (to_tuple(*it));
        }
        m_value = m_value + pybind11::make_tuple(new_tuple);
    }
}

template <class Iter, std::enable_if_t<!PickleType<typename std::iterator_traits<Iter>::value_type>::value, void*>>
void PickleObject::add_list(const std::string& name, Iter b, Iter e)
{
    if (m_status->is_ok()) {
        auto new_tuple = pybind11::tuple(0);
        for (auto it = b; it < e; ++it) {
            auto ctxt = PickleContext(m_status, m_flags);
            mio::serialize(ctxt, *it);
            if (m_status) {
                new_tuple = new_tuple + pybind11::make_tuple(std::move(ctxt).value());
            }
        }
        m_value = m_value + pybind11::make_tuple(new_tuple);
    }
}

template <class T, std::enable_if_t<PickleType<T>::value, void*>>
IOResult<T> PickleObject::expect_element(const std::string& name, Tag<T> /*tag*/)
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }
    if (pybind11::len(m_value) <= m_index) {
        return failure(StatusCode::OutOfRange, "Index of " + name + " out of range.");
    }

    return success(from_tuple_element<T>(m_value[m_index++]));
}

template <class T, std::enable_if_t<!PickleType<T>::value, void*>>
IOResult<T> PickleObject::expect_element(const std::string& name, Tag<T> tag)
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }

    if (pybind11::len(m_value) <= m_index) {
        return failure(StatusCode::OutOfRange, "Index of " + name + " out of range.");
    }

    auto ctxt = PickleContext(m_value[m_index++].template cast<pybind11::tuple>(), m_status, m_flags);
    return mio::deserialize(ctxt, tag);
}

template <class T>
IOResult<boost::optional<T>> PickleObject::expect_optional(const std::string& name, mio::Tag<T> tag)
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }

    BOOST_OUTCOME_TRY(auto&& r, expect_list(name, tag));

    if (r.size() == 0) {
        return success();
    }
    else if (r.size() == 1) {
        return success(r[0]);
    }

    return failure(StatusCode::OutOfRange, "Optional must be a tuple with only one element");
}

template <class T, std::enable_if_t<PickleType<T>::value, void*>>
IOResult<std::vector<T>> PickleObject::expect_list(const std::string& name, Tag<T> /*tag*/)
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }
    if (pybind11::len(m_value) <= m_index) {
        return failure(StatusCode::OutOfRange, "Index of " + name + " out of range.");
    }
    const auto& tuple = m_value[m_index++].template cast<pybind11::tuple>();
    std::vector<T> v;
    v.reserve(pybind11::len(tuple));
    for (size_t i = 0; i < pybind11::len(tuple); ++i) {
        v.emplace_back(from_tuple_element<T>(tuple[i]));
    }
    return success(std::move(v));
}

template <class T, std::enable_if_t<!PickleType<T>::value, void*>>
IOResult<std::vector<T>> PickleObject::expect_list(const std::string& name, Tag<T> tag)
{
    if (m_status->is_error()) {
        return failure(*m_status);
    }
    if (pybind11::len(m_value) <= m_index) {
        return failure(StatusCode::OutOfRange, "Index of " + name + " out of range.");
    }
    const auto& tuple = m_value[m_index++].template cast<pybind11::tuple>();
    std::vector<T> v;
    v.reserve(pybind11::len(tuple));
    for (size_t i = 0; i < pybind11::len(tuple); ++i) {
        auto ctxt = PickleContext(tuple[i], m_status, m_flags);
        auto r    = mio::deserialize(ctxt, tag);
        if (r) {
            v.emplace_back(std::move(r).value());
        }
        else {
            return failure(r.error());
        }
    }
    return success(std::move(v));
}

} // namespace mio

#endif
