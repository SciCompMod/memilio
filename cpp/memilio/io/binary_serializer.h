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
#ifndef MIO_IO_BINARY_SERIALIZER_H
#define MIO_IO_BINARY_SERIALIZER_H

#include "memilio/io/io.h"
#include "memilio/io/serializer_base.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/metaprogramming.h"
#include <memory>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>

namespace mio
{

/**
* In-memory stream of bytes.
*/
class ByteStream
{
public:
    /**
    * Create stream with n readable bytes initialized to 0.
    * The value of the readable bytes can be set through the pointer returned by
    * the data() member function.
    * @param n number of readable bytes.
    */
    ByteStream(size_t n)
        : m_buf(n, 0)
        , m_read_head(0)
    {
    }

    /**
    * Create empty stream.
    */
    ByteStream()
        : ByteStream(0)
    {
    }

    /**
    * Write bytes to the stream.
    * @param p pointer to bytes to be written.
    * @param s number of bytes to write.
    */
    void write(const unsigned char* const p, size_t s)
    {
        m_buf.insert(m_buf.end(), p, p + s);
    }

    /**
    * Read bytes from the stream.
    * @param p pointer to buffer to read bytes into.
    * @param s number of bytes to read.
    * @returns true if bytes can be read, false if there are not enough readable bytes in the stream.
    */
    bool read(unsigned char* p, size_t s)
    {
        if (s <= (m_buf.size() - m_read_head)) {
            auto read_begin = m_buf.begin() + m_read_head;
            auto read_end   = read_begin + s;
            std::copy(read_begin, read_end, p);
            m_read_head += s;
            return true;
        }
        return false;
    }

    /**
    * Reset the stream to n readable bytes of value 0.
    * @param n number of readable bytes.
    */
    void reset(size_t n)
    {
        m_buf.assign(n, 0);
        m_read_head = 0;
    }

    /**
    * Reset the stream to empty.
    */
    void reset()
    {
        reset(0);
    }

    /**
    * Get the pointer to the buffer of the stream.
    * @{
    */
    const unsigned char* data() const
    {
        return m_buf.data();
    }
    unsigned char* data()
    {
        return m_buf.data();
    }
    /**@}*/

    /**
    * Get the size of the buffer of the stream.
    */
    size_t data_size() const
    {
        return m_buf.size();
    }

private:
    std::vector<unsigned char> m_buf; ///< store of written/readable bytes
    size_t m_read_head; ///< index in the buffer where next byte is read/written
};

/**
* Stores a binary serialized object.
* See io.h for documentation of serialization.
*/
class BinarySerializerObject : public SerializerBase
{
public:
    BinarySerializerObject(ByteStream& stream, std::shared_ptr<IOStatus> status, int flags)
        : SerializerBase(status, flags)
        , m_stream(stream)
    {
    }

    /**
    * Add element of basic type to this object.
    * As an optimization, also handles trivial structs, Enums, etc (everything that
    * can be memcpy'd.)
    * @param name Name of the element.
    * @param value Value to be serialized.
    */
    template <class T, std::enable_if_t<std::is_trivial<T>::value, void*> = nullptr>
    void add_element(const std::string& name, const T& value);

    /**
    * Add string element to this object.
    * @param name Name of the element.
    * @param value Value to be serialized.
    */
    void add_element(const std::string& name, const std::string& value);

    /**
    * Add serializable value as an element to this object.
    * @param name Name of the element.
    * @param value Value to be serialized.
    */
    template <class T, std::enable_if_t<negation_v<std::is_trivial<T>>, void*> = nullptr>
    void add_element(const std::string& name, const T& value);

    /**
    * Get element of basic type from this object.
    * As an optimization, also handles trivial structs, Enums, etc (everything that
    * can be memcpy'd.)
    * @param name Name of the element.
    * @param tag Tag that determines the type of the element.
    * @returns Deserialized value.
    */
    template <class T, std::enable_if_t<std::is_trivial<T>::value, void*> = nullptr>
    IOResult<T> expect_element(const std::string& name, Tag<T> tag);

    /**
    * Get deserializable element from this object.
    * @param name Name of the element.
    * @param tag Tag that determines the type of the element.
    * @returns Deserialized value.
    */
    template <class T, std::enable_if_t<negation<std::is_trivial<T>>::value, void*> = nullptr>
    IOResult<T> expect_element(const std::string& name, Tag<T> tag);

    /**
    * Get serialized string element from this object.
    * @param name Name of the element.
    * @param tag Tag that determines the type of the element.
    * @returns Deserialized value.
    */
    IOResult<std::string> expect_element(const std::string& name, Tag<std::string> tag);

    /**
    * Add optional element to this object.
    * @param name Name of the element.
    * @param value Value to be serialized. Nullptr if empty.
    */
    template <class T>
    void add_optional(const std::string& name, const T* value);

    /**
    * Get optional element from this object.
    * @param name Name of the element.
    * @returns Deserialized optional value.
    */
    template <class T>
    IOResult<boost::optional<T>> expect_optional(const std::string& name, Tag<T>);

    /**
    * Add a list of elements to this object.
    * @param name Name of the list.
    * @param b Iterator to beginning of the list.
    * @param e Iterator to end of the list.
    */
    template <class Iter>
    void add_list(const std::string& name, Iter b, Iter e);

    /**
    * Get a list of elements from this object.
    * @param name Name of the list.
    * @returns List of deserialized elements.
    */
    template <class T>
    IOResult<std::vector<T>> expect_list(const std::string& name, Tag<T>);

private:
    ByteStream& m_stream; ///< Reference to a stream that stores the serialized bytes.
};

/**
* Serializes objects in binary format.
* See io.h for documentation of serialization.
*/
class BinarySerializerContext : public SerializerBase
{
public:
    BinarySerializerContext(ByteStream& stream, std::shared_ptr<IOStatus> status, int flags)
        : SerializerBase(status, flags)
        , m_stream(stream)
    {
    }

    /**
    * Begin serialization of a new object.
    * @param type Name of the type of the object. Only used if flag IOF_IncludeTypeInfo is set.
    */
    BinarySerializerObject create_object(const std::string& type)
    {
        auto obj = BinarySerializerObject(m_stream, m_status, m_flags);
        if (m_flags & IOF_IncludeTypeInfo)
        {
            obj.add_element("Type", type);
        }
        return obj;
    }

    /**
    * Begin deserialization of an object.
    * @param type Name of the type of the object. Only used if flag IOF_IncludeTypeInfo is set.
    */
    BinarySerializerObject expect_object(const std::string& type)
    {
        auto obj = BinarySerializerObject(m_stream, m_status, m_flags);
        if (m_flags & IOF_IncludeTypeInfo)
        {
            auto type_result = obj.expect_element("Type", Tag<std::string>{});
            if (!type_result) {
                *m_status = IOStatus(StatusCode::InvalidType, "No TypeInfo in stream.");
            } else if(type_result.value() != type) {
                *m_status = IOStatus(StatusCode::InvalidType, "Unexpected type in stream:" + type_result.value() + ". Expected " + type);
            }
        }
        return BinarySerializerObject(m_stream, m_status, m_flags);
    }

    /**
    * Serialize "naked" ints, etc. that are not contained in a complex object hierarchy,
    * because they don't have an object that they can be added to.
    * @param ctxt The serialization context.
    * @param t The value to be serialized.
    */
    template <class T, std::enable_if_t<std::is_trivial<T>::value, void*> = nullptr>
    friend void serialize_internal(BinarySerializerContext& ctxt, const T& t)
    {
        //add element to dummy object.
        //objects don't store anything in the stream, so this has no overhead.
        BinarySerializerObject obj(ctxt.m_stream, ctxt.m_status, ctxt.m_flags);
        obj.add_element("", t);
    }

    /**
    * Deserialize "naked" ints, etc. that are not contained in a complex object hierarchy.
    * @param ctxt The serialization context.
    * @param tag The value to be serialized.
    * @returns The deserialized value.
    */
    template <class T, std::enable_if_t<std::is_trivial<T>::value, void*> = nullptr>
    friend IOResult<T> deserialize_internal(BinarySerializerContext& ctxt, Tag<T> tag)
    {
        unused(tag);
        //read element from on dummy object.
        BinarySerializerObject obj(ctxt.m_stream, ctxt.m_status, ctxt.m_flags);
        return obj.expect_element("", Tag<T>{});
    }

private:
    ByteStream& m_stream; ///< Reference to a stream that stores the serialized bytes.
};

template <class T, std::enable_if_t<std::is_trivial<T>::value, void*>>
void BinarySerializerObject::add_element(const std::string& name, const T& value)
{
    mio::unused(name);
    auto p = reinterpret_cast<const unsigned char*>(std::addressof(value));
    m_stream.write(p, sizeof(value));
}

template <class T, std::enable_if_t<negation_v<std::is_trivial<T>>, void*>>
void BinarySerializerObject::add_element(const std::string& name, const T& value)
{
    mio::unused(name);
    auto ctxt = BinarySerializerContext(m_stream, m_status, m_flags);
    mio::serialize(ctxt, value);
}

inline void BinarySerializerObject::add_element(const std::string& name, const std::string& value)
{
    mio::unused(name);
    const auto size   = value.size();
    const auto p_size = reinterpret_cast<const unsigned char*>(std::addressof(size));
    const auto p_data = reinterpret_cast<const unsigned char*>(value.data());
    m_stream.write(p_size, sizeof(size));
    m_stream.write(p_data, size);
}

template <class T, std::enable_if_t<std::is_trivial<T>::value, void*>>
IOResult<T> BinarySerializerObject::expect_element(const std::string& name, Tag<T> /*tag*/)
{
    mio::unused(name);

    if (m_status->is_ok()) {
        T t;
        if (m_stream.read(reinterpret_cast<unsigned char*>(std::addressof(t)), sizeof(t))) {
            return mio::success(t);
        }
        else {
            *m_status =
                IOStatus(mio::StatusCode::UnknownError, "Unexpected EOF reading " + name + " from binary stream.");
        }
    }
    return failure(*m_status);
}

template <class T, std::enable_if_t<negation<std::is_trivial<T>>::value, void*>>
IOResult<T> BinarySerializerObject::expect_element(const std::string& name, Tag<T> tag)
{
    mio::unused(name);
    if (m_status->is_ok()) {
        auto ctxt = BinarySerializerContext(m_stream, m_status, m_flags);
        return mio::deserialize(ctxt, tag);
    }
    return failure(*m_status);
}

inline IOResult<std::string> BinarySerializerObject::expect_element(const std::string& name, Tag<std::string>)
{
    mio::unused(name);
    if (m_status->is_ok()) {
        size_t size;
        if (m_stream.read(reinterpret_cast<unsigned char*>(&size), sizeof(size))) {
            std::string t(size, 0);
            if (m_stream.read(reinterpret_cast<unsigned char*>(&t[0]), size)) {
                return success(t);
            }
        }
        *m_status = IOStatus(StatusCode::UnknownError, "Unexpected EOF reading " + name + " from binary stream.");
    }
    return failure(*m_status);
}

template <class Iter>
void BinarySerializerObject::add_list(const std::string& name, Iter b, Iter e)
{
    mio::unused(name);
    add_element("Size", size_t(e - b));
    for (; b != e; ++b) {
        add_element("Item", *b);
    }
}

template <class T>
IOResult<std::vector<T>> BinarySerializerObject::expect_list(const std::string& name, Tag<T>)
{
    mio::unused(name);
    BOOST_OUTCOME_TRY(size, expect_element("Size", Tag<size_t>{}));
    std::vector<T> v;
    v.reserve(size);
    for (auto i = size_t(0); i < size; ++i) {
        BOOST_OUTCOME_TRY(t, expect_element("Item", Tag<T>{}));
        v.emplace_back(std::move(t));
    }
    return success(v);
}

template <class T>
void BinarySerializerObject::add_optional(const std::string& name, const T* value)
{
    mio::unused(name);
    add_element("Exists", size_t(value ? 1 : 0));
    if (value) {
        add_element("Value", *value);
    }
}

template <class T>
IOResult<boost::optional<T>> BinarySerializerObject::expect_optional(const std::string& name, Tag<T>)
{
    mio::unused(name);
    BOOST_OUTCOME_TRY(size, expect_element("Exists", Tag<size_t>{}));
    if (size) {
        BOOST_OUTCOME_TRY(t, expect_element("Value", Tag<T>{}));
        return mio::success(t);
    }
    return mio::success(boost::optional<T>{});
}

/**
* Serialize an object into binary format.
* The format is not portable, object must be deserialized by the same (or identically compiled) program.
* @param t object to serialize.
* @return a ByteStream that contains the serialized bytes that represent the object.
* @tparam T the type of the object to be serialized.
*/
template <class T>
ByteStream serialize_binary(const T& t, int flags = 0)
{
    ByteStream stream;
    BinarySerializerContext ctxt(stream, std::make_shared<IOStatus>(), flags);
    mio::serialize(ctxt, t);
    return stream;
}


/**
* Deserialize an object from binary format.
* @param stream ByteStream that contains the bytes that represent the object.
* @param tag Tag that carries the type of the object to be deserialized.
* @return The deserialized object if succesful, an error otherwise.
* @tparam T the type of the object to be serialized.
*/
template <class T>
IOResult<T> deserialize_binary(ByteStream& stream, Tag<T>, int flags = 0)
{
    BinarySerializerContext ctxt(stream, std::make_shared<IOStatus>(), flags);
    return mio::deserialize(ctxt, Tag<T>{});
}

} // namespace mio

#endif
