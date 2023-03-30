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
    std::vector<unsigned char> m_buf;
    size_t m_read_head;

public:
    /**
    * Create stream with n readable bytes initialized to 0.
    * Set value of readable bytes by accessing buffer through member function data().
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
    * Write bytes to streams.
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
};

/**
* Stores a binary serialized object.
* See io.h for documentation of serialization.
*/
class BinarySerializerObject : public SerializerBase
{
public:
    BinarySerializerObject(std::shared_ptr<IOStatus> status, int flags, ByteStream& stream)
        : SerializerBase(status, flags)
        , m_stream(stream)
    {
    }

    template <class T, std::enable_if_t<std::is_trivial<T>::value, void*> = nullptr>
    void add_element(const std::string& name, const T& value);

    void add_element(const std::string& name, const std::string& value);

    template <class T, std::enable_if_t<negation_v<std::is_trivial<T>>, void*> = nullptr>
    void add_element(const std::string& name, const T& value);

    template <class T, std::enable_if_t<std::is_trivial<T>::value, void*> = nullptr>
    IOResult<T> expect_element(const std::string& name, Tag<T>);

    template <class T, std::enable_if_t<negation<std::is_trivial<T>>::value, void*> = nullptr>
    IOResult<T> expect_element(const std::string& name, Tag<T>);

    IOResult<std::string> expect_element(const std::string& name, Tag<std::string>);

    template <class T>
    void add_optional(const std::string& name, const T* value);

    template <class T>
    IOResult<boost::optional<T>> expect_optional(const std::string& name, Tag<T>);

    template <class Iter>
    void add_list(const std::string& name, Iter b, Iter e);

    template <class T>
    IOResult<std::vector<T>> expect_list(const std::string& name, Tag<T>);

private:
    ByteStream& m_stream;
};

/**
* Serializes objects in binary format.
* See io.h for documentation of serialization.
*/
class BinarySerializerContext : public SerializerBase
{
public:
    BinarySerializerContext(std::shared_ptr<IOStatus> status, int flags, ByteStream& stream)
        : SerializerBase(status, flags)
        , m_stream(stream)
    {
    }

    BinarySerializerObject create_object(const std::string& /*type*/)
    {
        return BinarySerializerObject(m_status, m_flags, m_stream);
    }

    BinarySerializerObject expect_object(const std::string& /*type*/)
    {
        //no way to do check if there really is an object here
        //maybe give option to add type info for debugging?
        return BinarySerializerObject(m_status, m_flags, m_stream);
    }

    template <class T, std::enable_if_t<std::is_trivial<T>::value, void*> = nullptr>
    friend void serialize_internal(BinarySerializerContext& ctxt, const T& t)
    {
        BinarySerializerObject obj(ctxt.m_status, ctxt.m_flags, ctxt.m_stream);
        obj.add_element("", t);
    }

    template <class T, std::enable_if_t<std::is_trivial<T>::value, void*> = nullptr>
    friend IOResult<T> deserialize_internal(BinarySerializerContext& ctxt, Tag<T>)
    {
        BinarySerializerObject obj(ctxt.m_status, ctxt.m_flags, ctxt.m_stream);
        return obj.expect_element("", Tag<T>{});
    }

private:
    ByteStream& m_stream;
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
    auto ctxt = BinarySerializerContext(m_status, m_flags, m_stream);
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
    T t;
    if (m_stream.read(reinterpret_cast<unsigned char*>(std::addressof(t)), sizeof(t))) {
        return mio::success(t);
    }
    else {
        return mio::failure(mio::StatusCode::UnknownError, "Unexpected EOF reading " + name + " from binary stream.");
    }
}

template <class T, std::enable_if_t<negation<std::is_trivial<T>>::value, void*>>
IOResult<T> BinarySerializerObject::expect_element(const std::string& name, Tag<T> tag)
{
    mio::unused(name);
    auto ctxt = BinarySerializerContext(m_status, m_flags, m_stream);
    return mio::deserialize(ctxt, tag);
}

inline IOResult<std::string> BinarySerializerObject::expect_element(const std::string& name, Tag<std::string>)
{
    mio::unused(name);
    size_t size;
    if (m_stream.read(reinterpret_cast<unsigned char*>(&size), sizeof(size))) {
        std::string t(size, 0);
        if (m_stream.read(reinterpret_cast<unsigned char*>(&t[0]), size)) {
            return success(t);
        }
    }
    return failure(StatusCode::UnknownError, "Unexpected EOF reading " + name + " from binary stream.");
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
        v.push_back(t);
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
ByteStream serialize_binary(const T& t)
{
    ByteStream stream;
    BinarySerializerContext ctxt(std::make_shared<IOStatus>(), 0, stream);
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
IOResult<T> deserialize_binary(ByteStream& stream, Tag<T>)
{
    BinarySerializerContext ctxt(std::make_shared<IOStatus>(), 0, stream);
    return mio::deserialize(ctxt, Tag<T>{});
}

} // namespace mio

#endif
