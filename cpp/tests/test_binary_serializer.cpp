#include "boost/none.hpp"
#include "boost/none_t.hpp"
#include "matchers.h"
#include "memilio/io/binary_serializer.h"
#include "memilio/io/io.h"
#include "gtest/gtest.h"

namespace {
    template<class T>
void check_roundtrip_serialize(const T& t)
{
    //binary serialized stream has no defined portable format
    //so no need to check contents of the stream after serialize
    //roundtrip check is sufficient 
    auto stream = mio::serialize_binary(t);
    auto result = mio::deserialize_binary(stream, mio::Tag<T>{});
    EXPECT_THAT(result, IsSuccess());
    EXPECT_EQ(result.value(), t);    
}
}

TEST(BinarySerializer, int)
{
    check_roundtrip_serialize(3);
}

TEST(BinarySerializer, float)
{
    check_roundtrip_serialize(3.1234f);
}

namespace {
enum class E
{
    E1 = 1,
    E2 = 2,
};
}

TEST(BinarySerializer, enum)
{
    check_roundtrip_serialize(E::E2);
}

namespace {

struct Foo {
    int i;
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Foo");
        obj.add_element("i", i);
    }

    template <class IOContext>
    static mio::IOResult<Foo> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Foo");
        auto i   = obj.expect_element("i", mio::Tag<int>{});
        return mio::apply(
            io,
            [](auto i_) {
                return Foo{i_};
            },
            i);
    }
    bool operator==(const Foo& other) const
    {
        return i == other.i;
    }
    bool operator!=(const Foo& other) const
    {
        return !(*this == other);
    }
};

}

TEST(BinarySerializer, simple_class)
{
    check_roundtrip_serialize(Foo{4});
}

struct Bar {
    std::string s;
    boost::optional<E> e;
    std::vector<Foo> v;

    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Bar");
        obj.add_element("s", s);
        obj.add_optional("e", e ? &e.value() : nullptr);
        obj.add_list("v", v.begin(), v.end());
    }
    template <class IOContext>
    static mio::IOResult<Bar> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Bar");
        auto s   = obj.expect_element("s", mio::Tag<std::string>{});
        auto e = obj.expect_optional("e", mio::Tag<E>{});
        auto v   = obj.expect_list("v", mio::Tag<Foo>{});
        return mio::apply(
            io,
            [](auto&& s_, auto&& e_, auto&& v_) {
                return Bar{s_, e_, v_};
            },
            s, e, v);
    }
    bool operator==(const Bar& other) const
    {
        return std::tie(s, e, v) == std::tie(other.s, other.e, other.v);
    }
    bool operator!=(const Bar& other) const
    {
        return !(*this == other);
    }
};

TEST(BinarySerializer, complex_class)
{
    check_roundtrip_serialize(Bar{"Hello", E::E1, {Foo{2}, Foo{3}}});
    check_roundtrip_serialize(Bar{"", boost::none, {}});
}

TEST(ByteStream, end)
{
    mio::ByteStream stream;
    int i = 3;
    auto p = (unsigned char*)&i;
    stream.write(p, sizeof(i));
    EXPECT_TRUE(stream.read(p, sizeof(i)));
    EXPECT_FALSE(stream.read(p, sizeof(i)));
}

TEST(ByteStream, data)
{
    mio::ByteStream stream(4);
    int i = 3;
    auto p = (unsigned char*)&i;
    std::copy(p, p + sizeof(i), stream.data());
    i = 4;
    EXPECT_TRUE(stream.read(p, sizeof(i)));
    EXPECT_EQ(i, 3);
}

TEST(ByteStream, reset)
{
    mio::ByteStream stream(0);
    stream.reset(4);
    int i;
    auto p = (unsigned char*)&i;
    EXPECT_TRUE(stream.read(p, sizeof(i)));
    EXPECT_EQ(i, 0);
}
