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
#include "boost/none.hpp"
#include "boost/none_t.hpp"
#include "matchers.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/epidemiology/damping_sampling.h"
#include "memilio/io/binary_serializer.h"
#include "memilio/io/io.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/uncertain_value.h"
#include "ode_secir/model.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/parameters.h"
#include "gtest/gtest.h"
#include <memory>

namespace
{
template <class T>
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
} // namespace

TEST(BinarySerializer, int)
{
    check_roundtrip_serialize(3);
}

TEST(BinarySerializer, float)
{
    check_roundtrip_serialize(3.1234f);
}

namespace
{
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

namespace
{

struct Foo {
    int i = 0;
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

} // namespace

TEST(BinarySerializer, simple_class)
{
    check_roundtrip_serialize(Foo{4});
}

namespace
{

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
        auto e   = obj.expect_optional("e", mio::Tag<E>{});
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

} // namespace

TEST(BinarySerializer, complex_class)
{
    check_roundtrip_serialize(Bar{"Hello", E::E1, {Foo{2}, Foo{3}}});
    check_roundtrip_serialize(Bar{"", boost::none, {}});
}

TEST(BinarySerializer, time_series)
{
    mio::TimeSeries<double> ts{48};
    for (int i = 0; i < 10; ++i) {
        ts.add_time_point(i, Eigen::VectorXd::Constant(48, double(i)));
    }
    auto stream = mio::serialize_binary(ts);
    auto result = mio::deserialize_binary(stream, mio::Tag<mio::TimeSeries<double>>{});
    EXPECT_THAT(result, IsSuccess());
    EXPECT_EQ(result.value().get_num_time_points(), Eigen::Index(10));
    EXPECT_EQ(result.value()[7][0], double(7));
}

TEST(BinarySerializer, model)
{
    //this test is only to make sure the correct number of bytes are serialized/deserialized
    //in a very complex object. correct serializing of single values is tested by other tests.
    mio::osecir::Model model{5};
    mio::set_log_level(mio::LogLevel::err);
    mio::osecir::set_params_distributions_normal(model, 0, 10, 0.01);
    mio::set_log_level(mio::LogLevel::warn);
    auto stream = mio::serialize_binary(model);
    auto result = mio::deserialize_binary(stream, mio::Tag<mio::osecir::Model>{});
    EXPECT_THAT(result, IsSuccess());
}

TEST(BinarySerializer, type_check)
{
    Foo foo;
    auto stream = mio::serialize_binary(foo, mio::IOF_IncludeTypeInfo);
    auto r = mio::deserialize_binary(stream, mio::Tag<Foo>{}, mio::IOF_IncludeTypeInfo);
    EXPECT_THAT(r, IsSuccess());
}

TEST(BinarySerializer, fail_type_check)
{
    Foo foo;
    auto stream = mio::serialize_binary(foo, mio::IOF_IncludeTypeInfo);
    auto r = mio::deserialize_binary(stream, mio::Tag<Bar>{}, mio::IOF_IncludeTypeInfo);
    EXPECT_THAT(r, IsFailure(mio::StatusCode::InvalidType));
}

TEST(ByteStream, end)
{
    mio::ByteStream stream;
    int i  = 3;
    auto p = (unsigned char*)&i;
    stream.write(p, sizeof(i));
    EXPECT_TRUE(stream.read(p, sizeof(i)));
    EXPECT_FALSE(stream.read(p, sizeof(i)));
}

TEST(ByteStream, data)
{
    mio::ByteStream stream(4);
    int i  = 3;
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
