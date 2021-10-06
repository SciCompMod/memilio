/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "memilio/io/json_serializer.h"
#include "memilio/utils/type_safe.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/uncertain_value.h"
#include "matchers.h"
#include "distributions_helpers.h"
#include "gtest/gtest.h"
#include <vector>

namespace jsontest
{
struct Foo {
    int i;
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Foo");
        obj.add_element("i", i);
    }

    template<class IOContext>
    static mio::IOResult<Foo> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Foo");
        auto i = obj.expect_element("i", mio::Tag<int>{});
        return mio::apply(io, [](auto i_) { return Foo{i_}; }, i);
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

struct Bar {
    std::string s;
    std::vector<Foo> v;

    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Bar");
        obj.add_element("s", s);
        obj.add_list("v", v.begin(), v.end());
    }
    template<class IOContext>
    static mio::IOResult<Bar> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Bar");
        auto s = obj.expect_element("s", mio::Tag<std::string>{});
        auto v = obj.expect_list("v", mio::Tag<Foo>{});
        return mio::apply(io, [](auto&& s_, auto&& v_) { return Bar {s_, v_}; }, s, v);
    }
    bool operator==(const Bar& other) const
    {
        return s == other.s && v == other.v;
    }
    bool operator!=(const Bar& other) const
    {
        return !(*this == other);
    }
};

struct DoubleList {
    std::vector<double> vd;
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("DoubleList");
        obj.add_list("vd", vd.begin(), vd.end());
    }
    template<class IOContext>
    static mio::IOResult<DoubleList> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("DoubleList");
        auto v = obj.expect_list("vd", mio::Tag<double>{});
        return mio::apply(io, [](auto&& v_) { return DoubleList{v_}; }, v);
    }

    bool operator==(const DoubleList& other) const
    {
        return vd == other.vd;
    }
    bool operator!=(const DoubleList& other) const
    {
        return !(*this == other);
    }
};

} // namespace jsontest

TEST(TestJsonSerializer, basic_type)
{
    mio::JsonSerializer js;
    serialize(js, 3);
    ASSERT_EQ(js.value(), Json::Value(3));

    mio::JsonSerializer js2{Json::Value(3)};
    auto r = deserialize(js2, mio::Tag<int>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    ASSERT_EQ(r.value(), 3);
}

TEST(TestJsonSerializer, basic_type_float)
{
    mio::JsonSerializer js;
    serialize(js, 3.0f);
    ASSERT_EQ(js.value(), Json::Value(3.0));

    mio::JsonSerializer js2{Json::Value(3.0)};
    auto r = deserialize(js2, mio::Tag<float>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    ASSERT_EQ(r.value(), 3.0f);

    mio::JsonSerializer js3{Json::Value(1e200)};
    ASSERT_THAT(print_wrap(deserialize(js3, mio::Tag<float>{})), testing::Not(IsSuccess()));
}

namespace jsontest
{
enum class E : int
{
    E0 = 0,
    E1,
    E2
};
} // namespace jsontest

TEST(TestJsonSerializer, enum)
{
    mio::JsonSerializer js;
    serialize(js, jsontest::E::E1);
    ASSERT_EQ(js.value(), Json::Value(1));

    mio::JsonSerializer js2{Json::Value(1)};
    auto r = deserialize(js2, mio::Tag<jsontest::E>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    EXPECT_EQ(r.value(), jsontest::E::E1);
}

TEST(TestJsonSerializer, tuple)
{
    mio::JsonSerializer js;
    auto tup = std::make_tuple(1, 2.0, std::string("Hello"));
    serialize(js, tup);
    auto expected_json        = Json::Value();
    expected_json["Element0"] = 1;
    expected_json["Element1"] = 2.0;
    expected_json["Element2"] = "Hello";
    EXPECT_EQ(js.value(), expected_json);

    mio::JsonSerializer js2{expected_json};
    auto r = deserialize(js2, mio::Tag<std::tuple<int, double, std::string>>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    EXPECT_EQ(r.value(), tup);
}

TEST(TestJsonSerializer, doublelist)
{
    mio::JsonSerializer js;
    jsontest::DoubleList dl{{0.1, 0.2, 0.3, 0.4}};
    serialize(js, dl);
    auto expected_json     = Json::Value();
    expected_json["vd"][0] = 0.1;
    expected_json["vd"][1] = 0.2;
    expected_json["vd"][2] = 0.3;
    expected_json["vd"][3] = 0.4;
    ASSERT_EQ(js.value(), expected_json);

    mio::JsonSerializer js2{expected_json};
    auto r = deserialize(js2, mio::Tag<jsontest::DoubleList>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    EXPECT_EQ(r.value(), dl);
}

TEST(TestJsonSerializer, aggregate)
{
    mio::JsonSerializer js;
    jsontest::Bar bar{"Hello", {{1}, {2}}};
    serialize(js, bar);
    auto expected_json         = Json::Value();
    expected_json["s"]         = "Hello";
    expected_json["v"][0]["i"] = 1;
    expected_json["v"][1]["i"] = 2;
    ASSERT_EQ(js.value(), expected_json);

    mio::JsonSerializer js2{expected_json};
    auto r = mio::deserialize(js2, mio::Tag<jsontest::Bar>());
    ASSERT_THAT(print_wrap(r), IsSuccess());
    EXPECT_EQ(r.value(), (jsontest::Bar{"Hello", {jsontest::Foo{1}, jsontest::Foo{2}}}));
}

namespace jsontest
{
DECL_TYPESAFE(int, TypeSafeInt);
}

TEST(TestJsonSerializer, typesafe)
{
    mio::JsonSerializer js;
    serialize(js, jsontest::TypeSafeInt(1));
    ASSERT_EQ(js.value(), Json::Value(1));

    mio::JsonSerializer js2{Json::Value(1)};
    auto r = deserialize(js2, mio::Tag<jsontest::TypeSafeInt>{});
    EXPECT_TRUE(r);
    EXPECT_EQ(r.value(), jsontest::TypeSafeInt(1));
}

namespace jsontest
{
struct Tag {
};
} // namespace jsontest

TEST(TestJsonSerializer, customindexarray)
{
    mio::JsonSerializer js;
    mio::CustomIndexArray<double, jsontest::Tag> a(mio::Index<jsontest::Tag>(2));
    a[mio::Index<jsontest::Tag>(0)] = 1.0;
    a[mio::Index<jsontest::Tag>(1)] = 2.0;
    serialize(js, a);
    Json::Value expected_value;
    expected_value["Dimensions"]  = Json::UInt64(2);
    expected_value["Elements"][0] = 1.0;
    expected_value["Elements"][1] = 2.0;
    EXPECT_EQ(js.value(), expected_value);

    mio::JsonSerializer js2{expected_value};
    auto r = deserialize(js2, mio::Tag<mio::CustomIndexArray<double, jsontest::Tag>>{});
    ASSERT_THAT(r, IsSuccess());
    EXPECT_EQ(r.value().size(), mio::Index<jsontest::Tag>(2));
    EXPECT_THAT(r.value(), testing::ElementsAre(1.0, 2.0));
}

TEST(TestJsonSerializer, normal_distribution)
{
    mio::JsonSerializer js;
    auto dist = mio::ParameterDistributionNormal(0.0, 1.0, 0.5, 0.1);
    Json::Value expected_value;
    expected_value["Type"] = "Normal";
    expected_value["Mean"] = 0.5;
    expected_value["StandardDev"] = 0.1;
    expected_value["LowerBound"] = 0.0;
    expected_value["UpperBound"] = 1.0;
    expected_value["PredefinedSamples"] = Json::Value(Json::arrayValue);
    serialize(js, dist);
    EXPECT_EQ(js.value(), expected_value);

    mio::JsonSerializer js2(expected_value);
    auto r = deserialize(js2, mio::Tag<std::shared_ptr<mio::ParameterDistribution>>{});
    EXPECT_THAT(print_wrap(r), IsSuccess());
    check_distribution(*r.value(), dist);
}

TEST(TestJsonSerializer, uniform_distribution)
{
    mio::JsonSerializer js;
    auto dist = mio::ParameterDistributionUniform(0.0, 1.0);
    Json::Value expected_value;
    expected_value["Type"] = "Uniform";
    expected_value["LowerBound"] = 0.0;
    expected_value["UpperBound"] = 1.0;
    expected_value["PredefinedSamples"] = Json::Value(Json::arrayValue);
    serialize(js, dist);
    EXPECT_EQ(js.value(), expected_value);

    mio::JsonSerializer js2(expected_value);
    auto r = deserialize(js2, mio::Tag<std::shared_ptr<mio::ParameterDistribution>>{});
    EXPECT_THAT(print_wrap(r), IsSuccess());
    check_distribution(*r.value(), dist);
}

TEST(TestJsonSerializer, serialize_uv)
{
    mio::JsonSerializer js;
    mio::UncertainValue uv(2.0);
    Json::Value expected_value;
    expected_value["Value"] = 2.0;
    serialize(js, uv);    
    EXPECT_EQ(js.value(), expected_value);
    
    uv.set_distribution(mio::ParameterDistributionNormal(-1.0, 2.0, 0.5, 0.1));
    serialize(js, uv);
    expected_value["Distribution"]["Type"] = "Normal";
    expected_value["Distribution"]["Mean"] = 0.5;
    expected_value["Distribution"]["StandardDev"] = 0.1;
    expected_value["Distribution"]["LowerBound"] = -1.0;
    expected_value["Distribution"]["UpperBound"] = 2.0;
    expected_value["Distribution"]["PredefinedSamples"] = Json::Value(Json::arrayValue);
    EXPECT_EQ(js.value(), expected_value);
}

TEST(TestJsonSerializer, deserialize_uv)
{    
    Json::Value json_uv;
    json_uv["Value"] = 2.0;
    {
        mio::JsonSerializer js{json_uv};
        auto r = mio::deserialize(js, mio::Tag<mio::UncertainValue>{});
        EXPECT_TRUE(r);
        EXPECT_EQ(double(r.value()), 2.0);
        EXPECT_EQ(r.value().get_distribution(), nullptr);
    }

    json_uv["Distribution"]["Type"] = "Normal";
    json_uv["Distribution"]["LowerBound"] = -1.0;
    json_uv["Distribution"]["UpperBound"] = 2.0;
    json_uv["Distribution"]["Mean"] = 0.5;
    json_uv["Distribution"]["StandardDev"] = 0.1;
    json_uv["Distribution"]["PredefinedSamples"] = Json::Value(Json::arrayValue);
    {
        mio::JsonSerializer js{json_uv};
        auto r = mio::deserialize(js, mio::Tag<mio::UncertainValue>{});
        EXPECT_TRUE(r);
        EXPECT_EQ(double(r.value()), 2.0);
        EXPECT_NE(r.value().get_distribution(), nullptr);
    }
}

namespace jsontest
{
    struct Param1
    {
        using Type = double;
        static constexpr Type get_default()
        {
            return 1.0;
        }
        static std::string name()
        {
            return "Param1";
        }
    };

    struct Param2
    {
        using Type = std::string;
        static Type get_default()
        {
            return "Hello";
        }
        static std::string name()
        {
            return "Param2";
        }
    };

    using ParameterSet = mio::ParameterSet<Param1, Param2>;
}

TEST(TestJsonSerializer, paramset)
{
    mio::JsonSerializer js;
    jsontest::ParameterSet params;
    serialize(js, params);
    Json::Value expected_value;
    expected_value["Param1"] = 1.0;
    expected_value["Param2"] = "Hello";
    EXPECT_EQ(js.value(), expected_value);

    mio::JsonSerializer js2{expected_value};
    auto r = deserialize(js2, mio::Tag<jsontest::ParameterSet>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    EXPECT_EQ(params, r.value());
}

TEST(TestJsonSerializer, matrix)
{    
    mio::JsonSerializer js;
    Eigen::MatrixXd m(2, 3);
    m << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
    serialize(js, m);
    Json::Value expected_value{Json::objectValue};
    expected_value["Rows"] = Json::Int64(2);
    expected_value["Columns"] = Json::Int64(3);
    expected_value["Elements"].append(1.0);
    expected_value["Elements"][1] = 2.0;
    expected_value["Elements"][2] = 3.0;
    expected_value["Elements"][3] = 4.0;
    expected_value["Elements"][4] = 5.0;
    expected_value["Elements"][5] = 6.0;
    EXPECT_EQ(js.value(), expected_value);

    mio::JsonSerializer js2(expected_value);
    auto r = mio::deserialize(js2, mio::Tag<Eigen::MatrixXd>{});
    EXPECT_TRUE(r);
    EXPECT_EQ(print_wrap(m), print_wrap(r.value()));
}

TEST(TestJsonSerializer, container)
{
    mio::JsonSerializer js;
    std::vector<int> v{1, 2, 3};
    mio::serialize(js, v);
    Json::Value expected_value;
    expected_value["Items"][0] = 1;
    expected_value["Items"][1] = 2;
    expected_value["Items"][2] = 3;
    EXPECT_EQ(js.value(), expected_value);

    mio::JsonSerializer js2{expected_value};
    auto r = mio::deserialize(js2, mio::Tag<std::vector<int>>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    EXPECT_THAT(r.value(), testing::ElementsAreArray(v.data(), 3));
}