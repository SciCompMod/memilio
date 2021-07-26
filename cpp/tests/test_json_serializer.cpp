#include "epidemiology_io/json_serializer.h"
#include "epidemiology/utils/type_safe.h"
#include "epidemiology/utils/custom_index_array.h"
#include "epidemiology/utils/parameter_set.h"
#include "epidemiology/utils/uncertain_value.h"
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
    static epi::IOResult<Foo> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Foo");
        auto i = obj.expect_element("i", epi::Tag<int>{});
        return epi::apply(io, [](auto i_) { return Foo{i_}; }, i);
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
    static epi::IOResult<Bar> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Bar");
        auto s = obj.expect_element("s", epi::Tag<std::string>{});
        auto v = obj.expect_list("v", epi::Tag<Foo>{});
        return epi::apply(io, [](auto&& s_, auto&& v_) { return Bar {s_, v_}; }, s, v);
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
    static epi::IOResult<DoubleList> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("DoubleList");
        auto v = obj.expect_list("vd", epi::Tag<double>{});
        return epi::apply(io, [](auto&& v_) { return DoubleList{v_}; }, v);
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
    epi::JsonSerializer js;
    serialize(js, 3);
    ASSERT_EQ(js.value(), Json::Value(3));

    epi::JsonSerializer js2{Json::Value(3)};
    auto r = deserialize(js2, epi::Tag<int>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    ASSERT_EQ(r.value(), 3);
}

TEST(TestJsonSerializer, basic_type_float)
{
    epi::JsonSerializer js;
    serialize(js, 3.0f);
    ASSERT_EQ(js.value(), Json::Value(3.0));

    epi::JsonSerializer js2{Json::Value(3.0)};
    auto r = deserialize(js2, epi::Tag<float>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    ASSERT_EQ(r.value(), 3.0f);

    epi::JsonSerializer js3{Json::Value(1e200)};
    ASSERT_THAT(print_wrap(deserialize(js3, epi::Tag<float>{})), testing::Not(IsSuccess()));
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
    epi::JsonSerializer js;
    serialize(js, jsontest::E::E1);
    ASSERT_EQ(js.value(), Json::Value(1));

    epi::JsonSerializer js2{Json::Value(1)};
    auto r = deserialize(js2, epi::Tag<jsontest::E>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    EXPECT_EQ(r.value(), jsontest::E::E1);
}

TEST(TestJsonSerializer, tuple)
{
    epi::JsonSerializer js;
    auto tup = std::make_tuple(1, 2.0, std::string("Hello"));
    serialize(js, tup);
    auto expected_json        = Json::Value();
    expected_json["Element0"] = 1;
    expected_json["Element1"] = 2.0;
    expected_json["Element2"] = "Hello";
    EXPECT_EQ(js.value(), expected_json);

    epi::JsonSerializer js2{expected_json};
    auto r = deserialize(js2, epi::Tag<std::tuple<int, double, std::string>>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    EXPECT_EQ(r.value(), tup);
}

TEST(TestJsonSerializer, doublelist)
{
    epi::JsonSerializer js;
    jsontest::DoubleList dl{{0.1, 0.2, 0.3, 0.4}};
    serialize(js, dl);
    auto expected_json     = Json::Value();
    expected_json["vd"][0] = 0.1;
    expected_json["vd"][1] = 0.2;
    expected_json["vd"][2] = 0.3;
    expected_json["vd"][3] = 0.4;
    ASSERT_EQ(js.value(), expected_json);

    epi::JsonSerializer js2{expected_json};
    auto r = deserialize(js2, epi::Tag<jsontest::DoubleList>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    EXPECT_EQ(r.value(), dl);
}

TEST(TestJsonSerializer, aggregate)
{
    epi::JsonSerializer js;
    jsontest::Bar bar{"Hello", {{1}, {2}}};
    serialize(js, bar);
    auto expected_json         = Json::Value();
    expected_json["s"]         = "Hello";
    expected_json["v"][0]["i"] = 1;
    expected_json["v"][1]["i"] = 2;
    ASSERT_EQ(js.value(), expected_json);

    epi::JsonSerializer js2{expected_json};
    auto r = epi::deserialize(js2, epi::Tag<jsontest::Bar>());
    ASSERT_THAT(print_wrap(r), IsSuccess());
    EXPECT_EQ(r.value(), (jsontest::Bar{"Hello", {jsontest::Foo{1}, jsontest::Foo{2}}}));
}

namespace jsontest
{
DECL_TYPESAFE(int, TypeSafeInt);
}

TEST(TestJsonSerializer, typesafe)
{
    epi::JsonSerializer js;
    serialize(js, jsontest::TypeSafeInt(1));
    ASSERT_EQ(js.value(), Json::Value(1));

    epi::JsonSerializer js2{Json::Value(1)};
    auto r = deserialize(js2, epi::Tag<jsontest::TypeSafeInt>{});
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
    epi::JsonSerializer js;
    epi::CustomIndexArray<double, jsontest::Tag> a(epi::Index<jsontest::Tag>(2));
    a[epi::Index<jsontest::Tag>(0)] = 1.0;
    a[epi::Index<jsontest::Tag>(1)] = 2.0;
    serialize(js, a);
    Json::Value expected_value;
    expected_value["Dimensions"]  = Json::UInt64(2);
    expected_value["Elements"][0] = 1.0;
    expected_value["Elements"][1] = 2.0;
    EXPECT_EQ(js.value(), expected_value);

    epi::JsonSerializer js2{expected_value};
    auto r = deserialize(js2, epi::Tag<epi::CustomIndexArray<double, jsontest::Tag>>{});
    ASSERT_THAT(r, IsSuccess());
    EXPECT_EQ(r.value().size(), epi::Index<jsontest::Tag>(2));
    EXPECT_THAT(r.value(), testing::ElementsAre(1.0, 2.0));
}

TEST(TestJsonSerializer, normal_distribution)
{
    epi::JsonSerializer js;
    auto dist = epi::ParameterDistributionNormal(0.0, 1.0, 0.5, 0.1);
    Json::Value expected_value;
    expected_value["Type"] = "Normal";
    expected_value["Mean"] = 0.5;
    expected_value["StandardDev"] = 0.1;
    expected_value["LowerBound"] = 0.0;
    expected_value["UpperBound"] = 1.0;
    expected_value["PredefinedSamples"] = Json::Value(Json::arrayValue);
    serialize(js, dist);
    EXPECT_EQ(js.value(), expected_value);

    epi::JsonSerializer js2(expected_value);
    auto r = deserialize(js2, epi::Tag<std::shared_ptr<epi::ParameterDistribution>>{});
    EXPECT_THAT(print_wrap(r), IsSuccess());
    check_distribution(*r.value(), dist);
}

TEST(TestJsonSerializer, uniform_distribution)
{
    epi::JsonSerializer js;
    auto dist = epi::ParameterDistributionUniform(0.0, 1.0);
    Json::Value expected_value;
    expected_value["Type"] = "Uniform";
    expected_value["LowerBound"] = 0.0;
    expected_value["UpperBound"] = 1.0;
    expected_value["PredefinedSamples"] = Json::Value(Json::arrayValue);
    serialize(js, dist);
    EXPECT_EQ(js.value(), expected_value);

    epi::JsonSerializer js2(expected_value);
    auto r = deserialize(js2, epi::Tag<std::shared_ptr<epi::ParameterDistribution>>{});
    EXPECT_THAT(print_wrap(r), IsSuccess());
    check_distribution(*r.value(), dist);
}

TEST(TestJsonSerializer, serialize_uv)
{
    epi::JsonSerializer js;
    epi::UncertainValue uv(2.0);
    Json::Value expected_value;
    expected_value["Value"] = 2.0;
    serialize(js, uv);    
    EXPECT_EQ(js.value(), expected_value);
    
    uv.set_distribution(epi::ParameterDistributionNormal(-1.0, 2.0, 0.5, 0.1));
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
        epi::JsonSerializer js{json_uv};
        auto r = epi::deserialize(js, epi::Tag<epi::UncertainValue>{});
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
        epi::JsonSerializer js{json_uv};
        auto r = epi::deserialize(js, epi::Tag<epi::UncertainValue>{});
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

    using ParameterSet = epi::ParameterSet<Param1, Param2>;
}

TEST(TestJsonSerializer, paramset)
{
    epi::JsonSerializer js;
    jsontest::ParameterSet params;
    serialize(js, params);
    Json::Value expected_value;
    expected_value["Param1"] = 1.0;
    expected_value["Param2"] = "Hello";
    EXPECT_EQ(js.value(), expected_value);

    epi::JsonSerializer js2{expected_value};
    auto r = deserialize(js2, epi::Tag<jsontest::ParameterSet>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    EXPECT_EQ(params, r.value());
}

TEST(TestJsonSerializer, matrix)
{    
    epi::JsonSerializer js;
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

    epi::JsonSerializer js2(expected_value);
    auto r = epi::deserialize(js2, epi::Tag<Eigen::MatrixXd>{});
    EXPECT_TRUE(r);
    EXPECT_EQ(print_wrap(m), print_wrap(r.value()));
}

TEST(TestJsonSerializer, container)
{
    epi::JsonSerializer js;
    std::vector<int> v{1, 2, 3};
    epi::serialize(js, v);
    Json::Value expected_value;
    expected_value["Items"][0] = 1;
    expected_value["Items"][1] = 2;
    expected_value["Items"][2] = 3;
    EXPECT_EQ(js.value(), expected_value);

    epi::JsonSerializer js2{expected_value};
    auto r = epi::deserialize(js2, epi::Tag<std::vector<int>>{});
    ASSERT_THAT(print_wrap(r), IsSuccess());
    EXPECT_THAT(r.value(), testing::ElementsAreArray(v.data(), 3));
}