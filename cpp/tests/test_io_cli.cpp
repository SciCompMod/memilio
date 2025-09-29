/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: René Schmieding
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
#include "matchers.h"
#include "memilio/io/cli.h"
#include "memilio/io/io.h"
#include "memilio/io/json_serializer.h"
#include "temp_file_register.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <list>

namespace
{

class TestCLI : public ::testing::Test
{
public:
    // allocates a char** argv, filled with the given values. arrays are freed during TearDown
    std::pair<int, char**> make_argv(const std::vector<std::string>& values)
    {
        const int argc = static_cast<int>(values.size());
        char** argv    = new char*[argc];
        // create and copy each string
        for (int i = 0; i < argc; i++) {
            const auto n = values[i].size() + 1;
            argv[i]      = new char[n];
            for (size_t j = 0; j < n; j++) {
                // make a manual strcpy, since msvc marks the function as unsafe
                // (this is not a concern as this function is only used for making test parameters)
                argv[i][j] = values[i][j];
            }
        }
        m_args.push_back({argc, argv});
        return m_args.back();
    }

    void TearDown() override
    {
        for (auto [argc, argv] : m_args) {
            free_argv(argc, argv);
        }
        m_args.clear();
    }

private:
    // deletes the char** created by make_argv
    void free_argv(int argc, char** argv)
    {
        for (int i = 0; i < argc; i++) {
            delete[] argv[i];
            argv[i] = nullptr;
        }
        delete[] argv;
    }

    std::vector<std::pair<int, char**>> m_args;
};

struct A {
    using Type = double;

    const static std::string name()
    {
        return "A";
    }
    const static std::string alias()
    {
        return "a";
    }
    const static std::string description()
    {
        return "This is A.";
    }
};

struct B {
    using Type = std::string;
    const static std::string name()
    {
        return "B";
    }
    // const static std::string alias()
    const static std::string description()
    {
        return "This is B. It has a needlessly long description that will most probably not fit into a single line on "
               "your terminal window, unless maybe on an ultrawide display.";
    }
};

struct C {
    using Type = std::vector<int>;
    const static std::string name()
    {
        return "C has a name that is too long";
    }
    const static std::string alias()
    {
        return "c";
    }
    // const static std::string description()
};

struct D {
    using Type = std::list<std::string>;
    const static std::string name()
    {
        return "D";
    }
    // const static std::string alias()
    // const static std::string description()
};

struct RequiredA {
    using Type = int;
    const static std::string name()
    {
        return "RequiredA";
    }
    // const static std::string alias()
    // const static std::string description()
    static bool is_required()
    {
        return true;
    }
};

struct RequiredB {
    using Type = int;
    const static std::string name()
    {
        return "RequiredB";
    }
    // const static std::string alias()
    // const static std::string description()
    static bool is_required()
    {
        return true;
    }
};

using Params = mio::ParameterSet<A, B, C, D>;

template <class T = char>
mio::cli::details::AbstractParameter test_param(const char* name, const char* alias, T* value = nullptr)
{
    return mio::cli::details::AbstractParameter(mio::Tag<T>{}, {name, alias}, std::shared_ptr<void>(value));
}

struct ErrorSerializable {
    void serialize(mio::JsonContext& c) const
    {
        c.set_error({mio::StatusCode::UnknownError, "Test"});
    }

    static mio::IOResult<ErrorSerializable> deserialize(mio::JsonContext&)
    {
        return mio::failure({mio::StatusCode::UnknownError, "Test"});
    }
};

} // namespace

TEST_F(TestCLI, test_identifier)
{
    // check that parsing, stripping and matching behave as intended
    using mio::cli::details::Identifier;
    using mio::cli::details::IdentifierType;
    mio::cli::details::DatalessParameter reference_parameter = {"A", "a"};
    // case: parse a name option; expect success, type set to Name, matching ref., strip to ref. name
    auto name_result = Identifier::parse("--A");
    ASSERT_THAT(name_result, IsSuccess());
    auto& name = name_result.value();
    EXPECT_EQ(name.type, IdentifierType::Name);
    EXPECT_EQ(name.strip(), reference_parameter.name);
    EXPECT_TRUE(name.matches_parameter(reference_parameter));
    // case: parse a alias option; expect success, type set to Alias, matching ref., strip to ref. alias
    auto alias_result = mio::cli::details::Identifier::parse("-a");
    auto& alias       = alias_result.value();
    EXPECT_EQ(alias.type, IdentifierType::Alias);
    EXPECT_EQ(alias.strip(), reference_parameter.alias);
    EXPECT_TRUE(alias.matches_parameter(reference_parameter));
    // case: parse raw name; expect failure
    auto raw_result = Identifier::parse("A");
    EXPECT_THAT(raw_result, IsFailure(mio::StatusCode::InvalidValue));
    // case: make raw name; expect type set to Raw, matching ref., strip to ref. name
    auto raw = Identifier::make_raw("A");
    EXPECT_EQ(raw.type, IdentifierType::Raw);
    EXPECT_EQ(raw.strip(), reference_parameter.name);
    EXPECT_TRUE(raw.matches_parameter(reference_parameter));
    // case: parse "obviously" invalid strings; expect failures
    EXPECT_THAT(Identifier::parse(""), IsFailure(mio::StatusCode::InvalidValue));
    EXPECT_THAT(Identifier::parse(" "), IsFailure(mio::StatusCode::InvalidValue));
    EXPECT_THAT(Identifier::parse("-"), IsFailure(mio::StatusCode::InvalidValue));
    EXPECT_THAT(Identifier::parse("--"), IsFailure(mio::StatusCode::InvalidValue));
}

TEST_F(TestCLI, test_option_verification)
{
    // check that the option verification, i.e. checking for duplicates and empty names, works correctly
    using mio::cli::details::AbstractSet;
    using mio::cli::details::Identifier;
    std::vector<mio::cli::details::AbstractParameter> params;
    // case: trivial list; expect success
    auto build_result = AbstractSet::build(params);
    EXPECT_THAT(print_wrap(build_result), IsSuccess());
    // case: good names and aliases; expect success and correctly filled maps
    params       = {test_param("A", "a"), test_param("B", ""), test_param("C", "c"), test_param("D", "")};
    build_result = AbstractSet::build(params);
    ASSERT_THAT(print_wrap(build_result), IsSuccess());
    auto set = build_result.value();
    EXPECT_TRUE(set.contains(Identifier::make_raw("A")));
    EXPECT_TRUE(set.contains(Identifier::make_raw("B")));
    EXPECT_TRUE(set.contains(Identifier::make_raw("C")));
    EXPECT_TRUE(set.contains(Identifier::make_raw("D")));
    EXPECT_TRUE(set.contains(Identifier::make_raw("a")));
    EXPECT_TRUE(set.contains(Identifier::make_raw("c")));
    EXPECT_FALSE(set.contains(Identifier::make_raw("E")));
    EXPECT_FALSE(set.contains(Identifier::make_raw("")));
    // case: alias collision; expect failure
    params.push_back(test_param("E", "a"));
    build_result = AbstractSet::build(params);
    EXPECT_THAT(print_wrap(build_result), IsFailure(mio::StatusCode::InvalidValue));
    params.pop_back();
    // case: name collision; expect failure
    params.push_back(test_param("A", "e"));
    build_result = AbstractSet::build(params);
    EXPECT_THAT(print_wrap(build_result), IsFailure(mio::StatusCode::InvalidValue));
    params.pop_back();
    // case: name empty; expect failure
    params.push_back(test_param("", "e"));
    build_result = AbstractSet::build(params);
    EXPECT_THAT(print_wrap(build_result), IsFailure(mio::StatusCode::InvalidValue));
}

TEST_F(TestCLI, test_get_param)
{
    // check that get_param works
    Params p;
    auto set = mio::cli::details::AbstractSet::build(p).value();
    const double AV{3.14159};
    const std::string BV{"BTestString"};
    // set params
    p.get<A>() = AV;
    p.get<B>() = BV;
    // test get by value
    // case: A with alias, expect success + value from name and alias
    auto a_by_name  = mio::cli::details::Identifier::parse("--A").value();
    auto a_by_alias = mio::cli::details::Identifier::parse("-a").value();
    ASSERT_TRUE(set.get_param(a_by_name).has_value());
    EXPECT_EQ(set.get_param(a_by_name).value(), AV);
    ASSERT_TRUE(set.get_param(a_by_alias).has_value());
    EXPECT_EQ(set.get_param(a_by_alias).value(), AV);
    // case: B without alias, expect success + value from name, error from alias
    auto b_by_name  = mio::cli::details::Identifier::make_raw("B");
    auto b_by_alias = mio::cli::details::Identifier::make_raw("b");
    ASSERT_TRUE(set.get_param(b_by_name).has_value());
    EXPECT_EQ(set.get_param(b_by_name).value(), BV);
    EXPECT_FALSE(set.get_param(b_by_alias).has_value());
    // case: D with default value, expect success + empty vector from name
    auto get_d_result = set.get_param(mio::cli::details::Identifier::make_raw("D"));
    ASSERT_TRUE(get_d_result.has_value());
    EXPECT_EQ(get_d_result.value().size(), Json::Value::ArrayIndex(0)); // empty vector
    // case: invalid name, expect error
    auto get_noname_result = set.get_param(mio::cli::details::Identifier::make_raw("NoName"));
    EXPECT_EQ(get_noname_result.error().formatted_message(),
              "Key not found: Could not get parameter: No such option \"NoName\".");
    EXPECT_EQ(get_noname_result.error().code(), mio::StatusCode::KeyNotFound);
    // case: error in serialize, expect error being forwarded

    EXPECT_THAT(print_wrap(test_param("", "", new ErrorSerializable).get()), IsFailure(mio::StatusCode::UnknownError));
}

TEST_F(TestCLI, test_set_param)
{
    Params p;
    auto set = mio::cli::details::AbstractSet::build(p).value();
    auto id  = mio::cli::details::Identifier::make_raw("A");
    // use normally
    EXPECT_THAT(print_wrap(set.set_param(id, std::string("5.2"))), IsSuccess());
    EXPECT_EQ(p.get<A>(), 5.2);
    // cause errors
    auto result = set.set_param(id, std::string("definitely a double"));
    EXPECT_THAT(print_wrap(result), IsFailure(mio::StatusCode::InvalidType));
    EXPECT_EQ(result.error().formatted_message(), "Invalid type: While setting \"A\": Json value is not a double.");
    result = set.set_param(mio::cli::details::Identifier::make_raw("NotAnOption"), std::string());
    EXPECT_THAT(print_wrap(result), IsFailure(mio::StatusCode::KeyNotFound));
    EXPECT_EQ(result.error().formatted_message(),
              "Key not found: Could not set parameter: No such option \"NotAnOption\".");
}

TEST_F(TestCLI, test_write_help)
{
    Params p;
    auto set = mio::cli::details::AbstractSet::build(p).value();
    std::stringstream ss;
    const std::string help_header = "Usage: TestSuite <option> <value> ...\n";
    const std::string help_header_with_defaults =
        "Usage: TestSuite [A D] <option> <value> ...\n"
        "Values of parameter options listed in [brackets] can be optionally entered in that order.\n";
    const std::string help_body =
        "All values must be entered as json values, i.e. the expression to the right of \"Name : \".\nNote that, when "
        "entering values, quotation marks may have to be escaped (\\\").\nOptions:\n  --help              -h     Show "
        "this dialogue and exit.\n  --print_option             Use with parameter option name(s) without \"--\" as "
        "value(s). Prints the current values of specified options in their correct json format, then exits.\n  "
        "--read_from_json           Takes a filepath as value. Reads and assigns parameter option values from the "
        "specified json file.\n  --write_to_json            Takes a filepath as value. Writes current values of all "
        "parameter options to the specified json file.\nParameter options:\n  --A                 -a     This is A.\n  "
        "--B                        This is B. It has a needlessly long description that will most probably not fit "
        "into a single line on your terminal window, unless maybe on an ultrawide display.\n  --C has a name that is "
        "too long\n                      -c\n  --D\n";
    // call write help directly
    mio::cli::details::write_help("TestSuite", set, {}, ss);
    EXPECT_EQ(ss.str(), help_header + help_body);
    ss.str("");
    ss.clear();
    // call write help directly, with defaults
    mio::cli::details::write_help("TestSuite", set, {"A", "D"}, ss);
    EXPECT_EQ(ss.str(), help_header_with_defaults + help_body);
    // call cli with -h option
    char c[]    = "-h", *argv[2];
    argv[1]     = c;
    auto result = mio::command_line_interface("TestSuite", 2, argv, p);
    ASSERT_THAT(print_wrap(result), IsFailure(mio::StatusCode::OK));
    EXPECT_EQ(result.error().code(), mio::StatusCode::OK);
    EXPECT_EQ(result.error().message(), help_header + help_body);
}

TEST_F(TestCLI, test_print_options)
{
    const std::vector<std::string> args{"", "--print_option", "a", "D"};
    auto [argc, argv] = make_argv(args);

    Params p;
    auto result = mio::command_line_interface("TestSuite", argc, argv, p);
    ASSERT_THAT(print_wrap(result), IsFailure(mio::StatusCode::OK));
    EXPECT_THAT(result.error().message(), "Option a:\n0.0\nOption D:\n[]\n");
}

TEST_F(TestCLI, test_import_export)
{
    TempFileRegister file_register;
    auto tmpfile                = file_register.get_unique_path("TestCLI-%%%%-%%%%.json");
    const std::string read_json = "{\"a\":5.4,\"D\":[\"d\",\"D\"]}";
    const std::string write_json =
        "{\n\t\"A\" : 5.4000000000000004,\n\t\"B\" : \"\",\n\t\"C has a name that is too long\" : "
        "[],\n\t\"D\" : \n\t[\n\t\t\"d\",\n\t\t\"D\"\n\t]\n}";
    Params p{};

    // write json input into tmpfile
    std::ofstream ofile(tmpfile);
    ofile << read_json;
    ofile.close();

    // test invalid file
    ASSERT_FALSE(mio::read_parameters_from_file(p, ""));

    // read tmpfile into p, then test its values
    ASSERT_TRUE(mio::read_parameters_from_file(p, tmpfile));
    EXPECT_EQ(p.get<A>(), 5.4);
    EXPECT_EQ(p.get<B>(), "");
    EXPECT_EQ(p.get<C>().size(), 0);
    EXPECT_EQ(p.get<D>().front(), "d");
    EXPECT_EQ(p.get<D>().back(), "D");

    // write the parameters back to tmpfile, check file contents
    ASSERT_TRUE(mio::write_parameters_to_file(p, tmpfile));
    std::ifstream ifile(tmpfile);
    std::string content{std::istreambuf_iterator<char>(ifile), std::istreambuf_iterator<char>()};
    ifile.close();
    EXPECT_EQ(content, write_json);

    // do read/write again using cli
    ofile.open(tmpfile);
    ofile << read_json;
    ofile.close();
    const std::vector<std::string> args{"", "--read_from_json", tmpfile, "--write_to_json", tmpfile};
    auto [argc, argv] = make_argv(args);
    ASSERT_TRUE(mio::command_line_interface("", argc, argv, p));
    ifile.open(tmpfile);
    content.assign(std::istreambuf_iterator<char>(ifile), std::istreambuf_iterator<char>());
    ifile.close();
    EXPECT_EQ(content, write_json);
}

TEST_F(TestCLI, test_default_options)
{
    // go through all possible routes default options can go through
    // generate argv
    const std::vector<std::string> args{"", "3.4", "\"TestValue\"", "[0, 8, 15]", "[\"Hello \",\"World!\"]"};
    auto [argc, argv] = make_argv(args);
    const std::vector<std::string> default_options{A::name(), B::name(), C::name(), D::name()};
    // case: set all params through defaults; expect success
    auto result = mio::command_line_interface<A, B, C, D>("TestSuite", argc, argv, default_options);
    ASSERT_THAT(print_wrap(result), IsSuccess());
    // also verify (on a subset of params) that values are being set
    EXPECT_EQ(result.value().get<A>(), 3.4);
    ASSERT_EQ(result.value().get<D>().size(), 2);
    // case: provide less params then defaults; expect success
    const int reduced_argc = argc - 2;
    result                 = mio::command_line_interface<A, B, C, D>("TestSuite", reduced_argc, argv, default_options);
    EXPECT_THAT(print_wrap(result), IsSuccess());
    EXPECT_EQ(result.value().get<A>(), 3.4);
    ASSERT_EQ(result.value().get<D>().size(), 0); // not set this time
    // case: use nonexistant param as default; expect failure (when verifying default options)
    const std::vector<std::string> wrong_options{"NoName"};
    EXPECT_THAT(print_wrap(mio::command_line_interface<A, B, C, D>("TestSuite", argc, argv, wrong_options)),
                IsFailure(mio::StatusCode::KeyNotFound));
    // case: provide more params then defaults; expect failure (when trying to parse identifier)
    const std::vector<std::string> reduced_options{"A", "B"};
    EXPECT_THAT(print_wrap(mio::command_line_interface<A, B, C, D>("TestSuite", argc, argv, reduced_options)),
                IsFailure(mio::StatusCode::InvalidValue));
    // missing case: pass invalid value; covered sufficiently through test_set_param
}

TEST_F(TestCLI, test_required_options)
{
    // go through all possible routes default options can go through
    // generate argv
    const std::vector<std::string> args_opt{"", "--RequiredA", "5", "--RequiredB", "2"};
    const std::vector<std::string> args_def{"", "5", "2"};
    auto [argc_opt, argv_opt] = make_argv(args_opt);
    auto [argc_def, argv_def] = make_argv(args_def);
    // case: set required params; expect success
    EXPECT_THAT(print_wrap(mio::command_line_interface<RequiredA, RequiredB>("TestSuite", argc_opt, argv_opt)),
                IsSuccess());
    // case: set required params through default option; expect success
    EXPECT_THAT(print_wrap(mio::command_line_interface<RequiredA, RequiredB>("TestSuite", argc_def, argv_def,
                                                                             {"RequiredA", "RequiredB"})),
                IsSuccess());
    // case: don't set required params; expect failure
    const int reduced_argc = argc_opt - 2;
    const int no_argc      = 1;
    EXPECT_THAT(print_wrap(mio::command_line_interface<RequiredA, RequiredB>("TestSuite", reduced_argc, argv_opt)),
                IsFailure(mio::StatusCode::InvalidValue));
    EXPECT_THAT(print_wrap(mio::command_line_interface<RequiredA, RequiredB>("TestSuite", no_argc, argv_opt)),
                IsFailure(mio::StatusCode::InvalidValue));
}

TEST_F(TestCLI, test_parameter_set_and_builder)
{
    // test that the builder and resulting set behave correctly
    // create a set using the builder
    auto parameters = mio::cli::ParameterSetBuilder().add<"A">(7).add<"B", double>(12, {"b", "desc", true}).build();
    // check that the DatalessParameter members are set correctly
    auto a = parameters.get<mio::TypeList<mio::cli::details::NamedType<"A">, int>>();
    EXPECT_EQ(a.name(), "A");
    EXPECT_EQ(a.alias(), "");
    EXPECT_EQ(a.description(), "");
    EXPECT_EQ(a.is_required(), false);
    auto b = parameters.get<mio::TypeList<mio::cli::details::NamedType<"B">, double>>();
    EXPECT_EQ(b.name(), "B");
    EXPECT_EQ(b.alias(), "b");
    EXPECT_EQ(b.description(), "desc");
    EXPECT_EQ(b.is_required(), true);
    // get and set values to test the parameter set
    EXPECT_EQ(parameters.get<"A">(), 7);
    parameters.get<"A">() = 2;
    EXPECT_EQ(parameters.get<"A">(), 2);
    EXPECT_EQ(parameters.get<"B">(), 12);
    parameters.get<"B">() = 0.125;
    EXPECT_EQ(parameters.get<"B">(), 0.125);
}

TEST_F(TestCLI, test_command_line_interface)
{
    // generate argv
    const std::vector<std::string> args{
        "",    "-a", "3.4", "--B", "\"TestValue\"",          "--C has a name that is too long",
        "[0,", "8,", "15]", "--D", "[\"Hello \",\"World!\"]"};
    auto [argc, argv] = make_argv(args);
    // trigger fail
    // by error in set_param
    EXPECT_FALSE(mio::command_line_interface<A>("TestSuite", argc, argv));
    // by missing value
    char* failv[3];
    failv[1]  = argv[1];
    failv[2]  = argv[3];
    auto fail = mio::command_line_interface<A, B, C, D>("TestSuite", 3, failv);
    ASSERT_FALSE(fail);
    EXPECT_EQ(fail.error().message(), "Missing value for option \"-a\".");
    // run cli normally
    auto result = mio::command_line_interface<A, B, C, D>("TestSuite", argc, argv);
    // check success
    ASSERT_THAT(print_wrap(result), IsSuccess());
    Params p = result.value();
    // test values
    // A
    EXPECT_EQ(p.get<A>(), 3.4);
    // B
    EXPECT_EQ(p.get<B>(), "TestValue");
    // C
    ASSERT_EQ(p.get<C>().size(), 3);
    EXPECT_EQ(p.get<C>()[0], 0);
    EXPECT_EQ(p.get<C>()[1], 8);
    EXPECT_EQ(p.get<C>()[2], 15);
    // D
    ASSERT_EQ(p.get<D>().size(), 2);
    EXPECT_EQ(p.get<D>().front(), "Hello ");
    EXPECT_EQ(p.get<D>().back(), "World!");
}

#endif // MEMILIO_HAS_JSONCPP
