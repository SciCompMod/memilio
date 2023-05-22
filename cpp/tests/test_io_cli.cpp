

/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include "memilio/config.h"
#include "memilio/io/cli.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <sstream>
#include <string>
#include <vector>
#include <list>

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

struct NameCollision {
    using Type = std::list<std::string>;
    const static std::string name()
    {
        return "A";
    }
    // const static std::string alias()
    // const static std::string description()
};

struct AliasCollision {
    using Type = std::list<std::string>;
    const static std::string name()
    {
        return "nonempty";
    }
    const static std::string alias()
    {
        return "a";
    }
    // const static std::string description()
};

struct NoName {
    using Type = std::list<std::string>;
    const static std::string name()
    {
        return "";
    }
    // const static std::string alias()
    // const static std::string description()
};

using Params = mio::ParameterSet<A, B, C, D>;
// using BadParams = mio::ParameterSet<A, CollisionA>;
TEST(TestCLI, test_option_verifier)
{
    EXPECT_DEBUG_DEATH(mio::details::cli::verify_options(mio::ParameterSet<A, NameCollision>()),
                       ".*Options may not have duplicate fields\\. \\(field required\\)");
    EXPECT_DEBUG_DEATH(mio::details::cli::verify_options(mio::ParameterSet<A, AliasCollision>()),
                       ".*Options may not have duplicate fields\\. \\(field optional\\)");
    EXPECT_DEBUG_DEATH(mio::details::cli::verify_options(mio::ParameterSet<NoName>()),
                       ".*Option is missing required field\\..*");
    EXPECT_DEBUG_DEATH(
        {
            mio::details::cli::verify_options(Params());
            assert(false && "DeathTestMessage");
        },
        "DeathTestMessage");
}

TEST(TestCLI, test_get_param)
{
    Params p;
    const double AV{3.14159};
    const std::string BV{"BTestString"};
    // set params
    p.get<A>() = AV;
    p.get<B>() = BV;
    // test get by value
    // A (with alias)
    ASSERT_TRUE(mio::details::cli::get_param(p, "A").has_value());
    EXPECT_EQ(mio::details::cli::get_param(p, "A").value(), AV);
    ASSERT_TRUE(mio::details::cli::get_param(p, "a").has_value());
    EXPECT_EQ(mio::details::cli::get_param(p, "a").value(), AV);
    // B (without alias)
    ASSERT_TRUE(mio::details::cli::get_param(p, "B").has_value());
    EXPECT_EQ(mio::details::cli::get_param(p, "B").value(), BV);
    EXPECT_FALSE(mio::details::cli::get_param(p, "b").has_value());
    // D (default value)
    ASSERT_TRUE(mio::details::cli::get_param(p, "D").has_value());
    EXPECT_EQ(mio::details::cli::get_param(p, "D").value().size(), Json::Value::ArrayIndex(0)); // empty vector
    // Errors
    EXPECT_EQ(mio::details::cli::get_param(p, "NoName").error().formatted_message(),
              "Key not found: No such option \"NoName\".");
    EXPECT_EQ(mio::details::cli::get_param(p, "NoName").error().code(), mio::StatusCode::KeyNotFound);
}

TEST(TestCLI, test_set_param)
{
    Params p;
    // set by name
    EXPECT_TRUE(mio::details::cli::set_param(p, "--A", "5.2"));
    EXPECT_EQ(p.get<A>(), 5.2);
    // set by alias
    EXPECT_TRUE(mio::details::cli::set_param(p, "-a", "2.7"));
    EXPECT_EQ(p.get<A>(), 2.7);
    // Errors
    EXPECT_EQ(mio::details::cli::set_param(p, "--A", "definitely a double").error().formatted_message(),
              "Invalid type: While setting \"A\": Json value is not a double.");
    EXPECT_EQ(mio::details::cli::set_param(p, "NoName", "").error().formatted_message(),
              "Key not found: Expected an option, got \"NoName\".");
    EXPECT_EQ(mio::details::cli::set_param(p, "--NoName", "").error().formatted_message(),
              "Key not found: No such option \"NoName\".");
}

TEST(TestCLI, test_write_help)
{
    std::stringstream ss;
    const std::string help =
        "Usage: TestSuite <option> <value> ...\nValues must be entered as json values, i.e. the expression to the "
        "right of \"Name : \".\nOptions:\n  --help              -h     Show this dialogue and exit.\n  --print_option  "
        "           Use with other option name(s) without \"--\" as value(s). Prints the current values of specified "
        "options in their correct json format, then exits. Note that quotation marks may have to be escaped (\\\")\n  "
        "--A                 -a     This is A.\n  --B                        This is B. It has a needlessly long "
        "description that will most probably not fit into a single line on your terminal window, unless maybe on an "
        "ultrawide display.\n  --C has a name that is too long\n                      -c     \n  --D                   "
        "     \n";
    // call write help directly
    mio::details::cli::write_help("TestSuite", Params(), ss);
    EXPECT_EQ(ss.str(), help);
    // call cli with -h option
    char c[] = "-h", *argv[2];
    argv[1]  = c;
    Params p;
    EXPECT_EXIT((void)mio::command_line_interface("TestSuite", 2, argv, p, std::cerr), testing::ExitedWithCode(0),
                testing::StrEq(help));
}

TEST(TestCLI, test_print_options)
{
    std::vector<std::string> args{"", "--print_option", "a", "D"};
    const int argc = 4;
    char* argv[argc];
    for (int i = 0; i < argc; i++) {
        const auto n = args[i].size() + 1;
        argv[i]      = new char[n];
        for (unsigned j = 0; j < n; j++) {
            argv[i][j] = args[i][j];
        }
    }
    Params p;
    EXPECT_EXIT((void)mio::command_line_interface("TestSuite", argc, argv, p, std::cerr), testing::ExitedWithCode(0),
                testing::StrEq("Option a:\n0.0\nOption D:\n[]\n"));
    for (int i = 0; i < argc; i++) {
        delete[] argv[i];
    }
}

TEST(TestCLI, test_command_line_interface)
{
    // generate argv
    std::vector<std::string> args{
        "",    "-a", "3.4", "--B", "\"TestValue\"",          "--C has a name that is too long",
        "[0,", "8,", "15]", "--D", "[\"Hello \",\"World!\"]"};
    const int argc = 11;
    char* argv[argc];
    for (int i = 0; i < argc; i++) {
        const auto n = args[i].size() + 1;
        argv[i]      = new char[n];
        for (unsigned j = 0; j < n; j++) {
            argv[i][j] = args[i][j];
        }
    }
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
    ASSERT_TRUE(result);
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
    // deallocate argv
    for (int i = 0; i < argc; i++) {
        delete[] argv[i];
    }
}

#endif // MEMILIO_HAS_JSONCPP
