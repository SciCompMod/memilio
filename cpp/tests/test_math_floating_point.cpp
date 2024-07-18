/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Rene Schmieding
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
#include "ad/ad.hpp"
#include "memilio/math/floating_point.h"

#include <gtest/gtest.h>

#include <array>
#include <string>

template <class FP>
class TestMathFloatingPoint : public ::testing::Test
{
public:
    using ParamType      = std::array<FP, 5>;
    using TruthTableType = std::array<std::array<bool, 3>, 4>;

    FP a = FP(1.1), b = FP(3.2); // arbitrary values s.t. a < b
    const ParamType params = {
        b - a, // abs_tol
        (b - a) / b, // rel_tol
        FP(0.0), // zero
        FP(10.0), // increased
        FP(0.1) // reduced
    };
};

using FpTypes = ::testing::Types<float, double, ad::gt1s<double>::type>;

TYPED_TEST_SUITE(TestMathFloatingPoint, FpTypes);

TYPED_TEST(TestMathFloatingPoint, abs_max)
{
    TypeParam v1 = this->a, v2 = this->b;
    // check for maximum
    EXPECT_TRUE(v2 == mio::abs_max(v1, v2));
    // check symmetries and signs
    EXPECT_TRUE(v2 == mio::abs_max(v2, v1));
    // specify type for signed comparisions, as -v may be an intermediate
    EXPECT_TRUE(v2 == mio::abs_max<TypeParam>(v1, -v2));
    EXPECT_TRUE(v2 == mio::abs_max<TypeParam>(-v2, v1));
    EXPECT_TRUE(v2 == mio::abs_max<TypeParam>(-v1, v2));
    EXPECT_TRUE(v2 == mio::abs_max<TypeParam>(v2, -v1));
    EXPECT_TRUE(v2 == mio::abs_max<TypeParam>(-v1, -v2));
    EXPECT_TRUE(v2 == mio::abs_max<TypeParam>(-v2, -v1));
}

/**
 * Check floating point comparisions of v1 and v2 against a truth table.
 * The table contents are described in the tests below.
 */
template <class FP, class Func>
void test_fp_compare(FP v1, FP v2, typename TestMathFloatingPoint<FP>::ParamType params, Func fp_compare,
                     typename TestMathFloatingPoint<FP>::TruthTableType truth_table)
{
    ASSERT_LT(v1, v2) << "This test is not set up correctly!";
    const auto [abs_tol, rel_tol, zero, increased, reduced] = params;
    // on error, print out the name of the test, which is the same as the tested function
    auto info = std::string("  With fp_compare as ") + testing::UnitTest::GetInstance()->current_test_info()->name();
    // check basics
    EXPECT_EQ(fp_compare(v1, v1, zero, zero), truth_table[0][0]) << info;
    EXPECT_EQ(fp_compare(v1, v2, zero, zero), truth_table[0][1]) << info;
    EXPECT_EQ(fp_compare(v2, v1, zero, zero), truth_table[0][2]) << info;
    // check equalities with absolute tolerances
    EXPECT_EQ(fp_compare(v1, v1, abs_tol, zero), truth_table[1][0]) << info;
    EXPECT_EQ(fp_compare(v1, v1, reduced * abs_tol, zero), truth_table[1][1]) << info;
    EXPECT_EQ(fp_compare(v1, v1, increased * abs_tol, zero), truth_table[1][2]) << info;
    // check equalities with relative tolerances
    EXPECT_EQ(fp_compare(v1, v1, zero, rel_tol), truth_table[1][0]) << info;
    EXPECT_EQ(fp_compare(v1, v1, zero, reduced * rel_tol), truth_table[1][1]) << info;
    EXPECT_EQ(fp_compare(v1, v1, zero, increased * rel_tol), truth_table[1][2]) << info;
    // check absolute tolerances
    EXPECT_EQ(fp_compare(v1, v2, abs_tol, zero), truth_table[2][0]) << info;
    EXPECT_EQ(fp_compare(v1, v2, reduced * abs_tol, zero), truth_table[2][1]) << info;
    EXPECT_EQ(fp_compare(v1, v2, increased * abs_tol, zero), truth_table[2][2]) << info;
    // check relative tolerances
    EXPECT_EQ(fp_compare(v1, v2, zero, rel_tol), truth_table[2][0]) << info;
    EXPECT_EQ(fp_compare(v1, v2, zero, reduced * rel_tol), truth_table[2][1]) << info;
    EXPECT_EQ(fp_compare(v1, v2, zero, increased * rel_tol), truth_table[2][2]) << info;
    // check absolute tolerances, with switched FPs
    EXPECT_EQ(fp_compare(v2, v1, abs_tol, zero), truth_table[3][0]) << info;
    EXPECT_EQ(fp_compare(v2, v1, reduced * abs_tol, zero), truth_table[3][1]) << info;
    EXPECT_EQ(fp_compare(v2, v1, increased * abs_tol, zero), truth_table[3][2]) << info;
    // check relative tolerances, with switched FPs
    EXPECT_EQ(fp_compare(v2, v1, zero, rel_tol), truth_table[3][0]) << info;
    EXPECT_EQ(fp_compare(v2, v1, zero, reduced * rel_tol), truth_table[3][1]) << info;
    EXPECT_EQ(fp_compare(v2, v1, zero, increased * rel_tol), truth_table[3][2]) << info;
}

TYPED_TEST(TestMathFloatingPoint, floating_point_equal)
{
    typename TestMathFloatingPoint<TypeParam>::TruthTableType truth_table = {{
        // {equal, strict less, strict greater}, both tolerances zero
        {true, false, false}, // basic fp operations without tolerances
        // {exact tolerance, reduced tolerance, increased tolerance} using either abs_tol or rel_tol
        {true, true, true}, // equality
        {true, false, true}, // less (or equal)
        {true, false, true} // greater (or equal)
    }};

    test_fp_compare(this->a, this->b, this->params, &mio::floating_point_equal<TypeParam>, truth_table);
}

TYPED_TEST(TestMathFloatingPoint, floating_point_less)
{
    typename TestMathFloatingPoint<TypeParam>::TruthTableType truth_table = {{
        // {equal, strict less, strict greater}, both tolerances zero
        {false, true, false}, // basic fp operations without tolerances
        // {exact tolerance, reduced tolerance, increased tolerance} using either abs_tol or rel_tol
        {false, false, false}, // equality
        {false, true, false}, // less (or equal)
        {false, false, false} // greater (or equal)
    }};
    test_fp_compare(this->a, this->b, this->params, &mio::floating_point_less<TypeParam>, truth_table);
}

TYPED_TEST(TestMathFloatingPoint, floating_point_less_equal)
{
    typename TestMathFloatingPoint<TypeParam>::TruthTableType truth_table = {{
        // {equal, strict less, strict greater}, both tolerances zero
        {true, true, false}, // basic fp operations without tolerances
        // {exact tolerance, reduced tolerance, increased tolerance} using either abs_tol or rel_tol
        {true, true, true}, // equality
        {true, true, true}, // less (or equal)
        {true, false, true} // greater (or equal)
    }};
    test_fp_compare(this->a, this->b, this->params, &mio::floating_point_less_equal<TypeParam>, truth_table);
}

TYPED_TEST(TestMathFloatingPoint, floating_point_greater)
{
    typename TestMathFloatingPoint<TypeParam>::TruthTableType truth_table = {{
        // {equal, strict less, strict greater}, both tolerances zero
        {false, false, true}, // basic fp operations without tolerances
        // {exact tolerance, reduced tolerance, increased tolerance} using either abs_tol or rel_tol
        {false, false, false}, // equality
        {false, false, false}, // less (or equal)
        {false, true, false} // greater (or equal)
    }};
    test_fp_compare(this->a, this->b, this->params, &mio::floating_point_greater<TypeParam>, truth_table);
}

TYPED_TEST(TestMathFloatingPoint, floating_point_greater_equal)
{
    typename TestMathFloatingPoint<TypeParam>::TruthTableType truth_table = {{
        // {equal, strict less, strict greater}, both tolerances zero
        {true, false, true}, // basic fp operations without tolerances
        // {exact tolerance, reduced tolerance, increased tolerance} using either abs_tol or rel_tol
        {true, true, true}, // equality
        {true, false, true}, // less (or equal)
        {true, true, true} // greater (or equal)
    }};
    test_fp_compare(this->a, this->b, this->params, &mio::floating_point_greater_equal<TypeParam>, truth_table);
}

TYPED_TEST(TestMathFloatingPoint, round_nth_decimal)
{
    using FP = double;
    FP value = static_cast<FP>(1.234567);
    EXPECT_EQ(mio::round_nth_decimal(value, 0), static_cast<FP>(1));
    EXPECT_EQ(mio::round_nth_decimal(value, 1), static_cast<FP>(1.2));
    EXPECT_EQ(mio::round_nth_decimal(value, 2), static_cast<FP>(1.23));
    EXPECT_EQ(mio::round_nth_decimal(value, 3), static_cast<FP>(1.235));
    EXPECT_EQ(mio::round_nth_decimal(value, 4), static_cast<FP>(1.2346));
    EXPECT_EQ(mio::round_nth_decimal(value, 5), static_cast<FP>(1.23457));
    EXPECT_EQ(mio::round_nth_decimal(value, 6), static_cast<FP>(1.234567));

    value = static_cast<FP>(-1.234567);
    EXPECT_EQ(mio::round_nth_decimal(value, 0), static_cast<FP>(-1));
    EXPECT_EQ(mio::round_nth_decimal(value, 1), static_cast<FP>(-1.2));
    EXPECT_EQ(mio::round_nth_decimal(value, 2), static_cast<FP>(-1.23));
    EXPECT_EQ(mio::round_nth_decimal(value, 3), static_cast<FP>(-1.235));
    EXPECT_EQ(mio::round_nth_decimal(value, 4), static_cast<FP>(-1.2346));
    EXPECT_EQ(mio::round_nth_decimal(value, 5), static_cast<FP>(-1.23457));
    EXPECT_EQ(mio::round_nth_decimal(value, 6), static_cast<FP>(-1.234567));

    value = static_cast<FP>(0.999999);
    EXPECT_EQ(mio::round_nth_decimal(value, 0), static_cast<FP>(1));
    EXPECT_EQ(mio::round_nth_decimal(value, 1), static_cast<FP>(1.0));
    EXPECT_EQ(mio::round_nth_decimal(value, 2), static_cast<FP>(1.0));
    EXPECT_EQ(mio::round_nth_decimal(value, 3), static_cast<FP>(1.0));
    EXPECT_EQ(mio::round_nth_decimal(value, 4), static_cast<FP>(1.0));
    EXPECT_EQ(mio::round_nth_decimal(value, 5), static_cast<FP>(1.0));
    EXPECT_EQ(mio::round_nth_decimal(value, 6), static_cast<FP>(0.999999));
}
