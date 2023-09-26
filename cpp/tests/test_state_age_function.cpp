/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Anna Wendler, Lena Ploetzke
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

#include "boost/fusion/functional/invocation/invoke.hpp"
#include "load_test_data.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/model.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/state_age_function.h"
#include <iostream>
#include <gtest/gtest.h>

TEST(TestStateAgeFunction, testSpecialMember)
{
    ScalarType dt = 0.5;
    mio::SmootherCosine smoothcos(1.0);

    // constructor
    EXPECT_EQ(smoothcos.get_parameter(), 1.0);
    EXPECT_EQ(smoothcos.get_support_max(dt), 1.0);

    // copy
    mio::SmootherCosine smoothcos2(smoothcos);
    EXPECT_EQ(smoothcos.get_state_age_function_type(), smoothcos2.get_state_age_function_type());
    EXPECT_EQ(smoothcos.get_parameter(), smoothcos2.get_parameter());
    EXPECT_EQ(smoothcos.get_support_max(dt), smoothcos2.get_support_max(dt));

    // check copy is true copy, not reference
    smoothcos.set_parameter(2.0);
    EXPECT_NE(smoothcos.get_parameter(), smoothcos2.get_parameter());
    smoothcos.set_parameter(1.0);

    // move
    mio::SmootherCosine smoothcos3(std::move(smoothcos2));
    EXPECT_EQ(smoothcos.get_state_age_function_type(), smoothcos3.get_state_age_function_type());
    EXPECT_EQ(smoothcos.get_parameter(), smoothcos3.get_parameter());
    EXPECT_EQ(smoothcos.get_support_max(dt), smoothcos3.get_support_max(dt));

    // copy assignment
    mio::SmootherCosine smoothcos4 = smoothcos3;
    EXPECT_EQ(smoothcos.get_state_age_function_type(), smoothcos4.get_state_age_function_type());
    EXPECT_EQ(smoothcos.get_parameter(), smoothcos4.get_parameter());
    EXPECT_EQ(smoothcos.get_support_max(dt), smoothcos4.get_support_max(dt));

    // check copy is true copy, not reference
    smoothcos.set_parameter(2.0);
    EXPECT_NE(smoothcos.get_parameter(), smoothcos4.get_parameter());
    smoothcos.set_parameter(1.0);

    // move assignment
    mio::SmootherCosine smoothcos5 = std::move(smoothcos4);
    EXPECT_EQ(smoothcos.get_state_age_function_type(), smoothcos5.get_state_age_function_type());
    EXPECT_EQ(smoothcos.get_parameter(), smoothcos5.get_parameter());
    EXPECT_EQ(smoothcos.get_support_max(dt), smoothcos5.get_support_max(dt));

    // also test the constructor of ExponentialDecay and ConstantFunction
    // copy and move (assignment) are defined in base class StateAgeFunction and are equal for all derived classes
    mio::ExponentialDecay expdecay(1.0);
    EXPECT_EQ(expdecay.get_parameter(), 1.0);

    mio::ConstantFunction constfunc(1.0);
    EXPECT_EQ(constfunc.get_parameter(), 1.0);
}

TEST(TestStateAgeFunction, testSettersAndGettersForParameter)
{
    ScalarType testvalue_before = 1.0;
    ScalarType testvalue_after  = 2.0;

    // test get and set for function parameter
    // only for SmootherCosine as set_parameter and get_parameter are equal for all derived classes
    mio::SmootherCosine smoothcos(testvalue_before);
    EXPECT_EQ(smoothcos.get_parameter(), testvalue_before);

    smoothcos.set_parameter(testvalue_after);
    EXPECT_EQ(smoothcos.get_parameter(), testvalue_after);
}

TEST(TestStateAgeFunction, testGetSupportMax)
{
    ScalarType dt = 0.5;

    // test get_support_max for all derived classes as this method can be overridden
    // Check that the maximum support is correct after setting the parameter object of a StateAgeFunction.
    mio::ExponentialDecay expdecay(1.0);
    EXPECT_NEAR(expdecay.get_support_max(dt), 23.5, 1e-14);
    expdecay.set_parameter(2.0);
    EXPECT_NEAR(expdecay.get_support_max(dt), 12.0, 1e-14);

    mio::SmootherCosine smoothcos(1.0);
    EXPECT_NEAR(smoothcos.get_support_max(dt), 1.0, 1e-14);
    smoothcos.set_parameter(2.0);
    EXPECT_NEAR(smoothcos.get_support_max(dt), 2.0, 1e-14);

    mio::ConstantFunction constfunc(1.0);
    EXPECT_NEAR(constfunc.get_support_max(dt), -2.0, 1e-14);
    constfunc.set_parameter(2.0);
    EXPECT_NEAR(constfunc.get_support_max(dt), -2.0, 1e-14);
}

TEST(TestStateAgeFunction, testSAFWrapperSpecialMember)
{
    ScalarType dt = 0.5;
    mio::SmootherCosine smoothcos(1.0);
    mio::StateAgeFunctionWrapper wrapper(smoothcos);

    // constructor
    EXPECT_EQ(wrapper.get_parameter(), 1.0);
    EXPECT_EQ(wrapper.get_support_max(dt), 1.0);

    // copy
    mio::StateAgeFunctionWrapper wrapper2(wrapper);
    EXPECT_EQ(wrapper.get_state_age_function_type(), wrapper2.get_state_age_function_type());
    EXPECT_EQ(wrapper.get_parameter(), wrapper2.get_parameter());
    EXPECT_EQ(wrapper.get_support_max(dt), wrapper2.get_support_max(dt));

    // test true copy, not reference
    wrapper.set_parameter(2.0);
    EXPECT_NE(wrapper.get_parameter(), wrapper2.get_parameter());
    wrapper.set_parameter(1.0);

    // move
    mio::StateAgeFunctionWrapper wrapper3(std::move(wrapper2));
    EXPECT_EQ(wrapper.get_state_age_function_type(), wrapper3.get_state_age_function_type());
    EXPECT_EQ(wrapper.get_parameter(), wrapper3.get_parameter());
    EXPECT_EQ(wrapper.get_support_max(dt), wrapper3.get_support_max(dt));

    // copy assignment
    mio::StateAgeFunctionWrapper wrapper4 = wrapper3;
    EXPECT_EQ(wrapper.get_state_age_function_type(), wrapper4.get_state_age_function_type());
    EXPECT_EQ(wrapper.get_parameter(), wrapper4.get_parameter());
    EXPECT_EQ(wrapper.get_support_max(dt), wrapper4.get_support_max(dt));

    // test true copy, not reference
    wrapper.set_parameter(2.0);
    EXPECT_NE(wrapper.get_parameter(), wrapper4.get_parameter());
    wrapper.set_parameter(1.0);

    // move assignment
    mio::StateAgeFunctionWrapper wrapper5 = std::move(wrapper4);
    EXPECT_EQ(wrapper.get_state_age_function_type(), wrapper5.get_state_age_function_type());
    EXPECT_EQ(wrapper.get_parameter(), wrapper5.get_parameter());
    EXPECT_EQ(wrapper.get_support_max(dt), wrapper5.get_support_max(dt));

    // Also test the constructor of StateAgeFunctionWrapper initialized with ExponentialDecay and ConstantFunction.
    // This way we can be sure that clone_impl works for all derived classes of StateAgeFunction.
    mio::ExponentialDecay expdecay(1.0);
    mio::StateAgeFunctionWrapper wrapper_exp(expdecay);
    EXPECT_NEAR(wrapper_exp.get_parameter(), 1.0, 1e-14);
    EXPECT_NEAR(wrapper_exp.get_support_max(dt), 23.5, 1e-14);

    mio::ConstantFunction constfunc(1.0);
    mio::StateAgeFunctionWrapper wrapper_const(constfunc);
    EXPECT_NEAR(wrapper_const.get_parameter(), 1.0, 1e-14);
    EXPECT_NEAR(wrapper_const.get_support_max(dt), -2.0, 1e-14);
}

TEST(TestStateAgeFunction, testSAFWrapperSettersAndGetters)
{
    // test only for SmootherCosine since we have checked set_parameter and get_parameter
    // as well as get_support_max for all derived classes in previous test

    ScalarType testvalue_before = 1.0;
    ScalarType testvalue_after  = 2.0;
    ScalarType dt               = 0.5;

    mio::SmootherCosine smoothcos(testvalue_before);
    mio::StateAgeFunctionWrapper wrapper(smoothcos);

    EXPECT_EQ(wrapper.get_parameter(), testvalue_before);
    EXPECT_EQ(wrapper.get_support_max(dt), testvalue_before);

    wrapper.set_parameter(testvalue_after);
    EXPECT_EQ(wrapper.get_parameter(), testvalue_after);
    EXPECT_EQ(wrapper.get_support_max(dt), testvalue_after);
}

TEST(TestStateAgeFunction, testComparisonOperator)
{
    // Check that StateAgeFunctions are only considered equal if they are of the same derived
    // class and have the same parameter
    mio::ExponentialDecay expdecay(0.5);
    mio::ExponentialDecay expdecay2(0.5);
    mio::ExponentialDecay expdecay3(1.0);
    mio::SmootherCosine smoothcos(0.5);

    EXPECT_TRUE(expdecay == expdecay2);
    EXPECT_FALSE(expdecay == expdecay3);
    EXPECT_FALSE(expdecay == smoothcos);

    // Check that it also holds when a StateAgeFunctionWrapper is set with the respective functions
    mio::StateAgeFunctionWrapper wrapper(expdecay);
    mio::StateAgeFunctionWrapper wrapper2(expdecay2);
    mio::StateAgeFunctionWrapper wrapper3(expdecay3);
    mio::StateAgeFunctionWrapper wrapper4(smoothcos);

    EXPECT_TRUE(wrapper == wrapper2);
    EXPECT_FALSE(wrapper == wrapper3);
    EXPECT_FALSE(wrapper == wrapper4);
}