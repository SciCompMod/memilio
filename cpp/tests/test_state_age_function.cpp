/*
* Copyright (C) 2020-2025 MEmilio
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

#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include <gtest/gtest.h>

TEST(TestStateAgeFunction, testSpecialMember)
{ // Copy and move (assignment) are defined in base class StateAgeFunction and are equal for all derived classes, therefore test for one Special Member.
    // Use GammaSurvivalFunction as all parameters are used for this StateAgeFunction.
    // Constructors of other members will be also tested.
    ScalarType dt = 0.5;
    mio::GammaSurvivalFunction gamma(1., 2., 3.);

    // constructor
    EXPECT_EQ(gamma.get_distribution_parameter(), 1.);
    EXPECT_EQ(gamma.get_location(), 2.);
    EXPECT_EQ(gamma.get_scale(), 3.);
    EXPECT_NEAR(gamma.get_support_max(dt), 71.5, 1e-14);
    EXPECT_NEAR(gamma.get_mean(), 5., 1e-14);

    // copy
    mio::GammaSurvivalFunction gamma2(gamma);
    EXPECT_EQ(gamma.get_state_age_function_type(), gamma2.get_state_age_function_type());
    EXPECT_EQ(gamma.get_distribution_parameter(), gamma2.get_distribution_parameter());
    EXPECT_EQ(gamma.get_location(), gamma2.get_location());
    EXPECT_EQ(gamma.get_scale(), gamma2.get_scale());
    EXPECT_EQ(gamma.get_support_max(dt), gamma2.get_support_max(dt));
    EXPECT_EQ(gamma.get_mean(dt), gamma2.get_mean(dt));

    // check copy is true copy, not reference.
    gamma.set_distribution_parameter(1.5);
    EXPECT_NE(gamma.get_distribution_parameter(), gamma2.get_distribution_parameter());
    gamma.set_distribution_parameter(1.0);

    // move
    mio::GammaSurvivalFunction gamma3(std::move(gamma2));
    EXPECT_EQ(gamma.get_state_age_function_type(), gamma3.get_state_age_function_type());
    EXPECT_EQ(gamma.get_distribution_parameter(), gamma3.get_distribution_parameter());
    EXPECT_EQ(gamma.get_location(), gamma3.get_location());
    EXPECT_EQ(gamma.get_scale(), gamma3.get_scale());
    EXPECT_EQ(gamma.get_support_max(dt), gamma3.get_support_max(dt));
    EXPECT_EQ(gamma.get_mean(dt), gamma3.get_mean(dt));

    // copy assignment
    mio::GammaSurvivalFunction gamma4 = gamma3;
    EXPECT_EQ(gamma.get_state_age_function_type(), gamma4.get_state_age_function_type());
    EXPECT_EQ(gamma.get_distribution_parameter(), gamma4.get_distribution_parameter());
    EXPECT_EQ(gamma.get_location(), gamma4.get_location());
    EXPECT_EQ(gamma.get_scale(), gamma4.get_scale());
    EXPECT_EQ(gamma.get_support_max(dt), gamma4.get_support_max(dt));
    EXPECT_EQ(gamma.get_mean(dt), gamma4.get_mean(dt));

    // check copy is true copy, not reference.
    gamma.set_scale(2.5);
    EXPECT_NE(gamma.get_scale(), gamma4.get_scale());
    gamma.set_scale(3.0);

    // move assignment
    mio::GammaSurvivalFunction gamma5 = std::move(gamma4);
    EXPECT_EQ(gamma.get_state_age_function_type(), gamma5.get_state_age_function_type());
    EXPECT_EQ(gamma.get_distribution_parameter(), gamma5.get_distribution_parameter());
    EXPECT_EQ(gamma.get_location(), gamma5.get_location());
    EXPECT_EQ(gamma.get_scale(), gamma5.get_scale());
    EXPECT_EQ(gamma.get_support_max(dt), gamma5.get_support_max(dt));
    EXPECT_EQ(gamma.get_mean(dt), gamma5.get_mean(dt));

    // Test constructors of the other derived classes.
    // Copy and move (assignment) are defined in base class StateAgeFunction and are equal for all derived classes.
    mio::ExponentialSurvivalFunction<double> exponential(1.5, 0.5);
    EXPECT_EQ(exponential.get_distribution_parameter(), 1.5);
    EXPECT_EQ(exponential.get_location(), 0.5);

    mio::SmootherCosine<double> smoothcos(2.0, 1.3);
    EXPECT_EQ(smoothcos.get_distribution_parameter(), 2.0);
    EXPECT_EQ(smoothcos.get_location(), 1.3);

    mio::LognormSurvivalFunction logn(0.1, 3.3, 0.5);
    EXPECT_EQ(logn.get_distribution_parameter(), 0.1);
    EXPECT_EQ(logn.get_location(), 3.3);
    EXPECT_EQ(logn.get_scale(), 0.5);

    mio::ConstantFunction<double> constfunc(1.4);
    EXPECT_EQ(constfunc.get_distribution_parameter(), 1.4);

    mio::ErlangDensity erl(2, 0.7);
    EXPECT_EQ(erl.get_distribution_parameter(), 2);
    EXPECT_EQ(erl.get_scale(), 0.7);
}

TEST(TestStateAgeFunction, testSettersAndGettersForParameters)
{
    // Deactivate temporarily log output for next tests. Errors are expected here.
    mio::set_log_level(mio::LogLevel::off);

    ScalarType testvalue_before = 0.5;
    ScalarType testvalue_after  = 2.1;

    // Test getter and setter for function parameter.
    // Only test for GammaSurvivalfunction as setter and getter are defined in the base class and thus equal for all derived classes.
    // Construct GammaSurvivalFunction with an invalid value for the scale parameter.
    mio::GammaSurvivalFunction gamma(testvalue_before, testvalue_before, -1);
    // Check if the invalid scale parameter is set to one.
    EXPECT_EQ(gamma.get_scale(), 1.);

    EXPECT_EQ(gamma.get_distribution_parameter(), testvalue_before);
    gamma.set_distribution_parameter(testvalue_after);
    EXPECT_EQ(gamma.get_distribution_parameter(), testvalue_after);

    EXPECT_EQ(gamma.get_location(), testvalue_before);
    gamma.set_location(testvalue_after);
    EXPECT_EQ(gamma.get_location(), testvalue_after);

    // Try to set scale to an invalid value and check if the value is corrected.
    gamma.set_scale(-1);
    EXPECT_EQ(gamma.get_scale(), 1);
    gamma.set_scale(testvalue_after);
    EXPECT_EQ(gamma.get_scale(), testvalue_after);

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestStateAgeFunction, testGetSupportMax)
{
    ScalarType dt = 0.5;

    // Test get_support_max for all derived classes as this method can be overridden.
    // Check that the maximum support is correct after setting the parameters of a StateAgeFunction.
    mio::ExponentialSurvivalFunction<double> exponential(1.0, 1.0);
    EXPECT_NEAR(exponential.get_support_max(dt), 24.5, 1e-14);
    exponential.set_distribution_parameter(2.0);
    EXPECT_NEAR(exponential.get_support_max(dt), 13.0, 1e-14);
    exponential.set_location(1.0);
    EXPECT_NEAR(exponential.get_support_max(dt), 13.0, 1e-14);
    exponential.set_scale(300.0);
    EXPECT_NEAR(exponential.get_support_max(dt), 13.0, 1e-14);

    mio::SmootherCosine<double> smoothcos(1.0, 1.);
    EXPECT_NEAR(smoothcos.get_support_max(dt), 2.0, 1e-14);
    smoothcos.set_distribution_parameter(2.0);
    EXPECT_NEAR(smoothcos.get_support_max(dt), 3.0, 1e-14);
    smoothcos.set_location(0);
    EXPECT_NEAR(smoothcos.get_support_max(dt), 2.0, 1e-14);
    smoothcos.set_scale(300.0);
    EXPECT_NEAR(smoothcos.get_support_max(dt), 2.0, 1e-14);

    mio::GammaSurvivalFunction gamma(1.0, 1, 0.5);
    EXPECT_NEAR(gamma.get_support_max(dt), 13, 1e-14);
    gamma.set_distribution_parameter(0.5);
    EXPECT_NEAR(gamma.get_support_max(dt), 11.5, 1e-14);
    gamma.set_location(0.0);
    EXPECT_NEAR(gamma.get_support_max(dt), 10.5, 1e-14);
    gamma.set_scale(1.0);
    EXPECT_NEAR(gamma.get_support_max(dt), 21, 1e-14);

    mio::LognormSurvivalFunction logn(0.1, 1.0, 0.5);
    EXPECT_NEAR(logn.get_support_max(dt), 2.0, 1e-14);
    logn.set_distribution_parameter(0.5);
    EXPECT_NEAR(logn.get_support_max(dt), 13.5, 1e-14);
    logn.set_location(0.0);
    EXPECT_NEAR(logn.get_support_max(dt), 12.5, 1e-14);
    logn.set_scale(0.1);
    EXPECT_NEAR(logn.get_support_max(dt), 2.5, 1e-14);

    // Deactivate temporarily log output for next tests.
    // Errors are expected here as the ConstantFunction is not suited to be a TransitionDistribution.
    // Support_max would be infinity here but is set to -2.0 in the get_support_max() method.
    mio::set_log_level(mio::LogLevel::off);

    mio::ConstantFunction<double> constfunc(1.0);
    EXPECT_NEAR(constfunc.get_support_max(dt), -2.0, 1e-14);
    constfunc.set_distribution_parameter(2.0);
    EXPECT_NEAR(constfunc.get_support_max(dt), -2.0, 1e-14);
    constfunc.set_location(2.0);
    EXPECT_NEAR(constfunc.get_support_max(dt), -2.0, 1e-14);
    constfunc.set_scale(2.0);
    EXPECT_NEAR(constfunc.get_support_max(dt), -2.0, 1e-14);

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);

    mio::ErlangDensity erl(2, 0.5);
    EXPECT_NEAR(erl.get_support_max(dt), 14, 1e-14);
    erl.set_distribution_parameter(1);
    EXPECT_NEAR(erl.get_support_max(dt), 12, 1e-14);
    erl.set_scale(0.1);
    EXPECT_NEAR(erl.get_support_max(dt), 3, 1e-14);
    erl.set_location(300);
    EXPECT_NEAR(erl.get_support_max(dt), 3, 1e-14);
}

TEST(TestStateAgeFunction, testSAFWrapperSpecialMember)
{ // Same as testSpecialMember for Wrapper.
    ScalarType dt = 0.5;
    mio::GammaSurvivalFunction gamma(1., 2., 3.);
    mio::StateAgeFunctionWrapper<double> wrapper(gamma);

    // constructor
    EXPECT_EQ(wrapper.get_distribution_parameter(), 1.);
    EXPECT_EQ(wrapper.get_location(), 2.);
    EXPECT_EQ(wrapper.get_scale(), 3.);
    EXPECT_NEAR(wrapper.get_support_max(dt), 71.5, 1e-14);
    EXPECT_NEAR(wrapper.get_mean(), 5., 1e-14);

    // copy
    mio::StateAgeFunctionWrapper<double> wrapper2(wrapper);
    EXPECT_EQ(wrapper.get_state_age_function_type(), wrapper2.get_state_age_function_type());
    EXPECT_EQ(wrapper.get_distribution_parameter(), wrapper2.get_distribution_parameter());
    EXPECT_EQ(wrapper.get_location(), wrapper2.get_location());
    EXPECT_EQ(wrapper.get_scale(), wrapper2.get_scale());
    EXPECT_EQ(wrapper.get_support_max(dt), wrapper2.get_support_max(dt));
    EXPECT_EQ(wrapper.get_mean(dt), wrapper2.get_mean(dt));

    // Check copy is true copy, not reference.
    wrapper.set_distribution_parameter(1.5);
    EXPECT_NE(wrapper.get_distribution_parameter(), wrapper2.get_distribution_parameter());
    wrapper.set_distribution_parameter(1.0);

    // move
    mio::StateAgeFunctionWrapper<double> wrapper3(std::move(wrapper2));
    EXPECT_EQ(wrapper.get_state_age_function_type(), wrapper3.get_state_age_function_type());
    EXPECT_EQ(wrapper.get_distribution_parameter(), wrapper3.get_distribution_parameter());
    EXPECT_EQ(wrapper.get_location(), wrapper3.get_location());
    EXPECT_EQ(wrapper.get_scale(), wrapper3.get_scale());
    EXPECT_EQ(wrapper.get_support_max(dt), wrapper3.get_support_max(dt));
    EXPECT_EQ(wrapper.get_mean(dt), wrapper3.get_mean(dt));

    // copy assignment
    mio::StateAgeFunctionWrapper<double> wrapper4 = wrapper3;
    EXPECT_EQ(wrapper.get_state_age_function_type(), wrapper4.get_state_age_function_type());
    EXPECT_EQ(wrapper.get_distribution_parameter(), wrapper4.get_distribution_parameter());
    EXPECT_EQ(wrapper.get_location(), wrapper4.get_location());
    EXPECT_EQ(wrapper.get_scale(), wrapper4.get_scale());
    EXPECT_EQ(wrapper.get_support_max(dt), wrapper4.get_support_max(dt));
    EXPECT_EQ(wrapper.get_mean(dt), wrapper4.get_mean(dt));

    // Check copy is true copy, not reference.
    wrapper.set_scale(2.3);
    EXPECT_NE(wrapper.get_scale(), wrapper4.get_scale());
    wrapper.set_scale(3.0);

    // move assignment
    mio::StateAgeFunctionWrapper<double> wrapper5 = std::move(wrapper4);
    EXPECT_EQ(wrapper.get_state_age_function_type(), wrapper5.get_state_age_function_type());
    EXPECT_EQ(wrapper.get_distribution_parameter(), wrapper5.get_distribution_parameter());
    EXPECT_EQ(wrapper.get_location(), wrapper5.get_location());
    EXPECT_EQ(wrapper.get_scale(), wrapper5.get_scale());
    EXPECT_EQ(wrapper.get_support_max(dt), wrapper5.get_support_max(dt));
    EXPECT_EQ(wrapper.get_mean(dt), wrapper5.get_mean(dt));

    // Also test the constructor of StateAgeFunctionWrapper initialized with other StateAgeFunction%s.
    // This way we can be sure that clone_impl works for all derived classes of StateAgeFunction.
    mio::ExponentialSurvivalFunction<double> exponential(1., 1.5);
    mio::StateAgeFunctionWrapper<double> wrapper_exp(exponential);
    EXPECT_EQ(wrapper_exp.get_distribution_parameter(), 1.0);
    EXPECT_EQ(wrapper_exp.get_location(), 1.5);
    EXPECT_NEAR(wrapper_exp.get_support_max(dt), 25, 1e-14);
    EXPECT_NEAR(wrapper_exp.get_mean(dt), 2.5, 1e-14);

    mio::SmootherCosine<double> smoothcos(2., 1.3);
    mio::StateAgeFunctionWrapper<double> wrapper_smoothcos(smoothcos);
    EXPECT_EQ(wrapper_smoothcos.get_distribution_parameter(), 2.);
    EXPECT_EQ(wrapper_smoothcos.get_location(), 1.3);
    EXPECT_NEAR(wrapper_smoothcos.get_support_max(dt), 2.0 + 1.3, 1e-14);
    EXPECT_NEAR(wrapper_smoothcos.get_mean(), 2.3, 1e-14);

    mio::LognormSurvivalFunction logn(0.1, 1, 0.5);
    mio::StateAgeFunctionWrapper<double> wrapper_logn(logn);
    EXPECT_EQ(wrapper_logn.get_distribution_parameter(), 0.1);
    EXPECT_EQ(wrapper_logn.get_location(), 1.0);
    EXPECT_EQ(wrapper_logn.get_scale(), 0.5);
    EXPECT_NEAR(wrapper_logn.get_support_max(dt), 2.0, 1e-14);
    // Mean value computed with python scipy-stats.
    EXPECT_NEAR(wrapper_logn.get_mean(0.1), 1.502506, 0.1);

    mio::ConstantFunction<double> constfunc(1.0);
    mio::StateAgeFunctionWrapper<double> wrapper_const(constfunc);
    EXPECT_EQ(wrapper_const.get_distribution_parameter(), 1.);
    // Deactivate temporarily log output as an error is expected here.
    mio::set_log_level(mio::LogLevel::off);

    EXPECT_NEAR(wrapper_const.get_support_max(dt), -2.0, 1e-14);
    EXPECT_NEAR(wrapper_const.get_mean(), 1., 1e-14);

    // Test if set_state_age_function works.
    mio::ErlangDensity erl(2, 0.5);
    wrapper_const.set_state_age_function(erl);
    EXPECT_EQ(wrapper_const.get_distribution_parameter(), 2);
    EXPECT_EQ(wrapper_const.get_scale(), 0.5);
    EXPECT_NEAR(wrapper_const.get_support_max(dt), 14, 1e-14);
    EXPECT_NEAR(wrapper_const.get_mean(dt), 1, 1e-14);

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestStateAgeFunction, testSAFWrapperSettersAndGetters)
{
    // Test getter and setter for Wrapper only for GammaSurvivalFunction since we have checked getter and setter
    // as well as get_support_max and get_mean for all derived classes in previous test.

    ScalarType dt = 0.5;

    mio::GammaSurvivalFunction gamma(1., 2., 3.);
    mio::StateAgeFunctionWrapper<double> wrapper(gamma);

    EXPECT_EQ(wrapper.get_distribution_parameter(), 1.);
    EXPECT_EQ(wrapper.get_location(), 2.);
    EXPECT_EQ(wrapper.get_scale(), 3.);
    EXPECT_NEAR(wrapper.get_support_max(dt), 71.5, 1e-14);
    EXPECT_NEAR(wrapper.get_mean(), 5.0, 1e-14);

    wrapper.set_distribution_parameter(4.5);
    wrapper.set_location(0.5);
    wrapper.set_scale(1.5);
    EXPECT_EQ(wrapper.get_distribution_parameter(), 4.5);
    EXPECT_EQ(wrapper.get_location(), 0.5);
    EXPECT_EQ(wrapper.get_scale(), 1.5);
    EXPECT_EQ(wrapper.get_support_max(dt), 50);
    EXPECT_NEAR(wrapper.get_mean(), 7.25, 1e-4);
}

TEST(TestStateAgeFunction, testComparisonOperator)
{
    // Check that StateAgeFunctions are only considered equal if they are of the same derived
    // class and have the same parameter.
    mio::GammaSurvivalFunction gamma(0.5, 0, 1);
    mio::GammaSurvivalFunction gamma2(0.5, 0, 1);
    mio::GammaSurvivalFunction gamma3(2, 0, 1);
    mio::GammaSurvivalFunction gamma4(0.5, 1.5, 1);
    mio::GammaSurvivalFunction gamma5(0.5, 0, 2);
    mio::LognormSurvivalFunction logn(0.5, 0, 1);

    EXPECT_TRUE(gamma == gamma2);
    EXPECT_FALSE(gamma == gamma3);
    EXPECT_FALSE(gamma == gamma4);
    EXPECT_FALSE(gamma == gamma5);
    EXPECT_FALSE(gamma == logn);

    // Check that it also holds when a StateAgeFunctionWrapper is set with the respective functions
    mio::StateAgeFunctionWrapper<double> wrapper(gamma);
    mio::StateAgeFunctionWrapper<double> wrapper2(gamma2);
    mio::StateAgeFunctionWrapper<double> wrapper3(gamma3);
    mio::StateAgeFunctionWrapper<double> wrapper4(gamma4);
    mio::StateAgeFunctionWrapper<double> wrapper5(gamma5);
    mio::StateAgeFunctionWrapper<double> wrapper6(logn);

    EXPECT_TRUE(wrapper == wrapper2);
    EXPECT_FALSE(wrapper == wrapper3);
    EXPECT_FALSE(wrapper == wrapper4);
    EXPECT_FALSE(wrapper == wrapper5);
    EXPECT_FALSE(wrapper == wrapper6);
}

TEST(TestStateAgeFunction, testGamma)
{
    // Check if the following connection between the GammaSurvivalFunction and ErlangDensity is fulfilled:
    // GammaSurvivalFunction(shape=s, rate=r)(tau)=(1/r) * \sum_{j=1}^{s} ErlangDensity(shape=j,rate=r)(tau),
    // with a natural number s and a real value r.

    ScalarType rate[]  = {0.5, 0.8, 1, 3};
    int shape[]        = {1, 3, 5, 10, 20};
    ScalarType times[] = {0, 0.5, 3, 5.555, 20, 70.34};

    for (ScalarType r : rate) {
        for (int s : shape) {
            mio::GammaSurvivalFunction survival(s, 0, 1. / r);
            // Value of the GammaSurvivalFunction at time 0 should be always 1.
            EXPECT_EQ(survival.eval(0), 1.0);
            mio::ErlangDensity density(1, 1. / r);
            // Value of ErlangDensity for negative values should be 0.
            EXPECT_EQ(density.eval(-1), 0.0);

            // Test for different times tau.
            for (ScalarType tau : times) {
                ScalarType f = 0;
                // Calculate sum of evaluations of densities as stated in the formula.
                for (int k = 1; k < s + 1; k++) {
                    density.set_distribution_parameter(k);
                    f += density.eval(tau);
                }
                // Multiply with (1/r) as stated in the formula and
                // compare with the evaluation of the GammaSurvivalFunction.
                EXPECT_NEAR(f / r, survival.eval(tau), 1e-10);
            }
        }
    }

    // Check if  GammaSurvivalFunction reduces to ExponentialSurvivalFunction with the
    // appropriate choice of parameters.
    for (ScalarType r : rate) {
        mio::GammaSurvivalFunction gamma(1, 0, 1. / r);
        mio::ExponentialSurvivalFunction<double> exp(r);
        for (ScalarType tau : times) {
            EXPECT_NEAR(gamma.eval(tau), exp.eval(tau), 1e-10);
        }
    }
}
