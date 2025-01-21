/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/utils/parameter_distributions.h"
#include <distributions_helpers.h>
#include <matchers.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

void check_distribution(const mio::ParameterDistribution& dist, const mio::ParameterDistribution& dist_read)
{

    struct CheckDistEqVisitor : public mio::ConstParameterDistributionVisitor {
        CheckDistEqVisitor(const mio::ParameterDistribution& other_distribution)
            : other(other_distribution)
        {
        }

        void visit(const mio::ParameterDistributionNormal& self) override
        {
            auto p_other_normal_distribution = dynamic_cast<const mio::ParameterDistributionNormal*>(&other);
            ASSERT_TRUE(p_other_normal_distribution != nullptr);

            EXPECT_THAT(self.get_mean(), FloatingPointEqual(p_other_normal_distribution->get_mean(), 1e-12, 1e-12));
            EXPECT_THAT(self.get_standard_dev(),
                        FloatingPointEqual(p_other_normal_distribution->get_standard_dev(), 1e-12, 1e-12));

            EXPECT_EQ(self.get_predefined_samples().size(),
                      p_other_normal_distribution->get_predefined_samples().size());
            for (size_t i = 0; i < self.get_predefined_samples().size(); i++) {
                EXPECT_THAT(self.get_predefined_samples()[i],
                            FloatingPointEqual(p_other_normal_distribution->get_predefined_samples()[i], 1e-12, 1e-12));
            }
        }
        void visit(const mio::ParameterDistributionUniform& self) override
        {
            auto p_other_uniform_distribution = dynamic_cast<const mio::ParameterDistributionUniform*>(&other);
            ASSERT_TRUE(p_other_uniform_distribution != nullptr);

            EXPECT_EQ(self.get_predefined_samples().size(),
                      p_other_uniform_distribution->get_predefined_samples().size());
            for (size_t i = 0; i < self.get_predefined_samples().size(); i++) {
                EXPECT_THAT(
                    self.get_predefined_samples()[i],
                    FloatingPointEqual(p_other_uniform_distribution->get_predefined_samples()[i], 1e-12, 1e-12));
            }
        }
        void visit(const mio::ParameterDistributionLogNormal& self) override
        {
            auto p_other_lognormal_distribution = dynamic_cast<const mio::ParameterDistributionLogNormal*>(&other);
            ASSERT_TRUE(p_other_lognormal_distribution != nullptr);

            EXPECT_EQ(self.get_predefined_samples().size(),
                      p_other_lognormal_distribution->get_predefined_samples().size());
            for (size_t i = 0; i < self.get_predefined_samples().size(); i++) {
                EXPECT_THAT(
                    self.get_predefined_samples()[i],
                    FloatingPointEqual(p_other_lognormal_distribution->get_predefined_samples()[i], 1e-12, 1e-12));
            }
        }
        void visit(const mio::ParameterDistributionExponential& self) override
        {
            auto p_other_exponential_distribution = dynamic_cast<const mio::ParameterDistributionExponential*>(&other);
            ASSERT_TRUE(p_other_exponential_distribution != nullptr);

            EXPECT_EQ(self.get_predefined_samples().size(),
                      p_other_exponential_distribution->get_predefined_samples().size());
            for (size_t i = 0; i < self.get_predefined_samples().size(); i++) {
                EXPECT_THAT(
                    self.get_predefined_samples()[i],
                    FloatingPointEqual(p_other_exponential_distribution->get_predefined_samples()[i], 1e-12, 1e-12));
            }
        }
        void visit(const mio::ParameterDistributionConstant& self) override
        {
            auto p_other_constant_distribution = dynamic_cast<const mio::ParameterDistributionConstant*>(&other);
            ASSERT_TRUE(p_other_constant_distribution != nullptr);

            EXPECT_EQ(self.get_predefined_samples().size(),
                      p_other_constant_distribution->get_predefined_samples().size());
            for (size_t i = 0; i < self.get_predefined_samples().size(); i++) {
                EXPECT_THAT(
                    self.get_predefined_samples()[i],
                    FloatingPointEqual(p_other_constant_distribution->get_predefined_samples()[i], 1e-12, 1e-12));
            }
        }
        const mio::ParameterDistribution& other;
    };

    CheckDistEqVisitor visitor(dist_read);
    dist.accept(visitor);
}
