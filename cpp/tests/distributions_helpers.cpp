#include <distributions_helpers.h>
#include <matchers.h>
#include <epidemiology/secir/parameter_studies.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

void check_distribution(const epi::ParameterDistribution& dist, const epi::ParameterDistribution& dist_read)
{

    struct CheckDistEqVisitor : public epi::ConstParameterDistributionVisitor {
        CheckDistEqVisitor(const epi::ParameterDistribution& other_distribution)
            : other(other_distribution)
        {
        }

        void visit(const epi::ParameterDistributionNormal& self) override
        {
            auto p_other_normal_distribution = dynamic_cast<const epi::ParameterDistributionNormal*>(&other);
            ASSERT_TRUE(p_other_normal_distribution != nullptr);

            EXPECT_THAT(self.get_mean(), FloatingPointEqual(p_other_normal_distribution->get_mean(), 1e-12, 1e-12));
            EXPECT_THAT(self.get_standard_dev(),
                        FloatingPointEqual(p_other_normal_distribution->get_standard_dev(), 1e-12, 1e-12));
            EXPECT_THAT(self.get_lower_bound(),
                        FloatingPointEqual(p_other_normal_distribution->get_lower_bound(), 1e-12, 1e-12));
            EXPECT_THAT(self.get_upper_bound(),
                        FloatingPointEqual(p_other_normal_distribution->get_upper_bound(), 1e-12, 1e-12));

            EXPECT_EQ(self.get_predefined_samples().size(),
                      p_other_normal_distribution->get_predefined_samples().size());
            for (size_t i = 0; i < self.get_predefined_samples().size(); i++) {
                EXPECT_THAT(self.get_predefined_samples()[i],
                            FloatingPointEqual(p_other_normal_distribution->get_predefined_samples()[i], 1e-12, 1e-12));
            }
        }
        void visit(const epi::ParameterDistributionUniform& self) override
        {
            auto p_other_uniform_distribution = dynamic_cast<const epi::ParameterDistributionUniform*>(&other);
            ASSERT_TRUE(p_other_uniform_distribution != nullptr);

            EXPECT_THAT(self.get_lower_bound(),
                        FloatingPointEqual(p_other_uniform_distribution->get_lower_bound(), 1e-12, 1e-12));
            EXPECT_THAT(self.get_upper_bound(),
                        FloatingPointEqual(p_other_uniform_distribution->get_upper_bound(), 1e-12, 1e-12));

            EXPECT_EQ(self.get_predefined_samples().size(),
                      p_other_uniform_distribution->get_predefined_samples().size());
            for (size_t i = 0; i < self.get_predefined_samples().size(); i++) {
                EXPECT_THAT(self.get_predefined_samples()[i],
                            FloatingPointEqual(p_other_uniform_distribution->get_predefined_samples()[i], 1e-12, 1e-12));
            }
        }
        const epi::ParameterDistribution& other;
    };

    CheckDistEqVisitor visitor(dist_read);
    dist.accept(visitor);
}