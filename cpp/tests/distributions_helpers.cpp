#include <distributions_helpers.h>
#include <epidemiology/parameter_studies/parameter_studies.h>
#include <gtest/gtest.h>

void check_dist(const epi::ParameterDistribution& dist, const epi::ParameterDistribution& dist_read)
{

    struct CheckDistEqVisitor : public epi::ConstParameterDistributionVisitor {
        CheckDistEqVisitor(const epi::ParameterDistribution& other_dist)
            : other(other_dist)
        {
        }

        void visit(const epi::ParameterDistributionNormal& self) override
        {
            auto p_other_normal_dist = dynamic_cast<const epi::ParameterDistributionNormal*>(&other);
            ASSERT_TRUE(p_other_normal_dist != nullptr);

            EXPECT_EQ(self.get_mean(), p_other_normal_dist->get_mean());
            EXPECT_EQ(self.get_standard_dev(), p_other_normal_dist->get_standard_dev());
            EXPECT_EQ(self.get_lower_bound(), p_other_normal_dist->get_lower_bound());
            EXPECT_EQ(self.get_upper_bound(), p_other_normal_dist->get_upper_bound());
        }
        void visit(const epi::ParameterDistributionUniform& self) override
        {
            auto p_other_uniform_dist = dynamic_cast<const epi::ParameterDistributionUniform*>(&other);
            ASSERT_TRUE(p_other_uniform_dist != nullptr);

            EXPECT_EQ(self.get_lower_bound(), p_other_uniform_dist->get_lower_bound());
            EXPECT_EQ(self.get_upper_bound(), p_other_uniform_dist->get_upper_bound());
        }
        const epi::ParameterDistribution& other;
    };

    CheckDistEqVisitor visitor(dist_read);
    dist.accept(visitor);
}