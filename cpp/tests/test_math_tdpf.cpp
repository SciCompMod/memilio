#include "memilio/math/time_dependent_parameter_functor.h"
#include "memilio/utils/random_number_generator.h"

#include <gtest/gtest.h>

#include <vector>

class TestMathTdpf : public ::testing::Test
{
public:
    const int num_evals = 1000;

    double fuzzy_val(double min, double max)
    {
        return mio::UniformDistribution<double>::get_instance()(m_rng, min, max);
    }

protected:
    void SetUp() override
    {
        log_rng_seeds(m_rng, mio::LogLevel::warn);
    }

private:
    mio::RandomNumberGenerator m_rng{};
};

TEST_F(TestMathTdpf, zero)
{
    // Test that the Zero-TDPF always returns zero, using a random evaluation point.

    // initialize
    mio::TimeDependentParameterFunctor tdpf;

    // verify output
    for (int i = 0; i < this->num_evals; i++) {
        auto random_t_eval = this->fuzzy_val(-1e+3, 1e+3);
        EXPECT_EQ(tdpf(random_t_eval), 0.0);
    }
}

TEST_F(TestMathTdpf, linearInterpolation)
{
    // Test that the LinearInterpolation-TDPF correctly reproduces a (piecewise) linear function, using random samples.
    // Since the initialization uses unsorted data, this also checks that the data gets sorted

    const double min = -1e+3, max = 1e+3; // reasonably large values for lin_fct height and slopes
    const double t_min = -1, t_max = 1, t_mid = this->fuzzy_val(t_min, t_max);
    const double slope1 = this->fuzzy_val(min, max), slope2 = this->fuzzy_val(min, max),
                 height = this->fuzzy_val(min, max);

    const auto pcw_lin_fct = [&](double t) {
        // continuous function with different slopes between t_min, t_mid and t_max, constant otherwise
        return height + slope1 * std::clamp(t - t_min, 0.0, t_mid - t_min) +
               slope2 * std::clamp(t - t_mid, 0.0, t_max - t_mid);
    };

    // initialize the data with the critical points
    std::vector<std::vector<double>> unsorted_data{
        {t_max, pcw_lin_fct(t_max)}, {t_min, pcw_lin_fct(t_min)}, {t_mid, pcw_lin_fct(t_mid)}};
    // randomly add a few more evaluations in between
    for (int i = 0; i < 10; i++) {
        const double t = this->fuzzy_val(-1.0, 1.0);
        unsorted_data.push_back({t, pcw_lin_fct(t)});
    }

    mio::TimeDependentParameterFunctor tdpf(mio::TimeDependentParameterFunctor::Type::LinearInterpolation,
                                            unsorted_data);

    // verify output
    for (int i = 0; i < this->num_evals; i++) {
        // sample in the interval [t_min - (t_max - t_min) / 4, t_max + (t_max - tmin) / 4]
        double random_t_eval = this->fuzzy_val(1.25 * t_min - 0.25 * t_max, 1.25 * t_max - 0.25 * t_min);
        EXPECT_NEAR(tdpf(random_t_eval), pcw_lin_fct(random_t_eval), 1e-10);
    }
}
