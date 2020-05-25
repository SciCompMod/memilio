#include "load_test_data.h"
#include <gtest/gtest.h>
#include <epidemiology/seir.h>

using real = double;

class TestCompareSeirWithJS : public testing::Test
{
protected:
    void SetUp() override
    {
        refData = load_test_data_csv<real>("data/seir-js-compare.csv");
        t0      = 0.;
        tmax    = 50.;
        dt      = 0.1002004008016032;

        params.populations.set_exposed_t0(10000);
        params.populations.set_infectious_t0(1000);
        params.populations.set_total_t0(1061000);
        params.populations.set_recovered_t0(1000);
        // suscetible now set with every other update
        // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
        params.times.set_incubation(5.2);
        params.times.set_cont_freq(2.7);
        params.times.set_infectious(2);

        // add two dampings
        params.dampings.add(epi::Damping(0., 1.0));
        params.dampings.add(epi::Damping(12.0, 0.4));
        
        params.dampings.set_smoothing(false);
    }

public:
    std::vector<std::vector<real>> refData;
    real t0;
    real tmax;
    real dt;
    epi::SeirParams params;
};

TEST_F(TestCompareSeirWithJS, integrate)
{
    std::vector<Eigen::VectorXd> result_x(0);

    auto result_t = simulate(t0, tmax, dt, params, result_x);

    ASSERT_EQ(refData.size(), result_x.size());
    ASSERT_EQ(refData.size(), result_t.size());

    for (size_t irow = 0; irow < result_x.size(); ++irow) {
        double t = refData[irow][0];
        ASSERT_NEAR(t, result_t[irow], 1e-12) << "at row " << irow;
        for (size_t icol = 0; icol < 4; ++icol) {
            double ref    = refData[irow][icol + 1];
            double actual = result_x[irow][icol];

            double tol = 1e-6 * ref;
            ASSERT_NEAR(ref, actual, tol) << "at row " << irow;
        }
    }
}

TEST(TestGroupSeir, noMigrationTheSameAsSingleSimulation)
{
    auto params = epi::SeirParams();
    params.populations.set_exposed_t0(100);
    params.populations.set_total_t0(10000);
    params.populations.set_recovered_t0(0);
    params.populations.set_infectious_t0(0);
    params.times.set_incubation(1);
    params.times.set_infectious(1);
    params.times.set_cont_freq(2.5);

    auto params_group1 = params;
    auto params_group2 = params;
    params_group1.times.set_cont_freq(0.5 * params_group2.times.get_cont_freq());

    auto t0   = 0;
    auto tmax = 10.0;
    auto dt   = 0.11;

    auto migration_f = [](auto&& arg1, auto&& arg2, auto&& matrix) {
        matrix = Eigen::MatrixXd::Identity(2, 2);
    }; //no migration

    std::vector<Eigen::VectorXd> result_groups;
    auto t_groups = epi::simulate_groups(0, tmax, dt, {params_group1, params_group2}, migration_f, result_groups);

    std::vector<Eigen::VectorXd> result_single1;
    auto t_single1 = epi::simulate(0, tmax, dt, params_group1, result_single1);

    std::vector<Eigen::VectorXd> result_single2;
    auto t_single2 = epi::simulate(0, tmax, dt, params_group2, result_single2);

    //without migration the groups should be the same when simulated together or apart
    EXPECT_DOUBLE_EQ(t_groups.back(), tmax);
    EXPECT_DOUBLE_EQ(t_groups.back(), t_single1.back());

    EXPECT_EQ(result_groups.back().size(), result_single1.back().size() + result_single2.back().size());
    for (int j = 0; j < result_single1.back().size(); j++) {
        //these only work if seir is modified to use rk45
        // EXPECT_NEAR(result_single1.back()[j], result_groups.back()[2 * j], result_groups.back()[2 * j] * 1e-5);
        // EXPECT_NEAR(result_single2.back()[j], result_groups.back()[2 * j + 1], result_groups.back()[2 * j + 1] * 1e-5);
    }
}
