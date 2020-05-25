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
        ASSERT_FLOAT_EQ(t, result_t[irow]) << "at row " << irow;
        for (size_t icol = 0; icol < 4; ++icol) {
            double ref    = refData[irow][icol + 1];
            double actual = result_x[irow][icol];

            double tol = 1e-6 * ref;
            ASSERT_NEAR(ref, actual, tol)  << "at row " << irow;
        }
    }
}
