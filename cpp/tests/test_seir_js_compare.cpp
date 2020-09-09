#include "load_test_data.h"
#include <epidemiology/secir/seir.h>
#include <gtest/gtest.h>

using real = double;

class TestCompareSeirWithJS : public testing::Test
{
protected:
    void SetUp() override
    {
        refData = load_test_data_csv<real>("seir-js-compare.csv");
        t0      = 0.;
        tmax    = 50.;
        dt      = 0.1002004008016032;

        double total_population = 1061000;
        params.populations.set({epi::SeirCompartments::E}, 10000);
        params.populations.set({epi::SeirCompartments::I}, 1000);
        params.populations.set({epi::SeirCompartments::R}, 1000);
        params.populations.set({epi::SeirCompartments::S}, total_population -
                                                               params.populations.get({epi::SeirCompartments::E}) -
                                                               params.populations.get({epi::SeirCompartments::I}) -
                                                               params.populations.get({epi::SeirCompartments::R}));
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
