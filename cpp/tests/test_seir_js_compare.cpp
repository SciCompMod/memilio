#include "load_test_data.h"
#include <epidemiology/secir/seir.h>
#include "epidemiology/math/euler.h"
#include "epidemiology/model/simulation.h"
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
        model.populations.set(10000, epi::SeirInfType::E);
        model.populations.set(1000, epi::SeirInfType::I);
        model.populations.set(1000, epi::SeirInfType::R);
        model.populations.set(total_population - model.populations.get(epi::SeirInfType::E) -
                                  model.populations.get(epi::SeirInfType::I) -
                                  model.populations.get(epi::SeirInfType::R),
                              epi::SeirInfType::S);
        // suscetible now set with every other update
        // model.nb_sus_t0   = model.nb_total_t0 - model.nb_exp_t0 - model.nb_inf_t0 - model.nb_rec_t0;
        model.parameters.times.set_incubation(5.2);
        model.parameters.times.set_cont_freq(2.7);
        model.parameters.times.set_infectious(2);

        // add two dampings
        model.parameters.dampings.add(epi::Damping(0., 1.0));
        model.parameters.dampings.add(epi::Damping(12.0, 0.4));

        model.parameters.dampings.set_smoothing(false);
    }

public:
    std::vector<std::vector<real>> refData;
    real t0;
    real tmax;
    real dt;
    epi::SeirModel model;
};

TEST_F(TestCompareSeirWithJS, integrate)
{
    auto result = epi::simulate<epi::SeirModel, epi::EulerIntegratorCore>(t0, tmax, dt, model);

    ASSERT_EQ(refData.size(), static_cast<size_t>(result.get_num_time_points()));

    for (Eigen::Index irow = 0; irow < result.get_num_time_points(); ++irow) {
        double t = refData[static_cast<size_t>(irow)][0];
        ASSERT_NEAR(t, result.get_times()[irow], 1e-12) << "at row " << irow;
        for (size_t icol = 0; icol < 4; ++icol) {
            double ref    = refData[static_cast<size_t>(irow)][icol + 1];
            double actual = result[irow][icol];

            double tol = 1e-6 * ref;
            ASSERT_NEAR(ref, actual, tol) << "at row " << irow;
        }
    }
}
