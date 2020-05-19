#include "epidemiology/migration.h"
#include "epidemiology/integrator.h"
#include "epidemiology/adapt_rk.h"
#include "gtest/gtest.h"

TEST(TestMigration, compareWithSingleIntegration)
{
    auto deriv_f = [] (auto&& y, auto&& t, auto&& ydt) { ydt[0] = cos(t); };
    auto migration_f = [] (auto&& idx, auto t, auto&& matrix) { matrix = Eigen::MatrixXd::Identity(1,1); };
    
    auto integrators = std::vector<std::unique_ptr<epi::IntegratorBase>>{};
    auto integrator = [&integrators, &deriv_f](){    
        auto tmp = std::make_unique<epi::RKIntegrator>(deriv_f, 1e-5, 1);
        auto p = tmp.get();
        integrators.push_back(std::move(tmp));
        return p;
    }();
    integrator->set_abs_tolerance(0);
    integrator->set_rel_tolerance(1e-5);

    auto t0 = 0.0;
    auto tmax = 10.0;
    auto dt = 0.1;

    auto v1 = std::vector<Eigen::VectorXd>(1, Eigen::VectorXd::Constant(1, 0.0));
    auto t1 = epi::ode_integrate(t0, tmax, dt, *integrator, v1);

    auto v2 = std::vector<Eigen::VectorXd>(1, Eigen::VectorXd::Constant(1, 0.0));
    auto t2 = epi::ode_integrate_with_migration(t0, tmax, dt, integrators, migration_f, v2);

    //without migration the groups should be the same when simulated together or apart
    EXPECT_FLOAT_EQ(t1.back(), tmax);
    EXPECT_FLOAT_EQ(t2.back(), tmax);

    //check end result
    for (size_t j = 0; j < v1.back().size(); j++)
    {
        EXPECT_NEAR(v1.back()[j], v2.back()[j], 1e-5);
    }
}