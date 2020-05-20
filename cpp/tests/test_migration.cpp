#define _USE_MATH_DEFINES

#include "epidemiology/migration.h"
#include "epidemiology/integrator.h"
#include "epidemiology/adapt_rk.h"
#include "gtest/gtest.h"

#include <cmath>

TEST(TestMigration, compareWithSingleIntegration)
{
    auto deriv_f1 = [](auto&& y, auto&& t, auto&& ydt) {
        ydt[0] = cos(t);
    };
    auto deriv_f2 = [](auto&& y, auto&& t, auto&& ydt) {
        ydt[0] = -sin(t);
    };
    auto migration_f = [](auto&& idx, auto&& t, auto&& matrix) {
        matrix = Eigen::MatrixXd::Identity(2, 2);
    };

    auto make_rk_integrator = [](auto&& f) {
        auto integrator = std::make_unique<epi::RKIntegrator>(f, 1e-5, 1);
        integrator->set_abs_tolerance(1e-5);
        integrator->set_rel_tolerance(1e-5);
        return std::move(integrator);
    };
    //two step init because vector constructors and initializer lists don't like move-only types...
    auto integrators = std::vector<std::unique_ptr<epi::IntegratorBase>>();
    integrators.emplace_back(make_rk_integrator(deriv_f1));
    integrators.emplace_back(make_rk_integrator(deriv_f2));

    auto t0   = 0.0;
    auto tmax = 2. * M_PI;
    auto dt   = 0.1;

    auto v_single1 = std::vector<Eigen::VectorXd>(1, Eigen::VectorXd::Constant(1, 0.0));
    auto t_single1 = epi::ode_integrate(t0, tmax, dt, *integrators[0], v_single1);

    auto v_single2 = std::vector<Eigen::VectorXd>(1, Eigen::VectorXd::Constant(1, 1.0));
    auto t_single2 = epi::ode_integrate(t0, tmax, dt, *integrators[1], v_single2);

    auto v_group = std::vector<Eigen::VectorXd>(1, (Eigen::VectorXd(2) << 0.0, 1.0).finished());
    auto t_group = epi::ode_integrate_with_migration(t0, tmax, dt, integrators, migration_f, v_group);

    EXPECT_DOUBLE_EQ(t_group.back(), tmax);

    //without migration the groups should be the same when simulated together or apart
    //check only at end t, intermediate steps will be different because of adaptive steps
    for (int j = 0; j < v_single1.back().size(); j++) {
        EXPECT_NEAR(v_single1.back()[j], v_group.back()[2 * j], 1e-5);
        EXPECT_NEAR(v_single2.back()[j], v_group.back()[2 * j + 1], 1e-5);
    }
}
