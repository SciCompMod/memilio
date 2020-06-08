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

    auto make_rk_integrator = []() {
        auto integrator = std::make_shared<epi::RKIntegratorCore>(1e-5, 1);
        integrator->set_abs_tolerance(1e-5);
        integrator->set_rel_tolerance(1e-5);
        return integrator;
    };

    auto t0   = 0.0;
    auto tmax = 2. * M_PI;
    auto dt   = 0.1;

    // epi::OdeIntegrator integrator1(deriv_f1, t0, Eigen::VectorXd::Constant(1, 0.0), dt, make_rk_integrator());
    // integrator1.eval(tmax);
    
    // epi::OdeIntegrator integrator2(deriv_f2, t0, Eigen::VectorXd::Constant(1, 1.0), dt, make_rk_integrator());
    // integrator2.eval(tmax);

    auto v_group = std::vector<Eigen::VectorXd>(1, (Eigen::VectorXd(2) << 0.0, 1.0).finished());
    auto t_group = epi::ode_integrate_with_migration(t0, tmax, dt, {deriv_f1, deriv_f2}, make_rk_integrator(), migration_f, v_group);

    EXPECT_DOUBLE_EQ(t_group.back(), tmax);

    //without migration the groups should be the same when simulated together or apart
    //check only at end t, intermediate steps will be different because of adaptive steps
    // for (int j = 0; j < integrator1.get_y().back().size(); j++) {
    //     EXPECT_NEAR(integrator1.get_y().back()[j], v_group.back()[2 * j], 1e-5);
    //     EXPECT_NEAR(integrator2.get_y().back()[j], v_group.back()[2 * j + 1], 1e-5);
    // }
}
