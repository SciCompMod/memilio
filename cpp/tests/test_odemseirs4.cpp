/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Henrik Zunker
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/time_series.h"
#include "ode_mseirs4/model.h"
#include "ode_mseirs4/infection_state.h"
#include "ode_mseirs4/parameters.h"

#include <gtest/gtest.h>

TEST(TestOdeMseirs4, simulateDefault)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.1;

    mio::omseirs4::Model<double> model;
    mio::TimeSeries<double> result = simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

class ModelTestOdeMseirs4 : public testing::Test
{
public:
    ModelTestOdeMseirs4()
        : model()
    {
    }

    double t0{};
    double tmax{};
    double dt{};
    double total_population{};
    mio::omseirs4::Model<double> model;

protected:
    void SetUp() override
    {
        t0   = 0.0;
        tmax = 1.0;
        dt   = 0.5;

        total_population = 1'000'000.0;

        // Initialize parameters (per day) and disable births/deaths for conservation test
        auto& params                                                  = model.parameters;
        params.get<mio::omseirs4::BaseTransmissionRate<double>>()     = 0.4; // b0
        params.get<mio::omseirs4::SeasonalAmplitude<double>>()        = 0.0; // no seasonality
        params.get<mio::omseirs4::SeasonalPhase<double>>()            = 0.0; // phase
        params.get<mio::omseirs4::NaturalBirthDeathRate<double>>()    = 0.0; // mu = 0
        params.get<mio::omseirs4::LossMaternalImmunityRate<double>>() = 1.0 / 90; // xi
        params.get<mio::omseirs4::ProgressionRate<double>>()          = 1.0 / 7; // sigma
        params.get<mio::omseirs4::RecoveryRate<double>>()             = 1.0 / 14; // nu
        params.get<mio::omseirs4::ImmunityWaningRate<double>>()       = 0.0; // gamma = 0
        // beta factors use defaults

        double M  = 5'000.0;
        double E1 = 300.0, E2 = 150.0, E3 = 80.0, E4 = 70.0;
        double I1 = 200.0, I2 = 100.0, I3 = 50.0, I4 = 50.0;
        double R1 = 40'000.0, R2 = 30'000.0, R3 = 20'000.0, R4 = 10'000.0;
        double S2 = 100'000.0, S3 = 50'000.0, S4 = 50'000.0;
        double assigned = M + (E1 + E2 + E3 + E4) + (I1 + I2 + I3 + I4) + (R1 + R2 + R3 + R4) + (S2 + S3 + S4);
        double S1       = std::max(0.0, total_population - assigned);

        using IS                                                = mio::omseirs4::InfectionState;
        model.populations[{mio::Index<IS>(IS::MaternalImmune)}] = M;
        model.populations[{mio::Index<IS>(IS::E1)}]             = E1;
        model.populations[{mio::Index<IS>(IS::E2)}]             = E2;
        model.populations[{mio::Index<IS>(IS::E3)}]             = E3;
        model.populations[{mio::Index<IS>(IS::E4)}]             = E4;
        model.populations[{mio::Index<IS>(IS::I1)}]             = I1;
        model.populations[{mio::Index<IS>(IS::I2)}]             = I2;
        model.populations[{mio::Index<IS>(IS::I3)}]             = I3;
        model.populations[{mio::Index<IS>(IS::I4)}]             = I4;
        model.populations[{mio::Index<IS>(IS::R1)}]             = R1;
        model.populations[{mio::Index<IS>(IS::R2)}]             = R2;
        model.populations[{mio::Index<IS>(IS::R3)}]             = R3;
        model.populations[{mio::Index<IS>(IS::R4)}]             = R4;
        model.populations[{mio::Index<IS>(IS::S1)}]             = S1;
        model.populations[{mio::Index<IS>(IS::S2)}]             = S2;
        model.populations[{mio::Index<IS>(IS::S3)}]             = S3;
        model.populations[{mio::Index<IS>(IS::S4)}]             = S4;

        model.check_constraints();
    }
};

TEST_F(ModelTestOdeMseirs4, checkPopulationConservation)
{
    auto result           = mio::simulate<double, mio::omseirs4::Model<double>>(t0, tmax, dt, model);
    double total_pop_last = result.get_last_value().sum();
    EXPECT_NEAR(total_pop_last, total_population, 1e-7);
}

TEST(TestOdeMseirs4, apply_constraints_parameters)
{
    mio::omseirs4::Model<double> model;

    auto& params = model.parameters;

    // valid default first
    EXPECT_EQ(params.apply_constraints(), false);

    // negative values should be clamped to 0
    params.get<mio::omseirs4::BaseTransmissionRate<double>>()     = -0.5;
    params.get<mio::omseirs4::ProgressionRate<double>>()          = -1.0;
    params.get<mio::omseirs4::RecoveryRate<double>>()             = -1.0;
    params.get<mio::omseirs4::ImmunityWaningRate<double>>()       = -1.0;
    params.get<mio::omseirs4::NaturalBirthDeathRate<double>>()    = -1.0;
    params.get<mio::omseirs4::LossMaternalImmunityRate<double>>() = -1.0;
    params.get<mio::omseirs4::Beta2Factor<double>>()              = -0.1;
    params.get<mio::omseirs4::Beta3Factor<double>>()              = -0.2;
    params.get<mio::omseirs4::Beta4Factor<double>>()              = -0.3;
    params.get<mio::omseirs4::SeasonalAmplitude<double>>()        = 1.5; // will be clamped to 1

    EXPECT_EQ(params.apply_constraints(), true);

    EXPECT_DOUBLE_EQ((double)params.get<mio::omseirs4::BaseTransmissionRate<double>>(), 0.0);
    EXPECT_DOUBLE_EQ((double)params.get<mio::omseirs4::ProgressionRate<double>>(), 0.0);
    EXPECT_DOUBLE_EQ((double)params.get<mio::omseirs4::RecoveryRate<double>>(), 0.0);
    EXPECT_DOUBLE_EQ((double)params.get<mio::omseirs4::ImmunityWaningRate<double>>(), 0.0);
    EXPECT_DOUBLE_EQ((double)params.get<mio::omseirs4::NaturalBirthDeathRate<double>>(), 0.0);
    EXPECT_DOUBLE_EQ((double)params.get<mio::omseirs4::LossMaternalImmunityRate<double>>(), 0.0);
    EXPECT_DOUBLE_EQ((double)params.get<mio::omseirs4::Beta2Factor<double>>(), 0.0);
    EXPECT_DOUBLE_EQ((double)params.get<mio::omseirs4::Beta3Factor<double>>(), 0.0);
    EXPECT_DOUBLE_EQ((double)params.get<mio::omseirs4::Beta4Factor<double>>(), 0.0);
    EXPECT_DOUBLE_EQ((double)params.get<mio::omseirs4::SeasonalAmplitude<double>>(), 1.0);
}

TEST(TestOdeMseirs4, check_constraints_parameters)
{
    mio::omseirs4::Model<double> model;
    auto& params = model.parameters;

    // default is valid
    ASSERT_EQ(params.check_constraints(), false);

    // amplitude must be in [0,1]
    params.get<mio::omseirs4::SeasonalAmplitude<double>>() = -0.1;
    ASSERT_EQ(params.check_constraints(), true);

    params.get<mio::omseirs4::SeasonalAmplitude<double>>() = 0.5;
    ASSERT_EQ(params.check_constraints(), false);

    params.get<mio::omseirs4::SeasonalAmplitude<double>>() = 1.2;
    ASSERT_EQ(params.check_constraints(), true);
}

TEST(TestOdeMseirs4, population_zero_no_nan)
{
    mio::omseirs4::Model<double> model;
    model.populations.set_total(0.0);

    auto dydt = Eigen::VectorXd((Eigen::Index)mio::omseirs4::InfectionState::Count);
    dydt.setZero();
    auto y0 = model.get_initial_values();
    model.get_derivatives(y0, y0, 0.0, dydt);

    for (int i = 0; i < dydt.size(); ++i) {
        EXPECT_FALSE(std::isnan(dydt[i]));
    }
}

TEST(TestOdeMseirs4, Simulation)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 1;

    mio::omseirs4::Model<double> model;

    using IS                                    = mio::omseirs4::InfectionState;
    model.populations[{mio::Index<IS>(IS::I1)}] = 10.0;
    model.populations[{mio::Index<IS>(IS::S1)}] = 990.0;

    model.parameters.get<mio::omseirs4::BaseTransmissionRate<double>>()  = 0.2;
    model.parameters.get<mio::omseirs4::ProgressionRate<double>>()       = 1.0 / 5.0;
    model.parameters.get<mio::omseirs4::RecoveryRate<double>>()          = 1.0 / 7.0;
    model.parameters.get<mio::omseirs4::SeasonalAmplitude<double>>()     = 0.0;
    model.parameters.get<mio::omseirs4::NaturalBirthDeathRate<double>>() = 0.0;

    model.check_constraints();

    std::unique_ptr<mio::OdeIntegratorCore<double>> integrator = std::make_unique<mio::EulerIntegratorCore<double>>();
    auto sim = mio::simulate(t0, tmax, dt, model, std::move(integrator));

    EXPECT_EQ(sim.get_num_time_points(), 2);
}

TEST(TestOdeMseirs4, normalized_transitions)
{
    // case: get derivative with isolated infection terms; expect derivative match computed values
    mio::omseirs4::Model<double> model;
    auto& params                                              = model.parameters;
    params.get<mio::omseirs4::BaseTransmissionRate<double>>() = 0.4;

    // disable other flows to isolate infection terms
    params.get<mio::omseirs4::NaturalBirthDeathRate<double>>()    = 0.0;
    params.get<mio::omseirs4::LossMaternalImmunityRate<double>>() = 0.0;
    params.get<mio::omseirs4::ProgressionRate<double>>()          = 0.0;
    params.get<mio::omseirs4::RecoveryRate<double>>()             = 0.0;
    params.get<mio::omseirs4::ImmunityWaningRate<double>>()       = 0.0;
    params.get<mio::omseirs4::SeasonalAmplitude<double>>()        = 0.0;

    using IS                                    = mio::omseirs4::InfectionState;
    model.populations[{mio::Index<IS>(IS::S1)}] = 500.0;
    model.populations[{mio::Index<IS>(IS::S2)}] = 300.0;
    model.populations[{mio::Index<IS>(IS::S3)}] = 200.0;
    model.populations[{mio::Index<IS>(IS::S4)}] = 100.0;
    model.populations[{mio::Index<IS>(IS::I1)}] = 10.0;
    model.populations[{mio::Index<IS>(IS::I2)}] = 5.0;
    model.populations[{mio::Index<IS>(IS::I3)}] = 2.0;
    model.populations[{mio::Index<IS>(IS::I4)}] = 3.0;

    auto y0   = model.get_initial_values();
    auto dydt = Eigen::VectorXd((Eigen::Index)IS::Count);
    model.get_derivatives(y0, y0, 0.0, dydt);

    const double N       = 500.0 + 300.0 + 200.0 + 100.0 + 10.0 + 5.0 + 2.0 + 3.0;
    const double I_total = 10.0 + 5.0 + 2.0 + 3.0;
    const double lambda1 = 0.4 * (I_total / N);
    const double lambda2 = 0.5 * lambda1;
    const double lambda3 = 0.35 * lambda1;
    const double lambda4 = 0.25 * lambda1;
    const double tol     = 1e-12;

    EXPECT_NEAR(dydt[(Eigen::Index)IS::S1], -lambda1 * 500.0, tol);
    EXPECT_NEAR(dydt[(Eigen::Index)IS::E1], lambda1 * 500.0, tol);
    EXPECT_NEAR(dydt[(Eigen::Index)IS::S2], -lambda2 * 300.0, tol);
    EXPECT_NEAR(dydt[(Eigen::Index)IS::E2], lambda2 * 300.0, tol);
    EXPECT_NEAR(dydt[(Eigen::Index)IS::S3], -lambda3 * 200.0, tol);
    EXPECT_NEAR(dydt[(Eigen::Index)IS::E3], lambda3 * 200.0, tol);
    EXPECT_NEAR(dydt[(Eigen::Index)IS::S4], -lambda4 * 100.0, tol);
    EXPECT_NEAR(dydt[(Eigen::Index)IS::E4], lambda4 * 100.0, tol);
}
