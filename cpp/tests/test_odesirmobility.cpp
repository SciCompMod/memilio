
#include "load_test_data.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "ode_sir_mobility/model.h"
#include "ode_sir_mobility/infection_state.h"
#include "ode_sir_mobility/parameters.h"
#include "memilio/math/euler.h"
#include "memilio/compartments/simulation.h"
#include <gtest/gtest.h>
#include <iomanip>
#include <vector>

TEST(TestOdeSirMobility, simulateDefault)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.1;

    size_t num_regions = 4;

    mio::osirmobility::Model model(num_regions);
    mio::TimeSeries<double> result = simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

TEST(TestOdeSirMobility, compareWithPreviousRun)
{
    // initialization
    double t0   = 0.;
    double tmax = 3.;
    double dt   = 0.1;

    size_t number_regions              = 4;
    size_t number_age_groups           = 1;
    size_t total_population_per_region = 5000;

    mio::osirmobility::Model model(number_regions, number_age_groups);

    for (size_t i = 0; i < number_regions; i++) {
        model.populations[{mio::Index<mio::osirmobility::Region, mio::AgeGroup, mio::osirmobility::InfectionState>(
            mio::osirmobility::Region(i), mio::AgeGroup(0), mio::osirmobility::InfectionState::Infected)}]  = 50;
        model.populations[{mio::Index<mio::osirmobility::Region, mio::AgeGroup, mio::osirmobility::InfectionState>(
            mio::osirmobility::Region(i), mio::AgeGroup(0), mio::osirmobility::InfectionState::Recovered)}] = 0;
        model.populations[{mio::Index<mio::osirmobility::Region, mio::AgeGroup, mio::osirmobility::InfectionState>(
            mio::osirmobility::Region(i), mio::AgeGroup(0), mio::osirmobility::InfectionState::Susceptible)}] =
            total_population_per_region -
            model.populations[{mio::Index<mio::osirmobility::Region, mio::AgeGroup, mio::osirmobility::InfectionState>(
                mio::osirmobility::Region(i), mio::AgeGroup(0), mio::osirmobility::InfectionState::Infected)}] -
            model.populations[{mio::Index<mio::osirmobility::Region, mio::AgeGroup, mio::osirmobility::InfectionState>(
                mio::osirmobility::Region(i), mio::AgeGroup(0), mio::osirmobility::InfectionState::Recovered)}];
    }
    model.parameters.set<mio::osirmobility::TransmissionProbabilityOnContact>(1.0);
    model.parameters.set<mio::osirmobility::TimeInfected>(2);

    model.parameters.get<mio::osirmobility::ContactPatterns>().get_baseline()(0, 0) = 2.7;
    model.parameters.get<mio::osirmobility::ContactPatterns>().add_damping(0.6, mio::SimulationTime(12.5));

    model.parameters.set<mio::osirmobility::ImpactTransmissionDuringCommuting>(1.);
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(0), 0.2});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(2), 0.6});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(2), mio::osirmobility::Region(0), 0.5});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(0), mio::osirmobility::Region(3), 1.0});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(3), 0.2});

    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(0),
                                                                  mio::osirmobility::Region(1)}] = {2};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(0),
                                                                  mio::osirmobility::Region(3)}] = {2};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(1),
                                                                  mio::osirmobility::Region(0)}] = {2};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(1),
                                                                  mio::osirmobility::Region(2)}] = {0};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(1),
                                                                  mio::osirmobility::Region(3)}] = {2};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(2),
                                                                  mio::osirmobility::Region(1)}] = {0};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(3),
                                                                  mio::osirmobility::Region(0)}] = {2};
    model.parameters.get<mio::osirmobility::PathIntersections>()[{mio::osirmobility::Region(3),
                                                                  mio::osirmobility::Region(1)}] = {2};

    std::vector<std::vector<double>> refData = load_test_data_csv<double>("ode-sir-mobility-compare.csv");
    auto integrator                          = std::make_shared<mio::EulerIntegratorCore>();
    auto result                              = mio::simulate<mio::osirmobility::Model>(t0, tmax, dt, model, integrator);

    ASSERT_EQ(refData.size(), static_cast<size_t>(result.get_num_time_points()));

    for (Eigen::Index irow = 0; irow < result.get_num_time_points(); ++irow) {
        double t     = refData[static_cast<size_t>(irow)][0];
        auto rel_tol = 1e-6;

        //test result diverges at damping because of changes, not worth fixing at the moment
        if (t > 11.0 && t < 13.0) {
            //strong divergence around damping
            rel_tol = 0.5;
        }
        else if (t > 13.0) {
            //minor divergence after damping
            rel_tol = 1e-2;
        }
        mio::unused(rel_tol);

        ASSERT_NEAR(t, result.get_times()[irow], 1e-12) << "at row " << irow;

        for (size_t icol = 0; icol < 12; ++icol) {
            double ref    = refData[static_cast<size_t>(irow)][icol + 1];
            double actual = result[irow][icol];

            double tol = rel_tol * ref;
            ASSERT_NEAR(ref, actual, tol) << "at row " << irow;
        }
    }
}

TEST(TestOdeSirMobility, checkPopulationConservation)
{
    // initialization
    double t0   = 0.;
    double tmax = 50.;
    double dt   = 0.1002004008016032;

    double population_per_region          = 1061000;
    mio::osirmobility::Region num_regions = 4;

    mio::osirmobility::Model model((size_t)num_regions);

    for (auto region = mio::osirmobility::Region(0); region < num_regions; ++region) {
        model.populations[{region, mio::AgeGroup(0), mio::osirmobility::InfectionState::Infected}]  = 1000;
        model.populations[{region, mio::AgeGroup(0), mio::osirmobility::InfectionState::Recovered}] = 1000;
        model.populations[{region, mio::AgeGroup(0), mio::osirmobility::InfectionState::Susceptible}] =
            population_per_region -
            model.populations[{region, mio::AgeGroup(0), mio::osirmobility::InfectionState::Infected}] -
            model.populations[{region, mio::AgeGroup(0), mio::osirmobility::InfectionState::Recovered}];
    }
    model.parameters.set<mio::osirmobility::TimeInfected>(2);
    model.parameters.set<mio::osirmobility::TransmissionProbabilityOnContact>(0.04);
    model.parameters.set<mio::osirmobility::ImpactTransmissionDuringCommuting>(1.);
    model.parameters.get<mio::osirmobility::ContactPatterns>().get_baseline()(0, 0) = 1.;
    model.parameters.get<mio::osirmobility::ContactPatterns>().add_damping(0.6, mio::SimulationTime(12.5));
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(0), 0.5});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(2), 0.8});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(2), mio::osirmobility::Region(0), 0.5});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(0), mio::osirmobility::Region(3), 1.0});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(3), 0.8});

    auto result        = mio::simulate<mio::osirmobility::Model>(t0, tmax, dt, model);
    double num_persons = 0.0;
    for (auto i = 0; i < result.get_last_value().size(); ++i) {
        num_persons += result.get_last_value()[i];
    }
    EXPECT_NEAR(num_persons, population_per_region * (size_t)num_regions, 1e-8);
}

TEST(TestOdeSirMobility, check_constraints_parameters)
{
    mio::osirmobility::Region num_regions = 2;

    mio::osirmobility::Model model((size_t)num_regions);
    model.parameters.set<mio::osirmobility::TimeInfected>(6);
    model.parameters.set<mio::osirmobility::TransmissionProbabilityOnContact>(0.04);
    model.parameters.set<mio::osirmobility::ImpactTransmissionDuringCommuting>(1.);
    model.parameters.get<mio::osirmobility::ContactPatterns>().get_baseline()(0, 0) = 10.;
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(0), 0.5});

    // model.check_constraints() combines the functions from population and parameters.
    // We only want to test the functions for the parameters defined in parameters.h
    ASSERT_EQ(model.parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    model.parameters.set<mio::osirmobility::TimeInfected>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osirmobility::TimeInfected>(6);
    model.parameters.set<mio::osirmobility::TransmissionProbabilityOnContact>(10.);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osirmobility::TransmissionProbabilityOnContact>(0.04);
    model.parameters.set<mio::osirmobility::ImpactTransmissionDuringCommuting>(10.);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::osirmobility::ImpactTransmissionDuringCommuting>(1.);
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(0), 10.5});
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.get<mio::osirmobility::CommutingRatio>().pop_back();
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(2), mio::osirmobility::Region(0), 0.5});
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSirMobility, apply_constraints_parameters)
{
    const double tol_times                = 1e-1;
    mio::osirmobility::Region num_regions = 2;

    mio::osirmobility::Model model((size_t)num_regions);
    model.parameters.set<mio::osirmobility::TimeInfected>(6);
    model.parameters.set<mio::osirmobility::TransmissionProbabilityOnContact>(0.04);
    model.parameters.set<mio::osirmobility::ImpactTransmissionDuringCommuting>(1.);
    model.parameters.get<mio::osirmobility::ContactPatterns>().get_baseline()(0, 0) = 10.;
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(0), 0.5});

    EXPECT_EQ(model.parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    model.parameters.set<mio::osirmobility::TimeInfected>(-2.5);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::osirmobility::TimeInfected>(), tol_times);

    model.parameters.set<mio::osirmobility::TransmissionProbabilityOnContact>(10.);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR(model.parameters.get<mio::osirmobility::TransmissionProbabilityOnContact>(), 0.0, 1e-14);

    model.parameters.set<mio::osirmobility::TransmissionProbabilityOnContact>(0.04);
    model.parameters.set<mio::osirmobility::ImpactTransmissionDuringCommuting>(10.);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR(model.parameters.get<mio::osirmobility::ImpactTransmissionDuringCommuting>(), 0.0, 1e-14);

    model.parameters.set<mio::osirmobility::ImpactTransmissionDuringCommuting>(1.);
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(0), 10.5});
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR(std::get<double>(model.parameters.get<mio::osirmobility::CommutingRatio>()[2]), 0.0, 1e-14);

    model.parameters.get<mio::osirmobility::CommutingRatio>().pop_back();
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(2), mio::osirmobility::Region(0), 0.5});
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    // EXPECT_EQ(model.parameters.get<mio::osirmobility::CommutingRatio>().size(), 2); // 1 by default + 1 added
    mio::set_log_level(mio::LogLevel::warn);
}