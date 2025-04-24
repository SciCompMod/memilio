
#include "load_test_data.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "ode_metapop/model.h"
#include "ode_metapop/infection_state.h"
#include "ode_metapop/parameters.h"
#include "memilio/math/euler.h"
#include "memilio/compartments/simulation.h"
#include <gtest/gtest.h>
#include <vector>
#include <Eigen/Dense>

class ModelTestOdeMetapop : public testing::Test
{
public:
    ModelTestOdeMetapop()
        : model(1, 1)
    {
    }
    double t0;
    double tmax;
    double dt;
    double total_population;
    mio::oseirmetapop::Model<double> model;

protected:
    void SetUp() override
    {
        t0   = 0.;
        tmax = 50.;
        dt   = 0.1;

        total_population = 1061000;

        model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
            mio::oseirmetapop::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Exposed)}]   = 10000;
        model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
            mio::oseirmetapop::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Infected)}]  = 1000;
        model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
            mio::oseirmetapop::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Recovered)}] = 1000;
        model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
            mio::oseirmetapop::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Susceptible)}] =
            total_population -
            model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
                mio::oseirmetapop::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Exposed)}] -
            model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
                mio::oseirmetapop::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Infected)}] -
            model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
                mio::oseirmetapop::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Recovered)}];

        model.set_commuting_strengths();
    }
};

TEST_F(ModelTestOdeMetapop, simulateDefault)
{
    mio::TimeSeries<double> result = simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

TEST_F(ModelTestOdeMetapop, checkPopulationConservation)
{
    auto result        = simulate(t0, tmax, dt, model);
    double num_persons = result.get_last_value().sum();
    EXPECT_NEAR(num_persons, total_population, 1e-8);
}

TEST_F(ModelTestOdeMetapop, check_constraints_parameters)
{
    model.parameters.set<mio::oseirmetapop::TimeExposed<double>>(5.2);
    model.parameters.set<mio::oseirmetapop::TimeInfected<double>>(6);
    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<double>>(0.04);

    model.parameters.get<mio::oseirmetapop::ContactPatterns<double>>().get_cont_freq_mat()[0].get_baseline()(0, 0) = 10;

    // model.check_constraints() combines the functions from population and parameters.
    // We only want to test the functions for the parameters defined in parameters.h
    ASSERT_EQ(model.parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::oseirmetapop::TimeExposed<double>>(-5.2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::oseirmetapop::TimeExposed<double>>(5.2);
    model.parameters.set<mio::oseirmetapop::TimeInfected<double>>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::oseirmetapop::TimeInfected<double>>(6);
    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<double>>(10.);
    ASSERT_EQ(model.parameters.check_constraints(), 1);
    mio::set_log_level(mio::LogLevel::warn);
}

TEST_F(ModelTestOdeMetapop, apply_constraints_parameters)
{
    const double tol_times = 1e-1;
    model.parameters.set<mio::oseirmetapop::TimeExposed<double>>(5.2);
    model.parameters.set<mio::oseirmetapop::TimeInfected<double>>(2);
    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<double>>(0.04);
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseirmetapop::ContactPatterns<double>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(10);

    EXPECT_EQ(model.parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::oseirmetapop::TimeExposed<double>>(-5.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseirmetapop::TimeExposed<double>>()[(mio::AgeGroup)0], tol_times);

    model.parameters.set<mio::oseirmetapop::TimeInfected<double>>(1e-5);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseirmetapop::TimeInfected<double>>()[(mio::AgeGroup)0], tol_times);

    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<double>>(10.);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR(model.parameters.get<mio::oseirmetapop::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0],
                0.0, 1e-14);
    mio::set_log_level(mio::LogLevel::warn);
}

// TEST_F(ModelTestOdeMetapop, compareSEIR)
// {
//     model.parameters.set<mio::oseirmetapop::TimeExposed<double>>(5.2);
//     model.parameters.set<mio::oseirmetapop::TimeInfected<double>>(2);
//     model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<double>>(1.);
//     model.parameters.get<mio::oseirmetapop::ContactPatterns<double>>().get_cont_freq_mat()[0].get_baseline()(0, 0) =
//         2.7;
//     model.parameters.get<mio::oseirmetapop::ContactPatterns<double>>().get_cont_freq_mat()[0].add_damping(
//         0.6, mio::SimulationTime(12.5));

//     model.check_constraints();

//     auto result                              = simulate(t0, tmax, dt, model);
//     std::vector<std::vector<double>> refData = load_test_data_csv<double>("seir-compare.csv");

//     ASSERT_EQ(refData.size(), static_cast<size_t>(result.get_num_time_points()));

//     for (Eigen::Index irow = 0; irow < result.get_num_time_points(); ++irow) {
//         double t     = refData[static_cast<size_t>(irow)][0];
//         auto rel_tol = 1e-6;

//         //test result diverges at damping because of changes, not worth fixing at the moment
//         if (t > 11.0 && t < 13.0) {
//             //strong divergence around damping
//             rel_tol = 0.5;
//         }
//         else if (t > 13.0) {
//             //minor divergence after damping
//             rel_tol = 1e-2;
//         }
//         mio::unused(rel_tol);

//         ASSERT_NEAR(t, result.get_times()[irow], 1e-12) << "at row " << irow;

//         for (size_t icol = 0; icol < refData[static_cast<size_t>(irow)].size() - 1; ++icol) {
//             double ref    = refData[static_cast<size_t>(irow)][icol + 1];
//             double actual = result[irow][icol];

//             double tol = rel_tol * ref;
//             ASSERT_NEAR(ref, actual, tol) << "at row " << irow;
//         }
//     }
// }

// TEST_F(ModelTestOdeMetapop, compareWithPreviousRun)
// {
//     // initialization
//     double t0   = 0.;
//     double tmax = 3.;
//     double dt   = 0.1;

//     size_t number_regions              = 4;
//     size_t number_age_groups           = 2;
//     size_t total_population_per_region = 5000;

//     mio::oseirmetapop::Model model(number_regions, number_age_groups);

//     for (size_t age = 0; age < number_age_groups; age++) {
//         for (size_t i = 0; i < number_regions; i++) {
//             model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
//                 mio::oseirmetapop::Region(i), mio::AgeGroup(age), mio::oseirmetapop::InfectionState::Infected)}]  = 50;
//             model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
//                 mio::oseirmetapop::Region(i), mio::AgeGroup(age), mio::oseirmetapop::InfectionState::Recovered)}] = 0;
//             model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
//                 mio::oseirmetapop::Region(i), mio::AgeGroup(age), mio::oseirmetapop::InfectionState::Susceptible)}] =
//                 total_population_per_region -
//                 model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup,
//                                               mio::oseirmetapop::InfectionState>(
//                     mio::oseirmetapop::Region(i), mio::AgeGroup(age), mio::oseirmetapop::InfectionState::Infected)}] -
//                 model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup,
//                                               mio::oseirmetapop::InfectionState>(
//                     mio::oseirmetapop::Region(i), mio::AgeGroup(age), mio::oseirmetapop::InfectionState::Recovered)}];
//         }
//     }
//     model.parameters.template set<mio::oseirmetapop::TransmissionProbabilityOnContact<>>(1.0);
//     model.parameters.set<mio::oseirmetapop::TimeInfected<>>(2);

//     model.parameters.get<mio::oseirmetapop::ContactPatterns<>>().get_cont_freq_mat()[0].get_baseline().setConstant(2.7);
//     model.parameters.get<mio::oseirmetapop::ContactPatterns<>>().get_cont_freq_mat()[0].add_damping(
//         0.6, mio::SimulationTime(12.5));

//     Eigen::MatrixXd mobility_data_commuter(4, 4);
//     mobility_data_commuter << 0., 0., 0., 1., 0.2, 0., 0.6, 0.2, 0.5, 0., 0.5, 0., 0., 0., 0., 1.;
//     model.parameters.template get<mio::oseirmetapop::CommutingStrengths<>>().get_cont_freq_mat()[0].get_baseline() =
//         mobility_data_commuter;

//     std::vector<std::vector<double>> refData = load_test_data_csv<double>("ode-seir-metapop-compare.csv");
//     std::shared_ptr<mio::IntegratorCore<double>> integrator = std::make_shared<mio::EulerIntegratorCore<double>>();
//     auto result = mio::simulate<double, mio::oseirmetapop::Model<>>(t0, tmax, dt, model, integrator);

//     ASSERT_EQ(refData.size(), static_cast<size_t>(result.get_num_time_points()));

//     for (Eigen::Index irow = 0; irow < result.get_num_time_points(); ++irow) {
//         double t     = refData[static_cast<size_t>(irow)][0];
//         auto rel_tol = 1e-6;

//         //test result diverges at damping because of changes, not worth fixing at the moment
//         if (t > 11.0 && t < 13.0) {
//             //strong divergence around damping
//             rel_tol = 0.5;
//         }
//         else if (t > 13.0) {
//             //minor divergence after damping
//             rel_tol = 1e-2;
//         }
//         mio::unused(rel_tol);

//         ASSERT_NEAR(t, result.get_times()[irow], 1e-12) << "at row " << irow;

//         for (size_t icol = 0; icol < 12; ++icol) {
//             double ref    = refData[static_cast<size_t>(irow)][icol + 1];
//             double actual = result[irow][icol];

//             double tol = rel_tol * ref;
//             ASSERT_NEAR(ref, actual, tol) << "at row " << irow;
//         }
//     }
// }
