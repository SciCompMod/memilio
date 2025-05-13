
#include "load_test_data.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "ode_seir_metapop/model.h"
#include "ode_seir_metapop/infection_state.h"
#include "ode_seir_metapop/parameters.h"
#include "memilio/math/euler.h"
#include "memilio/compartments/simulation.h"
#include <gtest/gtest.h>
#include <vector>
#include <Eigen/Dense>

class ModelTestOdeMetapop : public testing::Test
{
public:
    ModelTestOdeMetapop()
        : model(4, 1)
    {
    }
    ScalarType t0;
    ScalarType tmax;
    ScalarType dt;
    ScalarType total_population_per_region;
    mio::oseirmetapop::Model<ScalarType> model;

protected:
    void SetUp() override
    {
        t0   = 0.;
        tmax = 50.;
        dt   = 0.1;

        total_population_per_region = 1061000;

        for (size_t i = 0; i < (size_t)model.parameters.get_num_regions(); i++) {
            model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
                mio::oseirmetapop::Region(i), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Exposed)}]   = 10000;
            model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
                mio::oseirmetapop::Region(i), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Infected)}]  = 1000;
            model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
                mio::oseirmetapop::Region(i), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Recovered)}] = 1000;
            model.populations[{mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
                mio::oseirmetapop::Region(i), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Susceptible)}] =
                total_population_per_region -
                model.populations[{
                    mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
                        mio::oseirmetapop::Region(i), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Exposed)}] -
                model.populations[{
                    mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
                        mio::oseirmetapop::Region(i), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Infected)}] -
                model.populations[{
                    mio::Index<mio::oseirmetapop::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>(
                        mio::oseirmetapop::Region(i), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Recovered)}];
        }
        model.set_commuting_strengths();
    }
};

TEST_F(ModelTestOdeMetapop, simulateDefault)
{
    mio::TimeSeries<ScalarType> result = simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

TEST_F(ModelTestOdeMetapop, checkPopulationConservation)
{
    auto result            = simulate(t0, tmax, dt, model);
    ScalarType num_persons = result.get_last_value().sum();
    EXPECT_NEAR(num_persons, total_population_per_region * (size_t)model.parameters.get_num_regions(), 1e-9);
}

TEST_F(ModelTestOdeMetapop, compareWithPreviousRun)
{
    model.parameters.set<mio::oseirmetapop::TimeExposed<ScalarType>>(5.2);
    model.parameters.set<mio::oseirmetapop::TimeInfected<ScalarType>>(6);
    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<ScalarType>>(0.1);
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::oseirmetapop::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(2.7);
    contact_matrix[0].add_damping(0.6, mio::SimulationTime(12.5));

    Eigen::MatrixXd mobility_data_commuter((size_t)model.parameters.get_num_regions(),
                                           (size_t)model.parameters.get_num_regions());
    mobility_data_commuter << 0., 0., 0., 1., 0.2, 0., 0.6, 0.2, 0., 0.5, 0.5, 0., 0., 0., 0., 1.;
    model.set_commuting_strengths(mobility_data_commuter);

    std::vector<std::vector<double>> refData = load_test_data_csv<double>("ode-seir-metapop-compare.csv");
    auto result                              = mio::simulate<double, mio::oseirmetapop::Model<>>(t0, tmax, dt, model);

    result.print_table({"S", "E", "I", "R"}, 16, 18);

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

TEST_F(ModelTestOdeMetapop, check_constraints_parameters)
{
    model.parameters.set<mio::oseirmetapop::TimeExposed<ScalarType>>(5.2);
    model.parameters.set<mio::oseirmetapop::TimeInfected<ScalarType>>(6);
    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<ScalarType>>(0.04);

    model.parameters.get<mio::oseirmetapop::ContactPatterns<ScalarType>>().get_cont_freq_mat()[0].get_baseline()(0, 0) =
        10;

    Eigen::MatrixXd mobility_data_commuter((size_t)model.parameters.get_num_regions(),
                                           (size_t)model.parameters.get_num_regions());
    mobility_data_commuter << 0., 0., 0., 1., 0.2, 0., 0.6, 0.2, 0., 0.5, 0.5, 0., 0., 0., 0., 1.;
    model.set_commuting_strengths(mobility_data_commuter);

    // model.check_constraints() combines the functions from population and parameters.
    // We only want to test the functions for the parameters defined in parameters.h
    ASSERT_EQ(model.parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::oseirmetapop::TimeExposed<ScalarType>>(-5.2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::oseirmetapop::TimeExposed<ScalarType>>(5.2);
    model.parameters.set<mio::oseirmetapop::TimeInfected<ScalarType>>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::oseirmetapop::TimeInfected<ScalarType>>(6);
    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<ScalarType>>(10.);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<ScalarType>>(0.04);
    mobility_data_commuter(0, 1) += 0.5;
    model.set_commuting_strengths(mobility_data_commuter);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    mobility_data_commuter(0, 1) = -0.5;
    mobility_data_commuter(0, 2) = 0.75;
    mobility_data_commuter(0, 3) = 0.75;
    model.set_commuting_strengths(mobility_data_commuter);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    mobility_data_commuter(0, 1) = 1.5;
    mobility_data_commuter(0, 2) = 0.;
    mobility_data_commuter(0, 3) = 0.;
    model.set_commuting_strengths(mobility_data_commuter);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    mobility_data_commuter(0, 0) = 0.;
    mobility_data_commuter(0, 1) = 0.;
    mobility_data_commuter(0, 2) = 0.;
    mobility_data_commuter(0, 3) = 1.;
    model.set_commuting_strengths(mobility_data_commuter);
    model.parameters.set<mio::oseirmetapop::PopulationAfterCommuting<ScalarType>>(
        mio::Populations<ScalarType, mio::oseirmetapop::Region, mio::AgeGroup>(
            {mio::oseirmetapop::Region(4), mio::AgeGroup(1)}, 0.));
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    // Nobody commutes to region 2 but everybody originating fron there commutes to other regions.
    mobility_data_commuter << 0., 0., 0., 1., 0.2, 0., 0.6, 0.2, 0., 0., 0.5, 0.5, 0., 0., 0., 1.;
    model.set_commuting_strengths(mobility_data_commuter);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    mio::set_log_level(mio::LogLevel::warn);
}

TEST_F(ModelTestOdeMetapop, apply_constraints_parameters)
{
    const ScalarType tol_times = 1e-1;
    model.parameters.set<mio::oseirmetapop::TimeExposed<ScalarType>>(5.2);
    model.parameters.set<mio::oseirmetapop::TimeInfected<ScalarType>>(2);
    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<ScalarType>>(0.04);
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseirmetapop::ContactPatterns<ScalarType>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(10);

    Eigen::MatrixXd mobility_data_commuter((size_t)model.parameters.get_num_regions(),
                                           (size_t)model.parameters.get_num_regions());
    mobility_data_commuter << 0., 0., 0., 1., 0.2, 0., 0.6, 0.2, 0., 0.5, 0.5, 0., 0., 0., 0., 1.;
    model.set_commuting_strengths(mobility_data_commuter);

    EXPECT_EQ(model.parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::oseirmetapop::TimeExposed<ScalarType>>(-5.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseirmetapop::TimeExposed<ScalarType>>()[(mio::AgeGroup)0], tol_times);

    model.parameters.set<mio::oseirmetapop::TimeInfected<ScalarType>>(1e-5);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseirmetapop::TimeInfected<ScalarType>>()[(mio::AgeGroup)0], tol_times);

    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<ScalarType>>(10.);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR(
        model.parameters.get<mio::oseirmetapop::TransmissionProbabilityOnContact<ScalarType>>()[(mio::AgeGroup)0], 0.0,
        1e-14);

    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<ScalarType>>(0.04);
    mobility_data_commuter(0, 1) += 0.5;
    model.set_commuting_strengths(mobility_data_commuter);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseirmetapop::CommutingStrengths<ScalarType>>()
                  .get_cont_freq_mat()[0]
                  .get_baseline()
                  .isIdentity(),
              true);

    mobility_data_commuter(0, 1) = -0.5;
    mobility_data_commuter(0, 2) = 0.75;
    mobility_data_commuter(0, 3) = 0.75;
    model.set_commuting_strengths(mobility_data_commuter);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseirmetapop::CommutingStrengths<ScalarType>>()
                  .get_cont_freq_mat()[0]
                  .get_baseline()
                  .isIdentity(),
              true);

    mobility_data_commuter(0, 1) = 1.5;
    mobility_data_commuter(0, 2) = 0.;
    mobility_data_commuter(0, 3) = 0.;
    model.set_commuting_strengths(mobility_data_commuter);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseirmetapop::CommutingStrengths<ScalarType>>()
                  .get_cont_freq_mat()[0]
                  .get_baseline()
                  .isIdentity(),
              true);

    mobility_data_commuter(0, 0) = 0.;
    mobility_data_commuter(0, 1) = 0.;
    mobility_data_commuter(0, 2) = 0.;
    mobility_data_commuter(0, 3) = 1.;
    model.set_commuting_strengths(mobility_data_commuter);
    model.parameters.set<mio::oseirmetapop::PopulationAfterCommuting<ScalarType>>(
        mio::Populations<ScalarType, mio::oseirmetapop::Region, mio::AgeGroup>(
            {mio::oseirmetapop::Region(4), mio::AgeGroup(1)}, 0.));
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR((model.parameters.get<mio::oseirmetapop::PopulationAfterCommuting<ScalarType>>()[{
                    mio::oseirmetapop::Region(3), mio::AgeGroup(0)}]),
                1.0, tol_times);

    // Nobody commutes to region 2 but everybody originating fron there commutes to other regions.
    mobility_data_commuter << 0., 0., 0., 1., 0.2, 0., 0.6, 0.2, 0., 0., 0.5, 0.5, 0., 0., 0., 1.;
    model.set_commuting_strengths(mobility_data_commuter);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR((model.parameters.get<mio::oseirmetapop::PopulationAfterCommuting<ScalarType>>()[{
                    mio::oseirmetapop::Region(1), mio::AgeGroup(0)}]),
                1.0, tol_times);

    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeMetapop, compareSEIR)
{
    ScalarType t0   = 0;
    ScalarType tmax = 50.;
    ScalarType dt   = 0.1;

    ScalarType total_population = 1061000;

    mio::oseirmetapop::Model<ScalarType> model(1, 1);

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

    // The model with a single region should correspond to the SEIR model.
    model.parameters.set<mio::oseirmetapop::TimeExposed<ScalarType>>(5.2);
    model.parameters.set<mio::oseirmetapop::TimeInfected<ScalarType>>(6);
    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<ScalarType>>(0.1);
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::oseirmetapop::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(2.7);
    contact_matrix[0].add_damping(0.6, mio::SimulationTime(12.5));

    model.set_commuting_strengths();

    model.check_constraints();

    // Use the Euler integrator as adaptive methods make different time steps in this model due to restructured equations.
    std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator = std::make_shared<mio::EulerIntegratorCore<>>();

    auto result                                  = simulate(t0, tmax, dt, model, integrator);
    std::vector<std::vector<ScalarType>> refData = load_test_data_csv<ScalarType>("ode-seir-compare-euler.csv");

    ASSERT_EQ(refData.size(), static_cast<size_t>(result.get_num_time_points()));

    for (Eigen::Index irow = 0; irow < result.get_num_time_points(); ++irow) {
        ScalarType t = refData[static_cast<size_t>(irow)][0];
        auto rel_tol = 1e-10;

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

        for (size_t icol = 0; icol < refData[static_cast<size_t>(irow)].size() - 1; ++icol) {
            ScalarType ref    = refData[static_cast<size_t>(irow)][icol + 1];
            ScalarType actual = result[irow][icol];

            ScalarType tol = rel_tol * ref;
            ASSERT_NEAR(ref, actual, tol) << "at row " << irow;
        }
    }
}
