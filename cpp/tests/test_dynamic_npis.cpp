#include "epidemiology/secir/dynamic_npis.h"
#include "epidemiology/migration/migration.h"
#include "epidemiology/secir/secir.h"
#include "matchers.h"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

TEST(DynamicNPIs, init)
{
    epi::DynamicNPIs npis;
    EXPECT_EQ(npis.get_thresholds().size(), 0);
    EXPECT_EQ(npis.get_max_exceeded_threshold(0.0), npis.get_thresholds().end());
}

TEST(DynamicNPIs, set_threshold)
{
    epi::DynamicNPIs npis;
    npis.set_threshold(0.5, {epi::DampingSampling(123.0, epi::DampingLevel(0), epi::DampingType(0),
                                                  epi::SimulationTime(0.0), {}, Eigen::VectorXd(1))});
    npis.set_threshold(1.0, {epi::DampingSampling(543.0, epi::DampingLevel(0), epi::DampingType(0),
                                                  epi::SimulationTime(0.0), {}, Eigen::VectorXd(1))});

    EXPECT_EQ(npis.get_thresholds()[0].first, 1.0);
    EXPECT_EQ(npis.get_thresholds()[0].second[0].get_value().value(), 543.0);
    EXPECT_EQ(npis.get_thresholds()[1].first, 0.5);
    EXPECT_EQ(npis.get_thresholds()[1].second[0].get_value().value(), 123.0);
}

TEST(DynamicNPIs, get_threshold)
{
    epi::DynamicNPIs npis;
    npis.set_threshold(0.5, {epi::DampingSampling(0.5, epi::DampingLevel(0), epi::DampingType(0),
                                                  epi::SimulationTime(0.0), {}, Eigen::VectorXd(1))});
    npis.set_threshold(1.0, {epi::DampingSampling(0.5, epi::DampingLevel(0), epi::DampingType(0),
                                                  epi::SimulationTime(0.0), {}, Eigen::VectorXd(1))});

    EXPECT_EQ(npis.get_max_exceeded_threshold(2.0), npis.get_thresholds().begin());
    EXPECT_EQ(npis.get_max_exceeded_threshold(0.75), npis.get_thresholds().begin() + 1);
    EXPECT_EQ(npis.get_max_exceeded_threshold(1.0), npis.get_thresholds().begin() + 1);
    EXPECT_EQ(npis.get_max_exceeded_threshold(0.25), npis.get_thresholds().end());
    EXPECT_EQ(npis.get_max_exceeded_threshold(0.5), npis.get_thresholds().end());
}

TEST(DynamicNPIs, get_damping_indices)
{
    using Damping                 = epi::Damping<epi::RectMatrixShape>;
    using DampingMatrixExpression = epi::DampingMatrixExpression<epi::Dampings<Damping>>;
    DampingMatrixExpression dampexpr(1, 1);
    dampexpr.add_damping(0.0, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.3));
    dampexpr.add_damping(0.0, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.4));
    dampexpr.add_damping(0.0, epi::DampingLevel(1), epi::DampingType(1), epi::SimulationTime(0.6));
    dampexpr.add_damping(0.0, epi::DampingLevel(0), epi::DampingType(1), epi::SimulationTime(0.5));
    dampexpr.add_damping(0.0, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.7));
    dampexpr.add_damping(0.0, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.9));

    auto indices =
        epi::get_damping_indices(dampexpr, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.3),
                                 epi::SimulationTime(0.9)); //time period is open interval!

    EXPECT_THAT(indices, testing::ElementsAre(1, 4));
}

TEST(DynamicNPIs, get_active_damping)
{
    using Damping                 = epi::Damping<epi::RectMatrixShape>;
    using DampingMatrixExpression = epi::DampingMatrixExpression<epi::Dampings<Damping>>;
    DampingMatrixExpression dampexpr(1, 1);
    dampexpr.add_damping(0.3, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.3));
    dampexpr.add_damping(0.4, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.4));
    dampexpr.add_damping(0.6, epi::DampingLevel(1), epi::DampingType(1), epi::SimulationTime(0.6));
    dampexpr.add_damping(0.5, epi::DampingLevel(0), epi::DampingType(1), epi::SimulationTime(0.5));
    dampexpr.add_damping(0.7, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.7));
    dampexpr.add_damping(0.9, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.9));

    auto a = epi::get_active_damping(dampexpr, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.6));
    EXPECT_EQ(print_wrap(a), print_wrap(Eigen::MatrixXd::Constant(1, 1, 0.4)));

    auto b = epi::get_active_damping(dampexpr, epi::DampingLevel(1), epi::DampingType(1), epi::SimulationTime(0.6));
    EXPECT_EQ(print_wrap(b), print_wrap(Eigen::MatrixXd::Constant(1, 1, 0.6)));

    auto c = epi::get_active_damping(dampexpr, epi::DampingLevel(1), epi::DampingType(1), epi::SimulationTime(0.4));
    EXPECT_EQ(print_wrap(c), print_wrap(Eigen::MatrixXd::Constant(1, 1, 0.0)));
}

TEST(DynamicNPIs, get_active_damping_empty)
{
    using Damping                 = epi::Damping<epi::RectMatrixShape>;
    using DampingMatrixExpression = epi::DampingMatrixExpression<epi::Dampings<Damping>>;
    DampingMatrixExpression dampexpr(1, 1);

    auto d = epi::get_active_damping(dampexpr, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.6));
    EXPECT_EQ(print_wrap(d), print_wrap(Eigen::MatrixXd::Constant(1, 1, 0.0)));
}

TEST(DynamicNPIs, implement_empty)
{
    using Damping                 = epi::Damping<epi::ColumnVectorShape>;
    using DampingMatrixExpression = epi::DampingMatrixExpression<epi::Dampings<Damping>>;
    epi::DampingMatrixExpressionGroup<DampingMatrixExpression> dampexprs(2, 2);
    auto make_mask = [](auto& g) {
        return g;
    };

    auto dynamic_npis = std::vector<epi::DampingSampling>({epi::DampingSampling(
        0.8, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0), {0, 1}, Eigen::VectorXd::Ones(2))});
    epi::implement_dynamic_npis(dampexprs, dynamic_npis, epi::SimulationTime(0.45), epi::SimulationTime(0.6),
                                make_mask);

    EXPECT_EQ(dampexprs[0].get_dampings().size(), 2);
    EXPECT_DOUBLE_EQ(dampexprs[0].get_dampings()[0].get_time().get(), 0.45);
    EXPECT_THAT(dampexprs[0].get_dampings()[0].get_coeffs(), MatrixNear(Eigen::VectorXd::Constant(2, 0.8)));
    EXPECT_DOUBLE_EQ(dampexprs[0].get_dampings()[1].get_time().get(), 0.6);
    EXPECT_THAT(dampexprs[0].get_dampings()[1].get_coeffs(), MatrixNear(Eigen::VectorXd::Zero(2)));
    EXPECT_EQ(dampexprs[1].get_dampings().size(), 2);
    EXPECT_DOUBLE_EQ(dampexprs[1].get_dampings()[0].get_time().get(), 0.45);
    EXPECT_THAT(dampexprs[1].get_dampings()[0].get_coeffs(), MatrixNear(Eigen::VectorXd::Constant(2, 0.8)));
    EXPECT_DOUBLE_EQ(dampexprs[1].get_dampings()[1].get_time().get(), 0.6);
    EXPECT_THAT(dampexprs[1].get_dampings()[1].get_coeffs(), MatrixNear(Eigen::VectorXd::Zero(2)));
}

TEST(DynamicNPIs, implement)
{
    using Damping                 = epi::Damping<epi::RectMatrixShape>;
    using DampingMatrixExpression = epi::DampingMatrixExpression<epi::Dampings<Damping>>;
    epi::DampingMatrixExpressionGroup<DampingMatrixExpression> dampexprs(2, 3, 1);
    auto make_mask = [](auto& g) {
        return (Eigen::MatrixXd(3, 1) << g(0, 0), g(1, 0), g(2, 0)).finished();
    };

    dampexprs[0].add_damping(0.3, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.3));
    dampexprs[0].add_damping(0.4, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.4));
    dampexprs[0].add_damping(0.5, epi::DampingLevel(0), epi::DampingType(1), epi::SimulationTime(0.5));
    dampexprs[0].add_damping(0.7, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.7));

    dampexprs[1].add_damping(0.5, epi::DampingLevel(1), epi::DampingType(1), epi::SimulationTime(0.5));
    dampexprs[1].add_damping(0.9, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.9));

    {
        auto dynamic_npis = std::vector<epi::DampingSampling>({epi::DampingSampling(
            0.8, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0), {0, 1}, Eigen::MatrixXd::Ones(3, 1))});
        epi::implement_dynamic_npis(dampexprs, dynamic_npis, epi::SimulationTime(0.45), epi::SimulationTime(0.6),
                                    make_mask);
    }

    EXPECT_THAT(
        dampexprs[0].get_dampings(),
        testing::ElementsAre(
            Damping(0.3, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.3), 3, 1), //before npi
            Damping(0.4, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.4), 3, 1), //before npi
            Damping(0.8, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.45), 3, 1), //npi begins
            Damping(0.5, epi::DampingLevel(0), epi::DampingType(1), epi::SimulationTime(0.5), 3, 1), //other type
            Damping(0.4, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.6), 3, 1), //npi ends
            Damping(0.7, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.7), 3, 1))); //after npi

    EXPECT_THAT(
        dampexprs[1].get_dampings(),
        testing::ElementsAre(
            Damping(0.8, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.45), 3, 1), //npi begins
            Damping(0.5, epi::DampingLevel(1), epi::DampingType(1), epi::SimulationTime(0.5), 3, 1), //other type/level
            Damping(0.0, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.6), 3, 1), //npi ends
            Damping(0.9, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.9), 3, 1))); //after npi

    {
        auto dynamic_npis = std::vector<epi::DampingSampling>({epi::DampingSampling(
            0.6, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0), {0}, Eigen::MatrixXd::Ones(3, 1))});
        epi::implement_dynamic_npis(dampexprs, dynamic_npis, epi::SimulationTime(0.3), epi::SimulationTime(0.9),
                                    make_mask);
    }

    EXPECT_THAT(
        dampexprs[0].get_dampings(),
        testing::ElementsAre(
            Damping(0.6, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.3), 3, 1), //new npi begins
            Damping(0.8, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.45), 3,
                    1), //old npi begins, is kept because it's bigger
            Damping(0.5, epi::DampingLevel(0), epi::DampingType(1), epi::SimulationTime(0.5), 3, 1), //other type
            Damping(0.6, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.6), 3,
                    1), //old npi ends, down to value of new npi
            Damping(0.7, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.7), 3,
                    1))); //bigger than new npi, new npi ends at t = 0.9, but is already overwritten here by a higher value

    //second matrix not changed by the new npi
    EXPECT_THAT(
        dampexprs[1].get_dampings(),
        testing::ElementsAre(Damping(0.8, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.45), 3, 1),
                             Damping(0.5, epi::DampingLevel(1), epi::DampingType(1), epi::SimulationTime(0.5), 3, 1),
                             Damping(0.0, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.6), 3, 1),
                             Damping(0.9, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0.9), 3, 1)));
}

namespace epi_test
{

struct DummyModel {
    template<class V>
    DummyModel(V&& y0)
        : result(0.0, std::forward<V>(y0))
    {
        ON_CALL(*this, advance).WillByDefault([this](auto t) {  result.add_time_point(t, result.get_last_value()); });
        ON_CALL(*this, get_result).WillByDefault(testing::ReturnRef(result));
    }

    MOCK_METHOD(void, advance, (double), ());
    MOCK_METHOD(epi::TimeSeries<double>&, get_result, ());

    epi::TimeSeries<double> result;
};

//overload required for dynamic NPIs
template<class DummyModel>
double get_infections_relative(const DummyModel&, double, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    return y[0] / y.sum();
}

//overload required for because the mock is not a compartment model simulation
template<class DummyModel>
void calculate_migration_returns(Eigen::Ref<epi::TimeSeries<double>::Vector>, const DummyModel&,
                                 Eigen::Ref<const epi::TimeSeries<double>::Vector>, double, double)
{
}

} // namespace epi_test

TEST(DynamicNPIs, migration)
{
    epi::ModelNode<testing::NiceMock<epi_test::DummyModel>> node_from((Eigen::VectorXd(2) << 0.0, 1.0).finished());
    epi::ModelNode<testing::NiceMock<epi_test::DummyModel>> node_to((Eigen::VectorXd(2) << 0.0, 1.0).finished());

    auto last_state_safe = (Eigen::VectorXd(2) << 0.01, 0.99).finished();
    auto last_state_crit = (Eigen::VectorXd(2) << 0.02, 0.98).finished();

    epi::DynamicNPIs npis;
    npis.set_threshold(
        0.015 * 100'000,
        {epi::DampingSampling{
            1.0, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0), {0}, Eigen::VectorXd::Ones(2)}});
    npis.set_duration(epi::SimulationTime(5.0));
    npis.set_base_value(100'000);

    epi::MigrationCoefficientGroup coeffs(1, 2);
    epi::MigrationParameters parameters(coeffs);
    parameters.set_dynamic_npis_infected(npis);

    epi::MigrationEdge edge(parameters);

    ASSERT_EQ(edge.get_parameters().get_coefficients()[0].get_dampings().size(), 0); //initial

    edge.apply_migration(0.5, 0.5, node_from, node_to);

    ASSERT_EQ(edge.get_parameters().get_coefficients()[0].get_dampings().size(), 0); //not check at the beginning

    EXPECT_CALL(node_from.model, advance).Times(1).WillOnce([&](auto t) { node_from.model.result.add_time_point(t, last_state_safe); });
    node_from.evolve(3.0, 3.0);
    node_to.evolve(3.0, 3.0);
    edge.apply_migration(3.0, 3.0, node_from, node_to);

    ASSERT_EQ(edge.get_parameters().get_coefficients()[0].get_dampings().size(), 0); //threshold not exceeded

    EXPECT_CALL(node_from.model, advance).Times(1).WillOnce([&](auto t) { node_from.model.result.add_time_point(t, last_state_crit); });
    node_from.evolve(6.0, 3.0);
    node_to.evolve(6.0, 3.0);
    edge.apply_migration(6.0, 3.0, node_from, node_to);

    ASSERT_EQ(edge.get_parameters().get_coefficients()[0].get_dampings().size(), 2); //NPIs implemented
}

namespace epi_test
{

class MockSimulation
{
public:
    MockSimulation(epi::SecirModel m, double t0, double /*dt*/)
        : m_model(m)
        , m_result(t0, m.get_initial_values())
    {
    }
    auto& get_result() const
    {
        return m_result;
    }
    auto& get_result()
    {
        return m_result;
    }
    auto& get_model()
    {
        return m_model;
    }
    const auto& get_model() const
    {
        return m_model;
    }
    auto advance(double t)
    { 
        //simple simulation that is constant over time
        return m_result.add_time_point(t, m_result.get_last_value());
    }

    epi::SecirModel m_model;
    epi::TimeSeries<double> m_result;
};

} // namespace epi_test

TEST(DynamicNPIs, secir_threshold_safe)
{
    epi::SecirModel model(1);
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}] = 1.0;
    model.populations.set_difference_from_total({epi::AgeGroup(0), epi::InfectionState::Susceptible}, 100.0);
    
    epi::DynamicNPIs npis;
    npis.set_threshold(
        0.05 * 23'000,
        {epi::DampingSampling{
            1.0, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0), {0}, Eigen::VectorXd::Ones(1)}});
    npis.set_duration(epi::SimulationTime(5.0));
    npis.set_base_value(23'000);
    model.parameters.get<epi::DynamicNPIsInfected>() = npis;

    ASSERT_EQ(model.parameters.get<epi::ContactPatterns>().get_cont_freq_mat()[0].get_dampings().size(), 0);
    
    epi::SecirSimulation<epi_test::MockSimulation> sim(model);
    sim.advance(3.0);
    
    ASSERT_EQ(sim.get_model().parameters.get<epi::ContactPatterns>().get_cont_freq_mat()[0].get_dampings().size(), 0);
}

TEST(DynamicNPIs, secir_threshold_exceeded)
{
    epi::SecirModel model(1);
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}] = 10;
    model.populations.set_difference_from_total({epi::AgeGroup(0), epi::InfectionState::Susceptible}, 100);
    
    epi::DynamicNPIs npis;
    npis.set_threshold(
        0.05 * 50'000,
        {epi::DampingSampling{
            1.0, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0), {0}, Eigen::VectorXd::Ones(1)}});
    npis.set_duration(epi::SimulationTime(5.0));
    npis.set_base_value(50'000);
    model.parameters.get<epi::DynamicNPIsInfected>() = npis;

    ASSERT_EQ(model.parameters.get<epi::ContactPatterns>().get_cont_freq_mat()[0].get_dampings().size(), 0);
    
    epi::SecirSimulation<epi_test::MockSimulation> sim(model);
    sim.advance(3.0);
    
    ASSERT_EQ(sim.get_model().parameters.get<epi::ContactPatterns>().get_cont_freq_mat()[0].get_dampings().size(), 2);
}