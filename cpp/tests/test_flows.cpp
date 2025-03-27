/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding, Henrik Zunker 
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
#include "load_test_data.h"
#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/compartments/simulation.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/model.h"
#include "ode_seir/parameters.h"

#include "gtest/gtest.h"

using I     = mio::oseir::InfectionState;
using Flows = mio::TypeList<mio::Flow<I::Susceptible, I::Exposed>, mio::Flow<I::Exposed, I::Infected>,
                            mio::Flow<I::Infected, I::Recovered>>;

struct CatA : public mio::Index<CatA> {
    CatA(size_t i)
        : mio::Index<CatA>(i)
    {
    }
};
struct CatB : public mio::Index<CatB> {
    CatB(size_t i)
        : mio::Index<CatB>(i)
    {
    }
};
struct CatC : public mio::Index<CatC> {
    CatC(size_t i)
        : mio::Index<CatC>(i)
    {
    }
};

class TestModel : public mio::FlowModel<double, I, mio::Populations<double, I, CatA, CatB, CatC>,
                                        mio::oseir::Parameters<double>, Flows>
{
    using Base =
        mio::FlowModel<double, I, mio::Populations<double, I, CatA, CatB, CatC>, mio::oseir::Parameters<double>, Flows>;

public:
    TestModel(Populations::Index dimensions)
        : Base(Populations(dimensions, 0.), mio::oseir::Parameters(mio::AgeGroup(1)))
    {
    }
    void get_flows(Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> /*y*/, double /*t*/,
                   Eigen::Ref<Eigen::VectorXd> /*flows*/) const override
    {
    }
};

TEST(TestFlows, FlowChart)
{
    EXPECT_EQ(Flows().size(), 3);
    // Testing get (by index) function, verifying with source members.
    auto flow0        = mio::Flow<I::Susceptible, I::Exposed>().source;
    auto test_source0 = mio::type_at_index_t<0, Flows>::source;
    EXPECT_EQ(flow0, test_source0);

    auto flow1        = mio::Flow<I::Exposed, I::Infected>().source;
    auto test_source1 = mio::type_at_index_t<1, Flows>::source;
    EXPECT_EQ(flow1, test_source1);

    auto flow2        = mio::Flow<I::Infected, I::Recovered>().source;
    auto test_source2 = mio::type_at_index_t<2, Flows>::source;
    EXPECT_EQ(flow2, test_source2);

    // Testing get (by Flow) function.
    const size_t index0 = mio::index_of_type_v<mio::Flow<I::Susceptible, I::Exposed>, Flows>;
    EXPECT_EQ(index0, 0);

    const size_t index1 = mio::index_of_type_v<mio::Flow<I::Exposed, I::Infected>, Flows>;
    EXPECT_EQ(index1, 1);

    const size_t index2 = mio::index_of_type_v<mio::Flow<I::Infected, I::Recovered>, Flows>;
    EXPECT_EQ(index2, 2);
}

TEST(TestFlows, FlowSimulation)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.001;

    mio::oseir::Model<double> model(1);

    double total_population                                                      = 10000;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]   = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]  = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}] = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}];
    // suscetible now set with every other update
    // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
    model.parameters.set<mio::oseir::TimeExposed<double>>(5.2);
    model.parameters.set<mio::oseir::TimeInfected<double>>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(0.04);
    model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat()[0].get_baseline().setConstant(10);

    mio::set_log_level(mio::LogLevel::off); // Suppress log output of check_constraints and the Simulation.
    model.check_constraints();
    auto IC   = std::make_shared<mio::DefaultIntegratorCore<double>>();
    auto seir = mio::simulate_flows<double, mio::oseir::Model<double>>(t0, tmax, dt, model, IC);
    mio::set_log_level(mio::LogLevel::warn);

    // verify results (computed using flows)
    auto results = seir[0].get_last_value();
    EXPECT_NEAR(results[0], 9660.5835936179408, 1e-14);
    EXPECT_NEAR(results[1], 118.38410512338650, 1e-14);
    EXPECT_NEAR(results[2], 104.06636087558745, 1e-14);
    EXPECT_NEAR(results[3], 116.96594038308581, 1e-14);
    // test flow results
    auto flows_results = seir[1].get_last_value();
    EXPECT_NEAR(flows_results[0], 39.416406382059762, 1e-14);
    EXPECT_NEAR(flows_results[1], 21.032301258673261, 1e-14);
    EXPECT_NEAR(flows_results[2], 16.965940383085815, 1e-14);
}

TEST(TestFlows, CompareSimulations)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.001;

    mio::oseir::Model<double> model(1);

    double total_population                                                      = 10000;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]   = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]  = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}] = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}];
    // suscetible now set with every other update
    // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
    model.parameters.set<mio::oseir::TimeExposed<double>>(5.2);
    model.parameters.set<mio::oseir::TimeInfected<double>>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(0.04);
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(10);

    mio::set_log_level(mio::LogLevel::off); // Suppress log output of check_constraints and the Simulation.
    model.check_constraints();
    auto seir_sim_flows = simulate_flows(t0, tmax, dt, model);
    auto seir_sim       = simulate(t0, tmax, dt, model);
    mio::set_log_level(mio::LogLevel::warn);

    auto results_flows = seir_sim_flows[0].get_last_value();
    auto results       = seir_sim.get_last_value();

    EXPECT_NEAR(results[0], results_flows[0], 1e-10);
    EXPECT_NEAR(results[1], results_flows[1], 1e-10);
    EXPECT_NEAR(results[2], results_flows[2], 1e-10);
    EXPECT_NEAR(results[3], results_flows[3], 1e-10);
}

TEST(TestFlows, GetInitialFlows)
{
    mio::oseir::Model<double> m(1);
    EXPECT_EQ(m.get_initial_flows().size(), 3); // 3 == Flows().size()
    EXPECT_EQ(m.get_initial_flows().norm(), 0);
}

TEST(TestFlows, GetFlowIndex)
{
    // test get_flat_flow_index with some prime number products
    TestModel m({I::Count, CatA(11), CatB(5), CatC(7)});
    auto idx0 = m.get_flat_flow_index<I::Susceptible, I::Exposed>({CatA(0), CatB(0), CatC(1)});
    EXPECT_EQ(idx0, 3);

    auto idx1 = m.get_flat_flow_index<I::Susceptible, I::Exposed>({CatA(0), CatB(1), CatC(0)});
    EXPECT_EQ(idx1, 7 * 3);

    auto idx2 = m.get_flat_flow_index<I::Susceptible, I::Exposed>({CatA(1), CatB(0), CatC(0)});
    EXPECT_EQ(idx2, 5 * 7 * 3);

    auto idx3 = m.get_flat_flow_index<I::Susceptible, I::Exposed>({CatA(1), CatB(1), CatC(1)});
    EXPECT_EQ(idx3, (5 * 7 * 3) + (7 * 3) + (3));

    auto idx4 = m.get_flat_flow_index<I::Susceptible, I::Exposed>({CatA(10), CatB(4), CatC(6)});
    EXPECT_EQ(idx4, 10 * (5 * 7 * 3) + 4 * (7 * 3) + 6 * (3));
}
