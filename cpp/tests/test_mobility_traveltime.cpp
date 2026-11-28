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

#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/mobility/metapopulation_mobility_traveltime.h"
#include "memilio/mobility/traveltime_schedule.h"
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "memilio/math/floating_point.h"
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"

#include "gtest/gtest.h"
#include "utils.h"
#include <cmath>
#include <limits>

using FP    = double;
using Model = mio::oseir::Model<FP>;
using Sim   = mio::Simulation<FP, Model>;

static Model make_seir_model(double N, double E0)
{
    Model model(1);
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = static_cast<FP>(N - E0);
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = static_cast<FP>(E0);
    model.populations.set_total(static_cast<FP>(N));
    model.parameters.get<mio::oseir::ContactPatterns<FP>>().get_cont_freq_mat()[0].get_baseline().setConstant(10.0);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<FP>>(0.1);
    model.parameters.set<mio::oseir::TimeExposed<FP>>(4.0);
    model.parameters.set<mio::oseir::TimeInfected<FP>>(10.0);
    return model;
}

static mio::Graph<mio::TravelTimeNodeProperty<FP, Sim>, mio::TravelTimeEdge<FP>>
make_two_node_graph(double N0, double E0, double N1, double E1, FP travel_time, FP stay, FP mob_coeff)
{
    const int n_comp          = static_cast<int>(mio::oseir::InfectionState::Count);
    Eigen::VectorX<FP> coeffs = Eigen::VectorX<FP>::Constant(n_comp, mob_coeff);

    mio::Graph<mio::TravelTimeNodeProperty<FP, Sim>, mio::TravelTimeEdge<FP>> g;
    g.add_node(0, make_seir_model(N0, E0), FP{0}, stay);
    g.add_node(1, make_seir_model(N1, E1), FP{0}, stay);
    g.add_edge(0, 1, mio::MobilityParameters<FP>(coeffs), travel_time, std::vector<size_t>{0, 1});
    g.add_edge(1, 0, mio::MobilityParameters<FP>(coeffs), travel_time, std::vector<size_t>{1, 0});
    return g;
}

TEST(TravelTimeMobility, ScheduleBasicStructure)
{
    const FP travel_time = 0.05;
    const FP stay_dur    = 0.4;
    const size_t n_steps = 100;

    auto graph = make_two_node_graph(1000.0, 10.0, 800.0, 0.0, travel_time, stay_dur, 0.1);
    auto sched = mio::TravelTimeSchedule(graph, n_steps, mio::Limits<FP>::zero_tolerance());

    // Schedule must have one entry per edge.
    ASSERT_EQ(sched.num_edges(), graph.edges().size());

    // Each edge schedule must have n_steps entries.
    ASSERT_EQ(sched.num_steps(), n_steps);

    // Breakpoints must be non-empty for every node.
    ASSERT_EQ(sched.num_nodes(), graph.nodes().size());
    for (size_t ni = 0; ni < graph.nodes().size(); ++ni) {
        EXPECT_FALSE(sched.local_breakpoints(ni).empty()) << "local_breakpoints empty for node " << ni;
    }

    // first_mobility_step must be valid indices.
    for (size_t ei = 0; ei < graph.edges().size(); ++ei) {
        EXPECT_LT(sched.first_mobility_step(ei), n_steps) << "first_mobility_step out of range for edge " << ei;
    }
}

TEST(TravelTimeMobility, PopulationConservation)
{
    mio::LogLevelOverride llo(mio::LogLevel::off);
    const FP t0          = 0.0;
    const FP tmax        = 3.0;
    const FP travel_time = 0.04;
    const FP stay_dur    = 0.3;
    const double N0      = 5000.0;
    const double N1      = 3000.0;

    auto graph = make_two_node_graph(N0, 50.0, N1, 0.0, travel_time, stay_dur, 0.1);
    auto sim   = mio::make_traveltime_sim<FP>(t0, FP{1.0}, std::move(graph));
    sim.advance(tmax);

    const FP total_local0 = sim.get_graph().nodes()[0].property.local_sim.get_result().get_last_value().sum();
    const FP total_local1 = sim.get_graph().nodes()[1].property.local_sim.get_result().get_last_value().sum();
    const FP total_mob0   = sim.get_graph().nodes()[0].property.mobility_sim.get_result().get_last_value().sum();
    const FP total_mob1   = sim.get_graph().nodes()[1].property.mobility_sim.get_result().get_last_value().sum();

    // At the end of complete days, all commuters must be home.
    EXPECT_NEAR(total_mob0, 0.0, 5.0) << "Mobility node 0 not empty at end of simulation";
    EXPECT_NEAR(total_mob1, 0.0, 5.0) << "Mobility node 1 not empty at end of simulation";

    // Total system population is conserved (SEIR: no deaths).
    EXPECT_NEAR(total_local0 + total_local1, N0 + N1, 5.0) << "Total population not conserved";
}

TEST(TravelTimeMobility, ZeroMobilityCoefficient)
{
    const FP t0   = 0.0;
    const FP tmax = 5.0;

    // Reference plain single-node simulation.
    Model model0 = make_seir_model(2000.0, 20.0);
    auto single  = mio::Simulation<FP, Model>(model0, t0);
    single.advance(tmax);
    const Eigen::VectorX<FP> ref = single.get_result().get_last_value().eval();

    // TravelTime graph with zero mobility coefficient.
    auto graph = make_two_node_graph(2000.0, 20.0, 1500.0, 0.0, 0.05, 0.4, 0.0);
    auto sim   = mio::make_traveltime_sim<FP>(t0, FP{1.0}, std::move(graph));
    sim.advance(tmax);

    const auto& v_tt = sim.get_graph().nodes()[0].property.local_sim.get_result().get_last_value();

    // With zero mobility, node 0 should evolve as the isolated simulation.
    EXPECT_NEAR((v_tt - ref).norm(), 0.0, 1.0) << "Zero-mobility node deviates from isolated simulation";
}

TEST(TravelTimeMobility, LargeTravelTimeDoesNotCrash)
{
    mio::LogLevelOverride llo(mio::LogLevel::off); // suppress expected errors from correct_negative_compartments
    const FP t0 = 0.0;

    // Nearly all-day travel: 49% + 1% stay + 49% return = 99% of day in travel.
    auto graph = make_two_node_graph(1000.0, 10.0, 800.0, 0.0, 0.49, 0.01, 0.05);

    EXPECT_NO_THROW({
        auto sim = mio::make_traveltime_sim<FP>(t0, FP{1.0}, std::move(graph));
        sim.advance(FP{2.0});
    });
}

TEST(TravelTimeMobility, CorrectNegativeCompartments)
{
    // 2 age groups, 4 compartments each == size 8.
    // Negative values are corrected via map_to_nonnegative per age group slice.
    Eigen::VectorX<FP> v(8);
    v << 100.0, 50.0, -1e-5, 30.0, // group 0: negative at index 2
        200.0, 80.0, 40.0, -2e-5; // group 1: negative at index 7

    const FP sum0_before = v.head(4).sum();
    const FP sum1_before = v.tail(4).sum();

    mio::correct_negative_compartments<FP>(v, 2);

    for (int i = 0; i < v.size(); ++i) {
        EXPECT_GE(v[i], FP{0}) << "Negative value remains at index " << i;
    }

    // Population within each age group is conserved
    EXPECT_NEAR(v.head(4).sum(), sum0_before, 1e-12) << "Population not conserved in group 0";
    EXPECT_NEAR(v.tail(4).sum(), sum1_before, 1e-12) << "Population not conserved in group 1";
}

TEST(TravelTimeMobility, ScheduleNodeAssignmentTwoNode)
{
    const FP travel_time = 0.1;
    const FP stay_dur    = 0.3;
    const size_t n_steps = 100;

    auto graph = make_two_node_graph(1000.0, 0.0, 800.0, 0.0, travel_time, stay_dur, 0.1);
    auto sched = mio::TravelTimeSchedule(graph, n_steps, mio::Limits<FP>::zero_tolerance());

    // Edge 0->1 is edge index 0 (first edge added).
    const size_t ei = 0;
    ASSERT_LT(ei, sched.num_edges());

    const size_t first = sched.first_mobility_step(ei);
    EXPECT_LT(first, n_steps);

    // During the stay phase, in_mobility must be false and node must be dest (1).
    bool found_stay = false;
    for (size_t s = first; s < n_steps; ++s) {
        if (!sched.in_mobility_at(ei, s)) {
            found_stay = true;
            EXPECT_EQ(sched.node_at(ei, s), size_t{1}) << "Expected destination node during stay at step " << s;
        }
    }
    EXPECT_TRUE(found_stay) << "No stay phase found in edge schedule";
}
