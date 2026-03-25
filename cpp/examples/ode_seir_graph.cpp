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

/**
 * Two-patch SEIR metapopulation with stage-aligned mobility returns.
 *
 * The stage-aligned approach does:
 *   1. A SEIR-specific RK4 integrator (SeirRK4IntegratorCore) that caches the
 *      force-of-infection at each stage during node integration.
 *   2. Reuse the cached values to compute mobility returns with only apply_rhs calls - no
 *      additional compute_foi work per edge.
 *
 */

#include "ode_seir/model.h"
#include "ode_seir/mobility.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"

#include <iostream>
#include <iomanip>

int main()
{
    using FP = double;

    mio::set_log_level(mio::LogLevel::warn);

    const FP t0   = 0.0;
    const FP tmax = 30.0;
    const FP dt   = 0.5;

    // SEIR model with single age group
    const size_t num_groups = 1;
    mio::oseir::Model<FP> model(static_cast<int>(num_groups));

    model.parameters.set<mio::oseir::TimeExposed<FP>>(5.2);
    model.parameters.set<mio::oseir::TimeInfected<FP>>(6.0);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<FP>>(0.1);

    auto& cm = model.parameters.get<mio::oseir::ContactPatterns<FP>>();
    cm.get_cont_freq_mat()[0].get_baseline().setConstant(2.7);

    // Patch 1
    auto model1                                                                     = model;
    model1.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 9700;
    model1.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = 200;
    model1.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]    = 80;
    model1.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}]   = 20;

    // Patch 2
    auto model2                                                                     = model;
    model2.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 9990;
    model2.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = 5;
    model2.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]    = 3;
    model2.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}]   = 2;

    // Build mobility graph: 10% of each compartment commutes daily between the two patches
    const auto NC =
        static_cast<Eigen::Index>(mio::oseir::InfectionState::Count) * static_cast<Eigen::Index>(num_groups);

    mio::Graph<mio::SimulationNode<FP, mio::Simulation<FP, mio::oseir::Model<FP>>>, mio::MobilityEdge<FP>> g;
    g.add_node(1001, model1, t0);
    g.add_node(1002, model2, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant(NC, 0.1));
    g.add_edge(1, 0, Eigen::VectorXd::Constant(NC, 0.1));

    // make_seir_mobility_sim sets SeirRK4IntegratorCore on every node as integrator,
    // synchronizes node dt with graph dt, and returns a GraphSimulation.
    auto sim = mio::oseir::make_seir_mobility_sim<FP>(t0, dt, std::move(g));

    sim.advance(tmax);

    // Print results
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n=== Stage-Aligned SEIR Graph Simulation Results ===\n\n";

    for (size_t n = 0; n < sim.get_graph().nodes().size(); ++n) {
        auto& result = sim.get_graph().nodes()[n].property.get_result();
        auto last    = result.get_last_value();
        std::cout << "Patch " << n << " at t=" << result.get_last_time() << ":\n"
                  << "  S = " << last[0] << "\n"
                  << "  E = " << last[1] << "\n"
                  << "  I = " << last[2] << "\n"
                  << "  R = " << last[3] << "\n"
                  << "  Total = " << last.sum() << "\n\n";
    }

    return 0;
}
