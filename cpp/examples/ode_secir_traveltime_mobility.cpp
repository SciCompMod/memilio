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
 * Example for the travel-time-aware metapopulation mobility model with three patches. 
 *
 * Only A has initial infections. Individuals from A commute to B via C, spending time in C's mobility node during transit.
 * Due to the time spent in the mobility node, C residents get infected early on, and B receives infections both from
 * arriving A commuters and from C residents who commute directly to B.
 *
 * This example demonstrates the travel-time-aware Graph-ODE model from:
 *
 * H. Zunker et al., "Novel travel time aware metapopulation models ...", (2025), doi.org/10.1371/journal.pcbi.1012630
 */

#include "memilio/config.h"
#include "memilio/mobility/metapopulation_mobility_traveltime.h"
#include "ode_secir/model.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/parameters.h"
#include "memilio/compartments/simulation.h"

#include <iostream>
#include <iomanip>
#include <string>

using FP    = double;
using Model = mio::osecir::Model<FP>;
using Sim   = mio::osecir::Simulation<FP>;

static Model build_model(double N, double I0)
{
    Model model(1 /*age groups*/);

    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}]        = N - I0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = I0;

    model.parameters.set<mio::osecir::StartDay<FP>>(0);
    model.parameters.set<mio::osecir::Seasonality<FP>>(0.0);

    model.parameters.get<mio::osecir::TimeExposed<FP>>()[mio::AgeGroup(0)]            = 3.2;
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<FP>>()[mio::AgeGroup(0)] = 2.0;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms<FP>>()[mio::AgeGroup(0)]   = 5.8;
    model.parameters.get<mio::osecir::TimeInfectedSevere<FP>>()[mio::AgeGroup(0)]     = 9.5;
    model.parameters.get<mio::osecir::TimeInfectedCritical<FP>>()[mio::AgeGroup(0)]   = 7.1;

    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<FP>>()[mio::AgeGroup(0)]  = 0.08;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<FP>>()[mio::AgeGroup(0)]    = 0.7;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<FP>>()[mio::AgeGroup(0)]    = 0.09;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<FP>>()[mio::AgeGroup(0)]    = 0.25;
    model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<FP>>()[mio::AgeGroup(0)] = 0.45;
    model.parameters.get<mio::osecir::TestAndTraceCapacity<FP>>()                                = 35;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<FP>>()[mio::AgeGroup(0)]         = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere<FP>>()[mio::AgeGroup(0)]                 = 0.25;
    model.parameters.get<mio::osecir::DeathsPerCritical<FP>>()[mio::AgeGroup(0)]                 = 0.3;

    auto& cm                                 = model.parameters.get<mio::osecir::ContactPatterns<FP>>();
    cm.get_cont_freq_mat()[0].get_baseline() = Eigen::MatrixX<FP>::Constant(1, 1, 10.0);

    return model;
}

static void print_row(const std::string& label, const auto& node_prop)
{
    const auto& last = node_prop.local_sim.get_result().get_last_value();
    const int iS     = static_cast<int>(mio::osecir::InfectionState::Susceptible);
    const int iE     = static_cast<int>(mio::osecir::InfectionState::Exposed);
    const int iINS   = static_cast<int>(mio::osecir::InfectionState::InfectedNoSymptoms);
    std::cout << std::left << std::setw(8) << label << "  S=" << std::setw(10) << static_cast<long>(last[iS])
              << "  E=" << std::setw(8) << static_cast<long>(last[iE]) << "  INS=" << std::setw(8)
              << static_cast<long>(last[iINS]) << "  Total=" << static_cast<long>(last.sum()) << "\n";
}

int main()
{
    mio::set_log_level(mio::LogLevel::warn);

    const FP t0   = 0.0;
    const FP tmax = 30.0;
    const FP dt   = 0.5;

    // -- Population sizes -------------------------------------------------
    const double N_A  = 500'000.0;
    const double N_B  = 300'000.0;
    const double N_C  = 80'000.0;
    const double I0_A = N_A * 0.001; // 0.1 % infected at start
    const double I0_B = 0.0;
    const double I0_C = 0.0;

    // -- Mobility parameters -----------------------------------------------
    // A -> B trip passes through C
    // Choosing 1 h per patch means A traveler spend 1 h in C's mobility node
    const FP tt_patch_AC_CB = 1.0 / 24.0; // 1 h per patch
    const FP tt_direct_BC   = 0.5 / 24.0; // 30 min for direct B <-> C travel
    const FP stay_dur_AB    = 8.0 / 24.0;
    const FP stay_dur_C     = 6.0 / 24.0;

    const FP mob_coeff_AB = 0.05; // 5 % of A commute to B
    const FP mob_coeff_BC = 0.08; // 8 % of B/C commute to each other

    const int nc                       = static_cast<int>(mio::osecir::InfectionState::Count);
    const Eigen::VectorX<FP> coeffs_AB = Eigen::VectorX<FP>::Constant(nc, mob_coeff_AB);
    const Eigen::VectorX<FP> coeffs_BC = Eigen::VectorX<FP>::Constant(nc, mob_coeff_BC);

    // Nodes (index 0=A, 1=B, 2=C)
    mio::Graph<mio::TravelTimeNodeProperty<FP, Sim>, mio::TravelTimeEdge<FP>> graph;
    graph.add_node(0, build_model(N_A, I0_A), t0, stay_dur_AB);
    graph.add_node(1, build_model(N_B, I0_B), t0, stay_dur_AB);
    graph.add_node(2, build_model(N_C, I0_C), t0, stay_dur_C);

    // A -> B via C
    const FP tt_AB                    = 2.0 * tt_patch_AC_CB;
    const std::vector<size_t> path_AB = {0, 2, 1}; // A -> C -> B
    graph.add_edge(0, 1, mio::MobilityParameters<FP>(coeffs_AB), tt_AB, path_AB);

    // B <-> C
    graph.add_edge(1, 2, mio::MobilityParameters<FP>(coeffs_BC), tt_direct_BC, std::vector<size_t>{1, 2});
    graph.add_edge(2, 1, mio::MobilityParameters<FP>(coeffs_BC), tt_direct_BC, std::vector<size_t>{2, 1});

    // -- Run simulation ----------------------------------------------------
    auto sim = mio::make_traveltime_sim<FP>(t0, dt, std::move(graph));
    sim.advance(tmax);

    // -- Results -----------------------------------------------------------
    std::cout << "Patch A starts with " << I0_A << " infected. A traveler from A travels to B via C.\n"
              << "C is transit-only for A->B; C also has its own B<->C travelers.\n"
              << "B and C start healthy and infections reach them through transit mixing.\n\n";

    std::cout << std::string(62, '-') << "\n";
    std::cout << "\n Results after " << tmax << " days \n\n";
    print_row("A", sim.get_graph().nodes()[0].property);
    print_row("B", sim.get_graph().nodes()[1].property);
    print_row("C", sim.get_graph().nodes()[2].property);
    std::cout << std::string(62, '-') << "\n";

    std::cout << "\nNote: C imports infections early because A commuters\n"
              << "      spend ~1 h inside C's mobility node,\n"
              << "      infecting susceptible there. This is not captured by a\n"
              << "      standard metapopulation model without travel time,\n"
              << "      where transit is instantaneous.\n";

    return 0;
}
