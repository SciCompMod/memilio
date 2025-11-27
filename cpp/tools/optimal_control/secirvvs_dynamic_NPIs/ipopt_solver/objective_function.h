#pragma once

#include "../optimization_settings/optimization_settings.h"
#include "models/ode_secirvvs/model.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

#include "../helpers/integrator_selector.h"
#include "../helpers/make_time_grid.h"

#include "../constraints/infection_state_utils.h"

#include <vector>

template <typename FP>
FP objective_function(
    mio::Graph<mio::SimulationNode<FP, mio::osecirvvs::Simulation<FP>>, mio::MobilityEdge<FP>> graph_model,
    const SecirvvsOptimization& settings, const FP* ptr_parameters, size_t n)
{
    // ------------------------------------------------------------------ //
    // Evaluate the objective function of the model.                      //
    // Step 1. Define dampings based on 'const FP* parameters'.           //
    // Step 2. Evaluate the objective function based on                   //
    //         the parameters and the infection states in the simulation. //
    // ------------------------------------------------------------------ //
    // assert(n == 2 * settings.num_control_parameters());

    FP objective = 0.0;

    std::vector<FP> dynamic_NPI_strengths(settings.num_control_parameters());
    for (size_t control_index = 0; control_index < settings.num_control_parameters(); control_index++) {
        dynamic_NPI_strengths[control_index] = ptr_parameters[control_index];
    }

    for (auto& node : graph_model.nodes()) {
        set_dynamic_NPIs<FP>(settings, node.property.get_simulation().get_model(), dynamic_NPI_strengths);
    }

    for (auto& node : graph_model.nodes()) {
        node.property.get_simulation().set_integrator_core(
            std::move(make_integrator<FP>(settings.integrator_type(), settings.dt())));
    }

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());
    auto graph_simulation      = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), graph_model);

    FP total_population = 0.0;
    for (auto& node : graph_simulation.get_graph().nodes()) {
        total_population += node.property.get_simulation().get_model().populations.get_total();
    }

    for (size_t interval = 0; interval < settings.num_intervals(); interval++) {

        graph_simulation.advance(time_steps[interval + 1]);

        for (auto& node : graph_simulation.get_graph().nodes()) {
            auto& sim          = node.property.get_simulation();
            auto& model        = sim.get_model();
            auto& dynamic_npis = model.parameters.template get<mio::osecirvvs::DynamicNPIsInfectedSymptoms<FP>>();

            // --- Compute infections relative per 100k directly ---
            const auto& y = sim.get_result().get_last_value();

            FP sum_inf = 0.0;
            for (mio::AgeGroup i = 0; i < model.parameters.get_num_groups(); ++i) {
                sum_inf +=
                    sim.get_model().populations.get_from(y, {i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive});
                sum_inf += sim.get_model().populations.get_from(
                    y, {i, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed});
                sum_inf += sim.get_model().populations.get_from(
                    y, {i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity});
                sum_inf += sim.get_model().populations.get_from(
                    y, {i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity});
                sum_inf += sim.get_model().populations.get_from(
                    y, {i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed});
                sum_inf += sim.get_model().populations.get_from(
                    y, {i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed});
            }
            FP inf_rel = (sum_inf / model.populations.get_total()) * dynamic_npis.get_base_value();

            int index = 0;
            for (auto& cp : settings.control_parameters()) {
                FP threshold = cp.threshold();
                for (auto& d : cp.dampings()) {
                    FP strength = dynamic_NPI_strengths[index++];
                    if (inf_rel >= threshold) {
                        // Strength multiplied by the number of people affected
                        objective += strength * (model.populations.get_total() / total_population);
                    }
                }
            }
        }
    }

    return objective;
}
