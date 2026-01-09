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
FP objective_function(const SecirvvsOptimization& settings, const FP* ptr_parameters, size_t n)
{
    assert(n == settings.num_control_parameters() * settings.num_control_intervals());

    FP objective_sum = 0.0;

    // Loop over all simulation runs to account for uncertainty
    for (size_t run = 0; run < settings.num_simulation_runs(); run++) {

        // Set seed for this run
        size_t seed = settings.base_seed() + run;

        // Define graph model for this run
        auto graph_model = settings.optimization_model().get_graph_model<FP>(seed);

        // Extract control parameters
        std::vector<FP> parameters(settings.num_control_parameters() * settings.num_control_intervals());
        for (size_t control_interval = 0; control_interval < settings.num_control_intervals(); control_interval++) {
            for (size_t control_index = 0; control_index < settings.num_control_parameters(); control_index++) {
                size_t idx      = control_interval * settings.num_control_parameters() + control_index;
                parameters[idx] = settings.activation_function()(ptr_parameters[idx]);
            }
        }

        // Set control dampings for each node
        for (auto& node : graph_model.nodes()) {
            set_control_dampings<FP>(settings, node.property.get_simulation().get_model(), parameters);
        }

        // Set integrators for each node
        for (auto& node : graph_model.nodes()) {
            node.property.get_simulation().set_integrator_core(
                std::move(make_integrator<FP>(settings.integrator_type(), settings.dt())));
        }

        // Prepare time grid and simulation
        std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());
        auto graph_simulation      = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), graph_model);

        // Run simulation and compute objective for this run
        FP objective = 0.0;
        for (size_t interval = 0; interval < settings.num_intervals(); interval++) {

            graph_simulation.advance(time_steps[interval + 1]);

            size_t control_interval = interval / settings.pc_resolution();

            FP interval_cost = 0.0;
            for (size_t control_index = 0; control_index < settings.num_control_parameters(); control_index++) {
                interval_cost += settings.control_parameters()[control_index].cost() *
                                 parameters[control_interval * settings.num_control_parameters() + control_index];
            }
            objective += interval_cost / settings.num_intervals();
        }

        // Accumulate objective
        objective_sum += objective;
    }

    // Return average over all simulation runs
    return objective_sum / static_cast<FP>(settings.num_simulation_runs());
}
