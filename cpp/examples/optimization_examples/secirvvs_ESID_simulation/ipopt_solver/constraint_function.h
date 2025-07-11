#pragma once

#include <vector>

#include "../optimization_settings/secirvvs_optimization.h"
#include "models/ode_secirvvs/model.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

#include "../helpers/integrator_selector.h"
#include "../helpers/make_time_grid.h"

#include "../constraints/infection_state_utils.h"

template <typename FP>
void constraint_function(mio::osecirvvs::Model<FP> model, const SecirvvsOptimization& settings,
                         const FP* ptr_parameters, size_t n, FP* ptr_constraints, size_t m)
{
    // ----------------------------------------------------------------- //
    // Evaluate the constraints on the model.                            //
    // Step 1. Define dampings based on 'FP* ptr_parameters'.            //
    // Step 2. Gather information in 'path_values' in 'terminal_values'. //
    // Step 3. Fill information into 'FP* ptr_constraints'.              //
    // ----------------------------------------------------------------- //
    assert(n == settings.num_control_parameters() * settings.num_control_intervals());
    assert(m == settings.num_path_constraints() + settings.num_terminal_constraints());

    std::vector<FP> parameters(n);
    for (size_t i = 0; i < n; i++) {
        parameters[i] = settings.activation_function()(ptr_parameters[i]);
    }

    std::vector<FP> path_constraint_values(settings.num_path_constraints(), 0.0);
    std::vector<FP> terminal_constraint_values(settings.num_terminal_constraints(), 0.0);

    set_control_dampings<FP>(settings, model, parameters);

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());

    mio::osecirvvs::Simulation<FP> simulation(model, settings.t0(), settings.dt());

    simulation.set_integrator(make_integrator<FP>(settings.integrator_type(), settings.dt()));

    update_path_constraint<FP>(settings, simulation, path_constraint_values);
    for (size_t interval = 0; interval < settings.num_intervals(); interval++) {
        simulation.advance(time_steps[interval + 1]);
        simulation.get_dt() = settings.dt();
        update_path_constraint<FP>(settings, simulation, path_constraint_values);
    }
    update_terminal_constraint<FP>(settings, simulation, terminal_constraint_values);

    size_t idx = 0;
    for (size_t path_constraint_index = 0; path_constraint_index < settings.num_path_constraints();
         path_constraint_index++) {
        ptr_constraints[idx] = path_constraint_values[path_constraint_index];
        ptr_constraints[idx] = path_constraint_values[path_constraint_index];
        idx++;
    }
    for (size_t terminal_constraint_index = 0; terminal_constraint_index < settings.num_terminal_constraints();
         terminal_constraint_index++) {
        ptr_constraints[idx] = terminal_constraint_values[terminal_constraint_index];
        ptr_constraints[idx] = terminal_constraint_values[terminal_constraint_index];
        idx++;
    }
}
