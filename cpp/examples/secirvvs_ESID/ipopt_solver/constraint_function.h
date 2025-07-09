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
    std::vector<FP> terminal_constraint_values(settings.num_terminal_constraints(), 42.0);

    set_control_dampings<FP>(settings, model, parameters);

    auto integrator = make_integrator<FP>(settings.integrator_type(), settings.dt());

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());

    update_path_constraint<FP>(settings, model, path_constraint_values);
    for (size_t interval = 0; interval < settings.num_intervals(); interval++) {
        mio::TimeSeries<FP> result = mio::simulate<FP, mio::osecirvvs::Model<FP>>(
            time_steps[interval], time_steps[interval + 1], settings.dt(), model, integrator);
        const auto& final_state = result.get_last_value();

        for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
            for (size_t state_index = 0; state_index < num_infection_states(); state_index++) {
                size_t idx = age_group.get() * num_infection_states() + state_index;
                model.populations[{age_group, mio::osecirvvs::InfectionState(state_index)}] = final_state[idx];
            }
        }
        update_path_constraint<FP>(settings, model, path_constraint_values);
    }
    update_terminal_constraint<FP>(settings, model, terminal_constraint_values);

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
