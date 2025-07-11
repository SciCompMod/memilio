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
FP objective_function(mio::osecirvvs::Model<FP> model, const SecirvvsOptimization& settings, const FP* ptr_parameters,
                      size_t n)
{
    // ------------------------------------------------------------------ //
    // Evaluate the objective function of the model.                      //
    // Step 1. Define dampings based on 'const FP* parameters'.           //
    // Step 2. Evaluate the objective function based on                   //
    //         the parameters and the infection states in the simulation. //
    // ------------------------------------------------------------------ //
    assert(n == settings.num_control_parameters() * settings.num_control_intervals());

    FP objective = 0.0;

    std::vector<FP> parameters(n);
    for (size_t i = 0; i < n; i++) {
        parameters[i] = settings.activation_function()(ptr_parameters[i]);
    }

    set_control_dampings<FP>(settings, model, parameters);

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());

    // mio::osecirvvs::Simulation<FP> simulation(model, settings.t0(), settings.dt());
    // simulation.set_integrator(make_integrator<FP>(settings.integrator_type(), settings.dt()));

    for (size_t interval = 0; interval < settings.num_intervals(); interval++) {

        // simulation.advance(time_steps[interval + 1]);
        // simulation.get_dt() = settings.dt();

        size_t control_interval = interval / settings.pc_resolution();

        auto param_at = [&](const std::string& name) {
            size_t control_index = static_cast<size_t>(string_to_control(name));
            return parameters[control_index + control_interval * settings.num_control_parameters()];
        };

        auto cost = [&](const std::string& name) {
            size_t control_index = static_cast<size_t>(string_to_control(name));
            return settings.control_parameters()[control_index].cost();
        };

        FP interval_cost = 0.0;
        interval_cost += cost("SchoolClosure") * param_at("SchoolClosure");
        interval_cost += cost("HomeOffice") * param_at("HomeOffice");
        interval_cost += cost("PhysicalDistancingSchool") * param_at("PhysicalDistancingSchool");
        interval_cost += cost("PhysicalDistancingWork") * param_at("PhysicalDistancingWork");
        interval_cost += cost("PhysicalDistancingOther") * param_at("PhysicalDistancingOther");

        objective += interval_cost / settings.num_intervals();
    }

    return objective;
}
