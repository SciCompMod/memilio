#pragma once

#include <vector>

#include "tools/optimal_control/control_parameters/damping_controls.h"

#include "tools/optimal_control/helpers/integrator_selector.h"
#include "tools/optimal_control/helpers/make_time_grid.h"


template <typename FP, class OptimizationSettings>
FP objective_function(const OptimizationSettings& settings, const typename OptimizationSettings::template ModelTemplate<FP>& model, const FP* ptr_parameters,
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

    auto integrator = make_integrator<FP>(settings.integrator_type(), settings.dt());

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());

    OptimizationSettings::template SimulationTemplate<FP> sim(model, settings.t0(), settings.dt());
    sim.set_integrator(integrator);
    
    set_control_dampings<FP, OptimizationSettings>(settings, sim.get_model(), parameters);

    for (size_t interval = 0; interval < settings.num_intervals(); interval++) {

        size_t control_interval = interval / settings.pc_resolution();

        sim.get_dt() = settings.dt();
        sim.advance(time_steps[interval + 1]);
        const auto& final_state = sim.get_result().get_last_value();

        // Add costs of all control parameters for current interval
        FP interval_cost = 0.0;
        for(size_t idx_control_parameter = 0; idx_control_parameter < settings.num_control_parameters(); idx_control_parameter++)
        {
            interval_cost += settings.control_parameters()[idx_control_parameter].cost() * parameters[idx_control_parameter + control_interval * settings.num_control_parameters()]
        }
        objective += interval_cost / settings.num_intervals();
    }

    return objective;
}
