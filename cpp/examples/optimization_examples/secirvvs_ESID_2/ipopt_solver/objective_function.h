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
        parameters[i] = ptr_parameters[i];
    }

    set_control_dampings<FP>(settings, model, parameters);

    auto integrator = make_integrator<FP>(settings.integrator_type(), settings.dt());

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());

    for (size_t interval = 0; interval < settings.num_intervals(); interval++) {

        size_t control_interval = interval / settings.pc_resolution();
        auto param_at           = [&](const std::string& name) {
            size_t control_index = static_cast<size_t>(string_to_control(name));
            return parameters[control_index + control_interval * settings.num_control_parameters()];
        };
        FP school_closure             = param_at("SchoolClosure");
        FP home_office                = param_at("HomeOffice");
        FP physical_distancing_school = param_at("PhysicalDistancingSchool");
        FP physical_distancing_work   = param_at("PhysicalDistancingWork");
        FP physical_distancing_other  = param_at("PhysicalDistancingOther");

        mio::TimeSeries<FP> result = mio::simulate<FP, mio::osecirvvs::Model<FP>>(
            time_steps[interval], time_steps[interval + 1], settings.dt(), model, integrator);
        const auto& final_state = result.get_last_value();

        for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
            for (size_t state_index = 0; state_index < num_infection_states(); state_index++) {
                size_t idx = age_group.get() * num_infection_states() + state_index;
                model.populations[{age_group, mio::osecirvvs::InfectionState(state_index)}] = final_state[idx];
            }
        }
        objective += (school_closure + home_office + physical_distancing_school + physical_distancing_work +
                      physical_distancing_other) /
                     settings.num_intervals();
    }

    return objective;
}
