#pragma once

#include "control_parameters.h"

#include "models/ode_secirvvs/model.h"

#include "../optimization_settings/optimization_settings.h"
#include "../helpers/make_time_grid.h"

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cstddef>
#include <functional>

template <typename FP>
void set_control_dampings(const SecirvvsOptimization& settings, mio::osecirvvs::Model<FP>& model,
                          const std::vector<FP>& control_parameter_values)
{
    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());

    mio::UncertainContactMatrix<FP>& contacts = model.parameters.template get<mio::osecirvvs::ContactPatterns<FP>>();
    std::vector<mio::DampingSampling<FP>>& contact_dampings = contacts.get_dampings();

    for (size_t control_interval = 0; control_interval < settings.num_control_intervals(); control_interval++) {
        mio::SimulationTime<FP> time(time_steps[control_interval * settings.pc_resolution()]);

        for (size_t i = 0; i < settings.num_control_parameters(); i++) {
            FP control_value = control_parameter_values[control_interval * settings.num_control_parameters() + i];

            mio::DampingSampling<FP> damping =
                settings.control_parameters()[i].make_damping_sampling<FP>(control_value, time);

            contact_dampings.push_back(damping);
        }
    }
    contacts.make_matrix();
}
