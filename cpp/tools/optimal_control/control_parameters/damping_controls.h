#pragma once

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cstddef>
#include <functional>

#include "models/ode_secirvvs/parameters.h"

#include "tools/optimal_control/control_parameters/control_parameters.h"
#include "tools/optimal_control/helpers/make_time_grid.h"

template <typename FP>
mio::DampingSampling<FP> make_damping_from_control_parameter(mio::SimulationTime<FP> time, ControlParameter& control_parameter, FP value)
{

    return mio::DampingSampling<FP>(control_parameter.effectiveness() * value, time, control_parameter.m_damping().get_matrix_indices(), control_parameter.m_damping().get_group_weights());
}

template <typename FP, class OptimizationSettings>
void set_control_dampings(const OptimizationSettings& settings, const typename OptimizationSettings::template ModelTemplate<FP>& model,
                          const std::vector<FP>& control_parameters)
{
    assert(control_parameters.size() == settings.num_control_parameters() * settings.num_control_intervals());

    mio::UncertainContactMatrix<FP>& contacts = model.parameters.template get<ypename OptimizationSettings::template ContactPatterns<FP>>();
    std::vector<mio::DampingSampling<FP>>& contact_dampings = contacts.get_dampings();

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());

    for (size_t control_interval = 0; control_interval < settings.num_control_intervals(); control_interval++) {

        mio::SimulationTime<FP> time(time_steps[control_interval * settings.pc_resolution()]);
        
        for(size_t idx_control_parameter = 0; idx_control_parameter < settings.num_control_parameters(); idx_control_parameter++)
        {
            auto damping = make_damping_from_control_parameter<FP>(time, settings.control_parameters()[idx_control_parameter], control_parameters[idx_control_parameter + control_interval * settings.num_control_parameters()]);
            contact_dampings.push_back(damping)
        }
    }
    contacts.make_matrix();
}
