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
void set_dynamic_NPIs(const SecirvvsOptimization& settings, mio::osecirvvs::Model<FP>& model,
                      const std::vector<FP>& dynamic_NPI_strengths)
{
    mio::DynamicNPIs<FP>& dynamic_npis =
        model.parameters.template get<mio::osecirvvs::DynamicNPIsInfectedSymptoms<FP>>();

    int index = 0;
    for (const auto& cp : settings.control_parameters()) {
        FP threshold                                   = cp.threshold();
        std::vector<mio::DampingSampling<FP>> dampings = cp.get_dampings<FP>();
        for (auto& d : dampings) {
            // Scale the damping value by the control strength
            FP strength = dynamic_NPI_strengths[index++];
            d.set_value(strength * d.get_value());
        }

        dynamic_npis.set_threshold(threshold, dampings);
    }
}
