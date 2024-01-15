/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Martin J. Kuehn
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#include "ode_secir/parameter_space.h" // IWYU pragma: keep


<<<<<<< HEAD
=======
void set_params_distributions_normal(Model& model, double t0, double tmax, double dev_rel)
{
    auto set_distribution = [dev_rel](UncertainValue& v, double min_val = 0.001) {
        v.set_distribution(ParameterDistributionNormal(
            //add add limits for nonsense big values. Also mscv has a problem with a few doubles so this fixes it
            std::min(std::max(min_val, (1 - dev_rel * 2.6) * v), 0.1 * std::numeric_limits<double>::max()),
            std::min(std::max(min_val, (1 + dev_rel * 2.6) * v), 0.5 * std::numeric_limits<double>::max()),
            std::min(std::max(min_val, double(v)), 0.3 * std::numeric_limits<double>::max()),
            std::min(std::max(min_val, dev_rel * v), std::numeric_limits<double>::max())));
    };
>>>>>>> upstream/main

