/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn
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
#ifndef SEIR_PARAMETERS_H
#define SEIR_PARAMETERS_H

#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/parameter_set.h"

#include <vector>

namespace mio
{
namespace oseir
{

    /*******************************************
      * Define Parameters of the SEIR model *
    *******************************************/

    struct TransmissionRisk {
        using Type = double;
        static constexpr Type get_default()
        {
            return 1.0;
        }
    };
    struct StageTimeIncubationInv {
        using Type = double;
        static constexpr Type get_default()
        {
            return 1.0 / 5.2;
        }
    };
    struct StageTimeInfectiousInv {
        using Type = double;
        static constexpr Type get_default()
        {
            return 1.0 / 6.0;
        }
    };
    struct ContactFrequency {
        using Type = ContactMatrix;
        static Type get_default()
        {
            return ContactMatrix{1};
        }
    };

    using ParametersBase =
        ParameterSet<TransmissionRisk, StageTimeIncubationInv, StageTimeInfectiousInv, ContactFrequency>;

} // namespace oseir
} // namespace mio

#endif // SEIR_PARAMETERS_H
