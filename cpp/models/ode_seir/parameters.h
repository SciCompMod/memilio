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

#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/utils/parameter_set.h"

#include <vector>

namespace mio
{
namespace oseir
{

    /*******************************************
      * Define Parameters of the SEIR model *
    *******************************************/

    /**
     * @brief probability of getting infected from a contact
     */
    struct TransmissionProbabilityOnContact {
        using Type = UncertainValue;
        static Type get_default()
        {
            return Type(1.0);
        }
        static std::string name()
        {
            return "TransmissionProbabilityOnContact";
        }
    };

    /**
     * @brief the latent time in day unit
     */
    struct TimeExposed {
        using Type = UncertainValue;
        static Type get_default()
        {
            return Type(5.2);
        }
        static std::string name()
        {
            return "TimeExposed";
        }
    };

    /**
     * @brief the infectious time in day unit
     */
    struct TimeInfected {
        using Type = UncertainValue;
        static Type get_default()
        {
            return Type(6.0);
        }
        static std::string name()
        {
            return "TimeInfected";
        }
    };

    /**
     * @brief the contact patterns within the society are modelled using a ContactMatrix
     */
    struct ContactPatterns {
        using Type = ContactMatrix;
        static Type get_default()
        {
            return Type{1};
        }
        static std::string name()
        {
            return "ContactPatterns";
        }
    };

    using Parameters = ParameterSet<TransmissionProbabilityOnContact, TimeExposed, TimeInfected, ContactPatterns>;

} // namespace oseir
} // namespace mio

#endif // SEIR_PARAMETERS_H
