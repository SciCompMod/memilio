/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Lena Ploetzke
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

#ifndef IDE_SEIR_PARAMS_H
#define IDE_SEIR_PARAMS_H

#include "memilio/math/eigen.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/parameter_set.h"

namespace mio
{
namespace iseir
{
    /****************************************************
    * define some parameters used in the IDE SEIR model *
    *****************************************************/

    /**
    * @brief The time of latency used in the IDE SEIR model.
    * 
    * Latency time is the average time from infection to the onset of infectivity used in the model.
    * Latency time is measured in days.
    */
    struct LatencyTime {
        using Type = double;
        static constexpr Type get_default()
        {
            return 3.3;
        }
    };

    /**
    * @brief The infectious time used in the IDE SEIR model.
    * 
    * Infectious time is the average time from onset of infectivity to recovery used in the model.
    * Infectious time is measured in days.
    */
    struct InfectiousTime {
        using Type = double;
        static constexpr Type get_default()
        {
            return 8.2;
        }
    };

    /**
    * @brief The risk of transmission in the event of a contact used in the IDE SEIR model.
    * 
    * The transmission risk is the average risk to get infected in the event of a contact, 
    * given that the contact takes place between a susceptible and an infected person.
    */
    struct TransmissionRisk {
        using Type = double;
        static constexpr Type get_default()
        {
            return 0.1;
        }
    };

    /**
    * @brief The contact frequency is modeled using an UncertainContactMatrix.
    * 
    * The contact frequency is the average number of contact of an individual per day.
    * We use the type UncertainContactMatrix, because of the Randomness in this variable.
    * Via this parameter, dampings can be included to simulate non-pharmaceutical interventions.
    */
    struct ContactFrequency {
        using Type = UncertainContactMatrix;
        static Type get_default()
        {
            ContactMatrixGroup contact_matrix = ContactMatrixGroup(1, 1);
            contact_matrix[0]                 = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
            return Type(contact_matrix);
        }
    };

    // Define Parameterset for IDE SEIR model.
    using ParametersBase = ParameterSet<TransmissionRisk, LatencyTime, InfectiousTime, ContactFrequency>;

} // namespace iseir
} // namespace mio

#endif
