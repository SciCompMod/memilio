/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#ifndef LCT_SECIR_PARAMS_H
#define LCT_SECIR_PARAMS_H

#include "memilio/utils/parameter_set.h"
#include "lct_secir/infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/smoother.h"

#include <vector>

namespace mio
{
namespace lsecir
{

/**********************************************
* Define Parameters of the LCT-SECIHURD model *
**********************************************/

/**
 * @brief Defines Number of Subcompartments for each compartment in InfectionState.
 */
struct SubcompartmentNumbers {
    /* For consistency, also define Number for Compartments where are no subcompartments possible. 
    For these States set number on 1.*/
    using Type = std::vector<int>;
    static Type get_default()
    {
        return std::vector<int>((int)InfectionState::Count, 1);
    }

    static std::string name()
    {
        return "TransitionProbabilities";
    }
};

/**
 * @brief The latency time in the SECIR model in day unit.
 */
struct TimeLatency {
    using Type = ScalarType;
    static Type get_default()
    {
        return 2.0;
    }
    static std::string name()
    {
        return "TimeLatency";
    }
};

/**
 * @brief Average time spent in the TimeInfectedNoSymptoms before developing Symptoms or recover in the SECIR model in day unit.
 */
struct TimeInfectedNoSymptoms {
    using Type = ScalarType;
    static Type get_default()
    {
        return 1.0;
    }
    static std::string name()
    {
        return "TimeInfectedNoSymptoms";
    }
};

/**
 * @brief Average time spent in the TimeInfectedSymptoms before going to Hospital or recover in the SECIR model in day unit.
 */
struct TimeInfectedNoSymptoms {
    using Type = ScalarType;
    static Type get_default()
    {
        return 1.5;
    }
    static std::string name()
    {
        return "TimeInfectedNoSymptoms";
    }
};

/**
 * @brief Average time being in the Hospital before treated by ICU or recover in the SECIR model in day unit.
 */
struct TimeInfectedSevere {
    using Type = ScalarType;
    static Type get_default()
    {
        return 1.0;
    }
    static std::string name()
    {
        return "TimeInfectedSevere";
    }
};

/**
 * @brief Average time treated by ICU before dead or recover in the SECIR model in day unit.
 */
struct TimeInfectedCritical {
    using Type = ScalarType;
    static Type get_default()
    {
        return 1.0;
    }
    static std::string name()
    {
        return "TimeInfectedCritical";
    }
};

/**
* @brief Probability of getting infected from a contact.
*/
struct TransmissionProbabilityOnContact {
    using Type = ScalarType;
    static Type get_default()
    {
        return 1.0;
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContact";
    }
};

/**
 * @brief The contact patterns within the society are modelled using an UncertainContactMatrix.
 */
struct ContactPatterns {
    using Type = UncertainContactMatrix;

    static Type get_default()
    {
        ContactMatrixGroup contact_matrix = ContactMatrixGroup(1, 1);
        contact_matrix[0]                 = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
        return Type(contact_matrix);
    }
    static std::string name()
    {
        return "ContactPatterns";
    }
};

/**
* @brief The relative InfectedNoSymptoms infectability.
*/
struct RelativeTransmissionNoSymptoms {
    using Type = ScalarType;
    static Type get_default()
    {
        return 0.5;
    }
    static std::string name()
    {
        return "RelativeTransmissionNoSymptoms";
    }
};

/**
* @brief The risk of infection from symptomatic cases in the SECIR model.
*/
struct RiskOfInfectionFromSymptomatic {
    using Type = ScalarType;
    static Type get_default()
    {
        return 0.5;
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSymptomatic";
    }
};

/**
* @brief The percentage of asymptomatic cases in the SECIR model.
*/
struct RecoveredPerInfectedNoSymptoms {
    using Type = ScalarType;
    static Type get_default()
    {
        return 0.5;
    }
    static std::string name()
    {
        return "RecoveredPerInfectedNoSymptoms";
    }
};

/**
* @brief The percentage of hospitalized patients per infected patients in the SECIR model.
*/
struct SeverePerInfectedSymptoms {
    using Type = ScalarType;
    static Type get_default()
    {
        return 0.5;
    }
    static std::string name()
    {
        return "SeverePerInfectedSymptoms";
    }
};

/**
* @brief The percentage of ICU patients per hospitalized patients in the SECIR model
*/
struct CriticalPerSevere {
    using Type = ScalarType;
    static Type get_default()
    {
        return 0.5;
    }
    static std::string name()
    {
        return "CriticalPerSevere";
    }
};

/**
* @brief The percentage of dead patients per ICU patients in the SECIR model
*/
struct DeathsPerCritical {
    using Type = ScalarType;
    static Type get_default()
    {
        return 0.1;
    }
    static std::string name()
    {
        return "DeathsPerCritical";
    }
};

using ParametersBase =
    ParameterSet<SubcompartmentNumbers, TimeLatency, TimeInfectedNoSymptoms, TimeInfectedNoSymptoms, TimeInfectedSevere,
                 TimeInfectedCritical, TransmissionProbabilityOnContact, ContactPatterns,
                 RelativeTransmissionNoSymptoms, RiskOfInfectionFromSymptomatic, RecoveredPerInfectedNoSymptoms,
                 SeverePerInfectedSymptoms, CriticalPerSevere, DeathsPerCritical>;
                 // TODO: write constraint check for parameters 

} // namespace lsecir
} // namespace mio


#endif // LCT_SECIR_PARAMS_H