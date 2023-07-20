/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Anna Wendler, Lena Ploetzke
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
#ifndef IDE_SECIR_PARAMS_H
#define IDE_SECIR_PARAMS_H

#include "memilio/config.h"
#include "memilio/utils/parameter_set.h"
#include "ide_secir/infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/math/smoother.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/epidemiology/uncertain_matrix.h"

#include <memory>
#include <vector>

namespace mio
{
namespace isecir
{

/**********************************************
* Define Parameters of the IDE-SECIHURD model *
**********************************************/

/**
 * @brief Transition distribution for each transition in #InfectionTransition.
 *
 * As a default we use SmootherCosine functions for all transitions with m_parameter=2.
 */
struct TransitionDistributions {

    using Type = std::vector<StateAgeFunctionWrapper>;
    static Type get_default()
    {
        SmootherCosine smoothcos(2.0);
        StateAgeFunctionWrapper delaydistribution(smoothcos);
        return std::vector<StateAgeFunctionWrapper>((int)InfectionTransition::Count, delaydistribution);
    }

    static std::string name()
    {
        return "TransitionDistributions";
    }
};

/**
 * @brief Defines the probability for each possible transition to take this flow/transition.
 */
struct TransitionProbabilities {
    /*For consistency, also define TransitionProbabilities for each transition in #InfectionTransition. 
    Transition Probabilities should be set to 1 if there is no possible other flow from starting compartment.*/
    using Type = std::vector<ScalarType>;
    static Type get_default()
    {
        std::vector<ScalarType> probs((int)InfectionTransition::Count, 0.5);
        // Set the following probablities to 1 as there is no other option to go anywhere else.
        probs[Eigen::Index(InfectionTransition::SusceptibleToExposed)]        = 1;
        probs[Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
        return probs;
    }

    static std::string name()
    {
        return "TransitionProbabilities";
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
* @brief Probability of getting infected from a contact.
*/
struct TransmissionProbabilityOnContact {
    using Type = StateAgeFunctionWrapper;
    static Type get_default()
    {
        ConstantFunction constfunc(1.0);
        return StateAgeFunctionWrapper(constfunc);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContact";
    }
};

/**
* @brief The relative InfectedNoSymptoms infectability.
*/
struct RelativeTransmissionNoSymptoms {
    using Type = StateAgeFunctionWrapper;
    static Type get_default()
    {
        ConstantFunction constfunc(1.0);
        return StateAgeFunctionWrapper(constfunc);
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
    using Type = StateAgeFunctionWrapper;
    static Type get_default()
    {
        ConstantFunction constfunc(1.0);
        return StateAgeFunctionWrapper(constfunc);
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSymptomatic";
    }
};

// Define Parameterset for IDE SECIR model.
using ParametersBase =
    ParameterSet<TransitionDistributions, TransitionProbabilities, ContactPatterns, TransmissionProbabilityOnContact,
                 RelativeTransmissionNoSymptoms, RiskOfInfectionFromSymptomatic>;

/**
 * @brief Parameters of an age-resolved SECIR/SECIHURD model.
 */
class Parameters : public ParametersBase
{
public:
    Parameters()
        : ParametersBase()
    {
    }

    /**
     * @brief checks whether all Parameters satisfy their corresponding constraints and throws errors, if they do not
     * @return Returns 1 if one constraint is not satisfied, otherwise 0.
     */
    int check_constraints() const
    {
        for (int i = 0; i < 20; i++) {
            std::cout << this->get<TransmissionProbabilityOnContact>().eval(i) << "\n";
            if (this->get<TransmissionProbabilityOnContact>().eval(i) < 0.0 ||
                this->get<TransmissionProbabilityOnContact>().eval(i) > 1.0) {
                log_error("Constraint check: TransmissionProbabilityOnContact(i) smaller {:d} or larger {:d} at some "
                          "time {:d}",
                          0, 1, i);
                return 1;
            }
        }

        for (int i = 0; i < 20; i++) {
            if (this->get<RelativeTransmissionNoSymptoms>().eval(i) < 0.0 ||
                this->get<RelativeTransmissionNoSymptoms>().eval(i) > 1.0) {
                log_error("Constraint check: TransmissionProbabilityOnContact(i) smaller {:d} or larger {:d} at some "
                          "time {:d}",
                          0, 1, i);
                return 1;
            }
        }

        for (int i = 0; i < 20; i++) {
            if (this->get<RiskOfInfectionFromSymptomatic>().eval(i) < 0.0 ||
                this->get<RiskOfInfectionFromSymptomatic>().eval(i) > 1.0) {
                log_error("Constraint check: TransmissionProbabilityOnContact(i) smaller {:d} or larger {:d} at some "
                          "time {:d}",
                          0, 1, i);
                return 1;
            }
        }

        for (int i = 0; i < (int)InfectionTransition::Count; i++) {
            if (this->get<TransitionProbabilities>()[i] < 0.0 || this->get<TransitionProbabilities>()[i] > 1.0) {
                log_error("Constraint check: One parameter TransitionProbabilities smaller {:d} or larger {:d}", 0, 1);
                return 1;
            }
        }

        if (this->get<TransitionProbabilities>()[(int)InfectionTransition::SusceptibleToExposed] != 1.0) {
            log_error("Constraint check: Parameter transitiion probability for SusceptibleToExposed unequal to {:d}",
                      1);
            return 1;
        }

        if (this->get<TransitionProbabilities>()[(int)InfectionTransition::ExposedToInfectedNoSymptoms] != 1.0) {
            log_error(
                "Constraint check: Parameter transitiion probability for ExposedToInfectedNoSymptoms unequal to {:d}",
                1);
            return 1;
        }

        if (this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] +
                this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedNoSymptomsToRecovered] !=
            1.0) {
            log_error("Constraint check: Sum of transitiion probability for InfectedNoSymptomsToInfectedSymptoms and "
                      "InfectedNoSymptomsToRecovered unequal to {:d}",
                      1);
            return 1;
        }

        if (this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedSymptomsToInfectedSevere] +
                this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedSymptomsToRecovered] !=
            1.0) {
            log_error("Constraint check: Sum of transitiion probability for InfectedSymptomsToInfectedSevere and "
                      "InfectedSymptomsToRecovered unequal to {:d}",
                      1);
            return 1;
        }

        if (this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedSevereToInfectedCritical] +
                this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedSevereToRecovered] !=
            1.0) {
            log_error("Constraint check: Sum of transitiion probability for InfectedSevereToInfectedCritical and "
                      "InfectedSevereToRecovered unequal to {:d}",
                      1);
            return 1;
        }

        if (this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedCriticalToDead] +
                this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedCriticalToRecovered] !=
            1.0) {
            log_error("Constraint check: Sum of transitiion probability for InfectedCriticalToDead and "
                      "InfectedCriticalToRecovered unequal to {:d}",
                      1);
            return 1;
        }

        return 0;
    }

private:
    Parameters(ParametersBase&& base)
        : ParametersBase(std::move(base))
    {
    }
};

} // namespace isecir
} // namespace mio

#endif // IDE_SECIR_PARAMS_H