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

#include "memilio/utils/parameter_set.h"
#include "ide_secir/infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/smoother.h"

#include <vector>

namespace mio
{
namespace isecir
{

/**********************************************
* Define Parameters of the IDE-SECIHURD model *
**********************************************/

/**
* @brief Function describing the time spent in a compartment before transiting to next compartment.

* This class defines a function that specifies which proportion of individuals is still in the compartment after a certain infection_age 
* (i.e. time after entering the compartment) and has not yet progressed to the next.
* See also function \gamma in Overleaf.
* Currently, you can only use a smoother_cosine() function with different parameters for this purpose.
* TransitionDistribution implements a vector with DelayDistribution%s for each compartment.
*/
struct DelayDistribution {
    /**
    * @brief Default constructor of the class. Default for m_max_support is 2.0 (relatively random.)
    */
    DelayDistribution()
        : m_max_support{2.0}
    {
    }

    /**
     * @brief Construct a new DelayDistribution object
     * 
     * @param init_max_support specifies the right bound and therefore the support of the function.
     */
    DelayDistribution(ScalarType init_max_support)
        : m_max_support{init_max_support}
    {
    }
    /**
     * @brief DelayDistribution is currently defined as a smoothed cosine function.
     * 
     * Used function goes through points (0,1) and (m_max_suppor,0) and is interpolated in between using a smoothed cosine function.
     * 
     * @param infection_age time at which the function should be evaluated
     * @return ScalarType evaluation of the smoother cosine function
     */
    ScalarType Distribution(ScalarType infection_age)
    {
        return smoother_cosine(infection_age, 0.0, m_max_support, 1.0, 0.0);
    }

    /**
     * @brief Get the m_max_support object
     * 
     * Can be used to access the m_max_support object, which specifies the right bound of the support of the function.
     * 
     * @return ScalarType 
     */
    ScalarType get_max_support() const
    {
        return m_max_support;
    }

    /**
     * @brief Set the m_max_support object
     * 
     * Can be used to set the m_max_support object, which specifies the right bound of the support of the function.
     * 
     * @return ScalarType 
     */
    void set_max_support(ScalarType new_max_support)
    {
        m_max_support = new_max_support;
    }

private:
    ScalarType m_max_support; ///< specifies the right bound of the support of the DelayDistribution.
};

/**
 * @brief Transition distribution for each transition in InfectionTransition.
 * 
 * Note that Distribution from S->E is just a dummy.
 * This transition is calculated in a different way.
 * 
 */
struct TransitionDistributions {
    using Type = std::vector<DelayDistribution>;
    static Type get_default()
    {
        return std::vector<DelayDistribution>((int)InfectionTransition::Count, DelayDistribution());
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
    /*For consistency, also define TransitionProbabilities for each transition in InfectionTransition. 
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
// TODO: Dependeny on infection age
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

struct ExponentialDecay {
    ExponentialDecay()
        : funcparam{1.0}
    {
    }

    ExponentialDecay(ScalarType init_funcparam)
        : funcparam{init_funcparam}
    {
    }

    ScalarType Function(ScalarType infection_age)
    {
        return std::exp(-funcparam * infection_age);
    }

    ScalarType get_funcparam()
    {
        return funcparam;
    }

private:
    ScalarType funcparam{};
};
/**
* @brief Probability of getting infected from a contact.
*/
template <class TransmissionProbabilityDecayFunction>
struct TransmissionProbabilityOnContact {
    // corresponds to rho, depends on infection_age
    using Type = TransmissionProbabilityDecayFunction;
    static Type get_default()
    {
        return TransmissionProbabilityDecayFunction();
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContact";
    }
};

/**
* @brief The relative InfectedNoSymptoms infectability.
*/
template <class TransmissionProbabilityDecayFunction>
struct RelativeTransmissionNoSymptoms {
    // correspond to xi_C, depends on infection_age
    using Type = TransmissionProbabilityDecayFunction;
    static Type get_default()
    {
        return TransmissionProbabilityDecayFunction();
    }
    static std::string name()
    {
        return "RelativeTransmissionNoSymptoms";
    }
};

/**
* @brief The risk of infection from symptomatic cases in the SECIR model.
*/
template <class TransmissionProbabilityDecayFunction>
struct RiskOfInfectionFromSymptomatic {
    // corresponds to xi_I, depends on infection_age
    using Type = TransmissionProbabilityDecayFunction;
    static Type get_default()
    {
        return TransmissionProbabilityDecayFunction();
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSymptomatic";
    }
};

// Define Parameterset for IDE SECIR model.
using ParametersBase =
    ParameterSet<TransitionDistributions, TransitionProbabilities, ContactPatterns, TransmissionProbabilityOnContact<mio::isecir::ExponentialDecay>,
                 RelativeTransmissionNoSymptoms<mio::isecir::ExponentialDecay>, RiskOfInfectionFromSymptomatic<mio::isecir::ExponentialDecay>>;

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
    // Define Parameterset for IDE SEIR model.
    /*using ParametersBase =
        ParameterSet<TransitionDistributions, TransitionProbabilities, ContactPatterns,
                     TransmissionProbabilityOnContact<mio::isecir::ExponentialDecay>,
                     RelativeTransmissionNoSymptoms<mio::isecir::ExponentialDecay>,
                     RiskOfInfectionFromSymptomatic<mio::isecir::ExponentialDecay>>;*/

} // namespace isecir
} // namespace mio

#endif // IDE_SECIR_PARAMS_H