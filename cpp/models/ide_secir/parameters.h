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
* @brief Distribution of the time spent in a compartment before transiting to next compartment.
*/
struct DelayDistribution {
    DelayDistribution()
        : xright{2.0}
    {
    }
    DelayDistribution(ScalarType init_x_right)
        : xright{init_x_right}
    {
    }

    ScalarType Distribution(ScalarType infection_age)
    {
        return smoother_cosine(infection_age, 0.0, xright, 1.0, 0.0);
    }

    ScalarType get_xright()
    {
        return xright;
    }

private:
    ScalarType xright;
};

/**
 * @brief Transition distribution for each Transition in InfectionTransition.
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
 * @brief Parameters needed for TransitionDistributions.
 * 
 * Parameters stored in a vector for initialisation of the TransitionDistributions.
 * Currently, for each TransitionDistribution is only one paramter used (eg. xright).
 * For each possible Transition defined in InfectionTransition, there is exactly one parameter.
 * Note that for transition S -> E, this is just a dummy.
 */
struct TransitionParameters {
    using Type = std::vector<ScalarType>;
    static Type get_default()
    {
        return std::vector<ScalarType>((int)InfectionTransition::Count, 1.0);
    }

    static std::string name()
    {
        return "TransitionParameters";
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
    using Type = ScalarType;
    static Type get_default()
    {
        return 0.5;
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

// Define Parameterset for IDE SECIR model.
using ParametersBase =
    ParameterSet<TransitionDistributions, TransitionParameters, TransitionProbabilities, ContactPatterns,
                 TransmissionProbabilityOnContact, RelativeTransmissionNoSymptoms, RiskOfInfectionFromSymptomatic>;

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
     */
    void check_constraints() const
    {
        if (this->get<RiskOfInfectionFromSymptomatic>() < 0.0) {
            log_warning("Constraint check: Parameter RiskOfInfectionFromSymptomatic smaller {:d}", 0);
        }

        // TODO
    }

private:
    Parameters(ParametersBase&& base)
        : ParametersBase(std::move(base))
    {
    }

public:
    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Parameters> deserialize(IOContext& io)
    {
        BOOST_OUTCOME_TRY(base, ParametersBase::deserialize(io));
        return success(Parameters(std::move(base)));
    }
};

} // namespace isecir
} // namespace mio

#endif // IDE_SECIR_PARAMS_H