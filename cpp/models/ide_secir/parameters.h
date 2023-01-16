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

/*******************************************
    * Define Parameters of the IDE-SECIHURD model *
    *******************************************/

/**
* @brief Distribution of the time spent in a compartment before transiting to next compartment.
*/
struct DelayDistribution {
    DelayDistribution()
        : xright{6.0}
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

    //private:
    ScalarType xright;
};

struct TransitionDistributions {
    /*Transition distribution for InfectionTransitions. Note that Distribution from S->E is just a dummy.
    This transition is calculated in a different way. */
    using Type = std::vector<DelayDistribution>;
    static Type get_default()
    {
        return std::vector<DelayDistribution>((int)InfectionTransitions::Count, DelayDistribution());
    }

    static std::string name()
    {
        return "TransitionDistributions";
    }
};

struct TransitionParameters {
    // we need to initialize transition distributions with some parameters (for now just use one parameter per distribution, i.e. xright)
    // to be able to define different transition distributions for each transition
    // Here also transition S-> E is just a dummy
    using Type = std::vector<ScalarType>;
    static Type get_default()
    {
        return std::vector<ScalarType>((int)InfectionTransitions::Count, 1.0);
    }

    static std::string name()
    {
        return "TransitionParameters";
    }
};

struct TransitionProbabilities {
    /*For consistency, also define TransitionProbabilities for each transition in InfectionTransitions.*/
    using Type = std::vector<ScalarType>;
    static Type get_default()
    {
        return std::vector<ScalarType>((int)InfectionTransitions::Count, 0.5);
    }

    static std::string name()
    {
        return "TransitionProbabilities";
    }
};

/**
 * @brief the contact patterns within the society are modelled using an UncertainContactMatrix
 */
 // TODO: Dependeny on infection age
struct ContactPatterns {
    using Type = UncertainContactMatrix;

    static Type get_default()
    {
        return Type(1, static_cast<Eigen::Index>((size_t)1));
    }
    static std::string name()
    {
        return "ContactPatterns";
    }
};

struct ExponentialDecay{
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
        return std::exp(- funcparam * infection_age);
    }

    ScalarType get_funcparam()
    {
        return funcparam;
    }
    
    private:
    ScalarType funcparam{};

 };
/**
* @brief probability of getting infected from a contact
*/
template <class TransmissionProbabilityDecayFunction>
struct TransmissionProbabilityOnContact {
    // TODO: Abhaengigkeit von tau (und t), entspricht rho
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
* @brief the relative InfectedNoSymptoms infectability
*/
template <class TransmissionProbabilityDecayFunction>
struct RelativeTransmissionNoSymptoms {
    // TODO: Abhaengigkeit von tau (und t), entspricht xi_C
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
* @brief the risk of infection from symptomatic cases in the SECIR model
*/
template <class TransmissionProbabilityDecayFunction>
struct RiskOfInfectionFromSymptomatic {
    // TODO: Abhaengigkeit von tau (und t), entspricht xi_I
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

/* @brief risk of infection from symptomatic cases increases as test and trace capacity is exceeded.
*/
/*Martin sagt das sollen wir aktuell mal weglassen
struct MaxRiskOfInfectionFromSymptomatic {
    // Dies könnte man irgendwie benutzen für abhaengigkeit von RiskOfInfectionFromSymptomatic
    //von t in Abhaengigkeit der Inzidenz wie im ODE-Modell, akteull nutzlos
    // evtl benoetigen wir noch die Parameter : TestAndTraceCapacity ,DynamicNPIsInfectedSymptoms
    using Type = ScalarType;
    static Type get_default()
    {
        return 0.0;
    }
    static std::string name()
    {
        return "MaxRiskOfInfectionFromSymptomatic";
    }
};*/

// Define Parameterset for IDE SEIR model.
using ParametersBase =
    ParameterSet<TransitionDistributions, TransitionParameters, TransitionProbabilities, ContactPatterns,
                 TransmissionProbabilityOnContact<mio::isecir::ExponentialDecay>, 
                 RelativeTransmissionNoSymptoms<mio::isecir::ExponentialDecay>, 
                 RiskOfInfectionFromSymptomatic<mio::isecir::ExponentialDecay>>;

} // namespace isecir
} // namespace mio

#endif // IDE_SECIR_PARAMS_H