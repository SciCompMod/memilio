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
        : xright{2.0}
    {
    }
    DelayDistribution(ScalarType init_x_right)
        : xright{init_x_right}
    {
    }

    ScalarType Distribution(ScalarType tau)
    {
        return smoother_cosine(tau, 0.0, xright, 1.0, 0.0);
    }

    ScalarType get_xright()
    {
        return xright;
    }

    //private:
    ScalarType xright;
};

struct TransitionDistributions {
    //we need transition distributions for E->C, C->I, C->R, I->H, I->R, H->U,H->R, U->D,U->R (so 9 distributions)
    using Type = std::vector<DelayDistribution>;
    static Type get_default()
    {
        return std::vector<DelayDistribution>(9, DelayDistribution());
    }

    static std::string name()
    {
        return "TransitionDistributions";
    }
};

struct TransitionParameters {
    // we need to initialize transition distributions with some parameters (for now just use one parameter per distribution, i.e. xright)
    // to be able to define different transition distributions for each transition
    using Type = std::vector<ScalarType>;
    static Type get_default()
    {
        return std::vector<ScalarType>(9, 1.0);
    }

    static std::string name()
    {
        return "TransitionParameters";
    }
};

struct TransitionProbabilities {
    // we need probabilities for S->E, E->C, C->I, I->H, H->U, U->D (so 5 values) (from this we can deduce prob for c->R etc.)
    // for consistency we also define mu_S^E = \mu_E^C = 1 to be able to apply general formula to get population from flows
    // in total we need 6 values
    using Type = std::vector<ScalarType>;
    static Type get_default()
    {
        return std::vector<ScalarType>(6, 0.5);
    }

    static std::string name()
    {
        return "TransitionProbabilities";
    }
};

/**
 * @brief the contact patterns within the society are modelled using an UncertainContactMatrix
 */
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

/**
* @brief probability of getting infected from a contact
*/
struct TransmissionProbabilityOnContact {
    // TODO: Abhaengigkeit von tau (und t), entspricht rho
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
* @brief the relative InfectedNoSymptoms infectability
*/
struct RelativeTransmissionNoSymptoms {
    // TODO: Abhaengigkeit von tau (und t), entspricht xi_C
    using Type = ScalarType;
    static Type get_default()
    {
        return 1.0;
    }
    static std::string name()
    {
        return "RelativeTransmissionNoSymptoms";
    }
};

/**
* @brief the risk of infection from symptomatic cases in the SECIR model
*/
struct RiskOfInfectionFromSymptomatic {
    // TODO: Abhaengigkeit von tau (und t), entspricht xi_I
    using Type = ScalarType;
    static Type get_default()
    {
        return 1.0;
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSymptomatic";
    }
};

/**
* @brief risk of infection from symptomatic cases increases as test and trace capacity is exceeded.
*/
struct MaxRiskOfInfectionFromSymptomatic {
    // TODO: Dies irgendwie benutzen für abhaengigkeit von RiskOfInfectionFromSymptomatic
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
};

// Define Parameterset for IDE SEIR model.
using ParametersBase =
    ParameterSet<TransitionDistributions, TransitionParameters, TransitionProbabilities, ContactPatterns,
                 TransmissionProbabilityOnContact, RelativeTransmissionNoSymptoms, RiskOfInfectionFromSymptomatic>;

} // namespace isecir
} // namespace mio

#endif // IDE_SECIR_PARAMS_H