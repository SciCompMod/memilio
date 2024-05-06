/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Ralf Hannemann-Tamas
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

#ifndef ODESEAIR_PARAMETERS_H
#define ODESEAIR_PARAMETERS_H

#include "memilio/utils/parameter_set.h"

namespace mio
{
namespace oseair
{

/*******************************************
* Define Parameters of the SEAIR model     *
*******************************************/

/* This model is a extented SEIR type model of the COVID-19 pandemic in the US
 * that als includes asymptomatic and dead people.
 * A detailed description of the model can be found in the publication
 * Tsay et al. (2020), Modeling, state estimation, and optimal control for the US COVID-19 outbreak */

/**
 * @brief Social distancing.
 */
template <typename FP = double>
struct SocialDistancing {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.2);
    }
    static std::string name()
    {
        return "SocialDistancing";
    }
};

/**
 * @brief Quarantining.
 */
template <typename FP = double>
struct Quarantined {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.2);
    }
    static std::string name()
    {
        return "Quarantined";
    }
};

/**
 * @brief Rate of testing.
 */
template <typename FP = double>
struct TestingRate {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.2);
    }
    static std::string name()
    {
        return "TestingRate";
    }
};

/**
 * @brief Recovery rate.
 */
template <typename FP = double>
struct RecoveryRate {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.0067);
    }
    static std::string name()
    {
        return "RecoveryRate";
    }
};

/**
 * @brief Death Rate.
 */
template <typename FP = double>
struct DeathRate {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.0041);
    }
    static std::string name()
    {
        return "DeathRate";
    }
};

/**
 * @brief Inverse of the latent period of the virus.
 */
template <typename FP = double>
struct TimeExposed {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.5);
    }
    static std::string name()
    {
        return "TimeExposed";
    }
};

/**
 * @brief Infectious period for unconfirmed infected people.
 */
template <typename FP = double>
struct RecoveryRateFromAsymptomatic {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.1);
    }
    static std::string name()
    {
        return "RecoveryRateFromAsymptomatic";
    }
};

/**
 * @brief Rate recovered people become susceptible again.
 */
template <typename FP = double>
struct TimeRecoveredInv {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "TimeRecoveredInv";
    }
};

template <typename FP = double>
using ParametersBase =
    ParameterSet<SocialDistancing<FP>, Quarantined<FP>, TestingRate<FP>, RecoveryRate<FP>, DeathRate<FP>,
                 TimeExposed<FP>, RecoveryRateFromAsymptomatic<FP>, TimeRecoveredInv<FP>>;

/**
 * @brief Parameters of an SEAIR model.
 */
template <typename FP = double>
class Parameters : public ParametersBase<FP>
{
public:
    Parameters()
        : ParametersBase<FP>()
    {
    }

    /**
     * @brief Checks whether all Parameters satisfy their corresponding constraints and logs an error 
     * if constraints are not satisfied.
     * @return Returns 1 if one constraint is not satisfied, otherwise 0.   
     */
    int check_constraints() const
    {
        return 0;
    }

private:
    Parameters(ParametersBase<FP>&& base)
        : ParametersBase<FP>(std::move(base))
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
        BOOST_OUTCOME_TRY(auto&& base, ParametersBase<FP>::deserialize(io));
        return success(Parameters(std::move(base)));
    }
};

} // namespace oseair
} // namespace mio

#endif // ODESEAIR_PARAMETERS_H
