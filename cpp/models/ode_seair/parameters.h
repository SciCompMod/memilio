/*
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/utils/logging.h"

namespace mio
{
namespace oseair
{

/****************************************
 * Define Parameters of the SEAIR model *
 ****************************************/

/**
 * @brief Social distancing.
 */
template <typename FP>
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
template <typename FP>
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
template <typename FP>
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
template <typename FP>
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
template <typename FP>
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
template <typename FP>
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
template <typename FP>
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
template <typename FP>
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

template <typename FP>
using ParametersBase =
    ParameterSet<SocialDistancing<FP>, Quarantined<FP>, TestingRate<FP>, RecoveryRate<FP>, DeathRate<FP>,
                 TimeExposed<FP>, RecoveryRateFromAsymptomatic<FP>, TimeRecoveredInv<FP>>;

/**
 * @brief Parameters of an SEAIR model.
 */
template <typename FP>
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
     * @return Returns true if one constraint is not satisfied, otherwise false. 
     */
    bool check_constraints() const
    {
        if (this->template get<SocialDistancing<FP>>() < 0.0) {
            log_error("Constraint check: Parameter SocialDistancing smaller {}", 0);
            return true;
        }

        if (this->template get<Quarantined<FP>>() < 0.0) {
            log_error("Constraint check: Parameter Quarantined smaller {}", 0);
            return true;
        }

        const FP tol_times = 1e-1; // accepted tolerance for compartment stays
        if (this->template get<TimeExposed<FP>>() < tol_times) {
            log_error("Constraint check: Parameter TimeExposed {} smaller {}. Please "
                      "note that unreasonably small compartment stays lead to massively increased run time. "
                      "Consider to cancel and reset parameters.",
                      this->template get<TimeExposed<FP>>(), tol_times);
            return true;
        }

        if (this->template get<RecoveryRateFromAsymptomatic<FP>>() < 0.0) {
            log_error("Constraint check: Parameter RecoveryRateFromAsymptomatic smaller {}", 0);
            return true;
        }

        if (this->template get<TestingRate<FP>>() < 0.0) {
            log_error("Constraint check: Parameter TestingRate smaller {}", 0);
            return true;
        }

        if (this->template get<RecoveryRate<FP>>() < 0.0) {
            log_error("Constraint check: Parameter RecoveryRate smaller {}", 0);
            return true;
        }

        if (this->template get<DeathRate<FP>>() < 0.0) {
            log_error("Constraint check: Parameter DeathRate smaller {}", 0);
            return true;
        }

        if (this->template get<TimeRecoveredInv<FP>>() < 0.0) {
            log_error("Constraint check: Parameter TimeRecoveredInv smaller {}", 0);
            return true;
        }

        return false;
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
