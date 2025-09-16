/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#ifndef ODE_MSEIRS4_PARAMETERS_H
#define ODE_MSEIRS4_PARAMETERS_H

#include "memilio/utils/parameter_set.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/utils/logging.h"
#include <string>

namespace mio
{
namespace omseirs4
{
template <typename FP = ScalarType>
struct BaseTransmissionRate { // b0
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "BaseTransmissionRate";
    }
};

template <typename FP = ScalarType>
struct SeasonalAmplitude { // b1 in [0,1]
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "SeasonalAmplitude";
    }
};

template <typename FP = ScalarType>
struct SeasonalPhase { // phi in radians
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "SeasonalPhase";
    }
};

template <typename FP = ScalarType>
struct NaturalBirthDeathRate { // mu
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "NaturalBirthDeathRate";
    }
};

template <typename FP = ScalarType>
struct LossMaternalImmunityRate { // xi
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "LossMaternalImmunityRate";
    }
};

template <typename FP = ScalarType>
struct ProgressionRate { // sigma, E->I
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "ProgressionRate";
    }
};

template <typename FP = ScalarType>
struct RecoveryRate { // nu, I->R
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "RecoveryRate";
    }
};

template <typename FP = ScalarType>
struct ImmunityWaningRate { // gamma, R -> S stages
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "ImmunityWaningRate";
    }
};

template <typename FP = ScalarType>
struct Beta2Factor { // f2 multiplier for beta1
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(0.5);
    }
    static std::string name()
    {
        return "Beta2Factor";
    }
};

template <typename FP = ScalarType>
struct Beta3Factor { // f3 multiplier for beta1
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(0.35);
    }
    static std::string name()
    {
        return "Beta3Factor";
    }
};

template <typename FP = ScalarType>
struct Beta4Factor { // f4 multiplier for beta1
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(0.25);
    }
    static std::string name()
    {
        return "Beta4Factor";
    }
};

template <typename FP = ScalarType>
using ParametersBase =
    ParameterSet<BaseTransmissionRate<FP>, SeasonalAmplitude<FP>, SeasonalPhase<FP>, NaturalBirthDeathRate<FP>,
                 LossMaternalImmunityRate<FP>, ProgressionRate<FP>, RecoveryRate<FP>, ImmunityWaningRate<FP>,
                 Beta2Factor<FP>, Beta3Factor<FP>, Beta4Factor<FP>>;

template <typename FP = ScalarType>
class Parameters : public ParametersBase<FP>
{
public:
    Parameters()
        : ParametersBase<FP>()
    {
    }

    bool apply_constraints()
    {
        bool corrected    = false;
        auto clamp_nonneg = [&](auto& v) {
            if (v < 0) {
                v         = 0;
                corrected = true;
            }
        };
        clamp_nonneg(this->template get<BaseTransmissionRate<FP>>());
        clamp_nonneg(this->template get<SeasonalAmplitude<FP>>());
        if (this->template get<SeasonalAmplitude<FP>>() > 1) {
            this->template get<SeasonalAmplitude<FP>>() = 1;
            corrected                                   = true;
        }
        clamp_nonneg(this->template get<NaturalBirthDeathRate<FP>>());
        clamp_nonneg(this->template get<LossMaternalImmunityRate<FP>>());
        clamp_nonneg(this->template get<ProgressionRate<FP>>());
        clamp_nonneg(this->template get<RecoveryRate<FP>>());
        clamp_nonneg(this->template get<ImmunityWaningRate<FP>>());
        clamp_nonneg(this->template get<Beta2Factor<FP>>());
        clamp_nonneg(this->template get<Beta3Factor<FP>>());
        clamp_nonneg(this->template get<Beta4Factor<FP>>());
        return corrected;
    }

    bool check_constraints() const
    {
        if (this->template get<SeasonalAmplitude<FP>>() < 0 || this->template get<SeasonalAmplitude<FP>>() > 1)
            return true;
        return false;
    }
};

} // namespace omseirs4
} // namespace mio

#endif // ODE_MSEIRS4_PARAMETERS_H
