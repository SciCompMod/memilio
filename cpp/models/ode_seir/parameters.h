/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn
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
#ifndef SEIR_PARAMETERS_H
#define SEIR_PARAMETERS_H

#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/utils/parameter_set.h"

#include <vector>

namespace mio
{
namespace oseir
{

/*******************************************
      * Define Parameters of the SEIR model *
    *******************************************/

/**
     * @brief probability of getting infected from a contact
     */
template<typename FP=double>
struct TransmissionProbabilityOnContact {
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContact";
    }
};

/**
     * @brief the latent time in day unit
     */
template<typename FP=double>
struct TimeExposed {
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(5.2);
    }
    static std::string name()
    {
        return "TimeExposed";
    }
};

/**
     * @brief the infectious time in day unit
     */
template<typename FP=double>
struct TimeInfected {
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(6.0);
    }
    static std::string name()
    {
        return "TimeInfected";
    }
};

/**
     * @brief the contact patterns within the society are modelled using a ContactMatrix
     */
struct ContactPatterns {
    using Type = ContactMatrix;
    static Type get_default()
    {
        return Type{1};
    }
    static std::string name()
    {
        return "ContactPatterns";
    }
};

template<typename FP=double>
using ParametersBase = ParameterSet<TransmissionProbabilityOnContact<FP>, TimeExposed<FP>,
                                    TimeInfected<FP>, ContactPatterns>;

/**
 * @brief Parameters of an age-resolved SECIR/SECIHURD model.
 */
template<typename FP=double>
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
        if (this->template get<TimeExposed<FP>>() <= 0.0) {
            log_error("Constraint check: Parameter TimeExposed {:.4f} smaller or equal {:.4f}",
                      this->template get<TimeExposed<FP>>(), 0.0);
            return 1;
        }
        if (this->template get<TimeInfected<FP>>() <= 0.0) {
            log_error("Constraint check: Parameter TimeInfected {:.4f} smaller or equal {:.4f}",
                      this->template get<TimeInfected<FP>>(), 0.0);
            return 1;
        }
        if (this->template get<TransmissionProbabilityOnContact<FP>>() < 0.0 ||
            this->template get<TransmissionProbabilityOnContact<FP>>() > 1.0) {
            log_error(
                "Constraint check: Parameter TransmissionProbabilityOnContact {:.4f} smaller {:.4f} or greater {:.4f}",
                this->template get<TransmissionProbabilityOnContact<FP>>(), 0.0, 1.0);
            return 1;
        }
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
        BOOST_OUTCOME_TRY(base, ParametersBase<FP>::deserialize(io));
        return success(Parameters(std::move(base)));
    }
};

} // namespace oseir
} // namespace mio

#endif // SEIR_PARAMETERS_H
