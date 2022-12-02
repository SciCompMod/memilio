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
struct TransmissionProbabilityOnContact {
    using Type = UncertainValue;
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
struct TimeExposed {
    using Type = UncertainValue;
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
struct TimeInfected {
    using Type = UncertainValue;
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

using ParametersBase = ParameterSet<TransmissionProbabilityOnContact, TimeExposed, TimeInfected, ContactPatterns>;

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
        if (this->get<TimeExposed>() <= 0.0) {
            log_error("Constraint check: Parameter TimeExposed {:.4f} smaller or equal {:.4f}",
                      this->get<TimeExposed>(), 0.0);
            exit(EXIT_FAILURE);
        }
        if (this->get<TimeInfected>() <= 0.0) {
            log_error("Constraint check: Parameter TimeInfected {:.4f} smaller or equal {:.4f}",
                      this->get<TimeInfected>(), 0.0);
            exit(EXIT_FAILURE);
        }
        if (this->get<TransmissionProbabilityOnContact>() < 0.0 ||
            this->get<TransmissionProbabilityOnContact>() > 1.0) {
            log_error(
                "Constraint check: Parameter TransmissionProbabilityOnContact {:.4f} smaller {:.4f} or greater {:.4f}",
                this->get<TransmissionProbabilityOnContact>(), 0.0, 1.0);
            exit(EXIT_FAILURE);
        }
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

} // namespace oseir
} // namespace mio

#endif // SEIR_PARAMETERS_H
