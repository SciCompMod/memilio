/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Nils Wassmuth, Rene Schmieding, Martin J. Kuehn
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

#ifndef MIO_SDE_SEIRVV_PARAMETERS_H
#define MIO_SDE_SEIRVV_PARAMETERS_H

#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/utils/parameter_set.h"

namespace mio
{
namespace sseirvv
{

/******************************************
 * Define Parameters of the SSEIRVV model *
 ******************************************/

/**
 * @brief Probability of getting infected from a contact 
 * with variant 1.
 */
struct TransmissionProbabilityOnContactV1 {
    using Type = UncertainValue<>;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactV1";
    }
};

/**
 * @brief Probability of getting infected from a contact 
 * with variant 2.
 */
struct TransmissionProbabilityOnContactV2 {
    using Type = UncertainValue<>;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContactV2";
    }
};

/**
 * @brief The latent time of variant 1 in days.
 */
struct TimeExposedV1 {
    using Type = UncertainValue<ScalarType>;
    static Type get_default()
    {
        return Type(6.0);
    }
    static std::string name()
    {
        return "TimeExposedV1";
    }
};

/**
 * @brief The latent time of variant 2 in days.
 */
struct TimeExposedV2 {
    using Type = UncertainValue<ScalarType>;
    static Type get_default()
    {
        return Type(6.0);
    }
    static std::string name()
    {
        return "TimeExposedV2";
    }
};

/**
 * @brief The infectious time of variant 1 in days.
 */
struct TimeInfectedV1 {
    using Type = UncertainValue<ScalarType>;
    static Type get_default()
    {
        return Type(6.0);
    }
    static std::string name()
    {
        return "TimeInfectedV1";
    }
};

/**
 * @brief The infectious time of variant 2 in days.
 */
struct TimeInfectedV2 {
    using Type = UncertainValue<ScalarType>;
    static Type get_default()
    {
        return Type(6.0);
    }
    static std::string name()
    {
        return "TimeInfectedV2";
    }
};

/**
 * @brief The contact patterns within the society are modelled using a ContactMatrix.
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

using ParametersBase = ParameterSet<TransmissionProbabilityOnContactV1, TransmissionProbabilityOnContactV2, 
    TimeExposedV1, TimeExposedV2, TimeInfectedV1, TimeInfectedV2, ContactPatterns>;

/**
 * @brief Parameters of stochastic SEIRVV model.
 */
class Parameters : public ParametersBase
{
public:
    Parameters()
        : ParametersBase()
    {
    }

    /**
     * @brief Checks whether all Parameters satisfy their corresponding constraints and applies them, if they do not.
     * Time spans cannot be negative and probabilities can only take values between [0,1]. 
     *
     * Attention: This function should be used with care. It is necessary for some test problems to run through quickly,
     *            but in a manual execution of an example, check_constraints() may be preferred. Note that the apply_constraints()
     *            function can and will not set Parameters to meaningful values in an epidemiological or virological context,
     *            as all models are designed to be transferable to multiple diseases. Consequently, only acceptable
     *            (like 0 or 1 for probabilities or small positive values for time spans) values are set here and a manual adaptation
     *            may often be necessary to have set meaningful values.
     *
     * @return Returns true if one ore more constraint were corrected, false otherwise.  
     */
    bool apply_constraints()
    {
        double tol_times = 1e-1;

        int corrected = false;
        if (this->get<TimeExposedV1>() < tol_times) {
            log_warning("Constraint check: Parameter TimeExposedV1 changed from {:.4f} to {:.4f}. Please note that "
                        "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                        "and reset parameters.",
                        this->get<TimeExposedV1>(), tol_times);
            this->get<TimeExposedV1>() = tol_times;
            corrected                  = true;
        }
        if (this->get<TimeExposedV2>() < tol_times) {
            log_warning("Constraint check: Parameter TimeExposedV2 changed from {:.4f} to {:.4f}. Please note that "
                        "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                        "and reset parameters.",
                        this->get<TimeExposedV2>(), tol_times);
            this->get<TimeExposedV2>() = tol_times;
            corrected                  = true;
        }
        if (this->get<TimeInfectedV1>() < tol_times) {
            log_warning("Constraint check: Parameter TimeInfectedV1 changed from {:.4f} to {:.4f}. Please note that "
                        "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                        "and reset parameters.",
                        this->get<TimeInfectedV1>(), tol_times);
            this->get<TimeInfectedV1>() = tol_times;
            corrected                   = true;
        }
        if (this->get<TimeInfectedV2>() < tol_times) {
            log_warning("Constraint check: Parameter TimeInfectedV2 changed from {:.4f} to {:.4f}. Please note that "
                        "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                        "and reset parameters.",
                        this->get<TimeInfectedV2>(), tol_times);
            this->get<TimeInfectedV2>() = tol_times;
            corrected                   = true;
        }
        if (this->get<TransmissionProbabilityOnContactV1>() < 0.0 ||
            this->get<TransmissionProbabilityOnContactV1>() > 1.0) {
            log_warning("Constraint check: Parameter TransmissionProbabilityOnContactV1 changed from {:0.4f} to {:d} ",
                        this->get<TransmissionProbabilityOnContactV1>(), 0.0);
            this->get<TransmissionProbabilityOnContactV1>() = 0.0;
            corrected                                       = true;
        }
        if (this->get<TransmissionProbabilityOnContactV2>() < 0.0 ||
            this->get<TransmissionProbabilityOnContactV2>() > 1.0) {
            log_warning("Constraint check: Parameter TransmissionProbabilityOnContactV2 changed from {:0.4f} to {:d} ",
                        this->get<TransmissionProbabilityOnContactV2>(), 0.0);
            this->get<TransmissionProbabilityOnContactV2>() = 0.0;
            corrected                                       = true;
        }
        return corrected;
    }

    /**
     * @brief Checks whether all Parameters satisfy their corresponding constraints and logs an error 
     * if constraints are not satisfied.
     * @return Returns true if one constraint is not satisfied, otherwise false.   
     */
    bool check_constraints() const
    {
        ScalarType tol_times = 1e-1;

        if (this->get<TimeExposedV1>() < tol_times) {
            log_error("Constraint check: Parameter TimeExposedV1 {:.4f} smaller or equal {:.4f}. Please note that "
                      "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                      "and reset parameters.",
                      this->get<TimeExposedV1>(), 0.0);
            return true;
        }
        if (this->get<TimeExposedV2>() < tol_times) {
            log_error("Constraint check: Parameter TimeExposedV2 {:.4f} smaller or equal {:.4f}. Please note that "
                      "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                      "and reset parameters.",
                      this->get<TimeExposedV2>(), 0.0);
            return true;
        }
        if (this->get<TimeInfectedV1>() < tol_times) {
            log_error("Constraint check: Parameter TimeInfectedV1 {:.4f} smaller or equal {:.4f}. Please note that "
                      "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                      "and reset parameters.",
                      this->get<TimeInfectedV1>(), 0.0);
            return true;
        }
        if (this->get<TimeInfectedV2>() < tol_times) {
            log_error("Constraint check: Parameter TimeInfectedV2 {:.4f} smaller or equal {:.4f}. Please note that "
                      "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                      "and reset parameters.",
                      this->get<TimeInfectedV2>(), 0.0);
            return true;
        }
        if (this->get<TransmissionProbabilityOnContactV1>() < 0.0 ||
            this->get<TransmissionProbabilityOnContactV1>() > 1.0) {
            log_error(
                "Constraint check: Parameter TransmissionProbabilityOnContactV1 {:.4f} smaller {:.4f} or greater {:.4f}",
                this->get<TransmissionProbabilityOnContactV1>(), 0.0, 1.0);
            return true;
        }
        if (this->get<TransmissionProbabilityOnContactV2>() < 0.0 ||
            this->get<TransmissionProbabilityOnContactV2>() > 1.0) {
            log_error(
                "Constraint check: Parameter TransmissionProbabilityOnContactV2 {:.4f} smaller {:.4f} or greater {:.4f}",
                this->get<TransmissionProbabilityOnContactV2>(), 0.0, 1.0);
            return true;
        }
        return false;
    }

private:
    Parameters(ParametersBase&& base)
        : ParametersBase(std::move(base))
    {
    }

public:
    /**
     * Deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Parameters> deserialize(IOContext& io)
    {
        BOOST_OUTCOME_TRY(auto&& base, ParametersBase::deserialize(io));
        return success(Parameters(std::move(base)));
    }
};

} // namespace sseirvv
} // namespace mio

#endif // MIO_SDE_SEIRVV_PARAMETERS_H
