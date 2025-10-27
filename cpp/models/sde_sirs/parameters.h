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

#ifndef MIO_SDE_SIRS_PARAMETERS_H
#define MIO_SDE_SIRS_PARAMETERS_H

#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/utils/parameter_set.h"

namespace mio
{
namespace ssirs
{

/***************************************
 * Define Parameters of the SIRS model *
 ***************************************/

/**
 * @brief probability of getting infected from a contact
 */
template <typename FP>
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
 * @brief the infectious time in day unit
 */
template <typename FP>
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
 * @brief the infectious time in day unit
 */
template <typename FP>
struct TimeImmune {
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(6.0);
    }
    static std::string name()
    {
        return "TimeImmune";
    }
};

/**
 * @brief the contact patterns within the society are modelled using a ContactMatrix
 */
template <typename FP>
struct ContactPatterns {
    using Type = ContactMatrix<FP>;
    static Type get_default()
    {
        return Type{1};
    }
    static std::string name()
    {
        return "ContactPatterns";
    }
};

/**
 * @brief The start day in the SIRS model
 * The start day defines in which season the simulation can be started
 * If the start day is 180 and simulation takes place from t0=0 to
 * tmax=100 the days 180 to 280 of the year are simulated
 */
template <typename FP>
struct StartDay {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "StartDay";
    }
};

/**
 * @brief The seasonality in the SIRS model.
 * The seasonality is given as (1+k*sin()) where the sine
 * curve is below one in summer and above one in winter
 */
template <typename FP>
struct Seasonality {
    using Type = UncertainValue<FP>;
    static Type get_default()
    {
        return Type(0.);
    }
    static std::string name()
    {
        return "Seasonality";
    }
};

template <typename FP>
using ParametersBase =
    ParameterSet<TransmissionProbabilityOnContact<FP>, TimeInfected<FP>, ContactPatterns<FP>, TimeImmune<FP>, Seasonality<FP>, StartDay<FP>>;

/**
 * @brief Parameters of SIR model.
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
        FP tol_times = 1e-1;

        int corrected = false;
        if (this->template get<Seasonality<FP>>() < 0.0 || this->template get<Seasonality<FP>>() > 0.5) {
            log_warning("Constraint check: Parameter Seasonality changed from {:0.4f} to {:d}",
                        this->template get<Seasonality<FP>>(), 0);
            this->template set<Seasonality<FP>>(0);
            corrected = true;
        }
        if (this->template get<TimeInfected<FP>>() < tol_times) {
            log_warning("Constraint check: Parameter TimeInfected changed from {:.4f} to {:.4f}. Please note that "
                        "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                        "and reset parameters.",
                        this->template get<TimeInfected<FP>>(), tol_times);
            this->template get<TimeInfected<FP>>() = tol_times;
            corrected                              = true;
        }
        if (this->template get<TimeImmune<FP>>() < tol_times) {
            log_warning("Constraint check: Parameter TimeInfected changed from {:.4f} to {:.4f}. Please note that "
                        "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                        "and reset parameters.",
                        this->template get<TimeImmune<FP>>(), tol_times);
            this->template get<TimeImmune<FP>>() = tol_times;
            corrected                            = true;
        }
        if (this->template get<TransmissionProbabilityOnContact<FP>>() < 0.0 ||
            this->template get<TransmissionProbabilityOnContact<FP>>() > 1.0) {
            log_warning("Constraint check: Parameter TransmissionProbabilityOnContact changed from {:0.4f} to {:d} ",
                        this->template get<TransmissionProbabilityOnContact<FP>>(), 0.0);
            this->template get<TransmissionProbabilityOnContact<FP>>() = 0.0;
            corrected                                                  = true;
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
        FP tol_times = 1e-1;

        if (this->template get<Seasonality<FP>>() < 0.0 || this->template get<Seasonality<FP>>() > 0.5) {
            log_error("Constraint check: Parameter Seasonality smaller {:d} or larger {:d}", 0, 0.5);
            return true;
        }
        if (this->template get<TimeInfected<FP>>() < tol_times) {
            log_error("Constraint check: Parameter TimeInfected {:.4f} smaller or equal {:.4f}. Please note that "
                      "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                      "and reset parameters.",
                      this->template get<TimeInfected<FP>>(), 0.0);
            return true;
        }
        if (this->template get<TimeImmune<FP>>() < tol_times) {
            log_error("Constraint check: Parameter TimeInfected {:.4f} smaller or equal {:.4f}. Please note that "
                      "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                      "and reset parameters.",
                      this->template get<TimeImmune<FP>>(), 0.0);
            return true;
        }
        if (this->template get<TransmissionProbabilityOnContact<FP>>() < 0.0 ||
            this->template get<TransmissionProbabilityOnContact<FP>>() > 1.0) {
            log_error(
                "Constraint check: Parameter TransmissionProbabilityOnContact {:.4f} smaller {:.4f} or greater {:.4f}",
                this->template get<TransmissionProbabilityOnContact<FP>>(), 0.0, 1.0);
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

} // namespace ssirs
} // namespace mio

#endif // MIO_SDE_SIRS_PARAMETERS_H
