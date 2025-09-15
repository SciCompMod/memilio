/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Lena Ploetzke
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

#ifndef MIO_GLCT_SECIR_PARAMS_H
#define MIO_GLCT_SECIR_PARAMS_H

#include "memilio/config.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/math/eigen.h"
#include "memilio/math/floating_point.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/logging.h"

namespace mio
{
namespace glsecir
{

/***********************************************
* Define Parameters of the GLCT-SECIHURD model *
***********************************************/

/// @brief Vector with the probability to start in any of the subcompartments of the Exposed compartment.
template <typename FP>
struct StartingProbabilitiesExposed {
    using Type = Eigen::VectorX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in the Exposed compartment.
     * @param[in] numExposed Number of subcompartments of the Exposed compartment.
     */
    static Type get_default(size_t numExposed)
    {
        Eigen::VectorX<FP> def = Eigen::VectorX<FP>::Zero(numExposed);
        def[0]                 = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesExposed";
    }
};

/// @brief Transition matrix of the Exposed compartment.
template <typename FP>
struct TransitionMatrixExposedToInfectedNoSymptoms {
    using Type = Eigen::MatrixX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in the Exposed compartment.
     * @param[in] numExposed Number of subcompartments of the Exposed compartment.
     * @param[in] timeExposed Average time spent in Exposed compartment in day unit.
     */
    static Type get_default(size_t numExposed, FP timeExposed = 1.)
    {
        Eigen::MatrixX<FP> def = Eigen::VectorX<FP>::Constant(numExposed, -(FP)numExposed / timeExposed).asDiagonal();
        def.diagonal(1).setConstant((FP)numExposed / timeExposed);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixExposedToInfectedNoSymptoms";
    }
};

/// @brief Vector with the probability to start in any of the subcompartments of the InfectedNoSymptoms compartment.
template <typename FP>
struct StartingProbabilitiesInfectedNoSymptoms {
    using Type = Eigen::VectorX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in InfectedNoSymptoms compartment.
     * @param[in] numInfectedNoSymptoms Number of subcompartments of the InfectedNoSymptoms compartment.
     */
    static Type get_default(size_t numInfectedNoSymptoms)
    {
        Eigen::VectorX<FP> def = Eigen::VectorX<FP>::Zero(numInfectedNoSymptoms);
        def[0]                 = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesInfectedNoSymptoms";
    }
};

/**
 * @brief Transition matrix of the phase-type distribution describing the stay time in the InfectedNoSymptoms
 *      compartment before developing symptoms.
 */
template <typename FP>
struct TransitionMatrixInfectedNoSymptomsToInfectedSymptoms {
    using Type = Eigen::MatrixX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in InfectedNoSymptoms compartment
     *   before developing symptoms.
     * @param[in] dimension Number of rows/columns of the transition matrix.
     * @param[in] time Average time spent in InfectedNoSymptoms before developing symptoms in day unit.
     */
    static Type get_default(size_t dimension, FP time = 1.)
    {
        Eigen::MatrixX<FP> def = Eigen::VectorX<FP>::Constant(dimension, -(FP)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((FP)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedNoSymptomsToInfectedSymptoms";
    }
};

/**
 * @brief Transition matrix of the phase-type distribution describing the stay time in the InfectedNoSymptoms
 *      compartment before recovery.
 */
template <typename FP>
struct TransitionMatrixInfectedNoSymptomsToRecovered {
    using Type = Eigen::MatrixX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in InfectedNoSymptoms compartment
     *   before recovery.
     * @param[in] dimension Number of rows/columns of the transition matrix.
     * @param[in] time Average time spent in InfectedNoSymptoms before recovery in day unit.
     */
    static Type get_default(size_t dimension, FP time = 1.)
    {
        Eigen::MatrixX<FP> def = Eigen::VectorX<FP>::Constant(dimension, -(FP)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((FP)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedNoSymptomsToInfectedSymptomsToRecovered";
    }
};

/// @brief Vector with the probability to start in any of the subcompartments of the InfectedSymptoms compartment.
template <typename FP>
struct StartingProbabilitiesInfectedSymptoms {
    using Type = Eigen::VectorX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in InfectedSymptoms compartment.
     * @param[in] numInfectedSymptoms Number of subcompartments of the InfectedSymptoms compartment.
     */
    static Type get_default(size_t numInfectedSymptoms)
    {
        Eigen::VectorX<FP> def = Eigen::VectorX<FP>::Zero(numInfectedSymptoms);
        def[0]                 = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesInfectedSymptoms";
    }
};

/**
 * @brief Transition matrix of the phase-type distribution describing the stay time in the InfectedNoSymptoms
 *      compartment before going to hospital.
 */
template <typename FP>
struct TransitionMatrixInfectedSymptomsToInfectedSevere {
    using Type = Eigen::MatrixX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in the InfectedSymptoms compartment
     *   before going to hospital.
     * @param[in] dimension Number of rows/columns of the transition matrix.
     * @param[in] time Average time spent in InfectedSymptoms before going to hospital in day unit.
     */
    static Type get_default(size_t dimension, FP time = 1.)
    {
        Eigen::MatrixX<FP> def = Eigen::VectorX<FP>::Constant(dimension, -(FP)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((FP)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedSymptomsToInfectedSevere";
    }
};

/**
 * @brief Transition matrix of the phase-type distribution describing the stay time in the InfectedSymptoms
 *      compartment before recovery.
 */
template <typename FP>
struct TransitionMatrixInfectedSymptomsToRecovered {
    using Type = Eigen::MatrixX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in the InfectedSymptoms compartment
     *   before recovery.
     * @param[in] dimension Number of rows/columns of the transition matrix.
     * @param[in] time Average time spent in InfectedSymptoms before recovery in day unit.
     */
    static Type get_default(size_t dimension, FP time = 1.)
    {
        Eigen::MatrixX<FP> def = Eigen::VectorX<FP>::Constant(dimension, -(FP)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((FP)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedSymptomsToRecovered";
    }
};

/// @brief Vector with the probability to start in any of the subcompartments of the InfectedSevere compartment.
template <typename FP>
struct StartingProbabilitiesInfectedSevere {
    using Type = Eigen::VectorX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in InfectedSevere compartment.
     * @param[in] numInfectedSevere Number of subcompartments of the InfectedSevere compartment.
     */
    static Type get_default(size_t numInfectedSevere)
    {
        Eigen::VectorX<FP> def = Eigen::VectorX<FP>::Zero(numInfectedSevere);
        def[0]                 = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesInfectedSevere";
    }
};

/**
 * @brief Transition matrix of the phase-type distribution describing the stay time in the InfectedSevere
 *      compartment before treated by ICU.
 */
template <typename FP>
struct TransitionMatrixInfectedSevereToInfectedCritical {
    using Type = Eigen::MatrixX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in InfectedSevere compartment
     *   before treated by ICU.
     * @param[in] dimension Number of rows/columns of the transition matrix.
     * @param[in] time Average time spent in InfectedSevere before treated by ICU in day unit.
     */
    static Type get_default(size_t dimension, FP time = 1.)
    {
        Eigen::MatrixX<FP> def = Eigen::VectorX<FP>::Constant(dimension, -(FP)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((FP)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedSevereToInfectedCritical";
    }
};

/**
 * @brief Transition matrix of the phase-type distribution describing the stay time in the InfectedSevere
 *      compartment before recovery.
 */
template <typename FP>
struct TransitionMatrixInfectedSevereToRecovered {
    using Type = Eigen::MatrixX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in InfectedSevere compartment
     *   before recovery.
     * @param[in] dimension Number of rows/columns of the transition matrix.
     * @param[in] time Average time spent in InfectedSevere before recovery in day unit.
     */
    static Type get_default(size_t dimension, FP time = 1.)
    {
        Eigen::MatrixX<FP> def = Eigen::VectorX<FP>::Constant(dimension, -(FP)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((FP)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedSevereToRecovered";
    }
};

/// @brief Vector with the probability to start in any of the subcompartments of the InfectedCritical compartment.
template <typename FP>
struct StartingProbabilitiesInfectedCritical {
    using Type = Eigen::VectorX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in InfectedCritical compartment.
     * @param[in] numInfectedCritical Number of subcompartments of the InfectedCritical compartment.
     */
    static Type get_default(size_t numInfectedCritical)
    {
        Eigen::VectorX<FP> def = Eigen::VectorX<FP>::Zero(numInfectedCritical);
        def[0]                 = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesInfectedCritical";
    }
};

/**
 * @brief Transition matrix of the phase-type distribution describing the stay time in the InfectedCritical
 *      compartment before death.
 */
template <typename FP>
struct TransitionMatrixInfectedCriticalToDead {
    using Type = Eigen::MatrixX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in InfectedCritical compartment
     *   before death.
     * @param[in] dimension Number of rows/columns of the transition matrix.
     * @param[in] time Average time treated by ICU before dying in day unit.
     */
    static Type get_default(size_t dimension, FP time = 1.)
    {
        Eigen::MatrixX<FP> def = Eigen::VectorX<FP>::Constant(dimension, -(FP)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((FP)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedCriticalToDead";
    }
};

/**
 * @brief Transition matrix of the phase-type distribution describing the stay time in the InfectedCritical
 *      compartment before recovery.
 */
template <typename FP>
struct TransitionMatrixInfectedCriticalToRecovered {
    using Type = Eigen::MatrixX<FP>;
    /**
     * @brief Default parameters can be used to get an Erlang distributed stay time in InfectedCritical compartment
     *   before recovery.
     * @param[in] dimension Number of rows/columns of the transition matrix.
     * @param[in] time Average time treated by ICU before recovery in day unit.
     */
    static Type get_default(size_t dimension, FP time = 1.)
    {
        Eigen::MatrixX<FP> def = Eigen::VectorX<FP>::Constant(dimension, -(FP)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((FP)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedCriticalToRecovered";
    }
};

/// @brief Probability of getting infected from a contact.
template <typename FP>
struct TransmissionProbabilityOnContact {
    using Type = FP;
    static Type get_default()
    {
        return Type(1.0);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContact";
    }
};

/// @brief The contact patterns within the society are modelled using an UncertainContactMatrix.
template <typename FP>
struct ContactPatterns {
    using Type = UncertainContactMatrix<FP>;

    static Type get_default()
    {
        ContactMatrixGroup<FP> contact_matrix = ContactMatrixGroup<FP>(1, 1);
        contact_matrix[0]                     = mio::ContactMatrix<FP>(Eigen::MatrixX<FP>::Constant(1, 1, 10.));
        return Type(contact_matrix);
    }
    static std::string name()
    {
        return "ContactPatterns";
    }
};

/// @brief The relative InfectedNoSymptoms infectability.
template <typename FP>
struct RelativeTransmissionNoSymptoms {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.5);
    }
    static std::string name()
    {
        return "RelativeTransmissionNoSymptoms";
    }
};

/// @brief The risk of infection from symptomatic cases in the GLCT-SECIR model.
template <typename FP>
struct RiskOfInfectionFromSymptomatic {
    using Type = FP;
    static Type get_default()
    {
        return Type(0.5);
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSymptomatic";
    }
};

/**
 * @brief The start day in the GLCT-SECIR model.
 *  The start day defines in which season the simulation is started.
 *  If the start day is 180 and simulation takes place from t0=0 to
 *  tmax=100 the days 180 to 280 of the year are simulated.
 */
template <typename FP>
struct StartDay {
    using Type = FP;
    static Type get_default(size_t)
    {
        return Type(0.0);
    }
    static std::string name()
    {
        return "StartDay";
    }
};

/**
 * @brief The seasonality in the GLCT-SECIR model.
 *  The seasonality is given as (1+k*sin()) where the sine
 *  curve is below one in summer and above one in winter.
 */
template <typename FP>
struct Seasonality {
    using Type = FP;
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
    ParameterSet<StartingProbabilitiesExposed<FP>, TransitionMatrixExposedToInfectedNoSymptoms<FP>,
                 StartingProbabilitiesInfectedNoSymptoms<FP>, TransitionMatrixInfectedNoSymptomsToInfectedSymptoms<FP>,
                 TransitionMatrixInfectedNoSymptomsToRecovered<FP>, StartingProbabilitiesInfectedSymptoms<FP>,
                 TransitionMatrixInfectedSymptomsToInfectedSevere<FP>, TransitionMatrixInfectedSymptomsToRecovered<FP>,
                 StartingProbabilitiesInfectedSevere<FP>, TransitionMatrixInfectedSevereToInfectedCritical<FP>,
                 TransitionMatrixInfectedSevereToRecovered<FP>, StartingProbabilitiesInfectedCritical<FP>,
                 TransitionMatrixInfectedCriticalToDead<FP>, TransitionMatrixInfectedCriticalToRecovered<FP>,
                 TransmissionProbabilityOnContact<FP>, ContactPatterns<FP>, RelativeTransmissionNoSymptoms<FP>,
                 RiskOfInfectionFromSymptomatic<FP>, StartDay<FP>, Seasonality<FP>>;

/// @brief Parameters of an GLCT-SECIR model.
template <typename FP>
class Parameters : public ParametersBase<FP>
{
public:
    /// @brief Default constructor.
    Parameters()
        : ParametersBase<FP>()
    {
    }

    /**
     * @brief Checks that all parameters satisfy their corresponding constraints and logs an error
     *      if constraints are not satisfied.
     *
     * @return Returns true if one or more constraints are not satisfied, false otherwise.
     */
    bool check_constraints() const
    {
        // --- Parameters affecting the transmission of the virus. ---
        if (this->template get<TransmissionProbabilityOnContact<FP>>() < 0.0 ||
            this->template get<TransmissionProbabilityOnContact<FP>>() > 1.0) {
            log_error("Constraint check: Parameter TransmissionProbabilityOnContact smaller {:d} or larger {:d}", 0, 1);
            return true;
        }

        if (this->template get<RelativeTransmissionNoSymptoms<FP>>() < 0.0 ||
            this->template get<RelativeTransmissionNoSymptoms<FP>>() > 1.0) {
            log_error("Constraint check: Parameter RelativeTransmissionNoSymptoms smaller {:d} or larger {:d}", 0, 1);
            return true;
        }

        if (this->template get<RiskOfInfectionFromSymptomatic<FP>>() < 0.0 ||
            this->template get<RiskOfInfectionFromSymptomatic<FP>>() > 1.0) {
            log_error("Constraint check: Parameter RiskOfInfectionFromSymptomatic smaller {:d} or larger {:d}", 0, 1);
            return true;
        }

        if (this->template get<Seasonality<FP>>() < 0.0 || this->template get<Seasonality<FP>>() > 0.5) {
            log_warning("Constraint check: Parameter Seasonality should lie between {:0.4f} and {:.4f}", 0.0, 0.5);
            return true;
        }

        // --- Parameters affecting the phase-type distributions. ---
        // --- Check that the dimensions are consistent. ---
        if ((this->template get<TransitionMatrixExposedToInfectedNoSymptoms<FP>>().cols() !=
             this->template get<TransitionMatrixExposedToInfectedNoSymptoms<FP>>().rows()) ||
            (this->template get<TransitionMatrixInfectedNoSymptomsToInfectedSymptoms<FP>>().cols() !=
             this->template get<TransitionMatrixInfectedNoSymptomsToInfectedSymptoms<FP>>().rows()) ||
            (this->template get<TransitionMatrixInfectedNoSymptomsToRecovered<FP>>().cols() !=
             this->template get<TransitionMatrixInfectedNoSymptomsToRecovered<FP>>().rows()) ||
            (this->template get<TransitionMatrixInfectedSymptomsToInfectedSevere<FP>>().cols() !=
             this->template get<TransitionMatrixInfectedSymptomsToInfectedSevere<FP>>().rows()) ||
            (this->template get<TransitionMatrixInfectedSymptomsToRecovered<FP>>().cols() !=
             this->template get<TransitionMatrixInfectedSymptomsToRecovered<FP>>().rows()) ||
            (this->template get<TransitionMatrixInfectedSevereToInfectedCritical<FP>>().cols() !=
             this->template get<TransitionMatrixInfectedSevereToInfectedCritical<FP>>().rows()) ||
            (this->template get<TransitionMatrixInfectedSevereToRecovered<FP>>().cols() !=
             this->template get<TransitionMatrixInfectedSevereToRecovered<FP>>().rows()) ||
            (this->template get<TransitionMatrixInfectedCriticalToDead<FP>>().cols() !=
             this->template get<TransitionMatrixInfectedCriticalToDead<FP>>().rows()) ||
            (this->template get<TransitionMatrixInfectedCriticalToRecovered<FP>>().cols() !=
             this->template get<TransitionMatrixInfectedCriticalToRecovered<FP>>().rows())) {
            log_error("Constraint check: At least one of the matrices used for the TransitionMatrix parameters is not "
                      "quadratic.");
            return true;
        }

        if (this->template get<StartingProbabilitiesExposed<FP>>().rows() !=
            this->template get<TransitionMatrixExposedToInfectedNoSymptoms<FP>>().rows()) {
            log_error("Constraint check: Dimensions of StartingProbabilitiesExposed and "
                      "TransitionMatrixExposedToInfectedNoSymptoms are not matching.");
            return true;
        }

        if (this->template get<StartingProbabilitiesInfectedNoSymptoms<FP>>().rows() !=
            this->template get<TransitionMatrixInfectedNoSymptomsToInfectedSymptoms<FP>>().rows() +
                this->template get<TransitionMatrixInfectedNoSymptomsToRecovered<FP>>().rows()) {
            log_error("Constraint check: Dimensions of StartingProbabilitiesInfectedNoSymptoms and "
                      "TransitionMatrices of InfectedNoSymptoms compartment are not matching.");
            return true;
        }

        if (this->template get<StartingProbabilitiesInfectedSymptoms<FP>>().rows() !=
            this->template get<TransitionMatrixInfectedSymptomsToInfectedSevere<FP>>().rows() +
                this->template get<TransitionMatrixInfectedSymptomsToRecovered<FP>>().rows()) {
            log_error("Constraint check: Dimensions of StartingProbabilitiesInfectedSymptoms and "
                      "TransitionMatrices of InfectedSymptoms compartment are not matching.");
            return true;
        }

        if (this->template get<StartingProbabilitiesInfectedSevere<FP>>().rows() !=
            this->template get<TransitionMatrixInfectedSevereToInfectedCritical<FP>>().rows() +
                this->template get<TransitionMatrixInfectedSevereToRecovered<FP>>().rows()) {
            log_error("Constraint check: Dimensions of StartingProbabilitiesInfectedSevere and "
                      "TransitionMatrices of InfectedSevere compartment are not matching.");
            return true;
        }

        if (this->template get<StartingProbabilitiesInfectedCritical<FP>>().rows() !=
            this->template get<TransitionMatrixInfectedCriticalToDead<FP>>().rows() +
                this->template get<TransitionMatrixInfectedCriticalToRecovered<FP>>().rows()) {
            log_error("Constraint check: Dimensions of StartingProbabilitiesInfectedCritical and "
                      "TransitionMatrices of InfectedCritical compartment are not matching.");
            return true;
        }

        // --- Check constraints of the starting probability vectors. ---
        if ((!floating_point_equal<FP>(1., this->template get<StartingProbabilitiesExposed<FP>>().sum())) ||
            (!floating_point_equal<FP>(1., this->template get<StartingProbabilitiesInfectedNoSymptoms<FP>>().sum())) ||
            (!floating_point_equal<FP>(1., this->template get<StartingProbabilitiesInfectedSymptoms<FP>>().sum())) ||
            (!floating_point_equal<FP>(1., this->template get<StartingProbabilitiesInfectedSevere<FP>>().sum())) ||
            (!floating_point_equal<FP>(1., this->template get<StartingProbabilitiesInfectedCritical<FP>>().sum()))) {
            log_warning(
                "Constraint check: At least one of the vectors for the starting probabilities does not sum to one.");
            return true;
        }

        if ((this->template get<StartingProbabilitiesExposed<FP>>().array() < -1e-10).any() ||
            (this->template get<StartingProbabilitiesInfectedNoSymptoms<FP>>().array() < -1e-10).any() ||
            (this->template get<StartingProbabilitiesInfectedSymptoms<FP>>().array() < -1e-10).any() ||
            (this->template get<StartingProbabilitiesInfectedSevere<FP>>().array() < -1e-10).any() ||
            (this->template get<StartingProbabilitiesInfectedCritical<FP>>().array() < -1e-10).any()) {
            log_warning("Constraint check: At least one of the vectors for the starting probabilities has at least one "
                        "negative entry.");
            return true;
        }

        // --- Check that we have no flows back from one compartment to the previous one
        // (only in between of the subcompartments). ---
        if (((this->template get<TransitionMatrixExposedToInfectedNoSymptoms<FP>>() *
              Eigen::VectorX<FP>::Ones(this->template get<TransitionMatrixExposedToInfectedNoSymptoms<FP>>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixExposedToInfectedNoSymptoms lead to a negative "
                "flow ExposedToInfectedNoSymptoms.");
            return true;
        }
        if (((this->template get<TransitionMatrixInfectedNoSymptomsToInfectedSymptoms<FP>>() *
              Eigen::VectorX<FP>::Ones(
                  this->template get<TransitionMatrixInfectedNoSymptomsToInfectedSymptoms<FP>>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning("Constraint check: The entries of TransitionMatrixInfectedNoSymptomsToInfectedSymptoms lead to "
                        "a negative "
                        "flow InfectedNoSymptomsToInfectedSymptoms.");
            return true;
        }
        if (((this->template get<TransitionMatrixInfectedNoSymptomsToRecovered<FP>>() *
              Eigen::VectorX<FP>::Ones(this->template get<TransitionMatrixInfectedNoSymptomsToRecovered<FP>>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixInfectedNoSymptomsToRecovered lead to a negative "
                "flow InfectedNoSymptomsToRecovered.");
            return true;
        }
        if (((this->template get<TransitionMatrixInfectedSymptomsToInfectedSevere<FP>>() *
              Eigen::VectorX<FP>::Ones(
                  this->template get<TransitionMatrixInfectedSymptomsToInfectedSevere<FP>>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixInfectedSymptomsToInfectedSevere lead to a negative "
                "flow InfectedSymptomsToInfectedSevere.");
            return true;
        }
        if (((this->template get<TransitionMatrixInfectedSymptomsToRecovered<FP>>() *
              Eigen::VectorX<FP>::Ones(this->template get<TransitionMatrixInfectedSymptomsToRecovered<FP>>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixInfectedSymptomsToRecovered lead to a negative "
                "flow InfectedSymptomsToRecovered.");
            return true;
        }
        if (((this->template get<TransitionMatrixInfectedSevereToInfectedCritical<FP>>() *
              Eigen::VectorX<FP>::Ones(
                  this->template get<TransitionMatrixInfectedSevereToInfectedCritical<FP>>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixInfectedSevereToInfectedCritical lead to a negative "
                "flow InfectedSevereToInfectedCritical.");
            return true;
        }
        if (((this->template get<TransitionMatrixInfectedSevereToRecovered<FP>>() *
              Eigen::VectorX<FP>::Ones(this->template get<TransitionMatrixInfectedSevereToRecovered<FP>>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning("Constraint check: The entries of TransitionMatrixInfectedSevereToRecovered lead to a negative "
                        "flow InfectedSevereToRecovered.");
            return true;
        }
        if (((this->template get<TransitionMatrixInfectedCriticalToDead<FP>>() *
              Eigen::VectorX<FP>::Ones(this->template get<TransitionMatrixInfectedCriticalToDead<FP>>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning("Constraint check: The entries of TransitionMatrixInfectedCriticalToDead lead to a negative "
                        "flow InfectedCriticalToDead.");
            return true;
        }
        if (((this->template get<TransitionMatrixInfectedCriticalToRecovered<FP>>() *
              Eigen::VectorX<FP>::Ones(this->template get<TransitionMatrixInfectedCriticalToRecovered<FP>>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixInfectedCriticalToRecovered lead to a negative "
                "flow InfectedCriticalToRecovered.");
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

} // namespace glsecir
} // namespace mio

#endif // MIO_GLCT_SECIR_PARAMS_H
