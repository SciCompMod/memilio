/* 
* Copyright (C) 2020-2024 MEmilio
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

/**
 * @brief Vector with the probability to start in any of the subcompartments of the Exposed compartment.
 */
struct StartingProbabilitiesExposed {
    using Type = Eigen::VectorXd;
    /** 
     * @brief Default parameters can be used to get Erlang distributed stay time in Exposed compartment.
     * @param[in] NumExposed Number of subcompartments of the Exposed compartment.
     */
    static Type get_default(size_t NumExposed)
    {
        Eigen::VectorXd def = Eigen::VectorXd::Zero(NumExposed);
        def[0]              = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesExposed";
    }
};

/**
 * @brief Transition Matrix of the Exposed compartment.
 */
struct TransitionMatrixExposedToInfectedNoSymptoms {
    using Type = Eigen::MatrixXd;
    /**
     * @brief Default parameters can be used to get Erlang distributed stay time in Exposed compartment.
     * @param[in] NumExposed Number of subcompartments of the Exposed compartment.
     * @param[in] TimeExposed Average time spent in Exposed compartment in day unit.
     */
    static Type get_default(size_t NumExposed, ScalarType TimeExposed = 2.)
    {
        Eigen::MatrixXd def = Eigen::VectorXd::Constant(NumExposed, -(ScalarType)NumExposed / TimeExposed).asDiagonal();
        def.diagonal(1).setConstant((ScalarType)NumExposed / TimeExposed);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixExposedToInfectedNoSymptoms";
    }
};

/**
 * @brief Vector with the probability to start in any of the subcompartments of the InfectedNoSymptoms compartment.
 */
struct StartingProbabilitiesInfectedNoSymptoms {
    using Type = Eigen::VectorXd;
    /**
     * @brief Default parameters can be used to get Erlang distributed stay time in InfectedNoSymptoms compartment.
     * @param[in] num_InfectedNoSymptoms Number of subcompartments of the InfectedNoSymptoms compartment.
     */
    static Type get_default(size_t num_InfectedNoSymptoms)
    {
        Eigen::VectorXd def = Eigen::VectorXd::Zero(num_InfectedNoSymptoms);
        def[0]              = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesInfectedNoSymptoms";
    }
};

/**
 * @brief Transition Matrix of the Phase-type distribution describing the stay time in the InfectedNoSymptoms 
 *      compartment before developing symptoms.
 */
struct TransitionMatrixInfectedNoSymptomsToInfectedSymptoms {
    using Type = Eigen::MatrixXd;
    /**
     * @brief Default parameters can be used to get Erlang distributed stay time in InfectedNoSymptoms compartment
     *   before developing symptoms.
     * @param[in] dimension Number of rows/columns of the Transitionmatrix.
     * @param[in] time Average time spent in InfectedNoSymptoms before developing symptoms in day unit.
     */
    static Type get_default(size_t dimension, ScalarType time = 1.)
    {
        Eigen::MatrixXd def = Eigen::VectorXd::Constant(dimension, -(ScalarType)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((ScalarType)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedNoSymptomsToInfectedSymptoms";
    }
};

/**
 * @brief Transition Matrix of the Phase-type distribution describing the stay time in the InfectedNoSymptoms 
 *      compartment before recovery.
 */
struct TransitionMatrixInfectedNoSymptomsToRecovered {
    using Type = Eigen::MatrixXd;
    /**
     * @brief Default parameters can be used to get Erlang distributed stay time in InfectedNoSymptoms compartment
     *   before recovery.
     * @param[in] dimension Number of rows/columns of the Transitionmatrix.
     * @param[in] time Average time spent in InfectedNoSymptoms before recovery in day unit.
     */
    static Type get_default(size_t dimension, ScalarType time = 1.)
    {
        Eigen::MatrixXd def = Eigen::VectorXd::Constant(dimension, -(ScalarType)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((ScalarType)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedNoSymptomsToInfectedSymptomsToRecovered";
    }
};

/**
 * @brief Vector with the probability to start in any of the subcompartments of the InfectedSymptoms compartment.
 */
struct StartingProbabilitiesInfectedSymptoms {
    using Type = Eigen::VectorXd;
    /**
     * @brief Default parameters can be used to get Erlang distributed stay time in InfectedSymptoms compartment.
     * @param[in] num_InfectedSymptoms Number of subcompartments of the Infectedymptoms compartment.
     */
    static Type get_default(size_t num_InfectedSymptoms)
    {
        Eigen::VectorXd def = Eigen::VectorXd::Zero(num_InfectedSymptoms);
        def[0]              = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesInfectedSymptoms";
    }
};

/**
 * @brief Transition Matrix of the Phase-type distribution describing the stay time in the InfectedNoSymptoms 
 *      compartment before going to hospital.
 */
struct TransitionMatrixInfectedSymptomsToInfectedSevere {
    using Type = Eigen::MatrixXd;
    /**
     * @brief Default parameters can be used to get Erlang distributed stay time in InfectedSymptoms compartment
     *   before going to hospital.
     * @param[in] dimension Number of rows/columns of the Transitionmatrix.
     * @param[in] time Average time spent in InfectedSymptoms before going to hospital in day unit.
     */
    static Type get_default(size_t dimension, ScalarType time = 1.5)
    {
        Eigen::MatrixXd def = Eigen::VectorXd::Constant(dimension, -(ScalarType)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((ScalarType)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedSymptomsToInfectedSevere";
    }
};

/**
 * @brief Transition Matrix of the Phase-type distribution describing the stay time in the InfectedSymptoms 
 *      compartment before recovery.
 */
struct TransitionMatrixInfectedSymptomsToRecovered {
    using Type = Eigen::MatrixXd;
    /**
     * @brief Default parameters can be used to get Erlang distributed stay time in InfectedSymptoms compartment
     *   before recovery.
     * @param[in] dimension Number of rows/columns of the Transitionmatrix.
     * @param[in] time Average time spent in InfectedSymptoms before recovery in day unit.
     */
    static Type get_default(size_t dimension, ScalarType time = 1.5)
    {
        Eigen::MatrixXd def = Eigen::VectorXd::Constant(dimension, -(ScalarType)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((ScalarType)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedSymptomsToRecovered";
    }
};

/**
 * @brief Vector with the probability to start in any of the subcompartments of the InfectedSevere compartment.
 */
struct StartingProbabilitiesInfectedSevere {
    using Type = Eigen::VectorXd;
    /**
     * @param[in] num_InfectedSevere Number of subcompartments of the InfectedSevere compartment.
     */
    static Type get_default(size_t num_InfectedSevere)
    {
        Eigen::VectorXd def = Eigen::VectorXd::Zero(num_InfectedSevere);
        def[0]              = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesInfectedSevere";
    }
};

/**
 * @brief Transition Matrix of the Phase-type distribution describing the stay time in the InfectedSevere
 *      compartment before treated by ICU.
 */
struct TransitionMatrixInfectedSevereToInfectedCritical {
    using Type = Eigen::MatrixXd;
    /**
     * @brief Default parameters can be used to get Erlang distributed stay time in InfectedSevere compartment
     *   before treated by ICU.
     * @param[in] dimension Number of rows/columns of the Transitionmatrix.
     * @param[in] time Average time spent in InfectedSevere before treated by ICU in day unit.
     */
    static Type get_default(size_t dimension, ScalarType time = 1.)
    {
        Eigen::MatrixXd def = Eigen::VectorXd::Constant(dimension, -(ScalarType)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((ScalarType)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedSevereToInfectedCritical";
    }
};

/**
 * @brief Transition Matrix of the Phase-type distribution describing the stay time in the InfectedSevere
 *      compartment before recovery.
 */
struct TransitionMatrixInfectedSevereToRecovered {
    using Type = Eigen::MatrixXd;
    /**
     * @brief Default parameters can be used to get Erlang distributed stay time in InfectedSevere compartment
     *   before recovery.
     * @param[in] dimension Number of rows/columns of the Transitionmatrix.
     * @param[in] time Average time spent in InfectedSevere before recovery in day unit.
     */
    static Type get_default(size_t dimension, ScalarType time = 1.)
    {
        Eigen::MatrixXd def = Eigen::VectorXd::Constant(dimension, -(ScalarType)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((ScalarType)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedSevereToRecovered";
    }
};

/**
 * @brief Vector with the probability to start in any of the subcompartments of the InfectedCritical compartment.
 */
struct StartingProbabilitiesInfectedCritical {
    using Type = Eigen::VectorXd;
    /**
     * @param[in] num_InfectedCritical Number of subcompartments of the InfectedCritical compartment.
     */
    static Type get_default(size_t num_InfectedCritical)
    {
        Eigen::VectorXd def = Eigen::VectorXd::Zero(num_InfectedCritical);
        def[0]              = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesInfectedCritical";
    }
};

/**
 * @brief Transition Matrix of the Phase-type distribution describing the stay time in the InfectedCritical
 *      compartment before death.
 */
struct TransitionMatrixInfectedCriticalToDead {
    using Type = Eigen::MatrixXd;
    /**
     * @brief Default parameters can be used to get Erlang distributed stay time in InfectedCritical compartment
     *   before death.
     * @param[in] dimension Number of rows/columns of the Transitionmatrix.
     * @param[in] time Average time treated by ICU before dead in day unit.
     */
    static Type get_default(size_t dimension, ScalarType time = 1.)
    {
        Eigen::MatrixXd def = Eigen::VectorXd::Constant(dimension, -(ScalarType)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((ScalarType)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedCriticalToDead";
    }
};

/**
 * @brief Transition Matrix of the Phase-type distribution describing the stay time in the InfectedCritical
 *      compartment before recovery.
 */
struct TransitionMatrixInfectedCriticalToRecovered {
    using Type = Eigen::MatrixXd;
    /**
     * @brief Default parameters can be used to get Erlang distributed stay time in InfectedCritical compartment
     *   before recovery.
     * @param[in] dimension Number of rows/columns of the Transitionmatrix.
     * @param[in] time Average time treated by ICU before recovery in day unit.
     */
    static Type get_default(size_t dimension, ScalarType time = 1.)
    {
        Eigen::MatrixXd def = Eigen::VectorXd::Constant(dimension, -(ScalarType)dimension / time).asDiagonal();
        def.diagonal(1).setConstant((ScalarType)dimension / time);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedCriticalToRecovered";
    }
};

/**
 * @brief Probability of getting infected from a contact.
 */
struct TransmissionProbabilityOnContact {
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
 * @brief The contact patterns within the society are modelled using an UncertainContactMatrix.
 */
struct ContactPatterns {
    using Type = UncertainContactMatrix<ScalarType>;

    static Type get_default()
    {
        ContactMatrixGroup contact_matrix = ContactMatrixGroup(1, 1);
        contact_matrix[0]                 = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
        return Type(contact_matrix);
    }
    static std::string name()
    {
        return "ContactPatterns";
    }
};

/**
 * @brief The relative InfectedNoSymptoms infectability.
 */
struct RelativeTransmissionNoSymptoms {
    using Type = ScalarType;
    static Type get_default()
    {
        return 0.5;
    }
    static std::string name()
    {
        return "RelativeTransmissionNoSymptoms";
    }
};

/**
 * @brief The risk of infection from symptomatic cases in the GLCT-SECIR model.
 */
struct RiskOfInfectionFromSymptomatic {
    using Type = ScalarType;
    static Type get_default()
    {
        return 0.5;
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSymptomatic";
    }
};

/**
 * @brief The start day in the GLCT-SECIR model.
 * The start day defines in which season the simulation is started.
 * If the start day is 180 and simulation takes place from t0=0 to
 * tmax=100 the days 180 to 280 of the year are simulated.
 */
struct StartDay {
    using Type = ScalarType;
    static Type get_default()
    {
        return 0.;
    }
    static std::string name()
    {
        return "StartDay";
    }
};

/**
 * @brief The seasonality in the GLCT-SECIR model.
 * The seasonality is given as (1+k*sin()) where the sine
 * curve is below one in summer and above one in winter.
 */
struct Seasonality {
    using Type = ScalarType;
    static Type get_default()
    {
        return Type(0.);
    }
    static std::string name()
    {
        return "Seasonality";
    }
};

using ParametersBase =
    ParameterSet<StartingProbabilitiesExposed, TransitionMatrixExposedToInfectedNoSymptoms,
                 StartingProbabilitiesInfectedNoSymptoms, TransitionMatrixInfectedNoSymptomsToInfectedSymptoms,
                 TransitionMatrixInfectedNoSymptomsToRecovered, StartingProbabilitiesInfectedSymptoms,
                 TransitionMatrixInfectedSymptomsToInfectedSevere, TransitionMatrixInfectedSymptomsToRecovered,
                 StartingProbabilitiesInfectedSevere, TransitionMatrixInfectedSevereToInfectedCritical,
                 TransitionMatrixInfectedSevereToRecovered, StartingProbabilitiesInfectedCritical,
                 TransitionMatrixInfectedCriticalToDead, TransitionMatrixInfectedCriticalToRecovered,
                 TransmissionProbabilityOnContact, ContactPatterns, RelativeTransmissionNoSymptoms,
                 RiskOfInfectionFromSymptomatic, StartDay, Seasonality>;

/**
 * @brief Parameters of an GLCT-SECIR model.
 */
class Parameters : public ParametersBase
{
public:
    /**
     * @brief Default constructor.
     */
    Parameters()
        : ParametersBase()
    {
    }

    /**
     * @brief Checks whether all Parameters satisfy their corresponding constraints and logs an error
     *      if constraints are not satisfied.
     * @return Returns true if one (or more) constraint(s) are not satisfied, otherwise false. 
     */
    bool check_constraints() const
    {
        // --- Parameters affecting the transmission of the virus. ---
        if (this->get<TransmissionProbabilityOnContact>() < 0.0 ||
            this->get<TransmissionProbabilityOnContact>() > 1.0) {
            log_error("Constraint check: Parameter TransmissionProbabilityOnContact smaller {:d} or larger {:d}", 0, 1);
            return true;
        }

        if (this->get<RelativeTransmissionNoSymptoms>() < 0.0 || this->get<RelativeTransmissionNoSymptoms>() > 1.0) {
            log_error("Constraint check: Parameter RelativeTransmissionNoSymptoms smaller {:d} or larger {:d}", 0, 1);
            return true;
        }

        if (this->get<RiskOfInfectionFromSymptomatic>() < 0.0 || this->get<RiskOfInfectionFromSymptomatic>() > 1.0) {
            log_error("Constraint check: Parameter RiskOfInfectionFromSymptomatic smaller {:d} or larger {:d}", 0, 1);
            return true;
        }

        if (this->get<Seasonality>() < 0.0 || this->get<Seasonality>() > 0.5) {
            log_warning("Constraint check: Parameter Seasonality should lie between {:0.4f} and {:.4f}", 0.0, 0.5);
            return true;
        }

        // --- Parameters affecting the phase-type distributions. ---
        // --- Check that the dimensions are consistent. ---
        if ((this->get<TransitionMatrixExposedToInfectedNoSymptoms>().cols() !=
             this->get<TransitionMatrixExposedToInfectedNoSymptoms>().rows()) ||
            (this->get<TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>().cols() !=
             this->get<TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>().rows()) ||
            (this->get<TransitionMatrixInfectedNoSymptomsToRecovered>().cols() !=
             this->get<TransitionMatrixInfectedNoSymptomsToRecovered>().rows()) ||
            (this->get<TransitionMatrixInfectedSymptomsToInfectedSevere>().cols() !=
             this->get<TransitionMatrixInfectedSymptomsToInfectedSevere>().rows()) ||
            (this->get<TransitionMatrixInfectedSymptomsToRecovered>().cols() !=
             this->get<TransitionMatrixInfectedSymptomsToRecovered>().rows()) ||
            (this->get<TransitionMatrixInfectedSevereToInfectedCritical>().cols() !=
             this->get<TransitionMatrixInfectedSevereToInfectedCritical>().rows()) ||
            (this->get<TransitionMatrixInfectedSevereToRecovered>().cols() !=
             this->get<TransitionMatrixInfectedSevereToRecovered>().rows()) ||
            (this->get<TransitionMatrixInfectedCriticalToDead>().cols() !=
             this->get<TransitionMatrixInfectedCriticalToDead>().rows()) ||
            (this->get<TransitionMatrixInfectedCriticalToRecovered>().cols() !=
             this->get<TransitionMatrixInfectedCriticalToRecovered>().rows())) {
            log_error("Constraint check: At least one of the matrices used for the TransitionMatrix parameters is not "
                      "quadratic.");
            return true;
        }

        if (this->get<StartingProbabilitiesExposed>().rows() !=
            this->get<TransitionMatrixExposedToInfectedNoSymptoms>().rows()) {
            log_error("Constraint check: Dimensions of StartingProbabilitiesExposed and "
                      "TransitionMatrixExposedToInfectedNoSymptoms not matching.");
            return true;
        }

        if (this->get<StartingProbabilitiesInfectedNoSymptoms>().rows() !=
            this->get<TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>().rows() +
                this->get<TransitionMatrixInfectedNoSymptomsToRecovered>().rows()) {
            log_error("Constraint check: Dimensions of StartingProbabilitiesInfectedNoSymptoms and "
                      "TransitionMatrices of InfectedNoSymptoms compartment are not matching.");
            return true;
        }

        if (this->get<StartingProbabilitiesInfectedSymptoms>().rows() !=
            this->get<TransitionMatrixInfectedSymptomsToInfectedSevere>().rows() +
                this->get<TransitionMatrixInfectedSymptomsToRecovered>().rows()) {
            log_error("Constraint check: Dimensions of StartingProbabilitiesInfectedSymptoms and "
                      "TransitionMatrices of InfectedSymptoms compartment are not matching.");
            return true;
        }

        if (this->get<StartingProbabilitiesInfectedSevere>().rows() !=
            this->get<TransitionMatrixInfectedSevereToInfectedCritical>().rows() +
                this->get<TransitionMatrixInfectedSevereToRecovered>().rows()) {
            log_error("Constraint check: Dimensions of StartingProbabilitiesInfectedSevere and "
                      "TransitionMatrices of InfectedSevere compartment are not matching.");
            return true;
        }

        if (this->get<StartingProbabilitiesInfectedCritical>().rows() !=
            this->get<TransitionMatrixInfectedCriticalToDead>().rows() +
                this->get<TransitionMatrixInfectedCriticalToRecovered>().rows()) {
            log_error("Constraint check: Dimensions of StartingProbabilitiesInfectedCritical and "
                      "TransitionMatrices of InfectedCritical compartment are not matching.");
            return true;
        }

        // --- Check other restrictions. ---
        if (!floating_point_equal(1., this->get<StartingProbabilitiesExposed>().sum())) {
            log_warning("Constraint check: The vector StartingProbabilitiesExposed should add up to one.");
            return true;
        }
        if (!floating_point_equal(1., this->get<StartingProbabilitiesInfectedNoSymptoms>().sum())) {
            log_warning("Constraint check: The vector StartingProbabilitiesInfectedNoSymptoms should add up to one.");
            return true;
        }
        if (!floating_point_equal(1., this->get<StartingProbabilitiesInfectedSymptoms>().sum())) {
            log_warning("Constraint check: The vector StartingProbabilitiesInfectedSymptoms should add up to one.");
            return true;
        }
        if (!floating_point_equal(1., this->get<StartingProbabilitiesInfectedSevere>().sum())) {
            log_warning("Constraint check: The vector StartingProbabilitiesInfectedSevere should add up to one.");
            return true;
        }
        if (!floating_point_equal(1., this->get<StartingProbabilitiesInfectedCritical>().sum())) {
            log_warning("Constraint check: The vector StartingProbabilitiesInfectedCritical should add up to one.");
            return true;
        }

        if ((this->get<StartingProbabilitiesExposed>().array() < -1e-10).any()) {
            log_warning(
                "Constraint check: The vector StartingProbabilitiesExposed should just have non-negative entries.");
            return true;
        }
        if ((this->get<StartingProbabilitiesInfectedNoSymptoms>().array() < -1e-10).any()) {
            log_warning("Constraint check: The vector StartingProbabilitiesInfectedNoSymptoms should just have "
                        "non-negative entries.");
            return true;
        }
        if ((this->get<StartingProbabilitiesInfectedSymptoms>().array() < -1e-10).any()) {
            log_warning("Constraint check: The vector StartingProbabilitiesInfectedSymptoms should just have "
                        "non-negative entries.");
            return true;
        }
        if ((this->get<StartingProbabilitiesInfectedSevere>().array() < -1e-10).any()) {
            log_warning("Constraint check: The vector StartingProbabilitiesInfectedSevere should just have "
                        "non-negative entries.");
            return true;
        }
        if ((this->get<StartingProbabilitiesInfectedCritical>().array() < -1e-10).any()) {
            log_warning("Constraint check: The vector StartingProbabilitiesInfectedCritical should just have "
                        "non-negative entries.");
            return true;
        }

        // --- Check that we have no flows back from one compartment to the last (only in between of the subcompartments). ---
        if (((this->get<TransitionMatrixExposedToInfectedNoSymptoms>() *
              Eigen::VectorXd::Ones(this->get<TransitionMatrixExposedToInfectedNoSymptoms>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixExposedToInfectedNoSymptoms lead to a negative "
                "flow ExposedToInfectedNoSymptoms.");
            return true;
        }
        if (((this->get<TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>() *
              Eigen::VectorXd::Ones(this->get<TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning("Constraint check: The entries of TransitionMatrixInfectedNoSymptomsToInfectedSymptoms lead to "
                        "a negative "
                        "flow InfectedNoSymptomsToInfectedSymptoms.");
            return true;
        }
        if (((this->get<TransitionMatrixInfectedNoSymptomsToRecovered>() *
              Eigen::VectorXd::Ones(this->get<TransitionMatrixInfectedNoSymptomsToRecovered>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixInfectedNoSymptomsToRecovered lead to a negative "
                "flow InfectedNoSymptomsToRecovered.");
            return true;
        }
        if (((this->get<TransitionMatrixInfectedSymptomsToInfectedSevere>() *
              Eigen::VectorXd::Ones(this->get<TransitionMatrixInfectedSymptomsToInfectedSevere>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixInfectedSymptomsToInfectedSevere lead to a negative "
                "flow InfectedSymptomsToInfectedSevere.");
            return true;
        }
        if (((this->get<TransitionMatrixInfectedSymptomsToRecovered>() *
              Eigen::VectorXd::Ones(this->get<TransitionMatrixInfectedSymptomsToRecovered>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixInfectedSymptomsToRecovered lead to a negative "
                "flow InfectedSymptomsToRecovered.");
            return true;
        }
        if (((this->get<TransitionMatrixInfectedSevereToInfectedCritical>() *
              Eigen::VectorXd::Ones(this->get<TransitionMatrixInfectedSevereToInfectedCritical>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning(
                "Constraint check: The entries of TransitionMatrixInfectedSevereToInfectedCritical lead to a negative "
                "flow InfectedSevereToInfectedCritical.");
            return true;
        }
        if (((this->get<TransitionMatrixInfectedSevereToRecovered>() *
              Eigen::VectorXd::Ones(this->get<TransitionMatrixInfectedSevereToRecovered>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning("Constraint check: The entries of TransitionMatrixInfectedSevereToRecovered lead to a negative "
                        "flow InfectedSevereToRecovered.");
            return true;
        }
        if (((this->get<TransitionMatrixInfectedCriticalToDead>() *
              Eigen::VectorXd::Ones(this->get<TransitionMatrixInfectedCriticalToDead>().rows()))
                 .array() > 1e-10)
                .any()) {
            log_warning("Constraint check: The entries of TransitionMatrixInfectedCriticalToDead lead to a negative "
                        "flow InfectedCriticalToDead.");
            return true;
        }
        if (((this->get<TransitionMatrixInfectedCriticalToRecovered>() *
              Eigen::VectorXd::Ones(this->get<TransitionMatrixInfectedCriticalToRecovered>().rows()))
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
        BOOST_OUTCOME_TRY(auto&& base, ParametersBase::deserialize(io));
        return success(Parameters(std::move(base)));
    }
};

} // namespace glsecir
} // namespace mio

#endif // MIO_GLCT_SECIR_PARAMS_H