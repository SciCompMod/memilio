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

#ifndef GLCT_SECIR_PARAMS_H
#define GLCT_SECIR_PARAMS_H

#include "memilio/config.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/math/eigen.h"
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
     * @param[in] NumExposed Number of subcompartiments of the Exposed compartment.
     */
    static Type get_default(int NumExposed)
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
     * @param[in] NumExposed Number of subcompartiments of the Exposed compartment.
     * @param[in] TimeExposed Average time spent in Exposed in day unit.
     */
    static Type get_default(int NumExposed, ScalarType TimeExposed = 2.)
    {
        Eigen::MatrixXd def = Eigen::VectorXd::Constant(NumExposed, -NumExposed / TimeExposed).asDiagonal();
        def.diagonal(1).setConstant(NumExposed / TimeExposed);
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
     * @param[in] NumInfectedNoSymptoms Number of subcompartiments of the InfectedNoSymptoms compartment.
     */
    static Type get_default(int NumInfectedNoSymptoms)
    {
        Eigen::VectorXd def = Eigen::VectorXd::Zero(NumInfectedNoSymptoms);
        def[0]              = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesInfectedNoSymptoms";
    }
};

/**
 * @brief Transition Matrix of the InfectedNoSymptoms compartment.
 */
struct TransitionMatrixInfectedNoSymptomsToInfectedSymptoms {
    using Type = Eigen::MatrixXd;
    /**
     * @param[in] NumInfectedNoSymptoms Number of subcompartiments of the InfectedNoSymptoms compartment.
     * @param[in] TimeInfectedNoSymptoms Average time spent in InfectedNoSymptoms before developing Symptoms or recover
     *       in day unit.
     */
    static Type get_default(int NumInfectedNoSymptoms, ScalarType TimeInfectedNoSymptoms = 1.)
    {
        Eigen::MatrixXd def =
            Eigen::VectorXd::Constant(NumInfectedNoSymptoms, -NumInfectedNoSymptoms / TimeInfectedNoSymptoms)
                .asDiagonal();
        def.diagonal(1).setConstant(NumInfectedNoSymptoms / TimeInfectedNoSymptoms);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedNoSymptomsToInfectedSymptoms";
    }
};

/**
 * @brief Transition Matrix of the InfectedNoSymptoms compartment.
 */
struct TransitionMatrixInfectedNoSymptomsToRecovered {
    using Type = Eigen::MatrixXd;
    /**
     * @param[in] NumInfectedNoSymptoms Number of subcompartiments of the InfectedNoSymptoms compartment.
     * @param[in] TimeInfectedNoSymptoms Average time spent in InfectedNoSymptoms before developing Symptoms or recover
     *       in day unit.
     */
    static Type get_default(int NumInfectedNoSymptoms, ScalarType TimeInfectedNoSymptoms = 1.)
    {
        Eigen::MatrixXd def =
            Eigen::VectorXd::Constant(NumInfectedNoSymptoms, -NumInfectedNoSymptoms / TimeInfectedNoSymptoms)
                .asDiagonal();
        def.diagonal(1).setConstant(NumInfectedNoSymptoms / TimeInfectedNoSymptoms);
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
     * @param[in] NumInfectedSymptoms Number of subcompartiments of the InfectedSymptoms compartment.
     */
    static Type get_default(int NumInfectedSymptoms)
    {
        Eigen::VectorXd def = Eigen::VectorXd::Zero(NumInfectedSymptoms);
        def[0]              = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesInfectedSymptoms";
    }
};

/**
 * @brief Transition Matrix of the InfectedSymptoms compartment.
 */
struct TransitionMatrixInfectedSymptomsToInfectedSevere {
    using Type = Eigen::MatrixXd;
    /**
     * @param[in] NumInfectedSymptoms Number of subcompartiments of the InfectedSymptoms compartment.
     * @param[in] TimeInfectedSymptomsToInfectedSevere Average time spent in InfectedSymptoms before going to Hospital or recover 
     *      in day unit.
     */
    static Type get_default(int NumInfectedSymptoms, ScalarType TimeInfectedSymptomsToInfectedSevere = 1.5)
    {
        Eigen::MatrixXd def =
            Eigen::VectorXd::Constant(NumInfectedSymptoms, -NumInfectedSymptoms / TimeInfectedSymptomsToInfectedSevere)
                .asDiagonal();
        def.diagonal(1).setConstant(NumInfectedSymptoms / TimeInfectedSymptomsToInfectedSevere);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedSymptomsToInfectedSevere";
    }
};

struct TransitionMatrixInfectedSymptomsToRecovered {
    using Type = Eigen::MatrixXd;
    /**
     * @param[in] NumInfectedSymptoms Number of subcompartiments of the InfectedSymptoms compartment.
     * @param[in] TimeInfectedSymptomsToRecovered Average time spent in InfectedSymptoms before going to Hospital or recover 
     *      in day unit.
     */
    static Type get_default(int NumInfectedSymptoms, ScalarType TimeInfectedSymptomsToRecovered = 1.5)
    {
        Eigen::MatrixXd def =
            Eigen::VectorXd::Constant(NumInfectedSymptoms, -NumInfectedSymptoms / TimeInfectedSymptomsToRecovered)
                .asDiagonal();
        def.diagonal(1).setConstant(NumInfectedSymptoms / TimeInfectedSymptomsToRecovered);
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
     * @param[in] NumInfectedSevere Number of subcompartiments of the InfectedSevere compartment.
     */
    static Type get_default(int NumInfectedSevere)
    {
        Eigen::VectorXd def = Eigen::VectorXd::Zero(NumInfectedSevere);
        def[0]              = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesInfectedSevere";
    }
};

/**
 * @brief Transition Matrix of the InfectedSevere compartment.
 */
struct TransitionMatrixInfectedSevereToInfectedCritical {
    using Type = Eigen::MatrixXd;
    /**
     * @param[in] NumInfectedSevere Number of subcompartiments of the InfectedSevere compartment.
     * @param[in] TimeInfectedSevereToInfectedCritical Average time being in the Hospital before treated by ICU or recover in day unit.
     */
    static Type get_default(int NumInfectedSevere, ScalarType TimeInfectedSevereToInfectedCritical = 1.)
    {
        Eigen::MatrixXd def =
            Eigen::VectorXd::Constant(NumInfectedSevere, -NumInfectedSevere / TimeInfectedSevereToInfectedCritical)
                .asDiagonal();
        def.diagonal(1).setConstant(NumInfectedSevere / TimeInfectedSevereToInfectedCritical);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedSevereToInfectedCritical";
    }
};

/**
 * @brief Transition Matrix of the InfectedSevere compartment.
 */
struct TransitionMatrixInfectedSevereToRecovered {
    using Type = Eigen::MatrixXd;
    /**
     * @param[in] NumInfectedSevere Number of subcompartiments of the InfectedSevere compartment.
     * @param[in] TimeInfectedSevereToRecovered  Average time being in the Hospital before treated by ICU or recover in day unit.
     */
    static Type get_default(int NumInfectedSevere, ScalarType TimeInfectedSevereToRecovered = 1.)
    {
        Eigen::MatrixXd def =
            Eigen::VectorXd::Constant(NumInfectedSevere, -NumInfectedSevere / TimeInfectedSevereToRecovered)
                .asDiagonal();
        def.diagonal(1).setConstant(NumInfectedSevere / TimeInfectedSevereToRecovered);
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
     * @param[in] NumInfectedCritical Number of subcompartiments of the InfectedCritical compartment.
     */
    static Type get_default(int NumInfectedCritical)
    {
        Eigen::VectorXd def = Eigen::VectorXd::Zero(NumInfectedCritical);
        def[0]              = 1.;
        return def;
    }
    static std::string name()
    {
        return "StartingProbabilitiesInfectedCritical";
    }
};

/**
 * @brief Transition Matrix of the InfectedCritical compartment.
 */
struct TransitionMatrixInfectedCriticalToDead {
    using Type = Eigen::MatrixXd;
    /**
     * @param[in] NumInfectedCritical Number of subcompartiments of the InfectedCritical compartment.
     * @param[in] TimeInfectedCriticalToDead Average time treated by ICU before dead in day unit.
     */
    static Type get_default(int NumInfectedCritical, ScalarType TimeInfectedCriticalToDead = 1.)
    {
        Eigen::MatrixXd def =
            Eigen::VectorXd::Constant(NumInfectedCritical, -NumInfectedCritical / TimeInfectedCriticalToDead)
                .asDiagonal();
        def.diagonal(1).setConstant(NumInfectedCritical / TimeInfectedCriticalToDead);
        return def;
    }
    static std::string name()
    {
        return "TransitionMatrixInfectedCriticalToDead";
    }
};

/**
 * @brief Transition Matrix of the InfectedCritical compartment.
 */
struct TransitionMatrixInfectedCriticalToRecovered {
    using Type = Eigen::MatrixXd;
    /**
     * @param[in] NumInfectedCritical Number of subcompartiments of the InfectedCritical compartment.
     * @param[in] TimeInfectedCriticalToRecovered Average time treated by ICU before recovery in day unit.
     */
    static Type get_default(int NumInfectedCritical, ScalarType TimeInfectedCriticalToRecovered = 1.)
    {
        Eigen::MatrixXd def =
            Eigen::VectorXd::Constant(NumInfectedCritical, -NumInfectedCritical / TimeInfectedCriticalToRecovered)
                .asDiagonal();
        def.diagonal(1).setConstant(NumInfectedCritical / TimeInfectedCriticalToRecovered);
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
     * @brief checks whether all Parameters satisfy their corresponding constraints and throws errors, if they do not.
     * @return Returns true if one (or more) constraint(s) are not satisfied, otherwise false. 
     */
    bool check_constraints() const
    {
        log_warning("Tests for sum of alpha=1 and >0 missing and -A times vecotr(1) bigger than zero.");
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
            log_error("Constraint check: Parameter  RiskOfInfectionFromSymptomatic smaller {:d} or larger {:d}", 0, 1);
            return true;
        }

        if (this->get<Seasonality>() < 0.0 || this->get<Seasonality>() > 0.5) {
            log_warning("Constraint check: Parameter Seasonality should lie between {:0.4f} and {:.4f}", 0.0, 0.5);
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

#endif // LCT_SECIR_PARAMS_H