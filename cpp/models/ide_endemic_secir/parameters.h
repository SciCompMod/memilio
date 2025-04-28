#ifndef IDE_END_SECIR_PARAMS_H
#define IDE_END_SECIR_PARAMS_H

#include "memilio/config.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/parameter_set.h"
#include "ide_endemic_secir/infection_state.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/epidemiology/uncertain_matrix.h"

#include <memory>
#include <cstddef>
#include <vector>

namespace mio
{
namespace endisecir
{

/**************************************************
* Define Parameters of the IDE-END-SECIHURD model *
**************************************************/

/**
 * @brief Transitions distribution for each transition in #Infection Transition.
 *
 * we use as a default SmootherCosine functions for all transitions with m_paramter=2.
 */

struct TransitionDistributions {
    using Type = std::vector<StateAgeFunctionWrapper>;
    static Type get_default()
    {
        SmootherCosine smoothcos(2.0);
        StateAgeFunctionWrapper delaydistribution(smoothcos);
        return Type((int)InfectionTransition::Count, delaydistribution);
    }

    static std::string mame()
    {
        return "TransitionDistributions";
    }
};

/**
 * @brief Defines the probability for each possible transitions to take this transition-
 */

struct TransitionProbabilities {
    /*For consistency, also define TransitionProbabilities for each transition in #InfectionTransition. 
    Transition Probabilities should be set to 1 if there is no possible other flow from starting compartment.*/
    using Type = std::vector<ScalarType>;
    static Type get_default()
    {
        std::vector<ScalarType> probs((int)InfectionTransition::Count, 0.5);
        // Set the following probablities to 1 as there is no other option to go anywhere else.
        probs[Eigen::Index(InfectionTransition::SusceptibleToExposed)]        = 1;
        probs[Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
        return Type(probs);
    }
    static std::string name()
    {
        return "TransitionsProbabilities";
    }
};

/**
 * @brief The contact patterns within the society are modelled using an UncertainContactMatrix.
 */
struct ContactPatterns {
    using Type = UncertainContactMatrix<double>;
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
* @brief Probability of getting infected from a contact.
*/
struct TransmissionProbabilityOnContact {
    using Type = StateAgeFunctionWrapper;
    static Type get_default()
    {
        ConstantFunction constfunc(1.0);
        return Type(constfunc);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContact";
    }
};

/**
* @brief The relative InfectedNoSymptoms infectability.
*/
struct RelativeTransmissionNoSymptoms {
    using Type = StateAgeFunctionWrapper;
    static Type get_default()
    {
        ConstantFunction constfunc(1.0);
        return Type(constfunc);
    }
    static std::string name()
    {
        return "RelativeTrabsnissionNoSymptoms";
    }
};

/**
* @brief The risk of infection from symptomatic cases in the SECIR model.
*/
struct RiskOfInfectionFromSymptomatic {
    using Type = StateAgeFunctionWrapper;
    static Type get_default()
    {
        ConstantFunction constfunc(1.0);
        return Type(constfunc);
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSymptomatic";
    }
};

/** 
 * @brief The natural birth rate.
 */
struct NaturalBirthRate {
    using Type = ScalarType;
    static Type get_default()
    {
        return Type(5e-5);
    }
    static std::string name()
    {
        return "NaturalBirthRate";
    }
};

/**
 * @brief The natural death rate.
 */
struct NaturalDeathRate {
    using Type = ScalarType;
    static Type get_default()
    {
        return Type(3e-5);
    }
    static std::string name()
    {
        return "NaturalDeathRate";
    }
};

// Define Parameterset for IDE-END-SECIR model.
using ParametersBase =
    ParameterSet<TransitionDistributions, TransitionProbabilities, ContactPatterns, TransmissionProbabilityOnContact,
                 RelativeTransmissionNoSymptoms, RiskOfInfectionFromSymptomatic, NaturalBirthRate, NaturalDeathRate>;

/**
 * @brief Parameters of an endemic SECIR/SECIHURD model.
 */

class Parameters : public ParametersBase
{
public:
    /**
     * @brief Checks whether all Parameters satisfy their corresponding constraints and logs an error.
     */
    bool check_constraints() const
    {
        // For paramerers potentially depending on the infectious age, the values are checked equidistantly on a
        // realistic maximum window.
        size_t infectious_window_check = 50; // parameter defining minmal window on x-axis.

        //Check if all the probabilitys are between 0 and 1.
        for (size_t i = 0; i < infectious_window_check; i++) {
            if (this->get<TransmissionProbabilityOnContact>().eval((static_cast<ScalarType>(i))) < 0.0 ||
                this->get<TransmissionProbabilityOnContact>().eval(static_cast<ScalarType>(i)) > 1.0) {
                log_error(
                    "Constraint check: TransmissionProbabilityOnContact smaller {:d} or larger {:d} at some time {:d}",
                    0, 1, i);
                return true;
            }
        }

        for (size_t i = 0; i < infectious_window_check; i++) {
            if (this->get<RelativeTransmissionNoSymptoms>().eval((static_cast<ScalarType>(i))) < 0.0 ||
                this->get<RelativeTransmissionNoSymptoms>().eval(static_cast<ScalarType>(i)) > 1.0) {
                log_error(
                    "Constraint check: RelativeTransmissionNoSymptoms smaller {:d} or larger {:d} at some time {:d}", 0,
                    1, i);
                return true;
            }
        }

        for (size_t i = 0; i < infectious_window_check; i++) {
            if (this->get<RiskOfInfectionFromSymptomatic>().eval((static_cast<ScalarType>(i))) < 0.0 ||
                this->get<RiskOfInfectionFromSymptomatic>().eval(static_cast<ScalarType>(i)) > 1.0) {
                log_error(
                    "Constraint check: RiskOfInfectionFromSymptomatic smaller {:d} or larger {:d} at some time {:d}", 0,
                    1, i);
                return true;
            }
        }

        for (size_t i = 0; i < static_cast<int>(InfectionTransition::Count); i++) {
            if (this->get<TransitionProbabilities>()[i] < 0.0 || this->get<TransitionProbabilities>()[i] > 1.0) {
                log_error("Constraint check: One parameter in TransitionProbabilities smaller {:d} or larger {:d}", 0,
                          1);
                return true;
            }
        }

        // The TransitionProbabilities SusceptibleToExposed and ExposedToInfectedNoSymptoms should be 1.
        if (!floating_point_equal(
                this->get<TransitionProbabilities>()[static_cast<int>(InfectionTransition::SusceptibleToExposed)], 1.0,
                1e-14)) {
            log_error("Constraint check: Parameter transition probability for SusceptibleToExposed unequal to {:d}", 1);
            return true;
        }

        if (!floating_point_equal(this->get<TransitionProbabilities>()[static_cast<int>(
                                      InfectionTransition::ExposedToInfectedNoSymptoms)],
                                  1.0, 1e-14)) {
            log_error("Constraint check: Parameter transition probability for ExposedToInfectedNoSymptoms unequal "
                      "to {:d}",
                      1);
            return true;
        }

        // All TransitionProbabilities where the Transition starts in the same compartment should sum up to 1.

        if (!floating_point_equal(this->get<TransitionProbabilities>()[static_cast<int>(
                                      InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] +
                                      this->get<TransitionProbabilities>()[static_cast<int>(
                                          InfectionTransition::InfectedNoSymptomsToRecovered)],
                                  1.0, 1e-14)) {
            log_error("Constraint check: Sum of transition probability for InfectedNoSymptomsToInfectedSymptoms and "
                      "InfectedNoSymptomsToRecovered not equal to {:d}",
                      1);
            return true;
        }

        if (!floating_point_equal(this->get<TransitionProbabilities>()[static_cast<int>(
                                      InfectionTransition::InfectedSymptomsToInfectedSevere)] +
                                      this->get<TransitionProbabilities>()[static_cast<int>(
                                          InfectionTransition::InfectedSymptomsToRecovered)],
                                  1.0, 1e-14)) {
            log_error("Constraint check: Sum of transition probability for InfectedSymptomsToInfectedSevere and "
                      "InfectedSymptomsToRecovered not equal to {:d}",
                      1);
            return true;
        }

        if (!floating_point_equal(this->get<TransitionProbabilities>()[static_cast<int>(
                                      InfectionTransition::InfectedSevereToInfectedCritical)] +
                                      this->get<TransitionProbabilities>()[static_cast<int>(
                                          InfectionTransition::InfectedSevereToRecovered)],
                                  1.0, 1e-14)) {
            log_error("Constraint check: Sum of transition probability for InfectedSevereToInfectedCritical and "
                      "InfectedSevereToRecovered not equal to {:d}",
                      1);
            return true;
        }

        if (!floating_point_equal(
                this->get<TransitionProbabilities>()[static_cast<int>(InfectionTransition::InfectedCriticalToDead)] +
                    this->get<TransitionProbabilities>()[static_cast<int>(
                        InfectionTransition::InfectedCriticalToRecovered)],
                1.0, 1e-14)) {
            log_error("Constraint check: Sum of transition probability for InfectedCriticalToDead and "
                      "InfectedCriticalToRecovered not equal to {:d}",
                      1);
            return true;
        }

        /* The first entry of TransitionDistributions is not checked because the distribution S->E is never used 
            (and it makes no sense to use the distribution). The support does not need to be valid.*/
        for (size_t i = 1; i < static_cast<int>(InfectionTransition::Count); i++) {
            if (floating_point_less(this->get<TransitionDistributions>()[i].get_support_max(10), 0.0, 1e-14)) {
                log_error("Constraint check: One function in TransitionDistributions has invalid support.");
                return true;
            }
        }

        if (this->get<NaturalBirthRate>() < 0.0) {
            log_warning("Constraint check: Parameter NaturalBirthRate should be greater than {:d}", 0);
            return true;
        }

        if (this->get<NaturalDeathRate>() < 0.0) {
            log_warning("Constraint check: Parameter NaturalDeathRate should be greater than {:d}", 0);
            return true;
        }

        return false;
    }
};

} // namespace endisecir

} // namespace mio

#endif //IDE_END_SECIR_PARAMS_H