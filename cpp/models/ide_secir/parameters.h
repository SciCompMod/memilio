/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler, Lena Ploetzke
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
#ifndef IDE_SECIR_PARAMS_H
#define IDE_SECIR_PARAMS_H

#include "memilio/config.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/parameter_set.h"
#include "ide_secir/infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/math/smoother.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/age_group.h"

#include <memory>
#include <cstddef>
#include <vector>

namespace mio
{
namespace isecir
{

/**********************************************
* Define Parameters of the IDE-SECIHURD model *
**********************************************/

/**
 * @brief Transition distribution for each transition in #InfectionTransition.
 *
 * For each transition, the corresponding transition distribution can be chosen independently. 
 * The choice of distributions determines how many initial time points are required to initialize the model, see
 * get_global_support_max() in models/ide_secir/model.h.
 *
 * As a default we use SmootherCosine functions for all transitions with m_parameter=2.
 */
struct TransitionDistributions {

    using Type = CustomIndexArray<std::vector<StateAgeFunctionWrapper<ScalarType>>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        SmootherCosine smoothcos(2.0);
        StateAgeFunctionWrapper<ScalarType> delaydistribution(smoothcos);
        std::vector<StateAgeFunctionWrapper<ScalarType>> state_age_function_vector((int)InfectionTransition::Count,
                                                                                   delaydistribution);
        return Type(size, state_age_function_vector);
    }

    static std::string name()
    {
        return "TransitionDistributions";
    }
};

/**
 * @brief Defines the probability for each possible transition to take this flow/transition.
 */
struct TransitionProbabilities {
    /*For consistency, also define TransitionProbabilities for each transition in #InfectionTransition. 
    Transition Probabilities should be set to 1 if there is no possible other flow from starting compartment.*/
    using Type = CustomIndexArray<std::vector<ScalarType>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        std::vector<ScalarType> probs((int)InfectionTransition::Count, 0.5);
        // Set the following probablities to 1 as there is no other option to go anywhere else.
        probs[Eigen::Index(InfectionTransition::SusceptibleToExposed)]        = 1;
        probs[Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
        return Type(size, probs);
    }

    static std::string name()
    {
        return "TransitionProbabilities";
    }
};

/**
 * @brief The contact patterns within the society are modelled using an UncertainContactMatrix.
 */
struct ContactPatterns {
    using Type = UncertainContactMatrix<ScalarType>;

    static Type get_default(AgeGroup size)
    {
        ContactMatrixGroup<ScalarType> contact_matrix =
            ContactMatrixGroup<ScalarType>(1, static_cast<Eigen::Index>((size_t)size));
        contact_matrix[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(
            static_cast<Eigen::Index>((size_t)size), static_cast<Eigen::Index>((size_t)size), 10.));
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
    using Type = CustomIndexArray<StateAgeFunctionWrapper<ScalarType>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        ConstantFunction constfunc(1.0);
        return Type(size, constfunc);
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

    using Type = CustomIndexArray<StateAgeFunctionWrapper<ScalarType>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        ConstantFunction constfunc(1.0);
        return Type(size, constfunc);
    }
    static std::string name()
    {
        return "RelativeTransmissionNoSymptoms";
    }
};

/**
* @brief The risk of infection from symptomatic cases in the SECIR model.
*/
struct RiskOfInfectionFromSymptomatic {
    using Type = CustomIndexArray<StateAgeFunctionWrapper<ScalarType>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        ConstantFunction constfunc(1.0);
        return Type(size, constfunc);
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSymptomatic";
    }
};

/**
 * @brief Sets the day in a year at which a simulation with an IDE-SECIR model is started.
 *
 * The value 0.0 corresponds to the 1st of January at 0:00 am, 31.25 is the 1st of February at 6:00 am, and so on.
 * The start day defines in which season the simulation is started.
 * If the start day is 180 and simulation takes place from t0=0 to
 * tmax=100 the days 180 to 280 of the year are simulated.
 * The parameter is used in the formula of the seasonality in the model.
 */
struct StartDay {
    using Type = ScalarType;
    static Type get_default(AgeGroup)
    {
        return 0.;
    }
    static std::string name()
    {
        return "StartDay";
    }
};

/**
 * @brief The seasonality parameter k in the IDE-SECIR model.
 * The formula for the seasonality used in the model is given as (1+k*sin()) where the sine
 * curve is below one in summer and above one in winter.
 */
struct Seasonality {
    using Type = ScalarType;
    static Type get_default(AgeGroup)
    {
        return Type(0.);
    }
    static std::string name()
    {
        return "Seasonality";
    }
};

// Define Parameterset for IDE-SECIR model.
using ParametersBase =
    ParameterSet<TransitionDistributions, TransitionProbabilities, ContactPatterns, TransmissionProbabilityOnContact,
                 RelativeTransmissionNoSymptoms, RiskOfInfectionFromSymptomatic, StartDay, Seasonality>;

/**
 * @brief Parameters of an age-resolved SECIR/SECIHURD model.
 */
class Parameters : public ParametersBase
{
public:
    Parameters(AgeGroup num_agegroups)
        : ParametersBase(num_agegroups)
        , m_num_groups{num_agegroups}
    {
    }

    /**
     * @brief Checks whether all Parameters satisfy their corresponding constraints and logs an error.
     * @param dt Time step size which is used to get model specific StateAgeFunction%s support.
     * @return Returns true if one (or more) constraint(s) are not satisfied, otherwise false.
     */
    bool check_constraints() const
    {
        for (AgeGroup group = AgeGroup(0); group < m_num_groups; ++group) {
            // For parameters potentially depending on the infectious age, values are checked
            // equidistantly on a realistic maximum window.
            // Please note that this is an incomplete check on correctness.
            size_t infectious_window_check = 50; // parameter defining minimal window on x-axis
            for (size_t i = 0; i < infectious_window_check; i++) {
                if (this->get<TransmissionProbabilityOnContact>()[group].eval((ScalarType)i) < 0.0 ||
                    this->get<TransmissionProbabilityOnContact>()[group].eval((ScalarType)i) > 1.0) {
                    log_error("Constraint check: TransmissionProbabilityOnContact smaller {} or larger {} at some "
                              "time {}",
                              0, 1, i);
                    return true;
                }
            }

            for (size_t i = 0; i < infectious_window_check; i++) {
                if (this->get<RelativeTransmissionNoSymptoms>()[group].eval((ScalarType)i) < 0.0 ||
                    this->get<RelativeTransmissionNoSymptoms>()[group].eval((ScalarType)i) > 1.0) {
                    log_error("Constraint check: RelativeTransmissionNoSymptoms smaller {} or larger {} at some "
                              "time {}",
                              0, 1, i);
                    return true;
                }
            }

            for (size_t i = 0; i < infectious_window_check; i++) {
                if (this->get<RiskOfInfectionFromSymptomatic>()[group].eval((ScalarType)i) < 0.0 ||
                    this->get<RiskOfInfectionFromSymptomatic>()[group].eval((ScalarType)i) > 1.0) {
                    log_error("Constraint check: RiskOfInfectionFromSymptomatic smaller {} or larger {} at some "
                              "time {}",
                              0, 1, i);
                    return true;
                }
            }

            for (size_t i = 0; i < (int)InfectionTransition::Count; i++) {
                if (this->get<TransitionProbabilities>()[group][i] < 0.0 ||
                    this->get<TransitionProbabilities>()[group][i] > 1.0) {
                    log_error("Constraint check: One parameter in TransitionProbabilities smaller {} or larger {}",
                              0, 1);
                    return true;
                }
            }

            if (!floating_point_equal(
                    this->get<TransitionProbabilities>()[group][(int)InfectionTransition::SusceptibleToExposed], 1.0,
                    1e-14)) {
                log_error("Constraint check: Parameter transition probability for SusceptibleToExposed unequal to {}",
                          1);
                return true;
            }

            if (!floating_point_equal(
                    this->get<TransitionProbabilities>()[group][(int)InfectionTransition::ExposedToInfectedNoSymptoms],
                    1.0, 1e-14)) {
                log_error("Constraint check: Parameter transition probability for ExposedToInfectedNoSymptoms unequal "
                          "to {}",
                          1);
                return true;
            }

            if (!floating_point_equal(this->get<TransitionProbabilities>()[group][(
                                          int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] +
                                          this->get<TransitionProbabilities>()[group][(
                                              int)InfectionTransition::InfectedNoSymptomsToRecovered],
                                      1.0, 1e-14)) {
                log_error(
                    "Constraint check: Sum of transition probability for InfectedNoSymptomsToInfectedSymptoms and "
                    "InfectedNoSymptomsToRecovered not equal to {}",
                    1);
                return true;
            }

            if (!floating_point_equal(this->get<TransitionProbabilities>()[group][(
                                          int)InfectionTransition::InfectedSymptomsToInfectedSevere] +
                                          this->get<TransitionProbabilities>()[group][(
                                              int)InfectionTransition::InfectedSymptomsToRecovered],
                                      1.0, 1e-14)) {
                log_error("Constraint check: Sum of transition probability for InfectedSymptomsToInfectedSevere and "
                          "InfectedSymptomsToRecovered not equal to {}",
                          1);
                return true;
            }

            if (!floating_point_equal(this->get<TransitionProbabilities>()[group][(
                                          int)InfectionTransition::InfectedSevereToInfectedCritical] +
                                          this->get<TransitionProbabilities>()[group][(
                                              int)InfectionTransition::InfectedSevereToRecovered],
                                      1.0, 1e-14)) {
                log_error("Constraint check: Sum of transition probability for InfectedSevereToInfectedCritical and "
                          "InfectedSevereToRecovered not equal to {}",
                          1);
                return true;
            }

            if (!floating_point_equal(
                    this->get<TransitionProbabilities>()[group][(int)InfectionTransition::InfectedCriticalToDead] +
                        this->get<TransitionProbabilities>()[group]
                                                            [(int)InfectionTransition::InfectedCriticalToRecovered],
                    1.0, 1e-14)) {
                log_error("Constraint check: Sum of transition probability for InfectedCriticalToDead and "
                          "InfectedCriticalToRecovered not equal to {}",
                          1);
                return true;
            }

            /* The first entry of TransitionDistributions is not checked because the distribution S->E is never used 
            (and it makes no sense to use the distribution). The support does not need to be valid.*/
            for (size_t i = 1; i < (int)InfectionTransition::Count; i++) {
                if (floating_point_less(this->get<TransitionDistributions>()[group][i].get_support_max(10), 0.0,
                                        1e-14)) {
                    log_error("Constraint check: One function in TransitionDistributions has invalid support.");
                    return true;
                }
            }
        }

        if (this->get<Seasonality>() < 0.0 || this->get<Seasonality>() > 0.5) {
            log_warning("Constraint check: Parameter Seasonality should lie between {} and {}", 0.0, 0.5);
            return true;
        }

        return false;
    }

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

private:
    Parameters(ParametersBase&& base)
        : ParametersBase(std::move(base))
        , m_num_groups(get<ContactPatterns>().get_cont_freq_mat().get_num_groups())
    {
    }

private:
    AgeGroup m_num_groups;
};
} // namespace isecir

} // namespace mio

#endif // IDE_SECIR_PARAMS_H