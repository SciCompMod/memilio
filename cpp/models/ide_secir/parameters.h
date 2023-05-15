/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include "memilio/utils/parameter_set.h"
#include "ide_secir/infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/smoother.h"

#include <memory>
#include <vector>

namespace mio
{
namespace isecir
{

/**********************************************
* Define Parameters of the IDE-SECIHURD model *
**********************************************/

/**
* @brief Function describing the time spent in a compartment before transiting to next compartment.

* This class defines a function that specifies which proportion of individuals is still in the compartment after a certain state_age 
* (i.e. time after entering the compartment) and has not yet progressed to the next.
* See also function \gamma in Overleaf.
* Currently, you can only use a smoother_cosine() function with different parameters for this purpose.
* TransitionDistribution implements a vector with DelayDistribution%s for each compartment.
*/
struct DelayDistribution {
    /**
    * @brief Default constructor of the class. Default for m_max_support is 2.0 (relatively random.)
    */
    DelayDistribution()
        : m_max_support{2.0}
    {
    }

    /**
     * @brief Construct a new DelayDistribution object
     * 
     * @param init_max_support specifies the right bound and therefore the support of the function.
     */
    DelayDistribution(ScalarType init_max_support)
        : m_max_support{init_max_support}
    {
    }
    /**
     * @brief DelayDistribution is currently defined as a smoothed cosine function.
     * 
     * Used function goes through points (0,1) and (m_max_support,0) and is interpolated in between using a smoothed cosine function.
     * 
     * @param state_age time at which the function should be evaluated
     * @return ScalarType evaluation of the smoother cosine function
     */
    ScalarType Distribution(ScalarType state_age)
    {
        return smoother_cosine(state_age, 0.0, m_max_support, 1.0, 0.0);
    }

    /**
     * @brief Get the m_max_support object
     * 
     * Can be used to access the m_max_support object, which specifies the right bound of the support of the function.
     * 
     * @return ScalarType 
     */
    ScalarType get_max_support() const
    {
        return m_max_support;
    }

    /**
     * @brief Set the m_max_support object
     * 
     * Can be used to set the m_max_support object, which specifies the right bound of the support of the function.
     * 
     * @return ScalarType 
     */
    void set_max_support(ScalarType new_max_support)
    {
        m_max_support = new_max_support;
    }

private:
    ScalarType m_max_support; ///< specifies the right bound of the support of the DelayDistribution.
};

/**
 * @brief Transition distribution for each transition in InfectionTransition.
 * 
 * Note that Distribution from S->E is just a dummy.
 * This transition is calculated in a different way.
 * 
 */
struct TransitionDistributions {
    using Type = std::vector<DelayDistribution>;
    static Type get_default()
    {
        return std::vector<DelayDistribution>((int)InfectionTransition::Count, DelayDistribution());
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
    /*For consistency, also define TransitionProbabilities for each transition in InfectionTransition. 
    Transition Probabilities should be set to 1 if there is no possible other flow from starting compartment.*/
    using Type = std::vector<ScalarType>;
    static Type get_default()
    {
        std::vector<ScalarType> probs((int)InfectionTransition::Count, 0.5);
        // Set the following probablities to 1 as there is no other option to go anywhere else.
        probs[Eigen::Index(InfectionTransition::SusceptibleToExposed)]        = 1;
        probs[Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
        return probs;
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
    using Type = UncertainContactMatrix;

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
 * @brief Function depending on the state age, i.e. the time time already spent in some InfectionState.
 */
struct StateAgeFunction {
    /**
    * @brief Default constructor of the class. Default for m_funcparam is 1.0 (relatively random.)
    */
    StateAgeFunction()
        : m_funcparam{1.0}
    {
    }

    /**
     * @brief Constructs a new StateAgeFunction object
     * 
     * @param init_funcparam specifies the initial function parameter of the function.
     */
    StateAgeFunction(ScalarType init_funcparam)
        : m_funcparam{init_funcparam}
    {
    }

    /**
     * @brief Defines function depending on state_age. The default is an exponential function.
     *
     * The function also depends on the function parameter m_funcparam which allows to further specify the function.
     * (depending on the used function).
     * 
     * @param state_age time at which the function should be evaluated
     * @return ScalarType evaluation of the function at state_age. 
     */
    virtual ScalarType Function(ScalarType state_age)
    {
        return std::exp(-m_funcparam * state_age);
    }

    /**
     * @brief Get the m_funcparam object
     * 
     * Can be used to access the m_funcparam object, which specifies the used function.
     * 
     * @return ScalarType 
     */
    ScalarType get_funcparam()
    {
        return m_funcparam;
    }

    /**
     * @brief Set the m_funcparam object.
     * 
     * Can be used to set the m_funcparam object, which specifies the used function.
     */
    void set_funcparam(ScalarType new_funcparam)
    {
        m_funcparam = new_funcparam;
    }

    /**
     * @brief Clones unique pointer to a StateAgeFunction.
     * 
     * @return std::unique_ptr<StateAgeFunction> unique pointer to a StateAgeFunction
     */
    std::unique_ptr<StateAgeFunction> clone() const
    {
        return std::unique_ptr<StateAgeFunction>(clone_impl());
    }

protected:
    // pure virtual function that implements cloning in derived classes
    virtual StateAgeFunction* clone_impl() const = 0;

    // specifies Function
    ScalarType m_funcparam{};
};

/**
 * @brief Class that defines an exponential decay function depending on the state age.
 */
struct ExponentialDecay : public StateAgeFunction {

    /**
    * @brief Default constructor of the class.
    */
    ExponentialDecay()
        : StateAgeFunction{}
    {
    }

    /**
     * @brief Constructs a new ExponentialDecay object
     * 
     * @param init_funcparam specifies the initial function parameter of the function.
     */
    ExponentialDecay(ScalarType init_funcparam)
        : StateAgeFunction{init_funcparam}
    {
    }

    /**
     * @brief Defines exponential decay function depending on state_age.
     *
     * The parameter m_funcparam defines how fast the exponential function decays.
     * 
     * @param state_age time at which the function should be evaluated
     * @return ScalarType evaluation of the function at state_age. 
     */
    ScalarType Function(ScalarType state_age) override
    {
        return std::exp(-m_funcparam * state_age);
    }

protected:
    /**
     * @brief Implements clone for ExponentialDecay.
     * 
     * @return StateAgeFunction* 
     */
    StateAgeFunction* clone_impl() const override
    {
        return new ExponentialDecay(*this);
    }
};

/**
 * @brief Class that defines an smoother cosine function depending on the state age.
 */
struct SmootherCosine : public StateAgeFunction {

    /**
    * @brief Default constructor of the class.
    */
    SmootherCosine()
        : StateAgeFunction{}
    {
    }

    /**
     * @brief Constructs a new SmootherCosine object
     * 
     * @param init_funcparam specifies the initial function parameter of the function.
     */
    SmootherCosine(ScalarType init_funcparam)
        : StateAgeFunction{init_funcparam}
    {
    }

    /**
     * @brief Defines smoother cosine function depending on state_age.
     *
     * Used function goes through points (0,1) and (m_funcparam,0) and is interpolated in between using a smoothed cosine function.
     * 
     * @param state_age time at which the function should be evaluated
     * @return ScalarType evaluation of the function at state_age. 
     */
    ScalarType Function(ScalarType state_age) override
    {
        return smoother_cosine(state_age, 0.0, m_funcparam, 1.0, 0.0);
    }

protected:
    /**
     * @brief Clones unique pointer to a StateAgeFunction.
     * 
     * @return std::unique_ptr<StateAgeFunction> unique pointer to a StateAgeFunction
     */
    StateAgeFunction* clone_impl() const override
    {
        return new SmootherCosine(*this);
    }
};

/**
 * @brief Wrapper for StateAgeFunction that allows to set a StateAgeFunction from outside. 
 */
struct StateAgeFunctionWrapper {

    /**
    * @brief Default constructor of the class. Sets m_function to constant function 1 as default.
    */
    StateAgeFunctionWrapper()
        : m_function{}
    {
        // Set m_function to a default function, choose constant function 1
        ExponentialDecay expdecay(0);
        m_function = expdecay.clone();
    }

    /**
     * @brief Copy constructor for StateAgeFunctionWrapper. 
     */
    StateAgeFunctionWrapper(StateAgeFunctionWrapper const& other)
        : m_function(other.m_function->clone())
    {
    }

    /**
     * @brief Move constructor for StateAgeFunctionWrapper. 
     */
    StateAgeFunctionWrapper(StateAgeFunctionWrapper&& other) = default;

    /**
     * @brief Copy assignment for StateAgeFunctionWrapper. 
     */
    StateAgeFunctionWrapper& operator=(StateAgeFunctionWrapper const& other)
    {
        m_function = other.m_function->clone();
        return *this;
    }

    /**
     * @brief Move assignment for StateAgeFunctionWrapper. 
     */
    StateAgeFunctionWrapper& operator=(StateAgeFunctionWrapper&& other) = default;

    /**
     * @brief Destructor for StateAgeFunctionWrapper. 
     */
    ~StateAgeFunctionWrapper() = default;

    /**
     * @brief Set the StateAgeFunction object
     *
     * @return std::unique_ptr<StateAgeFunction>
     */
    void set_state_age_function(StateAgeFunction& new_function)
    {
        m_function = new_function.clone();
    }

    /**
     * @brief Accesses Function of m_function.
     *
     * @param state_age time at which the function should be evaluated
     * @return ScalarType evaluation of the function at state_age. 
     */
    ScalarType Function(ScalarType state_age)
    {
        return m_function->Function(state_age);
    }

    /**
     * @brief Get the m_funcparam object of m_function.
     * 
     * @return ScalarType 
     */
    ScalarType get_funcparam()
    {
        return m_function->get_funcparam();
    }

    /**
     * @brief Set the m_funcparam object of m_function. 
     */
    void set_funcparam(ScalarType new_funcparam)
    {
        m_function->set_funcparam(new_funcparam);
    }

private:
    std::unique_ptr<StateAgeFunction> m_function;
};

/**
* @brief Probability of getting infected from a contact.
*/
struct TransmissionProbabilityOnContact {
    using Type = StateAgeFunctionWrapper;
    static Type get_default()
    {
        return StateAgeFunctionWrapper();
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
        return StateAgeFunctionWrapper();
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
    using Type = StateAgeFunctionWrapper;
    static Type get_default()
    {
        return StateAgeFunctionWrapper();
    }
    static std::string name()
    {
        return "RiskOfInfectionFromSymptomatic";
    }
};

// Define Parameterset for IDE SECIR model.
using ParametersBase =
    ParameterSet<TransitionDistributions, TransitionProbabilities, ContactPatterns, TransmissionProbabilityOnContact,
                 RelativeTransmissionNoSymptoms, RiskOfInfectionFromSymptomatic>;

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
     * @return Returns 1 if one constraint is not satisfied, otherwise 0.
     */
    int check_constraints() const
    {
        for (int i = 0; i < (int)InfectionTransition::Count; i++) {
            if (this->get<TransitionProbabilities>()[i] < 0.0 || this->get<TransitionProbabilities>()[i] > 1.0) {
                log_error("Constraint check: One parameter TransitionProbabilities smaller {:d} or larger {:d}", 0, 1);
                return 1;
            }
        }

        if (this->get<TransitionProbabilities>()[(int)InfectionTransition::SusceptibleToExposed] != 1.0) {
            log_error("Constraint check: Parameter transitiion probability for SusceptibleToExposed unequal to {:d}",
                      1);
            return 1;
        }

        if (this->get<TransitionProbabilities>()[(int)InfectionTransition::ExposedToInfectedNoSymptoms] != 1.0) {
            log_error(
                "Constraint check: Parameter transitiion probability for ExposedToInfectedNoSymptoms unequal to {:d}",
                1);
            return 1;
        }

        if (this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] +
                this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedNoSymptomsToRecovered] !=
            1.0) {
            log_error("Constraint check: Sum of transitiion probability for InfectedNoSymptomsToInfectedSymptoms and "
                      "InfectedNoSymptomsToRecovered unequal to {:d}",
                      1);
            return 1;
        }

        if (this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedSymptomsToInfectedSevere] +
                this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedSymptomsToRecovered] !=
            1.0) {
            log_error("Constraint check: Sum of transitiion probability for InfectedSymptomsToInfectedSevere and "
                      "InfectedSymptomsToRecovered unequal to {:d}",
                      1);
            return 1;
        }

        if (this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedSevereToInfectedCritical] +
                this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedSevereToRecovered] !=
            1.0) {
            log_error("Constraint check: Sum of transitiion probability for InfectedSevereToInfectedCritical and "
                      "InfectedSevereToRecovered unequal to {:d}",
                      1);
            return 1;
        }

        if (this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedCriticalToDead] +
                this->get<TransitionProbabilities>()[(int)InfectionTransition::InfectedCriticalToRecovered] !=
            1.0) {
            log_error("Constraint check: Sum of transitiion probability for InfectedCriticalToDead and "
                      "InfectedCriticalToRecovered unequal to {:d}",
                      1);
            return 1;
        }

        return 0;
    }

private:
    Parameters(ParametersBase&& base)
        : ParametersBase(std::move(base))
    {
    }
};

} // namespace isecir
} // namespace mio

#endif // IDE_SECIR_PARAMS_H