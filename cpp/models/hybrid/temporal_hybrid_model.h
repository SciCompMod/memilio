/* 
* Copyright (C) 2020-2025 German Aerospace Center (DLR-SC)
*
* Authors: Julia Bicker, René Schmieding
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

#ifndef MIO_TEMPORAL_HYBRID_MODEL_H
#define MIO_TEMPORAL_HYBRID_MODEL_H

#include <algorithm>
#include <functional>
#include <type_traits>

#include "memilio/config.h"

namespace mio
{
namespace hybrid
{

/**
 * @brief Convert one model into another. Some template specilizations of this function can be found in conversion_functions.h
 * @tparam CurrentModel Simulation type of the model currently used.
 * @tparam TargetModel Simulation type of the model the current simulation should be converted to.
 * @param[in,out] CurrentModel Simulation of currently used model.
 * @param[in,out] TargetModel Simulation of target model.
 */
template <class CurrentModel, class TargetModel>
void convert_model(const CurrentModel&, TargetModel&) = delete;

/**
 * @brief A temporal-hybrid simulation. 
 * The temporal-hybrid simulation switches between two models during the course of time according to a given condition. This requires a specilization of the convert_model function for the two models used.
 * @tparam Model1 (Simulation) Type of the first model used.
 * @tparam Model2 (Simulation) Type of the second model used.
 * @tparam ResultType1 Result type of the first model. The results of both models are needed to evalute the switching condition.
 * @tparam ResultType2 Result type of the second model. The results of both models are needed to evalute the switching condition.
 */

template <class Model1, class Model2, class ResultType1, class ResultType2>
class TemporalHybridSimulation
{

public:
    //Functions returning the result/current state of both models
    using result1_function = std::function<ResultType1(const Model1&, ScalarType t)>;
    using result2_function = std::function<ResultType2(const Model2&, ScalarType t)>;

    //Should return true when the simulation should be continued with the model that is not used currently i.e. a switch needs to be applied
    using switching_condition =
        std::function<bool(const ResultType1& state_model1, const ResultType2& state_model2, bool model1_used)>;

    /**
     * @brief Create a temporal-hybrid simulation
     * @param[in] model1 First model/simulation used for the hybrid simulation.
     * @param[in] model2 Second model/simulation used for the hybrid simulation.
     * @param[in] result1 Function returning the result/current state of first model.
     * @param[in] result2 Function returning the result/current state of second model.
     * @param[in] initially_use_model1 Boolean specifying which model to use at simulation start.
     * @param[in] t0 Start time of the simulation.
     * @param[in] dt Timestep with which the switching is checked.
     */
    TemporalHybridSimulation(Model1&& model1, Model2&& model2, const result1_function& result1,
                             const result2_function& result2, bool initially_use_model1, ScalarType t0 = 0,
                             ScalarType dt = 0.1)
        : m_model1(std::move(model1))
        , m_model2(std::move(model2))
        , m_result1(result1)
        , m_result2(result2)
        , m_using_model1(initially_use_model1)
        , m_t(t0)
        , m_dt(dt)
    {
    }

    /**
     * @brief Advance simulation to tmax.
     * @param[in] tmax End time point of the simulation
     * @param[in] switch_model Switching condition that is checked every m_dt step.
     */
    void advance(ScalarType tmax, const switching_condition& switch_model)
    {
        while (m_t < tmax) {
            bool condition = switch_model(get_result_model1(), get_result_model2(), m_using_model1);
            if (m_using_model1 &&
                condition) { //currently model1 is used, but the condition to switch to model2 is fulfilled
                convert_model(m_model1, m_model2);
                m_using_model1 = false;
            }
            else if (
                !m_using_model1 &&
                condition) { //currently model2 is used, but the condition to switch is fulfilled i.e. we need to switch to model1
                convert_model(m_model2, m_model1);
                m_using_model1 = true;
            }
            //else{Switching condition is not fulfilled and currently used model is just advanced}
            ScalarType next_step = std::min(m_dt, tmax - m_t);
            if (m_using_model1) {
                m_model1.advance(m_t + next_step);
            }
            else {
                m_model2.advance(m_t + next_step);
            }
            m_t += next_step;
        }
    }

    /**
     * @brief Get the result of model 1.
     * @return Result of model 1 using the function m_result1.
     */
    ResultType1 get_result_model1() const
    {
        return m_result1(m_model1, m_t);
    }

    /**
     * @brief Get the result of model 2.
     * @return Result of model 2 using the function m_result1.
     */
    ResultType2 get_result_model2() const
    {
        return m_result2(m_model2, m_t);
    }

    /**
     * @brief Returns first model used for the simulation.
     */
    const auto& get_model1() const
    {
        return m_model1;
    }
    auto& get_model1()
    {
        return m_model1;
    }

    /**
     * @brief Returns second model used for the simulation.
     */
    const auto& get_model2() const
    {
        return m_model2;
    }
    auto& get_model2()
    {
        return m_model2;
    }

    /**
     * @brief Returns whether the first model is currently used for simulation.
     */
    const auto& using_model1() const
    {
        return m_using_model1;
    }
    auto& using_model1()
    {
        return m_using_model1;
    }

private:
    Model1 m_model1; ///< First model used for the simulation.
    Model2 m_model2; ///< Second model used for the simulation.
    result1_function m_result1; ///< Result function of first model.
    result2_function m_result2; ///< Result function of second model.
    bool m_using_model1; ///< Boolean specifying whether model 1 is currently used for simulation.
    ScalarType m_t; ///< Current time step.
    ScalarType m_dt; ///< Step size with which the switching condition is checked.
};

} //namespace hybrid

} //namespace mio

#endif //MIO_TEMPORAL_HYBRID_MODEL_H
