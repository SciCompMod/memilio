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
namespace mio
{
namespace hybrid
{

template <class CurrentModel, class TargetModel>
void convert_model(const CurrentModel&, TargetModel&) = delete;

template <class Model1, class Model2, class ResultType1, class ResultType2>
class TemporalHybridSimulation
{

public:
    using result1_function = ResultType1 (*)(const Model1&, double t);
    using result2_function = ResultType2 (*)(const Model2&, double t);

    //Should return true when the simulation should be continued with model2
    using switching_condition =
        std::function<bool(const ResultType1& state_model1, const ResultType2& state_model2, bool model1_used)>;

    TemporalHybridSimulation(Model1& model1, Model2& model2, result1_function result1, result2_function result2,
                             bool initially_use_model1, double t0 = 0, double dt = 0.1)
        : m_model1(model1)
        , m_model2(model2)
        , m_result1(result1)
        , m_result2(result2)
        , m_using_model1(initially_use_model1)
        , m_t(t0)
        , m_dt(dt)
    {
    }

    void advance(double tmax, const switching_condition& switch_model)
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
            double next_step = std::min(m_dt, tmax - m_t);
            if (m_using_model1) {
                m_model1.advance(m_t + next_step);
            }
            else {
                m_model2.advance(m_t + next_step);
            }
            m_t += next_step;
        }
    }

    ResultType1 get_result_model1()
    {
        return m_result1(m_model1, m_t);
    }

    ResultType2 get_result_model2()
    {
        return m_result2(m_model2, m_t);
    }

    const auto& get_model1() const
    {
        return m_model1;
    }

    auto& get_model1()
    {
        return m_model1;
    }

    const auto& get_model2() const
    {
        return m_model2;
    }

    auto& get_model2()
    {
        return m_model2;
    }

    const auto& using_model1() const
    {
        return m_using_model1;
    }

    auto& using_model1()
    {
        return m_using_model1;
    }

private:
    Model1 m_model1;
    Model2 m_model2;
    result1_function m_result1;
    result2_function m_result2;
    bool m_using_model1;
    double m_t;
    double m_dt;
};

} //namespace hybrid

} //namespace mio

#endif //MIO_TEMPORAL_HYBRID_MODEL_H
