/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Julia Bicker
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

#ifndef MIO_ABM_STATE_TRANSITION_DIST_H
#define MIO_ABM_STATE_TRANSITION_DIST_H

#include "abm/personal_rng.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/utils/logging.h"
#include <memory>
#include <utility>
#include <vector>

namespace mio
{
namespace abm
{

struct StateTransitionDist {

    /**
     * @brief Default constructor.
     */
    StateTransitionDist()
    {
    }

    /**
     * @brief Virtual destructor.
     */
    virtual ~StateTransitionDist() = default;

    /**
     * @brief Copy constructor.
     */
    StateTransitionDist(const StateTransitionDist& other) = default;

    /**
     * @brief Move constructor
     */
    StateTransitionDist(StateTransitionDist&& other) = default;

    /**
     * @brief Copy assignment operator.
     */
    StateTransitionDist& operator=(const StateTransitionDist& other) = default;

    /**
     * @brief Move assignment operator.
     */
    StateTransitionDist& operator=(StateTransitionDist&& other) = default;

    /**
     * @brief Returns the distribution parameters as vector.
     */
    virtual std::vector<double> params() const = 0;

    /**
     * @brief Returns name of the distribution.
     */
    virtual std::string name() const = 0;

    /**
     * @brief Returns a sampled value of the distribution.
     */
    virtual double get(PersonalRandomNumberGenerator& rng) = 0;

    /**
     * @brief Returns unique pointer to a StateTransitionDist.
     */
    std::unique_ptr<StateTransitionDist> clone() const
    {
        return std::unique_ptr<StateTransitionDist>(clone_impl());
    }

protected:
    /**
     * @brief Pure virtual method that implements cloning.
     */
    virtual StateTransitionDist* clone_impl() const = 0;
};

struct LogNormal : StateTransitionDist {

    LogNormal(double p1, double p2)
        : StateTransitionDist()
        , m_dist(p1, p2)
    {
    }

    LogNormal(std::vector<double> params)
        : StateTransitionDist()
        , m_dist(params[0], params[1])
    {
    }

    std::vector<double> params() const override
    {
        return {m_dist.params.m(), m_dist.params.s()};
    }

    std::string name() const override
    {
        return "LogNormal";
    }

    double get(PersonalRandomNumberGenerator& rng) override
    {
        return m_dist.get_distribution_instance()(rng, m_dist.params);
    }

    template <class IOContext, class IOObject>
    static IOResult<LogNormal> deserialize_elements(IOContext& io, IOObject& obj)
    {
        auto params = obj.expect_list("params", Tag<double>{});
        auto p      = apply(
            io,
            [](auto&& params_) {
                auto dist = LogNormal(params_);
                return dist;
            },
            params);
        if (p) {
            return success(p.value());
        }
        else {
            return p.as_failure();
        }
    }

    template <class IOContext>
    static IOResult<LogNormal> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("LogNormal");
        return deserialize_elements(io, obj);
    }

protected:
    /**
     * @brief Implements clone for LogNormal.
     * 
     * @return Pointer to StateTransitionDist.
     */
    StateTransitionDist* clone_impl() const override
    {
        return new LogNormal(*this);
    }

private:
    LogNormalDistribution<double>::ParamType m_dist;
};

struct Exponential : StateTransitionDist {

    Exponential(double p1)
        : StateTransitionDist()
        , m_dist(p1)
    {
    }

    Exponential(std::vector<double> params)
        : StateTransitionDist()
        , m_dist(params[0])
    {
    }

    std::vector<double> params() const override
    {
        return {m_dist.params.lambda()};
    }

    std::string name() const override
    {
        return "Exponential";
    }

    double get(PersonalRandomNumberGenerator& rng) override
    {
        return m_dist.get_distribution_instance()(rng, m_dist.params);
    }

    template <class IOContext, class IOObject>
    static IOResult<Exponential> deserialize_elements(IOContext& io, IOObject& obj)
    {
        auto params = obj.expect_list("params", Tag<double>{});
        auto p      = apply(
            io,
            [](auto&& params_) {
                auto dist = Exponential(params_);
                return dist;
            },
            params);
        if (p) {
            return success(p.value());
        }
        else {
            return p.as_failure();
        }
    }

    template <class IOContext>
    static IOResult<Exponential> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Exponential");
        return deserialize_elements(io, obj);
    }

protected:
    /**
     * @brief Implements clone for Exponential.
     * 
     * @return Pointer to StateTransitionDist.
     */
    StateTransitionDist* clone_impl() const override
    {
        return new Exponential(*this);
    }

private:
    ExponentialDistribution<double>::ParamType m_dist;
};

struct Constant : StateTransitionDist {

    Constant(double p1)
        : StateTransitionDist()
        , m_dist(p1)
    {
    }

    Constant(std::vector<double> params)
        : StateTransitionDist()
        , m_dist(params[0])
    {
    }

    std::vector<double> params() const override
    {
        return {m_dist};
    }

    std::string name() const override
    {
        return "Constant";
    }

    double get(PersonalRandomNumberGenerator& /*rng*/) override
    {
        return m_dist;
    }

    template <class IOContext, class IOObject>
    static IOResult<Constant> deserialize_elements(IOContext& io, IOObject& obj)
    {
        auto params = obj.expect_list("params", Tag<double>{});
        auto p      = apply(
            io,
            [](auto&& params_) {
                auto dist = Constant(params_);
                return dist;
            },
            params);
        if (p) {
            return success(p.value());
        }
        else {
            return p.as_failure();
        }
    }

    template <class IOContext>
    static IOResult<Exponential> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Constant");
        return deserialize_elements(io, obj);
    }

protected:
    /**
     * @brief Implements clone for Constant.
     * 
     * @return Pointer to StateTransitionDist.
     */
    StateTransitionDist* clone_impl() const override
    {
        return new Constant(*this);
    }

private:
    double m_dist;
};

struct StateTransitionDistWrapper {

    StateTransitionDistWrapper()
    {
    }

    StateTransitionDistWrapper(StateTransitionDist& dist)
        : m_dist(dist.clone())
    {
    }

    StateTransitionDistWrapper(const StateTransitionDistWrapper& other)
        : m_dist(other.m_dist->clone())
    {
    }

    StateTransitionDistWrapper(StateTransitionDistWrapper&& other)
        : m_dist(other.m_dist->clone())
    {
    }

    StateTransitionDistWrapper& operator=(StateTransitionDistWrapper const& other)
    {
        m_dist = other.m_dist->clone();
        return *this;
    }

    StateTransitionDistWrapper& operator=(StateTransitionDistWrapper&& other)
    {
        m_dist = other.m_dist->clone();
        return *this;
    };

    ~StateTransitionDistWrapper() = default;

    std::vector<double> params() const
    {
        return m_dist->params();
    }

    std::string name() const
    {
        return m_dist->name();
    }

    double get(PersonalRandomNumberGenerator& rng)
    {
        return m_dist->get(rng);
    }

    /**
     * serialize this.
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("StateTransitionDistWrapper");
        obj.add_element("Type", m_dist->name());
        obj.add_list("params", m_dist->params().begin(), m_dist->params().end());
    }

private:
    std::unique_ptr<StateTransitionDist> m_dist;
};

/**
 * deserialize a state transition dist as a shared_ptr.
 * @see mio::deserialize
 */
template <class IOContext>
IOResult<StateTransitionDistWrapper> deserialize_internal(IOContext& io, Tag<StateTransitionDistWrapper>)
{
    auto obj  = io.expect_object("StateTransitionDistWrapper");
    auto type = obj.expect_element("Type", Tag<std::string>{});
    if (type) {
        if (type.value() == "LogNormal") {
            BOOST_OUTCOME_TRY(auto&& r, LogNormal::deserialize_elements(io, obj));
            return StateTransitionDistWrapper(r);
        }
        else if (type.value() == "Exponential") {
            BOOST_OUTCOME_TRY(auto&& r, Exponential::deserialize_elements(io, obj));
            return StateTransitionDistWrapper(r);
        }
        else {
            return failure(StatusCode::InvalidValue,
                           "Type of StateTransitionDistWrapper " + type.value() + " not valid.");
        }
    }
    return failure(type.error());
}
} // namespace abm
} // namespace mio

#endif //MIO_ABM_STATE_TRANSITION_DIST_H
