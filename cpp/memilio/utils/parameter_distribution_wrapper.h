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
#ifndef PARAMETER_DISTRIBUTION_WRAPPER_H
#define PARAMETER_DISTRIBUTION_WRAPPER_H

#include "memilio/io/io.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"
#include "parameter_distributions.h"
#include <memory>
#include <string>

namespace mio
{

class AbstractParameterDistribution
{
public:
    template <class Impl>
    AbstractParameterDistribution(Impl&& dist)
        : m_dist(std::make_shared<Impl>(std::move(dist)))
        , sample_impl1([](void* dist, RandomNumberGenerator& rng) {
            return static_cast<Impl*>(dist)->get_sample(rng);
        })
        , sample_impl2([](void* dist, abm::PersonalRandomNumberGenerator& rng) {
            return static_cast<Impl*>(dist)->get_sample(rng);
        })
    {
    }

    AbstractParameterDistribution(AbstractParameterDistribution& other)
        : m_dist(other.m_dist)
        , sample_impl1(other.sample_impl1)
        , sample_impl2(other.sample_impl2)
    {
    }

    AbstractParameterDistribution(AbstractParameterDistribution&& other)
        : m_dist(other.m_dist)
        , sample_impl1(other.sample_impl1)
        , sample_impl2(other.sample_impl2)
    {
    }

    AbstractParameterDistribution(const AbstractParameterDistribution& other)
        : m_dist(other.m_dist)
        , sample_impl1(other.sample_impl1)
        , sample_impl2(other.sample_impl2)
    {
    }

    AbstractParameterDistribution()
        : m_dist(nullptr)
        , sample_impl1([](void* /*dist*/, RandomNumberGenerator& /*rng*/) {
            log_critical("AbstractParameterDistribution does not hold a distribution.");
            exit(static_cast<int>(StatusCode::UnknownError));
            return -1.;
        })
        , sample_impl2([](void* /*dist*/, abm::PersonalRandomNumberGenerator& /*rng*/) {
            log_critical("AbstractParameterDistribution does not hold a distribution.");
            exit(static_cast<int>(StatusCode::UnknownError));
            return -1.;
        })
    {
    }

    AbstractParameterDistribution& operator=(AbstractParameterDistribution&& other) = default;

    AbstractParameterDistribution& operator=(const AbstractParameterDistribution& other) = default;

    double get(RandomNumberGenerator& rng) const
    {
        return sample_impl1(m_dist.get(), rng);
    }

    double get(abm::PersonalRandomNumberGenerator& rng) const
    {
        return sample_impl2(m_dist.get(), rng);
    }

    std::vector<double> params() const
    {
        return static_cast<ParameterDistribution*>(m_dist.get())->params();
    }

    /**
     * serialize a AbstractParameterDistribution.
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        static_cast<ParameterDistribution*>(m_dist.get())->serialize(io);
    }

private:
    std::shared_ptr<void> m_dist;
    double (*sample_impl1)(void*, RandomNumberGenerator&);
    double (*sample_impl2)(void*, abm::PersonalRandomNumberGenerator&);
};

/**
 * deserialize a AbstractParameterDistribution.
 * @see mio::deserialize
 */
template <class IOContext>
IOResult<AbstractParameterDistribution> deserialize_internal(IOContext& io, Tag<AbstractParameterDistribution>)
{
    auto obj  = io.expect_object("ParameterDistribution");
    auto type = obj.expect_element("Type", Tag<std::string>{});
    auto a    = type.value();
    std::cout << a << std::endl;
    if (type) {
        if (type.value() == "Uniform") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionUniform::deserialize_elements(io, obj));
            return mio::success(AbstractParameterDistribution(std::move(r)));
        }
        else if (type.value() == "Normal") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionNormal::deserialize_elements(io, obj));
            return mio::success(AbstractParameterDistribution(std::move(r)));
        }
        else if (type.value() == "LogNormal") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionLogNormal::deserialize_elements(io, obj));
            return mio::success(AbstractParameterDistribution(std::move(r)));
        }
        else if (type.value() == "Exponential") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionExponential::deserialize_elements(io, obj));
            return mio::success(AbstractParameterDistribution(std::move(r)));
        }
        else if (type.value() == "Constant") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionConstant::deserialize_elements(io, obj));
            return mio::success(AbstractParameterDistribution(std::move(r)));
        }
        else {
            return failure(StatusCode::InvalidValue, "Type of ParameterDistribution in AbstractParameterDistribution" +
                                                         type.value() + " not valid.");
        }
    }
    return failure(type.error());
}

class ParameterDistributionWrapper
{
public:
    ParameterDistributionWrapper()
        : m_dist(nullptr)
    {
    }

    ParameterDistributionWrapper(ParameterDistribution& dist)
        : m_dist(std::unique_ptr<ParameterDistribution>(dist.clone()))
    {
    }

    ParameterDistributionWrapper(const ParameterDistributionWrapper& other)
    {
        m_dist = (other.m_dist == nullptr) ? nullptr : std::unique_ptr<ParameterDistribution>(other.m_dist->clone());
    }

    ParameterDistributionWrapper(ParameterDistributionWrapper&& other)
    {
        m_dist = (other.m_dist == nullptr) ? nullptr : std::unique_ptr<ParameterDistribution>(other.m_dist->clone());
    }

    ParameterDistributionWrapper& operator=(ParameterDistributionWrapper const& other)
    {
        m_dist = (other.m_dist == nullptr) ? nullptr : std::unique_ptr<ParameterDistribution>(other.m_dist->clone());
        return *this;
    }

    ParameterDistributionWrapper& operator=(ParameterDistributionWrapper&& other)
    {
        m_dist = (other.m_dist == nullptr) ? nullptr : std::unique_ptr<ParameterDistribution>(other.m_dist->clone());
        return *this;
    };

    ~ParameterDistributionWrapper() = default;

    std::vector<double> params() const
    {
        if (m_dist == nullptr) {
            log_error("Distribution is not defined. Parameters cannot be deduced.");
        }
        return m_dist->params();
    }

    template <class RNG>
    double get(RNG& rng)
    {
        if (m_dist == nullptr) {
            log_error("Distribution is not defined. Value cannot be sampled.");
        }
        return m_dist->get_sample(rng);
    }

    /**
     * serialize this.
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        m_dist->serialize(io);
    }

private:
    std::unique_ptr<ParameterDistribution> m_dist;
};

/**
 * deserialize a ParameterDistributionWrapper.
 * @see mio::deserialize
 */
template <class IOContext>
IOResult<ParameterDistributionWrapper> deserialize_internal(IOContext& io, Tag<ParameterDistributionWrapper>)
{

    auto obj  = io.expect_object("ParameterDistribution");
    auto type = obj.expect_element("Type", Tag<std::string>{});
    if (type) {
        if (type.value() == "Uniform") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionUniform::deserialize_elements(io, obj));
            return ParameterDistributionWrapper(r);
        }
        else if (type.value() == "Normal") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionNormal::deserialize_elements(io, obj));
            return ParameterDistributionWrapper(r);
        }
        else if (type.value() == "LogNormal") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionLogNormal::deserialize_elements(io, obj));
            return ParameterDistributionWrapper(r);
        }
        else if (type.value() == "Exponential") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionExponential::deserialize_elements(io, obj));
            return ParameterDistributionWrapper(r);
        }
        else if (type.value() == "Constant") {
            BOOST_OUTCOME_TRY(auto&& r, ParameterDistributionConstant::deserialize_elements(io, obj));
            return ParameterDistributionWrapper(r);
        }
        else {
            return failure(StatusCode::InvalidValue, "Type of ParameterDistribution in ParameterDistributionWrapper" +
                                                         type.value() + " not valid.");
        }
    }
    return failure(type.error());
}

} // namespace mio

#endif //PARAMETER_DISTRIBUTION_WRAPPER_H
