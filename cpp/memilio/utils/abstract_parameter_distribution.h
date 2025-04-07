/* 
* Copyright (C) 2020-2025 MEmilio
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
#ifndef ABSTRACT_PARAMETER_DISTRIBUTION_H
#define ABSTRACT_PARAMETER_DISTRIBUTION_H

#include "memilio/io/io.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"
#include "parameter_distributions.h"
#include <memory>
#include <string>

namespace mio
{

/**
 * @brief This class represents an arbitrary ParameterDistribution.
 * @see mio::ParameterDistribution
 * This class can for instance be used for model parameters that should have an arbitrary distribution.
 */
class AbstractParameterDistribution
{
public:
    /**
     * The implementation handed to the constructor should have get_sample function 
     * overloaded with mio::RandomNumberGenerator and mio::abm::PersonalRandomNumberGenerator as input arguments
     */
    template <class Impl>
    AbstractParameterDistribution(Impl&& dist)
        : m_dist(std::make_shared<Impl>(std::move(dist)))
        , sample_impl1([](void* d, RandomNumberGenerator& rng) {
            return static_cast<Impl*>(d)->get_sample(rng);
        })
        , sample_impl2([](void* d, abm::PersonalRandomNumberGenerator& rng) {
            return static_cast<Impl*>(d)->get_sample(rng);
        })
    {
    }

    AbstractParameterDistribution(AbstractParameterDistribution& other) = default;

    AbstractParameterDistribution(AbstractParameterDistribution&& other) = default;

    AbstractParameterDistribution(const AbstractParameterDistribution& other) = default;

    AbstractParameterDistribution()
        : m_dist(nullptr)
        , sample_impl1([](void* /*dist*/, RandomNumberGenerator& /*rng*/) {
            log_critical("AbstractParameterDistribution does not hold a distribution.");
            if (true) {
                exit(static_cast<int>(StatusCode::UnknownError));
            }
            else {
                return -1.;
            }
        })
        , sample_impl2([](void* /*dist*/, abm::PersonalRandomNumberGenerator& /*rng*/) {
            log_critical("AbstractParameterDistribution does not hold a distribution.");
            if (true) {
                exit(static_cast<int>(StatusCode::UnknownError));
            }
            else {
                return -1.;
            }
        })
    {
    }

    AbstractParameterDistribution& operator=(AbstractParameterDistribution&& other) = default;

    AbstractParameterDistribution& operator=(const AbstractParameterDistribution& other) = default;

    bool operator<(const AbstractParameterDistribution& other) const
    {
        return static_cast<ParameterDistribution*>(m_dist.get())
            ->smaller_impl(*static_cast<ParameterDistribution*>(other.m_dist.get()));
    }

    /**
     * @brief Returns a value sampled with the given distribution.
     * @param[in] rng RandomNumberGenerator used for sampling. 
     */
    double get(RandomNumberGenerator& rng) const
    {
        return sample_impl1(m_dist.get(), rng);
    }

    /**
     * @brief Returns a value sampled with the given distribution.
     * @param[in] rng abm::PersonalRandomNumberGenerator used for sampling. 
     */
    double get(abm::PersonalRandomNumberGenerator& rng) const
    {
        return sample_impl2(m_dist.get(), rng);
    }

    /**
     * @brief Get the parameters of the given distribution.
     */
    std::vector<double> params() const
    {
        return static_cast<ParameterDistribution*>(m_dist.get())->params();
    }

    /**
     * serialize an AbstractParameterDistribution.
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        static_cast<ParameterDistribution*>(m_dist.get())->serialize(io);
    }

private:
    std::shared_ptr<void> m_dist; ///< Underlying distribtuion.
    double (*sample_impl1)(
        void*,
        RandomNumberGenerator&); ///< Sample function of the distribution which gets a RandomNumberGenerator as rng.
    double (*sample_impl2)(
        void*,
        abm::
            PersonalRandomNumberGenerator&); ///< Sample function of the distribution which gets a abm::PersonalRandomNumberGenerator as rng.
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

} // namespace mio

#endif //ABSTRACT_PARAMETER_DISTRIBUTION_H
