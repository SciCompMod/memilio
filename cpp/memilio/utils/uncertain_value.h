/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Martin J. Kuehn, Martin Siggel, Daniel Abele
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
#ifndef UNCERTAINVALUE_H
#define UNCERTAINVALUE_H

#include "memilio/config.h"
#include "memilio/utils/memory.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/random_number_generator.h"
#include <memory>
#include <ostream>

namespace mio
{

/**
 * @brief The UncertainValue class consists of a 
 *        scalar value and a Distribution object
 * 
 * The UncertainValue class represents a model parameter that
 * can take a scalar value but that is subjected to a uncertainty.
 * The uncertainty is represented by a distribution object of kind
 * ParameterDistribution and the current scalar value can be 
 * replaced by drawing a new sample from the the distribution
 * @tparam FP underlying floating point type, e.g., double
 */
template <typename FP = ScalarType>
class UncertainValue
{
public:
    UncertainValue(FP v, const ParameterDistribution& dist)
        : m_value(v)
        , m_dist(dist.clone())
    {
    }

    UncertainValue(FP v = static_cast<FP>(0.))
        : m_value(v)
    {
    }

    UncertainValue(UncertainValue&& other) = default;

    /**
    * @brief Create an UncertainValue by cloning scalar value 
    *        and distribution of another UncertainValue
    */
    UncertainValue(const UncertainValue& other)
        : m_value(other.m_value)
    {
        if (other.m_dist) {
            m_dist.reset(other.m_dist->clone());
        }
    }

    /**
    * @brief Set an UncertainValue from another UncertainValue
    *        containing a scalar and a distribution
    */
    UncertainValue& operator=(const UncertainValue& other)
    {
        UncertainValue tmp(other);
        m_value = tmp.m_value;
        std::swap(m_dist, tmp.m_dist);
        return *this;
    }

    /**
     * @brief Conversion to scalar by returning the scalar contained in UncertainValue
     */
    FP value() const
    {
        return m_value;
    }

    /**
     * @brief Conversion to scalar reference by returning the scalar contained in UncertainValue
     */
    operator FP&()
    {
        return m_value;
    }
    operator const FP&() const
    {
        return m_value;
    }

    /**
     * @brief Set an UncertainValue from a scalar, distribution remains unchanged.
     */
    UncertainValue& operator=(FP v)
    {
        m_value = v;
        return *this;
    }

    /**
     * @brief Sets the distribution of the value.
     *
     * The function uses copy semantics, i.e. it copies
     * the distribution object.
     */
    void set_distribution(const ParameterDistribution& dist)
    {
        m_dist.reset(dist.clone());
    }

    /**
     * @brief Returns the parameter distribution.
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<const ParameterDistribution> get_distribution() const
    {
        return m_dist.get();
    }

    /**
     * @brief Returns the parameter distribution.
     *
     * If it is not set, a nullptr is returned.
     */
    observer_ptr<ParameterDistribution> get_distribution()
    {
        return m_dist.get();
    }

    /**
     * @brief Sets the value by sampling from the distribution
     *        and returns the new value
     *
     * If no distribution is set, the value is not changed.
     */
    FP draw_sample()
    {
        if (m_dist) {
            m_value = m_dist->get_sample(mio::thread_local_rng());
        }

        return m_value;
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("UncertainValue");
        if (!(io.flags() & IOF_OmitValues)) {
            obj.add_element("Value", value());
        }
        if (!(io.flags() & IOF_OmitDistributions)) {
            obj.add_optional("Distribution", get_distribution().get());
        }
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<UncertainValue> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("UncertainValue");
        if (!(io.flags() & IOF_OmitValues) && !(io.flags() & IOF_OmitDistributions)) {
            auto v = obj.expect_element("Value", Tag<FP>{});
            auto d = obj.expect_optional("Distribution", Tag<std::shared_ptr<ParameterDistribution>>{});
            return apply(
                io,
                [](auto&& v_, auto&& d_) {
                    auto uv = UncertainValue(v_);
                    if (d_) {
                        uv.set_distribution(**d_);
                    }
                    return uv;
                },
                v, d);
        }
        else if (!(io.flags() & IOF_OmitValues) && (io.flags() & IOF_OmitDistributions)) {
            auto v = obj.expect_element("Value", Tag<FP>{});
            return apply(
                io,
                [](auto&& v_) {
                    return UncertainValue(v_);
                },
                v);
        }
        else if ((io.flags() & IOF_OmitValues) && !(io.flags() & IOF_OmitDistributions)) {
            auto d = obj.expect_optional("Distribution", Tag<std::shared_ptr<ParameterDistribution>>{});
            return apply(
                io,
                [](auto&& d_) {
                    auto uv = UncertainValue();
                    if (d_) {
                        uv.set_distribution(**d_);
                    }
                    return uv;
                },
                d);
        }
        else {
            return failure(StatusCode::InvalidValue,
                           "Incompatible Flags in IO Context: IOF_OmitValues & IOF_OmitDistributions.");
        }
    }

private:
    FP m_value;
    std::unique_ptr<ParameterDistribution> m_dist;
};

/**
 * @brief Format UncertainValues using their value for logging with spdlog.
 */
template <class FP>
const FP& format_as(const UncertainValue<FP>& uv)
{
    // uses UncertainValue<FP>::operator const FP&() const
    return uv;
}

// gtest printer
// TODO: should be extended when UncertainValue gets operator== that compares distributions as well
template <typename FP = double>
inline void PrintTo(const UncertainValue<FP>& uv, std::ostream* os)
{
    (*os) << uv.value();
}

} // namespace mio

#endif // UNCERTAINVALUE_H
