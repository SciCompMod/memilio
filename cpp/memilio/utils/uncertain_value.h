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
 * @brief Represents a scalar value with an optional uncertainty distribution.
 *
 * The UncertainValue class wraps a floating-point value and an optional
 * ParameterDistribution. You can sample a new value from the distribution
 * via draw_sample(), or leave the distribution unset to keep a fixed value.
 *
 * @tparam FP Underlying floating-point type (e.g., double, AD types).
 */
template <typename FP>
class UncertainValue
{
public:
    /**
     * @brief Default-constructs with value = 0.0 and no distribution.
     */
    UncertainValue(FP v = FP(0.0))
        : m_value(v)
    {
    }

    /**
     * @brief Constructs with an initial value and a distribution.
     *
     * @param v    Initial scalar value.
     * @param dist Distribution to sample from (copied).
     */
    UncertainValue(FP v, const ParameterDistribution& dist)
        : m_value(v)
        , m_dist(dist.clone())
    {
    }

    /**
     * @brief Construct from any arithmetic type (e.g., int, double).
     *
     * Enables implicit conversion from arithmetic types to UncertainValue.
     */
    template <typename T, typename std::enable_if_t<std::is_convertible<T, FP>::value, int> = 0>
    UncertainValue(T v)
        : m_value(static_cast<FP>(v))
    {
    }

    UncertainValue(UncertainValue<FP>&& other) noexcept  = default;
    UncertainValue& operator=(UncertainValue&&) noexcept = default;

    /**
     * @brief Copy-construct from another UncertainValue.
     *
     * Clones both the scalar value and, if present, the distribution.
     */
    UncertainValue(const UncertainValue<FP>& other)
        : m_value(other.m_value)
    {
        if (other.m_dist) {
            m_dist.reset(other.m_dist->clone());
        }
    }

    /**
     * @brief Copy-assign from another UncertainValue.
     *
     * Uses copy-and-swap to copy value and distribution.
     */
    UncertainValue<FP>& operator=(const UncertainValue<FP>& other)
    {
        UncertainValue<FP> tmp(other);
        std::swap(m_value, tmp.m_value);
        std::swap(m_dist, tmp.m_dist);
        return *this;
    }

    /**
     * @brief Assign a new scalar value (distribution unchanged).
     */
    UncertainValue<FP>& operator=(const FP& v)
    {
        m_value = v;
        return *this;
    }

    /**
     * @brief Assign from any arithmetic type (distribution unchanged).
     */
    template <typename T, typename std::enable_if_t<std::is_convertible<T, FP>::value, int> = 0>
    UncertainValue<FP>& operator=(T v)
    {
        m_value = static_cast<FP>(v);
        return *this;
    }

    // Arithmetic compound operators (value-only)
    UncertainValue<FP>& operator+=(const FP& v)
    {
        m_value += v;
        return *this;
    }
    UncertainValue<FP>& operator-=(const FP& v)
    {
        m_value -= v;
        return *this;
    }
    UncertainValue<FP>& operator*=(const FP& v)
    {
        m_value *= v;
        return *this;
    }
    UncertainValue<FP>& operator/=(const FP& v)
    {
        m_value /= v;
        return *this;
    }

    /**
     * @brief Explicit access to the underlying scalar.
     */
    FP value() const
    {
        return m_value;
    }

    /**
     * @brief Implicit convert to FP& (modifiable).
     */
    operator FP&()
    {
        return m_value;
    }

    /**
     * @brief Implicit convert to const FP&.
     */
    operator const FP&() const
    {
        return m_value;
    }

    /**
     * @brief Set the uncertainty distribution.
     *
     * Copies the provided distribution.
     */
    void set_distribution(const ParameterDistribution& dist)
    {
        m_dist.reset(dist.clone());
    }

    /**
     * @brief Get a pointer to the distribution (const).
     *
     * @return A pointer to the parameter distribution, or nullptr if no distribution is set.
     */
    observer_ptr<const ParameterDistribution> get_distribution() const
    {
        return m_dist.get();
    }

    /**
     * @brief Get a pointer to the distribution (modifiable).
     *
     * @return A pointer to the parameter distribution, or nullptr if no distribution is set.
     */
    observer_ptr<ParameterDistribution> get_distribution()
    {
        return m_dist.get();
    }

    /**
     * @brief Sets the value by sampling from the distribution
     *        and returns the new value.
     *
     * If no distribution is set, the value remains unchanged.
     *
     * @return The new (or unchanged) scalar value.
     */
    FP draw_sample()
    {
        if (m_dist) {
            m_value = m_dist->get_sample(mio::thread_local_rng());
        }
        return m_value;
    }

    /**
     * @brief Serialize this object.
     *
     * @tparam IOContext Serialization context type.
     * @param io         IO context.
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
     * @brief Deserialize an UncertainValue.
     *
     * @tparam IOContext Serialization context type.
     * @param io         IO context.
     * @return Deserialized UncertainValue.
     */
    template <class IOContext>
    static IOResult<UncertainValue<FP>> deserialize(IOContext& io)
    {
        auto obj         = io.expect_object("UncertainValue");
        const bool omitV = io.flags() & IOF_OmitValues;
        const bool omitD = io.flags() & IOF_OmitDistributions;

        if (!omitV && !omitD) {
            auto v = obj.expect_element("Value", Tag<FP>{});
            auto d = obj.expect_optional("Distribution", Tag<std::shared_ptr<ParameterDistribution>>{});
            return apply(
                io,
                [](auto&& v_, auto&& d_) {
                    UncertainValue<FP> uv(v_);
                    if (d_) {
                        uv.set_distribution(**d_);
                    }
                    return uv;
                },
                v, d);
        }
        else if (!omitV && omitD) {
            auto v = obj.expect_element("Value", Tag<FP>{});
            return apply(
                io,
                [](auto&& v_) {
                    return UncertainValue<FP>(v_);
                },
                v);
        }
        else if (omitV && !omitD) {
            auto d = obj.expect_optional("Distribution", Tag<std::shared_ptr<ParameterDistribution>>{});
            return apply(
                io,
                [](auto&& d_) {
                    UncertainValue<FP> uv;
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
    std::unique_ptr<ParameterDistribution> m_dist = nullptr;
};

// Helper to format as scalar for spdlog
template <typename FP>
const FP& format_as(const UncertainValue<FP>& uv)
{
    return uv;
}

// gtest printer
template <typename FP>
inline void PrintTo(const UncertainValue<FP>& uv, std::ostream* os)
{
    (*os) << uv.value();
}

// Free-standing comparison and arithmetic operators

// Comparison with scalar on right
template <typename FP, typename T>
inline bool operator>(const UncertainValue<FP>& lhs, T rhs)
{
    return lhs.value() > rhs;
}
template <typename FP, typename T>
inline bool operator<(const UncertainValue<FP>& lhs, T rhs)
{
    return lhs.value() < rhs;
}
template <typename FP, typename T>
inline bool operator>=(const UncertainValue<FP>& lhs, T rhs)
{
    return lhs.value() >= rhs;
}
template <typename FP, typename T>
inline bool operator<=(const UncertainValue<FP>& lhs, T rhs)
{
    return lhs.value() <= rhs;
}
template <typename FP, typename T>
inline bool operator==(const UncertainValue<FP>& lhs, T rhs)
{
    return lhs.value() == rhs;
}
template <typename FP, typename T>
inline bool operator!=(const UncertainValue<FP>& lhs, T rhs)
{
    return lhs.value() != rhs;
}

// Comparison with scalar on left
template <typename FP, typename T>
inline bool operator>(T lhs, const UncertainValue<FP>& rhs)
{
    return lhs > rhs.value();
}
template <typename FP, typename T>
inline bool operator<(T lhs, const UncertainValue<FP>& rhs)
{
    return lhs < rhs.value();
}
template <typename FP, typename T>
inline bool operator>=(T lhs, const UncertainValue<FP>& rhs)
{
    return lhs >= rhs.value();
}
template <typename FP, typename T>
inline bool operator<=(T lhs, const UncertainValue<FP>& rhs)
{
    return lhs <= rhs.value();
}
template <typename FP, typename T>
inline bool operator==(T lhs, const UncertainValue<FP>& rhs)
{
    return lhs == rhs.value();
}
template <typename FP, typename T>
inline bool operator!=(T lhs, const UncertainValue<FP>& rhs)
{
    return lhs != rhs.value();
}

// Comparison between two UncertainValue
template <typename FP>
inline bool operator>(const UncertainValue<FP>& lhs, const UncertainValue<FP>& rhs)
{
    return lhs.value() > rhs.value();
}
template <typename FP>
inline bool operator<(const UncertainValue<FP>& lhs, const UncertainValue<FP>& rhs)
{
    return lhs.value() < rhs.value();
}
template <typename FP>
inline bool operator>=(const UncertainValue<FP>& lhs, const UncertainValue<FP>& rhs)
{
    return lhs.value() >= rhs.value();
}
template <typename FP>
inline bool operator<=(const UncertainValue<FP>& lhs, const UncertainValue<FP>& rhs)
{
    return lhs.value() <= rhs.value();
}
template <typename FP>
inline bool operator==(const UncertainValue<FP>& lhs, const UncertainValue<FP>& rhs)
{
    return lhs.value() == rhs.value();
}
template <typename FP>
inline bool operator!=(const UncertainValue<FP>& lhs, const UncertainValue<FP>& rhs)
{
    return lhs.value() != rhs.value();
}

// Arithmetic between two UncertainValue (returns FP)
template <typename FP>
inline FP operator*(const UncertainValue<FP>& lhs, const UncertainValue<FP>& rhs)
{
    return lhs.value() * rhs.value();
}
template <typename FP>
inline FP operator/(const UncertainValue<FP>& lhs, const UncertainValue<FP>& rhs)
{
    return lhs.value() / rhs.value();
}
template <typename FP>
inline FP operator+(const UncertainValue<FP>& lhs, const UncertainValue<FP>& rhs)
{
    return lhs.value() + rhs.value();
}
template <typename FP>
inline FP operator-(const UncertainValue<FP>& lhs, const UncertainValue<FP>& rhs)
{
    return lhs.value() - rhs.value();
}

// Arithmetic with scalar on right
template <typename FP, typename T>
inline FP operator*(const UncertainValue<FP>& lhs, const T& rhs)
{
    return lhs.value() * rhs;
}
template <typename FP, typename T>
inline FP operator/(const UncertainValue<FP>& lhs, const T& rhs)
{
    return lhs.value() / rhs;
}
template <typename FP, typename T>
inline FP operator+(const UncertainValue<FP>& lhs, const T& rhs)
{
    return lhs.value() + rhs;
}
template <typename FP, typename T>
inline FP operator-(const UncertainValue<FP>& lhs, const T& rhs)
{
    return lhs.value() - rhs;
}

// Arithmetic with scalar on left
template <typename FP, typename T>
inline FP operator*(const T& lhs, const UncertainValue<FP>& rhs)
{
    return lhs * rhs.value();
}
template <typename FP, typename T>
inline FP operator/(const T& lhs, const UncertainValue<FP>& rhs)
{
    return lhs / rhs.value();
}
template <typename FP, typename T>
inline FP operator+(const T& lhs, const UncertainValue<FP>& rhs)
{
    return lhs + rhs.value();
}
template <typename FP, typename T>
inline FP operator-(const T& lhs, const UncertainValue<FP>& rhs)
{
    return lhs - rhs.value();
}

// Compound assignment operators for FP and UncertainValue<FP>
template <typename FP>
FP& operator+=(FP& lhs, const UncertainValue<FP>& rhs)
{
    lhs += rhs.value();
    return lhs;
}

template <typename FP>
FP& operator-=(FP& lhs, const UncertainValue<FP>& rhs)
{
    lhs -= rhs.value();
    return lhs;
}

template <typename FP>
FP& operator*=(FP& lhs, const UncertainValue<FP>& rhs)
{
    lhs *= rhs.value();
    return lhs;
}

template <typename FP>
FP& operator/=(FP& lhs, const UncertainValue<FP>& rhs)
{
    lhs /= rhs.value();
    return lhs;
}

/**
 * @brief Stream insertion prints only the scalar value.
 */
template <typename FP>
inline std::ostream& operator<<(std::ostream& os, const UncertainValue<FP>& uv)
{
    return os << uv.value();
}

} // namespace mio

#endif // UNCERTAINVALUE_H
