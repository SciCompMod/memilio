/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele
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
#ifndef EPI_UTILS_TYPE_SAFE_H
#define EPI_UTILS_TYPE_SAFE_H

#include "memilio/io/io.h"
#include <ostream>

namespace mio
{

/**
 * typesafe wrapper around any type to make function arguments, tuple elements, etc. easily distinguishable.
 * e.g.
 *   class Foo : public Type<int, Foo>() {};
 *   class Bar : public Type<int, Bar>() {};
 *   void work(Foo f, Bar b);
 * instead of
 *   void work(int f, int b);
 * @tparam T type of the underlying value
 * @tparam Derived Concrete type derived from this type for CRTP.
 */
template <class T, class Derived>
class TypeSafe
{
public:
    using ValueType = T;

    /**
     * default constructor.
     */
    TypeSafe()
    {
    }

    /**
     * value constructor.
     */
    explicit TypeSafe(T t)
        : m_t(t)
    {
    }

    /**
     * conversion to underlying type.
     */
    explicit operator T() const
    {
        return m_t;
    }
    T get() const
    {
        return m_t;
    }

    /**
     * equality operators.
     */
    bool operator==(const Derived& other) const
    {
        return m_t == other.m_t;
    }
    bool operator!=(const Derived& other) const
    {
        return !(*this == other);
    }

    /**
     * stream operators.
     */
    friend std::ostream& operator<<(std::ostream& os, const Derived& ts)
    {
        return (os << ts.m_t);
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        mio::serialize(io, T(*this));
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Derived> deserialize(IOContext& io)
    {
        BOOST_OUTCOME_TRY(auto&& t, mio::deserialize(io, Tag<T>{}));
        return success(Derived(t));
    }

private:
    T m_t;
};

/**
 * base class to add operator ++, -- (pre- and post-) to a class derived from TypeSafe.
 * @tparam TS Concrete class derived from TypeSafe
 */
template <class TS>
class OperatorIncrementDecrement
{
public:
    TS& operator++()
    {
        return static_cast<TS&>(*this) = TS{static_cast<const TS&>(*this).get() + 1};
    }
    TS operator++(int)
    {
        auto tmp                = static_cast<TS&>(*this);
        static_cast<TS&>(*this) = TS{static_cast<const TS&>(*this).get() + 1};
        return tmp;
    }
    TS& operator--()
    {
        return static_cast<TS&>(*this) = TS{static_cast<const TS&>(*this).get() - 1};
    }
    TS operator--(int)
    {
        auto tmp                = static_cast<TS&>(*this);
        static_cast<TS&>(*this) = TS{static_cast<const TS&>(*this).get() - 1};
        return tmp;
    }
};

/**
 * base class to add default operator +, +=, -, -= to a class derived from TypeSafe.
 * @tparam TS Concrete class derived from TypeSafe
 */
template <class TS>
class OperatorAdditionSubtraction : public OperatorIncrementDecrement<TS>
{
public:
    TS operator+(const TS& other) const
    {
        return TS{static_cast<const TS&>(*this).get() + other.get()};
    }
    TS& operator+=(const TS& other)
    {
        return static_cast<TS&>(*this) = (*this + other);
    }
    TS operator-(const TS& other) const
    {
        return TS{static_cast<const TS&>(*this).get() - other.get()};
    }
    TS& operator-=(const TS& other)
    {
        return static_cast<TS&>(*this) = (*this - other);
    }
};

/**
 * base class to add operator *, *=, /, /= with a scalar to a class derived from TypeSafe.
 * @tparam TS Concrete class derived from TypeSafe
 * @tparam S Type of the scalar, default underlying type of TS
 */
template <class TS, class S> //S can't default to TS::ValueType because TS is imcomplete
class OperatorScalarMultiplicationDivision
{
public:
    TS operator*(const S& other) const
    {
        return TS{static_cast<const TS&>(*this).get() * other};
    }
    TS& operator*=(const S& other)
    {
        return static_cast<TS&>(*this) = (*this * other);
    }
    TS operator/(const S& other) const
    {
        return TS{static_cast<const TS&>(*this).get() / other};
    }
    TS& operator/=(const S& other)
    {
        return static_cast<TS&>(*this) = (*this / other);
    }
};

/**
 * base class to add operator <, <=, >, >= to a class derived from TypeSafe.
 * @tparam TS Concrete class derived from TypeSafe
 */
template <class TS>
class OperatorComparison
{
public:
    bool operator<(const TS& other) const
    {
        return static_cast<const TS&>(*this).get() < static_cast<const TS&>(other).get();
    }
    bool operator<=(const TS& other) const
    {
        return static_cast<const TS&>(*this).get() <= static_cast<const TS&>(other).get();
    }
    bool operator>(const TS& other) const
    {
        return static_cast<const TS&>(*this).get() > static_cast<const TS&>(other).get();
    }
    bool operator>=(const TS& other) const
    {
        return static_cast<const TS&>(*this).get() >= static_cast<const TS&>(other).get();
    }
};

/**
 * helper macro to declare a class that derives from TypeSafe.
 */
#define DECL_TYPESAFE(T, Name)                                                                                         \
    struct Name : public ::mio::TypeSafe<T, Name> {                                                                    \
        using TypeSafe<T, Name>::TypeSafe;                                                                             \
    }

} // namespace mio

#endif //EPI_UTILS_TYPE_SAFE_H
