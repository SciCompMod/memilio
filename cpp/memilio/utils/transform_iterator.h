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

#ifndef MEMILIO_UTILS_TRANSFORM_ITERATOR_H
#define MEMILIO_UTILS_TRANSFORM_ITERATOR_H

#include <utility>
#include <type_traits>

namespace mio
{

/**
 * Iterator that transforms the values produced by an underlying iterator.
 */
template <class UnderlyingIter, class Transform>
class TransformIterator
{
private:
    UnderlyingIter m_underlying_iter;
    Transform m_transform;

public:
    using reference  = decltype(std::declval<Transform>()(std::declval<typename UnderlyingIter::reference>()));
    using value_type = typename std::remove_reference<reference>::type;
    struct ProxyPointer {
        ProxyPointer(value_type&& v)
            : m_value(std::move(v))
        {
        }
        value_type m_value;
        value_type* operator->() const
        {
            return &m_value;
        }
    };
    using pointer           = std::conditional_t<std::is_lvalue_reference<reference>::value, value_type*, ProxyPointer>;
    using difference_type   = typename UnderlyingIter::difference_type;
    using iterator_category = typename UnderlyingIter::iterator_category;

    /**
     * Create a transforming iterator.
     * Each instance of the class owns a copy of the transform functor!
     * @param underlying_iter The underlying iterator.
     * @param transform Instance of the transform functor. Default value is default constructed if possible.
     * @see make_transform_iterator for a factory function with template argument deduction.
     * @{
     */
    template <class DTransform                                                          = Transform,
              std::enable_if_t<std::is_default_constructible<DTransform>::value, void*> = nullptr>
    explicit TransformIterator(UnderlyingIter underlying_iter, Transform transform = Transform{})
        : m_underlying_iter(underlying_iter)
        , m_transform(std::move(transform))
    {
    }
    template <class DTransform                                                           = Transform,
              std::enable_if_t<!std::is_default_constructible<DTransform>::value, void*> = nullptr>
    TransformIterator(UnderlyingIter underlying_iter, Transform transform)
        : m_underlying_iter(underlying_iter)
        , m_transform(std::move(transform))
    {
    }
    /**@}*/

    /**
     * Dereference the underlying iterator and transform the value.
     * Return value is a temporary value if the transform returns a temporary or rvalue reference.
     * Otherwise the return value is an lvalue reference.
     */
    reference operator*() const
    {
        return m_transform(*m_underlying_iter);
    }
    /**
     * Dereference the underlying iterator with offset and transform the value.
     * @return Result of the transformation. Lvalue reference if transform returns an lvalue, Rvalue otherwise.
     */
    reference operator[](difference_type d) const
    {
        return m_transform(m_underlying_iter[d]);
    }
    /**
     * Dereference the underlying iterator and transform the value.
     * @return Address of the result of the transformation if the transform returns an lvalue. Otherwise returns a proxy that stores the
     * transformation result internally.
     * @{
     */
    template <class Dummy = reference, std::enable_if_t<std::is_lvalue_reference<Dummy>::value, void*> = nullptr>
    pointer operator->() const
    {
        return &m_transform(*m_underlying_iter);
    }
    template <class Dummy = reference, std::enable_if_t<!std::is_lvalue_reference<Dummy>::value, void*> = nullptr>
    pointer operator->() const
    {
        return ProxyPointer{m_transform(*m_underlying_iter)};
    }
    /**@}*/

    /**
     * increment/decrement the underlying operator.
     */
    TransformIterator& operator++()
    {
        ++m_underlying_iter;
        return *this;
    }
    TransformIterator operator++(int)
    {
        return m_underlying_iter++;
    }
    TransformIterator& operator+=(difference_type d)
    {
        m_underlying_iter += d;
        return *this;
    }
    TransformIterator operator+(difference_type d) const
    {
        return TransformIterator(m_underlying_iter + d, m_transform);
    }
    friend TransformIterator operator+(difference_type d, const TransformIterator& i)
    {
        return i + d;
    }
    TransformIterator& operator--()
    {
        --m_underlying_iter;
        return *this;
    }
    TransformIterator operator--(int)
    {
        return m_underlying_iter--;
    }
    TransformIterator& operator-=(difference_type d)
    {
        m_underlying_iter -= d;
        return *this;
    }
    TransformIterator operator-(difference_type d) const
    {
        return {m_underlying_iter - d, m_transform};
    }
    difference_type operator-(TransformIterator d) const
    {
        return m_underlying_iter - d.m_underlying_iter;
    }
    /*@}*/

    /**
     * compare the underlying operator.
     * @{
     */
    bool operator<(const TransformIterator& other) const
    {
        return m_underlying_iter < other.m_underlying_iter;
    }
    bool operator<=(const TransformIterator& other) const
    {
        return m_underlying_iter <= other.m_underlying_iter;
    }
    bool operator>(const TransformIterator& other) const
    {
        return m_underlying_iter > other.m_underlying_iter;
    }
    bool operator>=(const TransformIterator& other) const
    {
        return m_underlying_iter >= other.m_underlying_iter;
    }
    bool operator!=(const TransformIterator& other) const
    {
        return m_underlying_iter != other.m_underlying_iter;
    }
    bool operator==(const TransformIterator& other) const
    {
        return m_underlying_iter == other.m_underlying_iter;
    }
    /**@}*/
};

/**
 * Create a TransformIterator.
 * @param iter The underlying iterator.
 * @param transform Instance of the transform functor. Default value is default constructed if possible.
 * @{
 */
template <class Iter, class Transform>
auto make_transform_iterator(Iter&& iter, Transform&& transform)
{
    return TransformIterator<std::decay_t<Iter>, std::decay_t<Transform>>(std::forward<Iter>(iter),
                                                                          std::forward<Transform>(transform));
}

} // namespace mio

#endif
