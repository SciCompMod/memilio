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
#ifndef MEMORY_H
#define MEMORY_H

#include <utility>

namespace mio
{

/**
 * @brief A non-owning pointer
 *
 * This class is identical to a raw pointer. However, the sematics
 * of the name makes the ownership more obvious.
 * One difference to the raw pointer is that it can't be deleted.
 */
template <typename T>
class observer_ptr
{
public:
    observer_ptr(T* p)
        : m_raw_ptr(p)
    {
    }

    observer_ptr(const observer_ptr& other) = default;

    /**
     * @brief Replaces the watched object
     */
    void reset(T* p = nullptr)
    {
        m_raw_ptr = p;
    }

    /**
     * @brief Returns a raw pointer to the watched object
     */
    T* get() const
    {
        return m_raw_ptr;
    }

    T* operator->() const
    {
        return m_raw_ptr;
    }

    T& operator*() const
    {
        return *m_raw_ptr;
    }

    operator bool() const
    {
        return m_raw_ptr != nullptr;
    }

    bool operator==(const observer_ptr& other) const
    {
        return m_raw_ptr == other.m_raw_ptr;
    }
    bool operator!=(const observer_ptr& other) const
    {
        return !(*this == other);
    }

private:
    T* m_raw_ptr;
};

template <typename T>
observer_ptr<T> make_observer(T* p)
{
    return observer_ptr<T>(p);
}

} // namespace mio

#endif // MEMORY_H
