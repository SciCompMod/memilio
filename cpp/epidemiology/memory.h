#ifndef MEMORY_H
#define MEMORY_H

#include <utility>

namespace epi
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

private:
    T* m_raw_ptr;
};

template <typename T>
observer_ptr<T> make_observer(T* p)
{
    return observer_ptr<T>(p);
}

} // namespace epi

#endif // MEMORY_H
