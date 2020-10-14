#ifndef EPI_SPAN_H
#define EPI_SPAN_H

#include <utility>
#include <cstddef>

namespace epi
{

/**
 * a reference to any contigiuous array of objects.
 */
template <class T>
class Span
{
public:
    /**
     * empty span.
     */
    Span() = default;

    /**
     * construct from a container.
     * e.g. std::vector or std::array
     */
    template <class Cont>
    Span(const Cont& c)
        : m_ptr(c.size() == 0 ? nullptr : std::addressof(*c.begin()))
        , m_size(c.size())
    {
    }

    /**
     * construct from a c array.
     */
    template <size_t N>
    Span(T (&p)[N])
        : m_ptr(p)
        , m_size(N)
    {
    }

    /**
     * construct from ptr and size.
     */
    Span(const T* ptr, size_t size)
        : m_ptr(ptr)
        , m_size(size)
    {
    }

    /**
     * get the adress of the beginning of array.
     */
    const T* get_ptr() const
    {
        return m_ptr;
    }

    /**
     * get an iterator to the first element.
     */
    const T* begin() const
    {
        return m_ptr;
    }

    /**
     * get an iterator to one past the last element.
     */
    const T* end() const
    {
        return m_ptr + m_size;
    }

    /**
     * get the number of elements in the array
     */
    size_t size() const
    {
        return m_size;
    }

private:
    const T* m_ptr = nullptr;
    size_t m_size  = 0;
};

} // namespace epi

#endif //EPI_SPAN_H