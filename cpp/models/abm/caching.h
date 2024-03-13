#ifndef CACHING_H_
#define CACHING_H_

namespace mio
{

template <class T>
struct Cache {
    // const access to the cached data
    const T& read() const
    {
        return data;
    }

    // access to cached data - only intended for writing, use read() to access data
    T& write()
    {
        return data;
    }

    bool is_valid() const
    {
        return m_is_valid;
    }

    void invalidate()
    {
        m_is_valid = false;
    }

    void validate()
    {
        m_is_valid = true;
    }

private:
    T data;
    bool m_is_valid = false;
};

} // namespace mio

#endif // CACHING_H_
