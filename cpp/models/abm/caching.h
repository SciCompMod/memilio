#ifndef CACHING_H_
#define CACHING_H_

#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/index.h"

// TODO: this implementation is not optimal

namespace mio
{

template <class T>
struct Cache {
    const T& read() const
    {
        return data;
    }

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

template <class CustomIndexArrayOfAtomics>
void setZero(CustomIndexArrayOfAtomics& x)
{
    for (Eigen::Index i = 0; i < x.array().size(); i++) {
        x.array()[i].store(0);
    }
}

template <class Atomic>
struct CopyableAtomic {
    Atomic value;
    CopyableAtomic() = default;
    CopyableAtomic(const typename Atomic::value_type& v)
        : value(v)
    {
    }
    CopyableAtomic(const Atomic& other)
        : value(other.load())
    {
    }
    CopyableAtomic(const CopyableAtomic& other)
        : CopyableAtomic(other.value)
    {
    }
    CopyableAtomic& operator=(const CopyableAtomic& other)
    {
        value.store(other.value.load());
        return *this;
    }
};

template <class Atomic, class... Categories>
class CustomAtomicIndexArray
{
    static_assert(Atomic::is_always_lock_free, "Atomic type must be lock free");

public:
    using Index = mio::Index<Categories...>;

    CustomAtomicIndexArray(const Index& dimensions, typename Atomic::value_type fill)
        : m_dimensions(dimensions)
        , m_data(product(m_dimensions), fill)
    {
        // TODO: why is product not in mio?
    }

    Atomic& operator[](const Index& i)
    {
        return m_data[mio::flatten_index(i, m_dimensions)].value;
    }

    const Atomic& operator[](const Index& i) const
    {
        return m_data[mio::flatten_index(i, m_dimensions)].value;
    }

    typename Atomic::value_type operator()(const Index& i) const
    {
        return m_data[mio::flatten_index(i, m_dimensions)].value.load();
    }

    template <class Category>
    mio::Index<Category> size() const
    {
        return get<Category>(m_dimensions);
    }

    void setZero() // TODO: this is only well defined for FP atomics
    {
        for (auto& ca : m_data) {
            ca.value.store(0.);
        }
    }

private:
    Index m_dimensions;

    std::vector<CopyableAtomic<Atomic>> m_data; // TODO: try using eigen::matrix without copyable
};

} // namespace mio

#endif // CACHING_H_