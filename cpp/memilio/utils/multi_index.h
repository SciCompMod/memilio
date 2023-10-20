#ifndef MIO_UTILS_MULTI_INDEX_H_
#define MIO_UTILS_MULTI_INDEX_H_

#include "memilio/utils/index.h"

namespace mio
{

template <class MultiIndex>
class MultiIndexRange
{
public:
    MultiIndexRange(MultiIndex dimensions)
        : m_dimensions(dimensions)
        , m_end(make_end_index())
    {
    }

    class MultiIndexIterator
    {
    public:
        using value_type = MultiIndex;
        using reference  = const MultiIndex&;

        explicit MultiIndexIterator(value_type index, reference dimensions)
            : m_index(index)
            , m_dims(dimensions)
        {
        }

        reference operator*() const
        {
            return m_index;
        }

        MultiIndexIterator& operator++()
        {
            if constexpr (MultiIndex::size > 0) {
                increment_index();
            }
            return *this;
        }

        MultiIndexIterator operator++(int)
        {
            auto retval = *this;
            ++(*this);
            return retval;
        }

        bool operator==(MultiIndexIterator& other) const
        {
            return m_index == other.m_index;
        }

        bool operator!=(MultiIndexIterator& other) const
        {
            return !(*this == other);
        }

    private:
        template <size_t I = MultiIndex::size - 1>
        inline void increment_index()
        {
            // TODO: should we allow dims[I] == 0? It makes sense not to (if any dims[.] is 0, the total "volume" of dims is 0, too), but it is easy to deal with (treat 0 as 1)
            assert(mio::get<I>(m_dims).get() > 0 && "Dimensions must be positive.");
            assert(mio::get<I>(m_index) < mio::get<I>(m_dims) && "Index out of bounds.");
            if constexpr (I > 0) {
                ++mio::get<I>(m_index);
                // usually, "==" would be enough. but we use ">=" in case m_dims[I] is 0
                if (mio::get<I>(m_index) >= mio::get<I>(m_dims)) {
                    mio::get<I>(m_index) = mio::get<I>(m_index).Zero();
                    increment_index<I - 1>();
                }
            }
            else {
                ++mio::get<0>(m_index);
            }
        }
        value_type m_index;
        reference m_dims;
    };

    MultiIndexIterator begin()
    {
        return MultiIndexIterator(MultiIndex::Zero(), m_dimensions);
    }
    MultiIndexIterator end()
    {
        return MultiIndexIterator(m_end, m_dimensions);
    }

private:
    MultiIndex make_end_index()
    {
        auto end         = MultiIndex::Zero();
        mio::get<0>(end) = mio::get<0>(m_dimensions);
        return end;
    }
    const MultiIndex m_dimensions; ///< Upper bound for each Index in MultiIndex.
    const MultiIndex m_end; ///< Index returned by
};

// TODO: rename, remove, ...
template <class... CategoryTags>
MultiIndexRange<Index<CategoryTags...>> make_range(Index<CategoryTags...> dimensions)
{
    return MultiIndexRange<Index<CategoryTags...>>(dimensions);
}

} // namespace mio

#endif