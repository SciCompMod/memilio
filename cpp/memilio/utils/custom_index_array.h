/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#ifndef CUSTOMINDEXARRAY_H
#define CUSTOMINDEXARRAY_H

#include "memilio/config.h"
#include "memilio/math/eigen.h"
#include "memilio/math/eigen_util.h"
#include "memilio/utils/index.h"
#include "memilio/utils/stl_util.h"

#include <vector>
#include <array>
#include <numeric>

namespace
{

//calculate the product of tuple elements
// std::apply or fold expression in C++17
template <int I, template <class...> class Index, class... Ts>
typename std::enable_if<(I == sizeof...(Ts)), size_t>::type product(Index<Ts...> const&)
{
    return 1;
}

template <int I, template <class...> class Index, class... Ts>
typename std::enable_if<(I < sizeof...(Ts)), size_t>::type product(Index<Ts...> const& t)
{
    return (size_t)mio::get<I>(t) * product<I + 1, Index, Ts...>(t);
}

template <template <class...> class Index, class... Ts>
size_t product(Index<Ts...> const& t)
{
    return product<0, Index, Ts...>(t);
}

} // namespace

namespace mio
{

namespace details
{

// calculate the Position of an element in a MultiIndex, given its type
template <class T, class Tuple>
struct IndexPosition;

template <class T, class... Types>
struct IndexPosition<T, Index<T, Types...>> {
    static const std::size_t value = 0;
};

template <class T, class U, class... Types>
struct IndexPosition<T, Index<U, Types...>> {
    static const std::size_t value = 1 + IndexPosition<T, Index<Types...>>::value;
};

// Internal implementation for flatten_index
template <size_t I, typename Index>
std::enable_if_t<(I == (Index::size - 1)), std::pair<size_t, size_t>> flatten_index(Index const& indices,
                                                                                    Index const& dimensions)
{
    assert(get<I>(indices) < get<I>(dimensions));
    return {(size_t)mio::get<I>(indices), (size_t)mio::get<I>(dimensions)};
}

template <size_t I, typename Index>
std::enable_if_t<(I < (Index::size - 1)), std::pair<size_t, size_t>> flatten_index(Index const& indices,
                                                                                   Index const& dimensions)
{
    assert(mio::get<I>(indices) < mio::get<I>(dimensions));

    size_t val, prod;
    std::tie(val, prod) = flatten_index<I + 1>(indices, dimensions);

    return {val + (size_t)mio::get<I>(indices) * prod, prod * (size_t)mio::get<I>(dimensions)};
}

template <typename T>
struct is_random_access_iterator
    : std::is_base_of<typename std::iterator_traits<T>::iterator_category, std::random_access_iterator_tag> {
};

} // namespace details

/**
 * @brief flatten_index takes a set of indices into a mutlidemsional array and calculates the flat index
 *
 * Given indices (i,j,k,...) of a tensor with dimensions (n,m,l,...), flatten_index calculates
 * the index of the corresponding element if the elements are sorted sequentially in a
 * row major fashion (that is right indices are incremented before left indices)
 *
 * @param indices a vector of indices of a hypothetical tensor
 * @param dimensions a vector of the dimension sizes of each dimension
 * @return the corresponding flat index
 */
template <typename MultiIndex>
size_t flatten_index(MultiIndex const& indices, MultiIndex const& dimensions)
{
    return details::flatten_index<0>(indices, dimensions).first;
}

/**
 * @brief A class template for an array with custom indices
 *
 * This class stores an array of elements that can be queried using
 * a MultiIndex. Each index in the MultiIndex is associated
 * with a category, or dimension into a multidimensional array.
 *
 *
 * Example:
 *
 * struct AgeGroup{};
 *
 * enum class Gender
 * {
 *    Female,
 *    Male,
 *    Diverse,
 *    Count, //number of values in the enum to allow iterating over the enum
 * };
 * struct Gender{};
 *
 * CustomIndexArray<size_t, AgeGroup, Gender> populations({Index<AgeGroup>(2), Index<Gender>(Gender::Count)});
 *
 * Here, populations represents a 2x3 size_t array (though the data is stored contigously).
 * An element can be accessed using a MultiIndex:
 *
 * auto x = populations[{Index<AgeGroup>(0), Index<Gender>(Gender::Female)}];
 * 
 * If the category is an enum class, the Index class can be omitted for convenience in almost all cases since
 * Index<E> is implicitly constructible from enum type E, e.g.:
 * 
 * CustomIndexArray<size_t, AgeGroup, Gender> populations({Index<AgeGroup>(2), Gender::Count});
 * auto x = populations[{Index<AgeGroup>(0), Gender::Female}];
 *
 * @tparam Typ the type stored in the array
 * @tparam Tags Types that tag the custom Index types, must be unique.
 */
template <class Typ, class... Tags>
class CustomIndexArray
{
public:
    using Type              = Typ;
    using Index             = ::mio::Index<Tags...>;
    using InternalArrayType = Eigen::Array<Type, Eigen::Dynamic, 1>;

    /**
     * @brief CustomIndexArray constructor, that initializes the array
     * to constant instances of `CustsomIndexArray::Type`.
     *
     * It forwards the arguments for initializer_list construction of the
     * contained objects.
     *
     * @tparam Ts The argument types of the constructor arguments of Type
     * @param args The constructor arguments of Type
     */
    template <class... Ts, typename std::enable_if_t<std::is_constructible<Type, Ts...>::value>* = nullptr>
    CustomIndexArray(Index const& dims, Ts&&... args)
        : m_dimensions{dims}
        , m_numel(product(dims))
        , m_y(InternalArrayType::Constant(m_numel, 1, Type{std::forward<Ts>(args)...}))
    {
    }

    /**
     * Initializes array with values from a range.
     * Each element of the array will be assigned the corresponding value from the range.
     * @tparam Iter Iterator class.
     * @param dims dimensions of the array.
     * @param b begin of the range of values.
     * @param e end of the range of values.
     */
    template <class Iter, typename std::enable_if_t<is_iterator<Iter>::value>* = nullptr>
    CustomIndexArray(Index const& dims, Iter b, Iter e)
        : m_dimensions(dims)
        , m_numel(product(dims))
        , m_y(m_numel, 1)
    {
        std::copy(b, e, begin());
    }

    /**
     * @brief numel returns the number of elements
     *
     * This corresponds to the product of the dimension sizes
     *
     * @return number of elements
     */
    size_t constexpr numel() const
    {
        return m_numel;
    }

    /**
     * @brief returns the size along the dimension provided as template parameter
     * @tparam Tag Tag of the queried dimension/category
     * @return size along a specified dimension
     */
    template <typename Tag>
    mio::Index<Tag> size() const
    {
        return get<Tag>(m_dimensions);
    }

    /**
     * @brief returns the size of the array along all dimensions.
     * @return multi-index with size of the array along all dimensions.
     */
    Index size() const
    {
        return m_dimensions;
    }

    /**
     * Resize all dimensions.
     * @param new_dims new dimensions.
     */
    void resize(Index new_dims)
    {
        m_dimensions = new_dims;
        m_numel      = product(m_dimensions);
        m_y.conservativeResize(m_numel);
    }

    /**
     * Resize a single dimension.
     * @param new dimension.
     */
    template <class Tag>
    void resize(mio::Index<Tag> new_dim)
    {
        std::get<mio::Index<Tag>>(m_dimensions.indices) = new_dim;
        m_numel                                         = product(m_dimensions);
        m_y.conservativeResize(m_numel);
    }

    /**
     * @brief array returns a reference to the internally stored flat array.
     * @return const reference to the CustomIndexArray::InternalArrayType instance
     */
    auto const& array() const
    {
        return m_y;
    }
    auto& array()
    {
        return m_y;
    }

    /**
     * @brief returns the entry of the array given a MultiIndex
     * @param MultiIndex
     * @return the value at the index
     */
    Type& operator[](Index const& index)
    {
        return m_y[get_flat_index(index)];
    }

    /**
     * @brief returns the entry of the array given a flat index index
     * @param index a flat index
     * @return the value at the index
     */
    Type const& operator[](Index const& index) const
    {
        return m_y[get_flat_index(index)];
    }

    /**
     * Assign the same value to each element of the array.
     * @param scalar scalar value.
     */
    CustomIndexArray& operator=(const Type& scalar)
    {
        m_y = scalar;
        return *this;
    }

    /**
     * Equality comparison.
     * @{
     */
    bool operator==(const CustomIndexArray& other) const
    {
        return this->m_dimensions == other.m_dimensions && (this->m_y.array() == other.m_y.array()).all();
    }
    bool operator!=(const CustomIndexArray& other) const
    {
        return !(*this == other);
    }
    /**@}*/

    /**
     * @brief get_flat_index returns the flat index into the stored array, given the
     * indices of each category
     * @param indices the custom indices for each category
     * @return a flat index into the data structure storing the compartment populations
     */
    size_t get_flat_index(Index const& index) const
    {
        return (Eigen::Index)flatten_index(index, m_dimensions);
    }

private:
    // Random Access Iterator for CustomIndexArray
    // To Do: As of Eigen 3.4, this is not needed anymore,
    // because the Eigen Matrix classes directly implement stl-compatible
    // iterators
    template <typename T>
    struct Iterator {
        using iterator_category = std::random_access_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = value_type*;
        using reference         = value_type&;

        Iterator(pointer ptr)
            : m_ptr(ptr)
        {
        }

        Iterator& operator=(pointer rhs)
        {
            m_ptr = rhs;
            return *this;
        }
        Iterator& operator=(const Iterator& rhs)
        {
            m_ptr = rhs.m_ptr;
            return *this;
        }
        Iterator& operator+=(const int& rhs)
        {
            m_ptr += rhs;
            return *this;
        }
        Iterator& operator-=(const int& rhs)
        {
            m_ptr -= rhs;
            return *this;
        }

        reference operator[](const difference_type& rhs)
        {
            return m_ptr[rhs];
        }
        value_type const& operator[](const difference_type& rhs) const
        {
            return m_ptr[rhs];
        }
        reference operator*()
        {
            return *(m_ptr);
        }
        value_type const& operator*() const
        {
            return *(m_ptr);
        }
        pointer operator->()
        {
            return m_ptr;
        }

        Iterator& operator++()
        {
            ++m_ptr;
            return *this;
        }
        Iterator operator++(int)
        {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }
        Iterator& operator--()
        {
            --m_ptr;
            return *this;
        }
        Iterator operator--(int)
        {
            Iterator tmp = *this;
            --(*this);
            return tmp;
        }

        Iterator operator+(const int& rhs) const
        {
            return Iterator(m_ptr + rhs);
        }
        Iterator operator-(const int& rhs) const
        {
            return Iterator(m_ptr - rhs);
        }
        difference_type operator-(const Iterator& rhs) const
        {
            return m_ptr - rhs.m_ptr;
        }

        friend bool operator==(const Iterator& a, const Iterator& b)
        {
            return a.m_ptr == b.m_ptr;
        }
        friend bool operator!=(const Iterator& a, const Iterator& b)
        {
            return !(a == b);
        }
        friend bool operator<(const Iterator& a, const Iterator& b)
        {
            return a.m_ptr < b.m_ptr;
        }
        friend bool operator<=(const Iterator& a, const Iterator& b)
        {
            return a.m_ptr <= b.m_ptr;
        }
        friend bool operator>(const Iterator& a, const Iterator& b)
        {
            return a.m_ptr > b.m_ptr;
        }
        friend bool operator>=(const Iterator& a, const Iterator& b)
        {
            return a.m_ptr >= b.m_ptr;
        }

    private:
        pointer m_ptr;
    };

    /**
     * @brief A Slice represents a slice of data along one dimension, given a start and
     * end index into that dimension. Its sole use is to provide an iterator for the data
     * along this slice.
     *
     * If dims = (d0,d1,...d_i,...,dN) are the dimension sizes of this array, and d_i is the
     * dimension size of the current slice, the indices into the data array are given by
     *
     *     i*d_i + j*d_i*d_right + k
     *
     * where
     *
     *       0 <= i <  d_right
     *   start <= j <= end
     *       0 <= k <  d_left
     *
     * and start and end are the start and end indices for this slice, d_left = d_0*d_1*...*d_(i-1)
     * and d_right = d_(i+1)*...*d_N. If d_i = d_0, d_left=1 and if d_i = d_N, d_right=1.
     *
     * A slice iterator will store the indices (i,j,k) and increment k before j before i, so that the
     * returned values are traversed in the flat index from front to back.
     *
     */
    template <typename Tag, typename iter_type,
              typename std::enable_if_t<details::is_random_access_iterator<iter_type>::value>* = nullptr>
    class Slice
    {
        using difference_type = typename iter_type::difference_type;

        template <typename T>
        class Iterator
        {
        public:
            using iterator_category = std::random_access_iterator_tag;
            using difference_type   = std::ptrdiff_t;
            using value_type        = T;
            using pointer           = value_type*;
            using reference         = value_type&;

            Iterator(iter_type begin_, size_t di_, size_t dr_, Seq<size_t> const& seq_, difference_type offset = 0)
                : data_begin(begin_)
                , di(di_)
                , dr(dr_)
                , seq(seq_)
                , inner_offset(offset)
            {
            }

            Iterator& operator=(size_t rhs)
            {
                inner_offset = rhs;
                return *this;
            }
            Iterator& operator=(const Iterator& rhs)
            {
                inner_offset = rhs.inner_offset;
                return *this;
            }
            Iterator& operator+=(const int& rhs)
            {
                inner_offset += rhs;
                return *this;
            }
            Iterator& operator-=(const int& rhs)
            {
                inner_offset -= rhs;
                return *this;
            }

            reference operator[](const difference_type& rhs)
            {
                return data_begin[outer_offset(inner_offset + rhs)];
            }
            value_type const& operator[](const difference_type& rhs) const
            {
                return data_begin[outer_offset(inner_offset + rhs)];
            }
            reference operator*()
            {
                return data_begin[outer_offset(inner_offset)];
            }
            value_type const& operator*() const
            {
                return *(data_begin[outer_offset(inner_offset)]);
            }
            pointer operator->()
            {
                return data_begin + outer_offset(inner_offset);
            }

            Iterator& operator++()
            {
                inner_offset++;
                return *this;
            }
            Iterator operator++(int)
            {
                Iterator tmp = *this;
                ++(*this);
                return tmp;
            }
            Iterator& operator--()
            {
                --inner_offset;
                return *this;
            }
            Iterator operator--(int)
            {
                Iterator tmp = *this;
                --(*this);
                return tmp;
            }

            Iterator operator+(const int& rhs) const
            {
                return Iterator(data_begin, di, dr, seq, inner_offset + rhs);
            }
            Iterator operator-(const int& rhs) const
            {
                return Iterator(data_begin, di, dr, seq, inner_offset - rhs);
            }

            friend bool operator==(const Iterator& a, const Iterator& b)
            {
                return a.inner_offset == b.inner_offset && a.data_begin == b.data_begin && a.di == b.di &&
                       a.dr == b.dr && a.seq.start == b.seq.start && a.seq.n == b.seq.n && a.seq.stride == b.seq.stride;
            }
            friend bool operator!=(const Iterator& a, const Iterator& b)
            {
                return !(a == b);
            }
            friend bool operator<(const Iterator& a, const Iterator& b)
            {
                return a.inner_offset < b.inner_offset;
            }
            friend bool operator<=(const Iterator& a, const Iterator& b)
            {
                return a.inner_offset <= b.inner_offset;
            }
            friend bool operator>(const Iterator& a, const Iterator& b)
            {
                return a.inner_offset > b.inner_offset;
            }
            friend bool operator>=(const Iterator& a, const Iterator& b)
            {
                return a.inner_offset >= b.inner_offset;
            }

        private:
            inline Slice::difference_type outer_offset(difference_type const& inner) const
            {

                // calculate the outer offset from the inner offset

                // first unravel the inner index into an index (i,j,k) for a 3-dim array with dims (dl, idx_sequence.n, dr)
                auto dv           = std::div(inner, seq.n * dr);
                difference_type i = dv.quot;
                dv                = std::div(dv.rem, dr);
                difference_type j = dv.quot * seq.stride + seq.start;
                difference_type k = dv.rem;

                // then flatten the index for a 3-dim array with dims (dl, di, dr)
                return i * di * dr + j * dr + k;
            }

            iter_type data_begin;
            size_t di, dr;
            Seq<size_t> seq;
            difference_type inner_offset;
        };

    public:
        using value_type     = Type;
        using reference      = Type&;
        using iterator       = Iterator<Type>;
        using const_iterator = Iterator<Type const>;

        /**
         * @brief Constructs a slice into the CustomIndexarray
         * @param dimensions the dimensions of the CustomIndexArray
         * @param start_iter An iterator to the first element of the data
         * @param idx_sequence_ A sequence of indices into the slice
         */
        Slice(Index const& dimensions, iter_type const& start_iter, Seq<size_t> idx_sequence_)
            : data_begin(start_iter)
            , idx_sequence(idx_sequence_)
            , m_dimensions(dimensions)
            , di(mio::get<Tag>(dimensions))
            , dr(product<details::IndexPosition<Tag, Index>::value>(dimensions) / di)
            , dl(product(dimensions) / (di * dr))
        {
            assert((size_t)idx_sequence.start + idx_sequence.n <= di);

            mio::get<Tag>(m_dimensions) = mio::Index<Tag>(idx_sequence.n);
        }

        // returns the number of elements in a slice
        size_t numel() const
        {
            return dl * dr * (idx_sequence.n);
        }

        // returns an stl-compatible random access iterator into the slice
        iterator begin()
        {
            return iterator(data_begin, di, dr, idx_sequence, 0);
        }

        // returns an stl-compatible random access iterator into the slice
        const_iterator begin() const
        {
            return const_iterator(data_begin, di, dr, idx_sequence, 0);
        }

        // returns an stl-compatible end random access iterator into the slice
        iterator end()
        {
            return iterator(data_begin, di, dr, idx_sequence, numel());
        }

        // returns an stl-compatible end random access iterator into the slice
        const_iterator end() const
        {
            return const_iterator(data_begin, di, dr, idx_sequence, numel());
        }

        // copies the slice elements into a CustomIndexArray of appropriate dimension
        CustomIndexArray<Type, Tags...> as_array()
        {
            CustomIndexArray<Type, Tags...> array(m_dimensions);
            Eigen::Index idx = 0;
            for (auto& v : *this) {
                array.array()[idx++] = v;
            };
            return array;
        }

        // comparison operators
        friend bool operator==(const Slice& a, const Slice& b)
        {
            return a.data_begin == b.data_begin && a.idx_sequence.start == b.idx_sequence.start &&
                   a.idx_sequence.n == b.idx_sequence.n && a.idx_sequence.stride == b.idx_sequence.stride &&
                   a.m_dimensions == b.m_dimensions;
        }

        friend bool operator!=(const Slice& a, const Slice& b)
        {
            return !(a == b);
        }

        /**
         * @brief Creates a subslice from the current slice
         * @tparam The Tag corresponding to the dimension of the slice
         * @param idx_seq An index sequence, consisting of the first index,
         *                the number of indices and a stride
         * @return The subslice
         */
        template <typename OtherTag>
        Slice<OtherTag, iterator> slice(Seq<size_t> idx_sequence_)
        {
            return Slice<OtherTag, iterator>(m_dimensions, begin(), idx_sequence_);
        }

        /**
         * @brief Creates a subslice from the current slice
         * @tparam The Tag corresponding to the dimension of the slice
         * @param idx_seq An index sequence, consisting of the first index,
         *                the number of indices and a stride
         * @return The subslice
         */
        template <typename OtherTag>
        Slice<OtherTag, const_iterator> slice(Seq<size_t> idx_sequence_) const
        {
            return Slice<OtherTag, const_iterator>(m_dimensions, begin(), idx_sequence_);
        }

        /**
         * Assign same value to each element of the slice.
         * @param scalar Scalar value.
         */
        Slice& operator=(const Type& scalar)
        {
            for (auto&& e : *this) {
                e = scalar;
            }
            return *this;
        }

    private:
        iter_type data_begin;
        Seq<size_t> idx_sequence;
        Index m_dimensions;
        size_t di, dr, dl;
    };

    // An array storying the size of each category
    Index m_dimensions;

    // number of elements stored
    size_t m_numel;

    // An array containing the elements
    InternalArrayType m_y{};

public:
    using value_type     = Type;
    using reference      = Type&;
    using iterator       = Iterator<Type>;
    using const_iterator = Iterator<Type const>;

    /**
     * @brief Get a start iterator for the elements
     * @return random access iterator
     */
    iterator begin()
    {
        return iterator(array().data());
    }

    /**
     * @brief Get a start iterator for the elements
     * @return random access iterator
     */
    const_iterator begin() const
    {
        return const_iterator(array().data());
    }

    /**
     * @brief Get an end iterator for the elements
     * @return random access iterator
     */
    iterator end()
    {
        return iterator(array().data() + array().size());
    }

    /**
     * @brief Get an end iterator for the elements
     * @return random access iterator
     */
    const_iterator end() const
    {
        return const_iterator(array().data() + array().size());
    }

    /**
     * @brief Creates a slice into the multidimensional array.
     * Selects a sequence of indices of one dimension.
     * @tparam The Tag corresponding to the dimension of the slice
     * @param idx_seq An index sequence, consisting of the first index,
     *                the number of indices and a stride
     * @return The slice
     * @{
     */
    template <typename Tag>
    Slice<Tag, iterator> slice(Seq<size_t> idx_seq)
    {
        return Slice<Tag, iterator>(m_dimensions, begin(), idx_seq);
    }
    template <typename Tag>
    Slice<Tag, const_iterator> slice(Seq<size_t> idx_seq) const
    {
        return Slice<Tag, const_iterator>(m_dimensions, begin(), idx_seq);
    }
    /**@}*/

    /**
     * @brief Creates a slice into the multidimensional array.
     * Selects a single index of one dimension.
     * @tparam The Tag corresponding to the dimension of the slice
     * @param idx index to be selected for the slice.
     * @{
     */
    template <typename Tag>
    Slice<Tag, iterator> slice(mio::Index<Tag> idx)
    {
        return slice<Tag>({(size_t)idx, 1});
    }
    template <typename Tag>
    Slice<Tag, const_iterator> slice(mio::Index<Tag> idx) const
    {
        return slice<Tag>({(size_t)idx, 1});
    }
    template <typename Tag, std::enable_if_t<std::is_enum<Tag>::value, void*> = nullptr>
    Slice<Tag, iterator> slice(Tag idx)
    {
        return slice<Tag>(mio::Index<Tag>(idx));
    }
    template <typename Tag, std::enable_if_t<std::is_enum<Tag>::value, void*> = nullptr>
    Slice<Tag, const_iterator> slice(Tag idx) const
    {
        return slice<Tag>(mio::Index<Tag>(idx));
    }
    /**@}*/

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Array");
        obj.add_element("Dimensions", m_dimensions);
        obj.add_list("Elements", begin(), end());
    }

protected:
    /**
     * deserialize an object of a class derived from this class.
     * @see mio::deserialize
     */
    template <class IOContext, class Derived>
    static IOResult<Derived> deserialize(IOContext& io, Tag<Derived>)
    {
        auto obj  = io.expect_object("Array");
        auto dims = obj.expect_element("Dimensions", Tag<Index>{});
        auto els  = obj.expect_list("Elements", Tag<Type>{});
        return apply(
            io,
            [](auto&& d, auto&& e) -> IOResult<Derived> {
                auto a = Derived(d);
                if (a.numel() != e.size())
                    return failure(StatusCode::OutOfRange, "Dimensions of Array don't match the number of elements.");
                std::copy(e.begin(), e.end(), mio::begin(a.m_y));
                return success(a);
            },
            dims, els);
    }

public:
    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<CustomIndexArray> deserialize(IOContext& io)
    {
        return deserialize(io, Tag<CustomIndexArray>{});
    }
};

} // namespace mio

#endif // CUSTOMINDEXARRAY_H
