#ifndef CUSTOMINDEXARRAY_H
#define CUSTOMINDEXARRAY_H

#include "epidemiology/utils/eigen.h"
#include "epidemiology/utils/eigen_util.h"
#include "epidemiology/utils/type_safe.h"

#include <numeric>

namespace
{

//calculate the product of tuple elements
// std::apply or fold expression in C++17
template <int I, template <class...> class MultiIndex, class... Ts>
typename std::enable_if<I == sizeof...(Ts), size_t>::type product(MultiIndex<Ts...> const&)
{
        return 1;
}

template <int I, template <class...> class MultiIndex, class... Ts>
typename std::enable_if<I < sizeof...(Ts), size_t>::type product(MultiIndex<Ts...> const& t)
{
    return (size_t)t.template get<I>()*product<I+1, MultiIndex, Ts...>(t);
}


template <template <class...> class MultiIndex, class... Ts>
size_t product(MultiIndex<Ts...> const& t)
{
    return product<0, MultiIndex, Ts...>(t);
}

} // namespace

namespace epi
{


/**
 * @brief An Index is a typesafe wrapper for size_t that is associated with a Tag.
 *
 * The Tag can be any type such as an enum or an empty struct. The Index is used
 * in the MultiIndex class, which in turn is used to index into a CustomIndexArray:
 *
 * CustomIndexArray<Tag1, Tag2> a;
 * a[{Index<Tag1>(0), Index<Tag2>(14)}]
 *
 * Will retrieve the element associated with the indices 0 and 14 for the Tags Tag1
 * and Tag2 respectively.
 *
 * Optionally, the tag can be derived from Index to shorten the notation:
 *
 * struct Tag1 : public Index<Tag1>{
 *    Tag1(size_t val) : Index<Tag1>(val) {}
 * };
 * struct Tag2 : public Index<Tag2>{
 *    Tag2(size_t val) : Index<Tag2>(val) {}
 * };
 * CustomIndexArray<Tag1, Tag2> a;
 * a[{Tag1(0), Tag2(14)}]
 *
 * @tparam CategoryTag A tag for the typesafe index
 *
 */
template <typename CategoryTag>
class Index : public TypeSafe<size_t, Index<CategoryTag>>
            , public OperatorComparison<Index<CategoryTag>>
            , public OperatorAdditionSubtraction<Index<CategoryTag>>
{
public:
    using TypeSafe<size_t, Index<CategoryTag>>::TypeSafe;

    /**
     * @brief Constructor from enum, if CategoryTag is an enum
     */
    template <typename Dummy = CategoryTag,
                  std::enable_if_t<std::is_enum<Dummy>::value, void>* = nullptr>
    Index(Dummy val) : TypeSafe<size_t, Index<CategoryTag>>((size_t)val) {}

    /**
     * @brief Constructor from size_t
     * @param val
     */
    explicit Index(size_t val) : TypeSafe<size_t, Index<CategoryTag>>(val) {}
};


/**
 * @brief A MultiIndex combines several Index objects. It is used to index into a multidimensional
 * CustomIndexArray
 *
 * @tparam CategoryTag Variadic template parameter for the Tags used in the MultiIndex
 */
template <typename... CategoryTag>
class MultiIndex
{
public:
    static size_t constexpr size = sizeof...(CategoryTag);

    // constructor from Indices
    MultiIndex(Index<CategoryTag> const&...indices) : m_indices{indices...} {}

    // retrieves the Index at the Ith position
    template <size_t I>
    constexpr typename std::tuple_element<I, std::tuple<Index<CategoryTag>...> >::type& get() noexcept
    {
        return std::get<I>(m_indices);
    }

    // retrieves the Index at the Ith position
    template <size_t I>
    constexpr typename std::tuple_element<I, std::tuple<Index<CategoryTag>...> >::type const& get() const noexcept
    {
        return std::get<I>(m_indices);
    }

    // retrieves the Index for the tag Tag
    template <typename Tag>
    constexpr Index<Tag>& get() noexcept
    {
        return std::get<Index<Tag>>(m_indices);
    }

    // retrieves the Index for the tag Tag
    template <typename Tag>
    constexpr Index<Tag> const& get() const noexcept
    {
        return std::get<Index<Tag>>(m_indices);
    }

    // comparison operators
    bool operator==(MultiIndex const& other) const
    {
        return m_indices == other.m_indices;
    }

    bool operator!=(MultiIndex const& other) const
    {
        return !(this == other);
    }

private:
    std::tuple<Index<CategoryTag>...> m_indices;
};


namespace details {

    // calculate the Position of an element in a MultiIndex, given its type
    template <class T, class Tuple>
    struct MultiIndexPosition;

    template <class T, class... Types>
    struct MultiIndexPosition<T, MultiIndex<T, Types...>> {
        static const std::size_t value = 0;
    };

    template <class T, class U, class... Types>
    struct MultiIndexPosition<T, MultiIndex<U, Types...>> {
        static const std::size_t value = 1 + MultiIndexPosition<T, MultiIndex<Types...>>::value;
    };

    // Internal implementation for flatten_index
    template <size_t I, typename MultiIndex>
    std::enable_if_t< (I == (MultiIndex::size - 1) ), std::pair<size_t, size_t>> flatten_index(MultiIndex const& indices, MultiIndex const& dimensions)
    {
        assert(indices.template get<I>() < dimensions.template get<I>());
        return {(size_t)indices.template get<I>(), (size_t)dimensions.template get<I>()};
    }

    template <size_t I, typename MultiIndex>
    std::enable_if_t< (I < (MultiIndex::size - 1) ), std::pair<size_t, size_t>> flatten_index(MultiIndex const& indices, MultiIndex const& dimensions)
    {
        assert(indices.template get<I>() < dimensions.template get<I>());

        size_t val, prod;
        std::tie(val, prod) = flatten_index<I+1>(indices, dimensions);

        return {val + (size_t)indices.template get<I>()*prod, prod*(size_t)dimensions.template get<I>()};
    }

    template <typename T>
    struct is_random_access_iterator : std::is_base_of<
        typename std::iterator_traits<T>::iterator_category
        , std::random_access_iterator_tag>
    {};

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
 * enum GenderE
 * {
 *    Female,
 *    Male,
 *    Diverse,
 *    NumGenders = 3
 * };
 * struct Gender{};
 *
 * CustomIndexArray<size_t, AgeGroup, Gender> populations({Index<AgeGroup>(2), Index<Gender>(NumGenders)});
 *
 * Here, populations represents a 2x3 size_t array (though the data is stored contigously).
 * An element can be accessed using a MultiIndex:
 *
 * auto x = populations[{Index<AgeGroup>(0), Index<Gender>(Female)}];
 *
 * @tparam Typ the type stored in the array
 * @tparam Categories The custom Index types
 *
 */
template <class Typ, class... Tags>
class CustomIndexArray
{
public:

    using Type              = Typ;
    using MultiIndex        = ::epi::MultiIndex<Tags...>;
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
    template <class... Ts,
              typename std::enable_if_t<std::is_constructible<Type, Ts...>::value>* = nullptr>
    CustomIndexArray(MultiIndex const& dims, Ts&&... args)
        : m_dimensions{dims}
        , m_numel(product(dims))
        , m_y(InternalArrayType::Constant(m_numel, 1, {std::forward<Ts>(args)...}))
    {}

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
    Index<Tag> size() const {
        return m_dimensions.template get<Tag>();
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
    Type& operator[](MultiIndex const& index) {
        return m_y[get_flat_index(index)];
    }

    /**
     * @brief returns the entry of the array given a flat index index
     * @param index a flat index
     * @return the value at the index
     */
    Type const& operator[](MultiIndex const& index) const {
        return m_y[get_flat_index(index)];
    }


    /**
     * @brief get_flat_index returns the flat index into the stored array, given the
     * indices of each category
     * @param indices the custom indices for each category
     * @return a flat index into the data structure storing the compartment populations
     */
    size_t get_flat_index(MultiIndex const& index) const
    {
        return (Eigen::Index)flatten_index(index, m_dimensions);
    }

private:

    // Random Access Iterator for CustomIndexArray
    // To Do: As of Eigen 3.4, this is not needed anymore,
    // because the Eigen Matrix classes directly implement stl-compatible
    // iterators
    template <typename T>
    struct Iterator
    {
        using iterator_category = std::random_access_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = value_type*;
        using reference         = value_type&;

        Iterator(pointer ptr) : m_ptr(ptr) {}

        Iterator& operator=(pointer rhs) { m_ptr = rhs; return *this;}
        Iterator& operator=(const Iterator &rhs) { m_ptr = rhs._ptr; return *this;}
        Iterator& operator+=(const int& rhs) { m_ptr += rhs; return *this;}
        Iterator& operator-=(const int& rhs) { m_ptr -= rhs; return *this;}

        reference operator[](const difference_type& rhs) {return m_ptr[rhs];}
        value_type const& operator[](const difference_type& rhs) const {return m_ptr[rhs];}
        reference operator*() { return *(m_ptr); }
        value_type const& operator*() const { return *(m_ptr); }
        pointer operator->() { return m_ptr; }

        Iterator& operator++() { ++m_ptr; return *this; }
        Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
        Iterator& operator--() { --m_ptr; return *this; }
        Iterator operator--(int) { Iterator tmp = *this; --(*this); return tmp; }

        Iterator operator+(const int& rhs) const { return Iterator(m_ptr + rhs);}
        Iterator operator-(const int& rhs) const { return Iterator(m_ptr - rhs);}

        friend bool operator== (const Iterator& a, const Iterator& b) { return a.m_ptr == b.m_ptr; }
        friend bool operator!= (const Iterator& a, const Iterator& b) { return !(a==b); }
        friend bool operator< (const Iterator&a, const Iterator& b) { return a.m_ptr <  b.m_ptr; }
        friend bool operator<=(const Iterator&a, const Iterator& b) { return a.m_ptr <= b.m_ptr; }
        friend bool operator> (const Iterator&a, const Iterator& b) { return a.m_ptr >  b.m_ptr; }
        friend bool operator>=(const Iterator&a, const Iterator& b) { return a.m_ptr >= b.m_ptr; }

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

            Iterator& operator=(size_t rhs) { inner_offset = rhs; return *this;}
            Iterator& operator=(const Iterator &rhs) { inner_offset = rhs.inner_offset; return *this;}
            Iterator& operator+=(const int& rhs) { inner_offset += rhs; return *this;}
            Iterator& operator-=(const int& rhs) { inner_offset -= rhs; return *this;}

            reference operator[](const difference_type& rhs) {return data_begin[outer_offset(inner_offset + rhs)];}
            value_type const& operator[](const difference_type& rhs) const {return data_begin[outer_offset(inner_offset + rhs)];}
            reference operator*() { return data_begin[outer_offset(inner_offset)]; }
            value_type const& operator*() const { return *(data_begin[outer_offset(inner_offset)]); }
            pointer operator->() { return data_begin + outer_offset(inner_offset); }

            Iterator& operator++() { inner_offset++; return *this; }
            Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
            Iterator& operator--() { --inner_offset; return *this; }
            Iterator operator--(int) { Iterator tmp = *this; --(*this); return tmp; }

            Iterator operator+(const int& rhs) const { return Iterator(data_begin, di, dr, seq, inner_offset + rhs);}
            Iterator operator-(const int& rhs) const { return Iterator(data_begin, di, dr, seq, inner_offset - rhs);}

            friend bool operator== (const Iterator& a, const Iterator& b) { return    a.inner_offset == b.inner_offset
                                                                                   && a.data_begin   == b.data_begin
                                                                                   && a.di           == b.di
                                                                                   && a.dr           == b.dr
                                                                                   && a.seq.start    == b.seq.start
                                                                                   && a.seq.n        == b.seq.n
                                                                                   && a.seq.stride   == b.seq.stride; }
            friend bool operator!= (const Iterator& a, const Iterator& b) { return !(a==b); }
            friend bool operator< (const Iterator&a, const Iterator& b) { return a.inner_offset <  b.inner_offset; }
            friend bool operator<=(const Iterator&a, const Iterator& b) { return a.inner_offset <= b.inner_offset; }
            friend bool operator> (const Iterator&a, const Iterator& b) { return a.inner_offset >  b.inner_offset; }
            friend bool operator>=(const Iterator&a, const Iterator& b) { return a.inner_offset >= b.inner_offset; }

        private:

            inline Slice::difference_type outer_offset(difference_type const& inner) const {

                // calculate the outer offset from the inner offset

                // first unravel the inner index into an index (i,j,k) for a 3-dim array with dims (dl, idx_sequence.n, dr)
                auto dv = std::div(inner, seq.n*dr);
                difference_type i = dv.quot;
                dv = std::div(dv.rem, dr);
                difference_type j = dv.quot*seq.stride + seq.start;
                difference_type k = dv.rem;

                // then flatten the index for a 3-dim array with dims (dl, di, dr)
                return i*di*dr + j*dr + k;
            }

            iter_type data_begin;
            size_t di, dr;
            Seq<size_t> seq;
            difference_type inner_offset;
        };

    public:

        using iterator        = Iterator<Type>;
        using const_iterator  = Iterator<Type const>;

        /**
         * @brief Slice represents a slice into the CustomIndexarray
         * @param dimensions the dimensions of the CustomIndexArray
         * @param start_iter An iterator to the first element of the data
         * @param idx_sequence_ A sequence of indices into the slice
         */
        Slice(MultiIndex const& dimensions,
              iter_type const& start_iter,
              Seq<size_t> idx_sequence_)
            : data_begin(start_iter)
            , idx_sequence(idx_sequence_)
            , m_dimensions(dimensions)
            , di(dimensions.template get<Tag>())
            , dr(product<details::MultiIndexPosition<Tag, MultiIndex>::value>(dimensions)/di)
            , dl(product(dimensions)/(di*dr))
        {
            assert( (size_t)idx_sequence.start + idx_sequence.n <= di );

            m_dimensions.template get<Tag>() = Index<Tag>(idx_sequence.n);
        }

        // returns the number of elements in a slice
        size_t numel() const {
            return dl*dr*(idx_sequence.n);
        }

        // returns an stl-compatible random access iterator into the slice
        iterator begin() {
            return iterator(data_begin, di, dr, idx_sequence, 0);
        }

         // returns an stl-compatible random access iterator into the slice
        const_iterator begin() const {
            return const_iterator(data_begin, di, dr, idx_sequence, 0);
        }

         // returns an stl-compatible end random access iterator into the slice
        iterator end() {
            return iterator(data_begin, di, dr, idx_sequence, numel());
        }

         // returns an stl-compatible end random access iterator into the slice
        const_iterator end() const {
            return const_iterator(data_begin, di, dr, idx_sequence, numel());
        }

        // copies the slice elements into a CustomIndexArray of appropriate dimension
        CustomIndexArray<Type, Tags...> as_array(){
            CustomIndexArray<Type, Tags...> array(m_dimensions);
            Eigen::Index idx = 0;
            for (auto& v : *this)
            {
                array.array()[idx++] = v;
            };
            return array;
        }

        // comparison operators
        friend bool operator== (const Slice& a, const Slice& b) {
                return    a.data_begin == b.data_begin
                       && a.idx_sequence.start == b.idx_sequence.start
                       && a.idx_sequence.n == b.idx_sequence.n
                       && a.idx_sequence.stride == b.idx_sequence.stride
                       && a.m_dimensions == b.m_dimensions;
        }

        friend bool operator!= (const Slice& a, const Slice& b) { return !(a==b); }

        /**
         * @brief slice creates a subslice from the current slice
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
         * @brief slice creates a subslice from the current slice
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

    private:

        iter_type data_begin;
        Seq<size_t> idx_sequence;
        MultiIndex m_dimensions;
        size_t di, dr, dl;
    };

    // An array storying the size of each category
    MultiIndex m_dimensions;

    // number of elements stored
    size_t m_numel;

    // An array containing the elements
    InternalArrayType m_y{};

public:

    using iterator          = Iterator<Type>;
    using const_iterator    = Iterator<Type const>;

    /**
     * @brief begin returns a start iterator for the elements
     * @return random access iterator
     */
    iterator begin()
    {
        return iterator(array().data());
    }

    /**
     * @brief begin returns a start iterator for the elements
     * @return random access iterator
     */
    const_iterator begin() const
    {
        return const_iterator(array().data());
    }

    /**
     * @brief begin returns an end iterator for the elements
     * @return random access iterator
     */
    iterator end()
    {
        return iterator(array().data() + array().size());
    }

    /**
     * @brief begin returns an end iterator for the elements
     * @return random access iterator
     */
    const_iterator end() const
    {
        return const_iterator(array().data() + array().size());
    }

    /**
     * @brief slice creates a slice into the multidimensional array
     * @tparam The Tag corresponding to the dimension of the slice
     * @param idx_seq An index sequence, consisting of the first index,
     *                the number of indices and a stride
     * @return The slice
     */
    template <typename Tag>
    Slice<Tag, iterator> slice(Seq<size_t> idx_seq)
    {
        return Slice<Tag, iterator>(m_dimensions, begin(), idx_seq);
    }

    /**
     * @brief slice creates a slice into the multidimensional array
     * @tparam The Tag corresponding to the dimension of the slice
     * @param idx_seq An index sequence, consisting of the first index,
     *                the number of indices and a stride
     * @return The slice
     */
    template <typename Tag>
    Slice<Tag, const_iterator> slice(Seq<size_t> idx_seq) const
    {
        return Slice<Tag, const_iterator>(m_dimensions, begin(), idx_seq);
    }

};

} // namespace epi

#endif // CUSTOMINDEXARRAY_H
