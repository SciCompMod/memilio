#ifndef EPI_TIME_SERIES_H
#define EPI_TIME_SERIES_H

#include <epidemiology/eigen_util.h>
#include <epidemiology/stl_util.h>
#include <Eigen/Core>
#include <vector>
#include <map>

namespace epi
{

namespace details
{
    /** round an integer to the nearest greater power of 2 */
    inline Eigen::Index next_pow2(Eigen::Index i);
} // namespace details

template <class FP, bool IsConstant>
class TimeSeriesValueIterator;
template <class FP, bool IsConstant>
class TimeSeriesTimeIterator;

/**
 * stores vectors of values at time points (or some other abstract variable)
 * the value at each time point is a vector.
 * size of the vector is determined at runtime but is equal for all time points.
 * grows efficiently (like std::vector) in time dimension.
 * @tparam FP any floating point like type accepted by Eigen
 */
template <class FP>
class TimeSeries
{
public:
    /** type that stores the data */
    using Matrix = Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    /** base type of expressions of vector values at a time point */
    using Vector = Eigen::Matrix<FP, Eigen::Dynamic, 1>;

    using size_type              = Eigen::Index;
    using iterator               = TimeSeriesValueIterator<FP, false>;
    using const_iterator         = TimeSeriesValueIterator<FP, true>;
    using time_iterator          = TimeSeriesTimeIterator<FP, false>;
    using const_time_iterator    = TimeSeriesTimeIterator<FP, true>;
    using value_type             = typename iterator::value_type;
    using reference              = typename iterator::reference;
    using difference_type        = typename iterator::difference_type;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    /**
     * initialize empty TimeSeries.
     * @param num_elements size of vector at each time point
     */
    TimeSeries(Eigen::Index num_elements)
        : m_data(num_elements + 1, 0)
        , m_num_time_points(0)
    {
        assert(num_elements >= 0);
    }

    /**
     * initialize TimeSeries with one time point.
     * size of vector at each time point determined from input.
     * @tparam Expr any type that can be assigned to TimeSeries::Vector
     * @param t time of initial time point
     * @param expr vector expression that evaluates to value at first time point
     */
    template <typename Expr>
    TimeSeries(FP t, Expr&& expr)
        : m_data(expr.rows() + 1, 1)
        , m_num_time_points(1)
    {
        auto col              = m_data.col(0);
        col(0)                = t;
        col.tail(expr.rows()) = expr;
    }

    /** copy ctor */
    TimeSeries(const TimeSeries& other)
        : m_data(other.get_num_elements() + 1, details::next_pow2(other.m_num_time_points))
        , m_num_time_points(other.m_num_time_points)
    {
        get_valid_block() = other.get_valid_block();
    }

    /** copy assignment */
    TimeSeries& operator=(const TimeSeries& other)
    {
        auto data = Matrix(other.get_num_elements() + 1, details::next_pow2(other.m_num_time_points));
        data.leftCols(other.m_num_time_points) = other.get_valid_block();
        m_data                                 = std::move(data);
        m_num_time_points                      = other.m_num_time_points;
        return *this;
    }

    /** move ctor and assignment */
    TimeSeries(TimeSeries&& other) = default;
    TimeSeries& operator=(TimeSeries&& other) = default;

    /**
     * number of time points in the series
     */
    Eigen::Index get_num_time_points() const
    {
        return m_num_time_points;
    }

    /**
     * number of elements of vector at each time point
     */
    Eigen::Index get_num_elements() const
    {
        return get_num_rows() - 1;
    }

    /**
     * number of rows in data storage (includes time)
     */
    Eigen::Index get_num_rows() const
    {
        return m_data.rows();
    }

    /**
     * add one uninitialized time point
     */
    Eigen::Ref<Vector> add_time_point()
    {
        add_time_point_noinit();
        return get_last_value();
    }

    /**
     * add one time point.
     * initialize time;
     */
    Eigen::Ref<Vector> add_time_point(FP t)
    {
        add_time_point_noinit();
        auto col = m_data.col(m_num_time_points - 1);
        col(0)   = t;
        return col.tail(get_num_elements());
    }

    /**
     * add one time point.
     * Initialize time and value.
     * Expr can be any vector expression assignable to TimeSeries::Vector.
     */
    template <class Expr>
    Eigen::Ref<Vector> add_time_point(FP t, Expr&& expr)
    {
        auto value = add_time_point(t);
        value      = expr;
        return value;
    }

    /**
     * time of time point at index i
     */
    FP& get_time(Eigen::Index i)
    {
        assert(i >= 0 && i < m_num_time_points);
        return m_data(0, i);
    }
    const FP& get_time(Eigen::Index i) const
    {
        assert(i >= 0 && i < m_num_time_points);
        return m_data(0, i);
    }

    /**
     * time of time point at index num_time_points - 1
     */
    FP& get_last_time()
    {
        return get_time(get_num_time_points() - 1);
    }
    const FP& get_last_time() const
    {
        return get_time(get_num_time_points() - 1);
    }

    /**
     * reference to value vector at time point i
     */
    Eigen::Ref<const Vector> get_value(Eigen::Index i) const
    {
        assert(i >= 0 && i < m_num_time_points);
        return m_data.col(i).segment(1, get_num_elements());
    }
    Eigen::Ref<Vector> get_value(Eigen::Index i)
    {
        assert(i >= 0 && i < m_num_time_points);
        return m_data.col(i).segment(1, get_num_elements());
    }
    Eigen::Ref<const Vector> operator[](Eigen::Index i) const
    {
        return get_value(i);
    }
    Eigen::Ref<Vector> operator[](Eigen::Index i)
    {
        return get_value(i);
    }

    /**
     * reference to value vector at time point (num_timepoints - 1)
     */
    Eigen::Ref<const Vector> get_last_value() const
    {
        return get_value(m_num_time_points - 1);
    }
    Eigen::Ref<Vector> get_last_value()
    {
        return get_value(m_num_time_points - 1);
    }

    /**
     * reserve capacity for n time points
     */
    void reserve(Eigen::Index n)
    {
        assert(n >= 0);
        if (n > get_capacity()) {
            m_data.conservativeResize(Eigen::NoChange, details::next_pow2(n));
        }
    }

    /**
     * current capacity
     */
    Eigen::Index get_capacity() const
    {
        return m_data.cols();
    }

    /**
     * raw data storage.
     * at least num_rows * num_time_points scalars.
     * packed t0, v0[0], v0[1], ..., t1, v1[0], ...
     */
    FP* data()
    {
        assert(m_data.innerStride() == 1);
        assert(m_data.outerStride() == m_data.rows());
        return m_data.data();
    }
    const FP* data() const
    {
        assert(m_data.innerStride() == 1);
        assert(m_data.outerStride() == m_data.rows());
        return m_data.data();
    }

    /*********************
     * 
     * Iterator interface to iterate over values.
     * vector at time point 0 -> vector at time point 1 -> ...
     * 
     *********************/
    iterator begin()
    {
        return {m_data, 0};
    }

    iterator end()
    {
        return {m_data, m_num_time_points};
    }

    const_iterator begin() const
    {
        return {m_data, 0};
    }

    const_iterator end() const
    {
        return {m_data, m_num_time_points};
    }

    const_iterator cbegin() const
    {
        return {m_data, 0};
    }

    const_iterator cend() const
    {
        return {m_data, m_num_time_points};
    }

    iterator rbegin()
    {
        return {m_data, 0};
    }

    iterator rend()
    {
        return {m_data, m_num_time_points};
    }

    const_iterator rbegin() const
    {
        return {m_data, 0};
    }

    const_iterator rend() const
    {
        return {m_data, m_num_time_points};
    }

    const_iterator crbegin() const
    {
        return {m_data, 0};
    }

    const_iterator crend() const
    {
        return {m_data, m_num_time_points};
    }

    /*********************
     * 
     * Iterator interface to iterate over times.
     * time at point 0 -> time at point 1 -> ...
     * 
     *********************/
    Range<std::pair<time_iterator, time_iterator>> get_times()
    {
        return make_range(time_iterator{m_data, 0}, time_iterator{m_data, m_num_time_points});
    }

    Range<std::pair<const_time_iterator, const_time_iterator>> get_times() const
    {
        return get_const_times();
    }

    Range<std::pair<const_time_iterator, const_time_iterator>> get_const_times() const
    {
        return make_range(const_time_iterator{m_data, 0}, const_time_iterator{m_data, m_num_time_points});
    }

    Range<std::pair<time_iterator, time_iterator>> get_reverse_times()
    {
        return make_range(const_time_iterator{m_data, 0}, const_time_iterator{m_data, m_num_time_points});
    }

    Range<std::pair<const_time_iterator, const_time_iterator>> get_reverse_times() const
    {
        return get_const_reverse_times();
    }

    Range<std::pair<const_time_iterator, const_time_iterator>> get_const_reverse_times() const
    {
        return make_range(const_time_iterator{m_data, 0}, const_time_iterator{m_data, m_num_time_points});
    }

private:
    void add_time_point_noinit()
    {
        reserve(m_num_time_points + 1);
        ++m_num_time_points;
    }
    /** currently occupied block of storage */
    auto get_valid_block()
    {
        return m_data.leftCols(get_num_time_points());
    }
    auto get_valid_block() const
    {
        return m_data.leftCols(get_num_time_points());
    }

    /** data storage */
    Matrix m_data;
    /** number of time points (i.e. occupied columns in m_data) */
    Eigen::Index m_num_time_points;
};

namespace details
{
    /** round integer up to the next higher power of 2 */
    inline Eigen::Index next_pow2(Eigen::Index i)
    {
        //https://stackoverflow.com/questions/1322510/given-an-integer-how-do-i-find-the-next-largest-power-of-two-using-bit-twiddlin
        --i;
        i |= i >> 1;
        i |= i >> 2;
        i |= i >> 4;
        i |= i >> 8;
        i |= i >> 16;
        if (sizeof(Eigen::Index) == 8) {
            i |= i >> 32;
        }
        ++i;
        return i;
    }

    /**
     * type traits for time series iterators
     */
    template <class FP, bool IsConst>
    struct TimeSeriesIterTraits {
        static bool is_const()
        {
            return IsConst;
        }
        using Matrix      = typename TimeSeries<FP>::Matrix;
        using MatrixRef   = std::conditional_t<IsConst, const Matrix, Matrix>&;
        using VectorValue = typename decltype(std::declval<MatrixRef>()
                                                  .col(std::declval<Eigen::Index>())
                                                  .tail(std::declval<Eigen::Index>()))::PlainObject;
        using VectorReference =
            decltype(std::declval<MatrixRef>().col(std::declval<Eigen::Index>()).tail(std::declval<Eigen::Index>()));
        using TimeValue     = FP;
        using TimeReference = std::conditional_t<IsConst, const FP&, FP&>;
    };

    /** base class for TimeSeries iterators that iterate by time point (i.e. column) 
     * @tparam Derived Iterator derived type, provides get_reference member function
     * @tparam FP floating point type of the TimeSeries
     * @tparam IsConstIter true for const_iterator
     * @tparam ValueType define iterator::value_type
     * @tparam ReferenceType define iterator::reference, must be the same type as returned by Derived::get_reference
    */
    template <class Derived, class FP, bool IsConstIter, class ValueType, class ReferenceType>
    class TimeSeriesIteratorBase
    {
    protected:
        using Traits    = details::TimeSeriesIterTraits<FP, IsConstIter>;
        using MatrixRef = typename Traits::MatrixRef;
        MatrixRef m_matrix;
        Eigen::Index m_col_idx;

    public:
        TimeSeriesIteratorBase(MatrixRef m, Eigen::Index col_idx = 0)
            : m_matrix(m)
            , m_col_idx(col_idx)
        {
        }

        using difference_type   = std::ptrdiff_t;
        using iterator_category = std::random_access_iterator_tag;
        using reference         = ReferenceType;
        using value_type        = ValueType;
        struct pointer {
            reference m_ref;
            reference* operator->()
            {
                return &m_ref;
            }
        };

        reference operator*()
        {
            return static_cast<Derived&>(*this).get_reference();
        }

        pointer operator->()
        {
            return *(*this);
        }

        reference operator[](difference_type i)
        {
            return *((*this) + i);
        }

        Derived& operator+=(difference_type i)
        {
            m_col_idx += i;
            return static_cast<Derived&>(*this);
        }

        Derived operator+(difference_type i) const
        {
            auto tmp = static_cast<const Derived&>(*this);
            tmp += i;
            return tmp;
        }

        friend Derived operator+(difference_type i, const TimeSeriesIteratorBase& b)
        {
            auto tmp = static_cast<const Derived&>(b);
            tmp += i;
            return tmp;
        }

        Derived& operator-=(difference_type i)
        {
            m_col_idx -= i;
            return static_cast<Derived&>(*this);
        }

        Derived operator-(difference_type i) const
        {
            auto tmp = static_cast<const Derived&>(*this);
            tmp -= i;
            return tmp;
        }

        difference_type operator-(const TimeSeriesIteratorBase& other)
        {
            return m_col_idx - other.m_col_idx;
        }

        Derived& operator++()
        {
            ++m_col_idx;
            return static_cast<Derived&>(*this);
        }

        Derived operator++(int)
        {
            auto tmp = static_cast<Derived&>(*this);
            ++(*this);
            return tmp;
        }

        Derived& operator--()
        {
            --m_col_idx;
            return static_cast<Derived&>(*this);
        }

        Derived operator--(int)
        {
            auto tmp = static_cast<Derived&>(*this);
            --(*this);
            return tmp;
        }

        bool operator==(const TimeSeriesIteratorBase& other) const
        {
            assert(&m_matrix == &other.m_matrix);
            return m_col_idx == other.m_col_idx;
        }

        bool operator!=(const TimeSeriesIteratorBase& other) const
        {
            assert(&m_matrix == &other.m_matrix);
            return !(*this == other);
        }

        bool operator<(const TimeSeriesIteratorBase& other) const
        {
            assert(&m_matrix == &other.m_matrix);
            return m_col_idx < other.m_col_idx;
        }

        bool operator>(const TimeSeriesIteratorBase& other) const
        {
            assert(&m_matrix == &other.m_matrix);
            return m_col_idx > other.m_col_idx;
        }

        bool operator<=(const TimeSeriesIteratorBase& other) const
        {
            assert(&m_matrix == &other.m_matrix);
            return m_col_idx <= other.m_col_idx;
        }

        bool operator>=(const TimeSeriesIteratorBase& other) const
        {
            assert(&m_matrix == &other.m_matrix);
            return m_col_idx >= other.m_col_idx;
        }
    };
} // namespace details

/** Iterate over vector values of a time series by time point */
template <class FP, bool IsConstIter>
class TimeSeriesValueIterator
    : public details::TimeSeriesIteratorBase<TimeSeriesValueIterator<FP, IsConstIter>, FP, IsConstIter,
                                             typename details::TimeSeriesIterTraits<FP, IsConstIter>::VectorValue,
                                             typename details::TimeSeriesIterTraits<FP, IsConstIter>::VectorReference>
{
private:
    using Base =
        details::TimeSeriesIteratorBase<TimeSeriesValueIterator<FP, IsConstIter>, FP, IsConstIter,
                                        typename details::TimeSeriesIterTraits<FP, IsConstIter>::VectorValue,
                                        typename details::TimeSeriesIterTraits<FP, IsConstIter>::VectorReference>;

    using Base::m_col_idx;
    using Base::m_matrix;

public:
    /** ctor from matrix ref and col idx*/
    using Base::Base;

    using iterator_category = typename Base::iterator_category;
    using difference_type   = typename Base::difference_type;
    using reference         = typename Base::reference;
    using value_type        = typename Base::value_type;
    using pointer           = typename Base::pointer;

    reference get_reference()
    {
        return m_matrix.col(m_col_idx).tail(m_matrix.rows() - 1);
    }
};

/** Iterate over vector values of a time series by time point */
template <class FP, bool IsConstIter>
class TimeSeriesTimeIterator
    : public details::TimeSeriesIteratorBase<TimeSeriesTimeIterator<FP, IsConstIter>, FP, IsConstIter,
                                             typename details::TimeSeriesIterTraits<FP, IsConstIter>::TimeValue,
                                             typename details::TimeSeriesIterTraits<FP, IsConstIter>::TimeReference>
{
private:
    using Base =
        details::TimeSeriesIteratorBase<TimeSeriesTimeIterator<FP, IsConstIter>, FP, IsConstIter,
                                        typename details::TimeSeriesIterTraits<FP, IsConstIter>::TimeValue,
                                        typename details::TimeSeriesIterTraits<FP, IsConstIter>::TimeReference>;

    using Base::m_col_idx;
    using Base::m_matrix;

public:
    /** ctor from matrix ref and col idx*/
    using Base::Base;

    using iterator_category = typename Base::iterator_category;
    using difference_type   = typename Base::difference_type;
    using reference         = typename Base::reference;
    using value_type        = typename Base::value_type;
    using pointer           = typename Base::pointer;

    reference get_reference()
    {
        return m_matrix.coeffRef(0, m_col_idx);
    }
};

} // namespace epi

#endif