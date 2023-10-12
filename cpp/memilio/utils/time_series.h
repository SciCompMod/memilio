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
#ifndef EPI_TIME_SERIES_H
#define EPI_TIME_SERIES_H

#include "memilio/io/io.h"
#include "memilio/math/eigen.h"
#include "memilio/math/eigen_util.h"
#include "memilio/utils/stl_util.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/math/floating_point.h"

#include <iterator>
#include <vector>
#include <map>
#include <ostream>

namespace mio
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
 * Time and values of a single point are stored together in memory: 
 * {t0, v0[0], v0[1], ...}, {t1, v1[0], v1[1], ...}, {t2, v20, ...
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

    using size_type                   = Eigen::Index;
    using iterator                    = TimeSeriesValueIterator<FP, false>;
    using const_iterator              = TimeSeriesValueIterator<FP, true>;
    using time_iterator               = TimeSeriesTimeIterator<FP, false>;
    using const_time_iterator         = TimeSeriesTimeIterator<FP, true>;
    using value_type                  = typename iterator::value_type;
    using reference                   = typename iterator::reference;
    using difference_type             = typename iterator::difference_type;
    using reverse_iterator            = std::reverse_iterator<iterator>;
    using const_reverse_iterator      = std::reverse_iterator<const_iterator>;
    using reverse_time_iterator       = std::reverse_iterator<time_iterator>;
    using const_reverse_time_iterator = std::reverse_iterator<const_time_iterator>;

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

    /**
     * @brief constructs TimeSeries instance and initializes it with zeros
     * @param num_time_points number of time steps
     * @param num_elements number of compartiments * number of groups
     * @return
     */
    static TimeSeries zero(Eigen::Index num_time_points, Eigen::Index num_elements)
    {
        TimeSeries value_matrix(num_elements);
        value_matrix.m_data            = Matrix::Zero(num_elements + 1, num_time_points);
        value_matrix.m_num_time_points = num_time_points;

        return value_matrix;
    }

    /** copy assignment */
    TimeSeries& operator=(const TimeSeries& other)
    {
        if (get_num_elements() == other.get_num_elements()) {
            //assign with preserved storage if possible
            reserve(other.m_num_time_points);
            m_num_time_points = other.m_num_time_points;
            get_valid_block() = other.get_valid_block();
        }
        else {
            //assign with new storage
            auto data = Matrix(other.get_num_elements() + 1, details::next_pow2(other.m_num_time_points));
            data.leftCols(other.m_num_time_points) = other.get_valid_block();
            m_data                                 = std::move(data);
            m_num_time_points                      = other.m_num_time_points;
        }
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
        get_last_time() = t;
        return get_last_value();
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
     * remove time point.
     * @param i index to remove
     */
    void remove_time_point(Eigen::Index i)
    {
        assert(i >= 0 && i < m_num_time_points);
        m_num_time_points -= 1;
        for (auto j = i; j < m_num_time_points; ++j) {
            m_data.col(j) = m_data.col(j + 1);
        }
    }

    void remove_last_time_point()
    {
        remove_time_point(m_num_time_points - 1);
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
     * packed {t0, v0[0], v0[1], ...}, {t1, v1[0], v1[1], ...}, ...
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

    /**
     * Matrix expression that contains one time point per column.
     * Each column contains the corresponding time at index 0.
     * @{
     */
    auto matrix()
    {
        return get_valid_block();
    }
    auto matrix() const
    {
        return get_valid_block();
    }
    /** @} */

    /*********************
     * 
     * Iterator interface to iterate over values.
     * vector at time point 0 -> vector at time point 1 -> ...
     * 
     *********************/
    iterator begin()
    {
        return {&m_data, 0};
    }

    iterator end()
    {
        return {&m_data, m_num_time_points};
    }

    const_iterator begin() const
    {
        return {&m_data, 0};
    }

    const_iterator end() const
    {
        return {&m_data, m_num_time_points};
    }

    const_iterator cbegin() const
    {
        return {&m_data, 0};
    }

    const_iterator cend() const
    {
        return {&m_data, m_num_time_points};
    }

    reverse_iterator rbegin()
    {
        return reverse_iterator{end()};
    }

    reverse_iterator rend()
    {
        return reverse_iterator{begin()};
    }

    const_reverse_iterator rbegin() const
    {
        return crbegin();
    }

    const_reverse_iterator rend() const
    {
        return crend();
    }

    const_reverse_iterator crbegin() const
    {
        return const_reverse_iterator{cend()};
    }

    const_reverse_iterator crend() const
    {
        return const_reverse_iterator{cbegin()};
    }

    /*********************
     * 
     * Iterator interface to iterate over times.
     * time at point 0 -> time at point 1 -> ...
     * 
     *********************/
    Range<std::pair<time_iterator, time_iterator>> get_times()
    {
        return make_range(time_iterator{&m_data, 0}, time_iterator{&m_data, m_num_time_points});
    }

    Range<std::pair<const_time_iterator, const_time_iterator>> get_times() const
    {
        return get_const_times();
    }

    Range<std::pair<const_time_iterator, const_time_iterator>> get_const_times() const
    {
        return make_range(const_time_iterator{&m_data, 0}, const_time_iterator{&m_data, m_num_time_points});
    }

    Range<std::pair<reverse_time_iterator, reverse_time_iterator>> get_reverse_times()
    {
        auto time_range = get_times();
        return make_range(reverse_time_iterator{time_range.end()}, reverse_time_iterator{time_range.begin()});
    }

    Range<std::pair<const_reverse_time_iterator, const_reverse_time_iterator>> get_reverse_times() const
    {
        return get_const_reverse_times();
    }

    Range<std::pair<const_reverse_time_iterator, const_reverse_time_iterator>> get_const_reverse_times() const
    {
        auto time_range = get_const_times();
        return make_range(const_reverse_time_iterator{time_range.end()},
                          const_reverse_time_iterator{time_range.begin()});
    }

    /**
     * print this object (googletest)
     */
    friend void PrintTo(const TimeSeries& self, std::ostream* os)
    {
        *os << '\n' << self.get_valid_block();
    }

    template<class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("TimeSeries");
        obj.add_element("NumTimePoints", m_num_time_points);
        obj.add_element("Data", get_valid_block());
    }

    template<class IOContext>
    static IOResult<TimeSeries> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("TimeSeries");
        auto nt = obj.expect_element("NumTimePoints", Tag<Eigen::Index>{});
        auto data = obj.expect_element("Data", Tag<Matrix>{});
        return apply(io, [](auto&& nt_, auto&& data_) {
            auto ts = TimeSeries(data_.rows() - 1);
            ts.m_data.resize(data_.rows(), data_.cols());
            ts.m_data = data_;
            ts.m_num_time_points = nt_;
            return ts;
        }, nt, data);
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
    if constexpr (sizeof(Eigen::Index) == 8)
    {
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
    using MatrixPtr   = std::conditional_t<IsConst, const Matrix, Matrix>*;
    using VectorValue = typename decltype(
        std::declval<MatrixPtr>()->col(std::declval<Eigen::Index>()).tail(std::declval<Eigen::Index>()))::PlainObject;
    using VectorReference =
        decltype(std::declval<MatrixPtr>()->col(std::declval<Eigen::Index>()).tail(std::declval<Eigen::Index>()));
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
    using MatrixPtr = typename Traits::MatrixPtr;
    MatrixPtr m_matrix;
    Eigen::Index m_col_idx = -1;

public:
    TimeSeriesIteratorBase(MatrixPtr m, Eigen::Index col_idx = 0)
        : m_matrix(m)
        , m_col_idx(col_idx)
    {
        assert(m_matrix != nullptr);
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

    reference operator*() const
    {
        assert(m_col_idx >= 0 && m_col_idx < m_matrix->cols());
        return static_cast<const Derived&>(*this).get_reference();
    }

    pointer operator->() const
    {
        return *(*this);
    }

    reference operator[](difference_type i) const
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

    difference_type operator-(const TimeSeriesIteratorBase& other) const
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
        assert(m_matrix == other.m_matrix);
        return m_col_idx == other.m_col_idx;
    }

    bool operator!=(const TimeSeriesIteratorBase& other) const
    {
        assert(m_matrix == other.m_matrix);
        return !(*this == other);
    }

    bool operator<(const TimeSeriesIteratorBase& other) const
    {
        assert(m_matrix == other.m_matrix);
        return m_col_idx < other.m_col_idx;
    }

    bool operator>(const TimeSeriesIteratorBase& other) const
    {
        assert(m_matrix == other.m_matrix);
        return m_col_idx > other.m_col_idx;
    }

    bool operator<=(const TimeSeriesIteratorBase& other) const
    {
        assert(m_matrix == other.m_matrix);
        return m_col_idx <= other.m_col_idx;
    }

    bool operator>=(const TimeSeriesIteratorBase& other) const
    {
        assert(m_matrix == other.m_matrix);
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

    reference get_reference() const
    {
        return m_matrix->col(m_col_idx).tail(m_matrix->rows() - 1);
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

    reference get_reference() const
    {
        return m_matrix->coeffRef(0, m_col_idx);
    }
};

/**
 * find the value in the time series at time t_search starting from the end.
 * @param ts TimeSeries to seach
 * @param t_search a time point
 * @param abs_tol absolute floating point tolerance for equality of time values
 * @param rel_tol relative floating point tolerance for equality of time values
 * @return TimeSeries::reverse_iterator that points to ts[t_search] or ts.rend()
 */
template <class TS, class FP>
decltype(std::declval<TS>().rend()) find_value_reverse(TS&& ts, FP t_search, FP abs_tol = 0, FP rel_tol = 0)
{
    auto iter_t = find_if(ts.get_reverse_times().begin(), ts.get_reverse_times().end(), [=](auto t) {
        return floating_point_equal(t, t_search, abs_tol, rel_tol);
    });
    if (iter_t != ts.get_reverse_times().end()) {
        return ts.rbegin() + (iter_t - ts.get_reverse_times().begin());
    }
    return ts.rend();
}

} // namespace mio

#endif
