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
/**
 * A collection of classes to simplify handling of matrix shapes in meta programming.
 * 
 * Matrix shape types follow this model:
 * - subtype Matrix that is an alias of an Eigen matrix type
 * - at least one constructor that sets the shape from dimensions
 * - static member function `get_shape_of` that takes a compatible matrix expression and returns it's shape
 * - const member functions rows() and cols()
 * - trivially copyable, moveable, assignable, move assignable
 * - equality comparable
 */

#ifndef EPI_UTILS_MATRIX_SHAPE_H
#define EPI_UTILS_MATRIX_SHAPE_H

#include "memilio/math/eigen.h"
#include "memilio/math/eigen_util.h"
#include "memilio/io/io.h"

namespace mio
{
/**
 * shape of a rectangular matrix.
 * variable rows and cols.
 */
template <typename FP = double>
class RectMatrixShape
{
public:
    using Matrix = Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic>;

    /**
     * construct the shape of a rectangular matrix.
     * @param r number of rows.
     * @param c number of columns.
     */
    RectMatrixShape(Eigen::Index r, Eigen::Index c)
        : m_rows(r)
        , m_cols(c)
    {
    }

    /**
     * extract the shape of a rectangular matrix.
     * @param m matrix to take the shape of.
     * @tparam ME matrix expression.
     */
    template <class ME>
    static RectMatrixShape get_shape_of(const Eigen::MatrixBase<ME>& m)
    {
        return {m.rows(), m.cols()};
    }

    /**
     * number of rows.
     */
    Eigen::Index rows() const
    {
        return m_rows;
    }
    /**
     * number of columns.
     */
    Eigen::Index cols() const
    {
        return m_cols;
    }

    /**
     * equality comparators.
     */
    bool operator==(const RectMatrixShape& other) const
    {
        return other.m_rows == m_rows && other.m_cols == m_cols;
    }
    bool operator!=(const RectMatrixShape& other) const
    {
        return !(*this == other);
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("RectMatrixShape");
        obj.add_element("Rows", rows());
        obj.add_element("Columns", cols());
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<RectMatrixShape> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("RectMatrixShape");
        auto r   = obj.expect_element("Rows", Tag<Eigen::Index>{});
        auto c   = obj.expect_element("Columns", Tag<Eigen::Index>{});
        return apply(
            io,
            [](auto&& r_, auto&& c_) -> IOResult<RectMatrixShape> {
                if (r_ > 0 && c_ > 0)
                    return success(RectMatrixShape{r_, c_});
                else
                    return failure(StatusCode::OutOfRange, "Rows and Columns of RectMatrixShape must be positive.");
            },
            r, c);
    }

private:
    Eigen::Index m_rows;
    Eigen::Index m_cols;
};

/**
 * shape of a square matrix.
 * rows() == cols()
 */
template <typename FP = double>
class SquareMatrixShape
{
public:
    using Matrix = Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic>;

    /**
     * construct a square matrix of dimensions r
     * @param r number of rows and columns
     */
    explicit SquareMatrixShape(Eigen::Index r)
        : m_rows(r)
    {
    }

    /**
     * extract the shape of a square matrix.
     * @param m matrix to take the shape of.
     * @tparam ME matrix expression.
     */
    template <class ME>
    static SquareMatrixShape get_shape_of(const Eigen::MatrixBase<ME>& m)
    {
        assert(m.rows() == m.cols());
        return SquareMatrixShape{m.rows()};
    }

    /**
     * number of rows.
     * equal to number of columns.
     */
    Eigen::Index rows() const
    {
        return m_rows;
    }
    /**
     * number of columns.
     * equal to number of rows.
     */
    Eigen::Index cols() const
    {
        return m_rows;
    }
    /**
     * number of rows or columns.
     */
    Eigen::Index size() const
    {
        return m_rows;
    }

    /**
     * equality comparators.
     */
    bool operator==(const SquareMatrixShape& other) const
    {
        return other.m_rows == m_rows;
    }
    bool operator!=(const SquareMatrixShape& other) const
    {
        return !(*this == other);
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("SquareMatrixShape");
        obj.add_element("Rows", rows());
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<SquareMatrixShape> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("SquareMatrixShape");
        auto r   = obj.expect_element("Rows", Tag<Eigen::Index>{});
        return apply(
            io,
            [](auto&& r_) -> IOResult<SquareMatrixShape> {
                if (r_ > 0)
                    return success(SquareMatrixShape{r_});
                else
                    return failure(StatusCode::OutOfRange, "Rows of SquareMatrixShape must be positive.");
            },
            r);
    }

private:
    Eigen::Index m_rows;
};

/**
 * shape of a column vector.
 * cols() == 1. 
 */
class ColumnVectorShape
{

public:
    using Matrix = Eigen::VectorXd;

    /**
     * construct the shape of a column vector.
     * @param r number of rows.
     */
    explicit ColumnVectorShape(Eigen::Index r)
        : m_rows(r)
    {
    }

    /**
     * extract the shape of a column vector.
     * @param m vector to take the shape of.
     * @tparam ME matrix expression.
     */
    template <class ME>
    static ColumnVectorShape get_shape_of(const Eigen::MatrixBase<ME>& m)
    {
        assert(m.cols() == 1);
        return ColumnVectorShape{m.rows()};
    }

    /**
     * number of rows.
     */
    Eigen::Index rows() const
    {
        return m_rows;
    }
    /**
     * number of columns.
     * equal to 1.
     */
    Eigen::Index cols() const
    {
        return 1;
    }
    /**
     * number of rows.
     */
    Eigen::Index size() const
    {
        return m_rows;
    }

    /**
     * equality comparators.
     */
    bool operator==(const ColumnVectorShape& other) const
    {
        return other.m_rows == m_rows;
    }
    bool operator!=(const ColumnVectorShape& other) const
    {
        return !(*this == other);
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("ColumnVectorShape");
        obj.add_element("Rows", rows());
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<ColumnVectorShape> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("ColumnVectorShape");
        auto r   = obj.expect_element("Rows", Tag<Eigen::Index>{});
        return apply(
            io,
            [](auto&& r_) -> IOResult<ColumnVectorShape> {
                if (r_ > 0)
                    return success(ColumnVectorShape{r_});
                else
                    return failure(StatusCode::OutOfRange, "Rows of ColumnVectorShape must be positive");
            },
            r);
    }

private:
    Eigen::Index m_rows;
};

} // namespace mio

#endif //EPI_UTILS_MATRIX_SHAPE_H
