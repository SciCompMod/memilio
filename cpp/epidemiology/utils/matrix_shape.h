/**
 * A collection of classes to simplify handling of matrix shapes in meta programming.
 * 
 * Matrix shape types follow this model:
 * - subtype Matrix that is an alias of an Eigen matrix type
 * - at least one constructor that sets the shape from dimensions
 * - static member function `get_shape_of` that takes a compatible matrix expression and returns it's shape
 * - const member functions rows() and cols()
 * - trivially copyable, moveable, assignable, move assignable
 */

#ifndef EPI_UTILS_MATRIX_SHAPE_H
#define EPI_UTILS_MATRIX_SHAPE_H

#include "epidemiology/utils/eigen.h"
#include "epidemiology/utils/eigen_util.h"

namespace epi
{
/**
 * shape of a rectangular matrix.
 * variable rows and cols.
 */
class RectMatrixShape
{
public:
    using Matrix = Eigen::MatrixXd;

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
    template <class ME, class = std::enable_if_t<is_matrix_expression<ME>::value, void>>
    static RectMatrixShape get_shape_of(const ME& m)
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
    bool operator==(const RectMatrixShape& other)
    {
        return other.m_rows == m_rows && other.m_cols == m_cols;
    }
    bool operator!=(const RectMatrixShape& other)
    {
        return !(*this == other);
    }

private:
    Eigen::Index m_rows;
    Eigen::Index m_cols;
};

/**
 * shape of a square matrix.
 * rows() == cols()
 */
class SquareMatrixShape
{
public:
    using Matrix = Eigen::MatrixXd;

    /**
     * construct a square matrix of dimensions r
     * @param r number of rows and columns
     */
    SquareMatrixShape(Eigen::Index r)
        : m_rows(r)
    {
    }

    /**
     * extract the shape of a square matrix.
     * @param m matrix to take the shape of.
     * @tparam ME matrix expression.
     */
    template <class ME, class = std::enable_if_t<is_matrix_expression<ME>::value, void>>
    static SquareMatrixShape get_shape_of(const ME& m)
    {
        assert(m.rows() == m.cols());
        return {m.rows()};
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
    bool operator==(const SquareMatrixShape& other)
    {
        return other.m_rows == m_rows;
    }
    bool operator!=(const SquareMatrixShape& other)
    {
        return !(*this == other);
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
    ColumnVectorShape(Eigen::Index r)
        : m_rows(r)
    {
    }

    /**
     * extract the shape of a column vector.
     * @param m vector to take the shape of.
     * @tparam ME matrix expression.
     */
    template <class ME, class = std::enable_if_t<is_matrix_expression<ME>::value, void>>
    static ColumnVectorShape get_shape_of(const ME& m)
    {
        assert(m.cols() == 1);
        return {m.rows()};
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
    bool operator==(const ColumnVectorShape& other)
    {
        return other.m_rows == m_rows;
    }
    bool operator!=(const ColumnVectorShape& other)
    {
        return !(*this == other);
    }

private:
    Eigen::Index m_rows;
};

} // namespace epi

#endif //EPI_UTILS_MATRIX_SHAPE_H