#include <Eigen/Core>
#include <utility>
#include <iostream>

namespace epi
{

/**
 * @brief sequence of indices
 */
template <typename T>
struct Seq {
    Seq(T start, T n, T stride = 1)
        : start(start)
        , n(n)
        , stride(stride)
    {
        assert(start >= 0);
        assert(n >= 0);
        assert(stride >= 1);
    }
    T start, n, stride = 1;
};

/**
 * @brief check if Eigen::Matrix type M is a dynamic vector type.
 */
template <class M>
struct is_dynamic_vector {
    static constexpr bool value =
        (std::remove_reference_t<M>::RowsAtCompileTime == Eigen::Dynamic && std::remove_reference_t<M>::ColsAtCompileTime == 1) ||
        (std::remove_reference_t<M>::RowsAtCompileTime == 1 && std::remove_reference_t<M>::ColsAtCompileTime == Eigen::Dynamic);
};

/**
 * @brief check if Eigen::Matrix type M is a dynamic matrix type.
 */
template <class M>
struct is_dynamic_matrix {
    static constexpr bool value =
        std::remove_reference_t<M>::RowsAtCompileTime == Eigen::Dynamic && std::remove_reference_t<M>::ColsAtCompileTime == Eigen::Dynamic;
};

/**
 * @brief number of rows (columns) of a row (column) major matrix.
 */
template <class M>
Eigen::Index major_size(M&& m)
{
    return std::remove_reference_t<M>::IsRowMajor ? m.rows() : m.cols();
}

/**
 * @brief number of columns (rows) of a row (column) major matrix.
 */
template <class M>
Eigen::Index minor_size(M&& m)
{
    return std::remove_reference_t<M>::IsRowMajor ? m.cols() : m.rows();
}

/**
 * @brief take a regular slice of a row or column vector.
 * The slices shares the same memory as the original vector, no copying is performed,
 * changes to the slice are also made to the original vector.
 * Assign to a different vector of compatible size if you need a copy.
 * @param v Row or column vector to take a slice of
 * @param elems sequence of row or column indices
 * @returns vector expression with selected entries from the input vector
 */
template <class V, std::enable_if_t<is_dynamic_vector<V>::value, int> = 0>
auto slice(V&& v, Seq<Eigen::Index> elems)
{
    return Eigen::Map<std::remove_reference_t<V>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
        v.data() + elems.start, elems.n, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{1, elems.stride});
}

/**
 * @brief take a regular slice of a matrix.
 * The slices shares the same memory as the original matrix, no copying is performed,
 * changes to the slice are also made to the original matrix.
 * Assign to a different matrix of compatible size if you need a copy.
 * @param m Matrix to take a slice of
 * @param rows sequence of row indices
 * @param cols sequence of column indices
 * @returns matrix expression with selected entries from the input matrix
 */
template <class M, std::enable_if_t<is_dynamic_matrix<M>::value, int> = 0>
auto slice(M&& m, Seq<Eigen::Index> rows, Seq<Eigen::Index> cols)
{
    assert(rows.start + rows.stride * rows.n <= m.rows());
    assert(cols.start + cols.stride * cols.n <= m.rows());

    auto majSpec   = std::remove_reference_t<M>::IsRowMajor ? rows : cols;
    auto minSpec   = std::remove_reference_t<M>::IsRowMajor ? cols : rows;
    auto majStride = majSpec.stride * minor_size(m);
    auto minStride = minSpec.stride;
    auto data      = m.data() + majSpec.start * minor_size(m) + minSpec.start;

    return Eigen::Map<std::remove_reference_t<M>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
        data, rows.n, cols.n, {majStride, minStride});
}

/**
 * @brief reshape the matrix.
 * Total number of entries before and after the reshape must be the same.
 * The new matrix shares memory with the input matrix, no copying is performed, changes to the new matrix
 * will be made to the input matrix as well.
 * Assign to another matrix of compatible size if you need a copy.
 * @param m matrix to reshape
 * @param rows number of rows of the new matrix
 * @param cols number of cols of the new matrix
 * @returns matrix expression with the same entries as the input matrix but new shape
 */
template <typename M>
auto reshape(M&& m, Eigen::Index rows, Eigen::Index cols)
{
    assert(rows * cols == m.rows() * m.cols());
    assert(rows >= 1);
    assert(cols >= 1);

    return Eigen::Map<std::remove_reference_t<M>>(m.data(), rows, cols);
}

/**
 * @brief overload gtest printer function for eigen matrices.
 * @note see https://stackoverflow.com/questions/25146997/teach-google-test-how-to-print-eigen-matrix
 */
template <class M>
struct MatrixPrintWrap : public M {
    friend void PrintTo(const MatrixPrintWrap& m, std::ostream* os)
    {
        (*os) << '\n' << m;
    }
};

/**
 * @brief wrap m for gtest printing
 * returns a reference to the original object, no copying or moving, mind the lifetime!
 */
template <class M>
const MatrixPrintWrap<M>& print_wrap(const M& m)
{
    return static_cast<const MatrixPrintWrap<M>&>(m);
}

} // namespace epi