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
#pragma once

#include "memilio/math/eigen.h"
#include "memilio/utils/metaprogramming.h"
#include <utility>
#include <iostream>
#include <iterator>

namespace mio
{

/**
 * @brief sequence of indices
 */
template <typename T>
struct Seq {
    Seq(T start_, T n_, T stride_ = 1)
        : start(start_)
        , n(n_)
        , stride(stride_)
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
    static constexpr bool value = (std::remove_reference_t<M>::RowsAtCompileTime == Eigen::Dynamic &&
                                   std::remove_reference_t<M>::ColsAtCompileTime == 1) ||
                                  (std::remove_reference_t<M>::RowsAtCompileTime == 1 &&
                                   std::remove_reference_t<M>::ColsAtCompileTime == Eigen::Dynamic);
};

/**
 * @brief check if Eigen::Matrix type M is a dynamic matrix type.
 */
template <class M>
struct is_dynamic_matrix {
    static constexpr bool value = std::remove_reference_t<M>::RowsAtCompileTime == Eigen::Dynamic &&
                                  std::remove_reference_t<M>::ColsAtCompileTime == Eigen::Dynamic;
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
 * helper to get the matrix type from an eigen expression 
 * with correct const volatile qualitfications.
 */
template <class M>
struct CVPlainMatrix {
    using Type = typename M::PlainMatrix;
};
template <class M>
struct CVPlainMatrix<Eigen::Ref<const M>> {
    using Type = const M;
};
template <class M>
using CVPlainMatrixT = typename CVPlainMatrix<M>::Type;

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
    return Eigen::Map<CVPlainMatrixT<std::decay_t<V>>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
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
    assert(cols.start + cols.stride * cols.n <= m.cols());

    auto majSpec   = std::remove_reference_t<M>::IsRowMajor ? rows : cols;
    auto minSpec   = std::remove_reference_t<M>::IsRowMajor ? cols : rows;
    auto majStride = majSpec.stride * minor_size(m);
    auto minStride = minSpec.stride;
    auto data      = m.data() + majSpec.start * minor_size(m) + minSpec.start;

    return Eigen::Map<CVPlainMatrixT<std::remove_reference_t<M>>, Eigen::Unaligned,
                      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(data, rows.n, cols.n, {majStride, minStride});
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

    return Eigen::Map<CVPlainMatrixT<std::remove_reference_t<M>>>(m.data(), rows, cols);
}

/**
 * template utility.
 * Defines value = true if M is an Eigen matrix expression.
 * Defines value = false, otherwise.
 */
template <class M>
using is_matrix_expression = std::is_base_of<Eigen::EigenBase<M>, M>;

/**
 * coefficient wise maximum of two matrices.
 * @param a a matrix expression
 * @param b a matrix expression of the same shape as a
 * @return a matrix expression the shape of a with each coefficient the maximum of the coefficients of a and b.
 */
template <class A, class B>
auto max(const Eigen::MatrixBase<A>& a, B&& b)
{
    return a.binaryExpr(std::forward<B>(b), [](auto a_i, auto b_i) {
        return std::max(a_i, b_i);
    });
}

/**
 * Maps a random access range (i.e. anything with size() and operator[], e.g. std::vector) onto a 
 * Eigen array expression. Returns a column array expression ´a´ where a[i] = f(v[i]).
 * The returned expression stores a reference to the range, lifetime of the range must exceed
 * lifetime of the return.
 * @param v a random access range.
 * @param f a function that returns a numeric scalar for each element of v. 
 * @return an array expression ´a´ the same size as v where a[i] = f(v[i]).
 */
template <class Rng, class F>
auto map(const Rng& v, F f)
{
    using Result = std::invoke_result_t<F, const typename Rng::value_type&>;
    using Scalar = std::decay_t<Result>;
    using Array  = Eigen::Array<Scalar, Eigen::Dynamic, 1>;
    return Array::NullaryExpr(v.size(), [f, &v](Eigen::Index i) -> Scalar {
        return f(v[size_t(i)]);
    });
}

namespace details
{
//true if elements returned by matrix(i, j) are references where matrix is of type M;
//false if the elements are temporaries, e.g. for expressions like Eigen::MatrixXd::Constant(r, c, v).
template <class M>
using IsElementReference =
    std::is_reference<decltype(std::declval<M>()(std::declval<Eigen::Index>(), std::declval<Eigen::Index>()))>;
} // namespace details

/**
 * iterate over elements of eigen matrix expressions in row major order.
 * @tparam M an eigen matrix expression type, can be in memory or lazy.
 * @tparam IsConst true if the elements of the matrix are not modifiable.
 */
template <class M, bool IsConst = false>
class RowMajorIterator
{
public:
    using MatrixRef = std::conditional_t<IsConst, const M&, M&>;
    using MatrixPtr = std::conditional_t<IsConst, const M*, M*>;
    static_assert(IsConst || details::IsElementReference<M>::value,
                  "Iterator must be const if matrix is not in memory.");

    using iterator_category = std::random_access_iterator_tag;
    using value_type        = typename M::Scalar;
    using reference         = std::conditional_t<IsConst, const value_type&, value_type&>;
    using difference_type   = Eigen::Index;

    // Operator-> returns a pointer to an element if the element can be referenced.
    // If the matrix returns temporaries, a proxy is returned that holds a local copy
    // and forwards the address of that copy. Adress is pointer-to-const because modifying
    // the local inside the proxy is almost certainly not what is intended.
    struct Proxy {
        value_type value;
        const value_type* operator->() const
        {
            return &value;
        }
    };
    using pointer = std::conditional_t<details::IsElementReference<MatrixRef>::value,
                                       std::conditional_t<IsConst, const value_type*, value_type*>, Proxy>;

    /**
     * Create an iterator that points to a specific element.
     * Indexing is done by flat indices i = r * nc + c where r is the row index, c is the column index and nc is the number of columns.
     * @param m reference of a matrix expression. Only a reference is stored, mind the lifetime of the object.
     * @param i flat index of the element pointed to.
     */
    RowMajorIterator(MatrixRef m, Eigen::Index i)
        : m_matrix(&m)
        , m_i(i)
    {
    }

    /**
     * pre increment operator.
     */
    RowMajorIterator& operator++()
    {
        ++m_i;
        return *this;
    }
    /**
     * post increment operator.
     */
    RowMajorIterator operator++(int)
    {
        auto cpy = *this;
        ++m_i;
        return cpy;
    }
    /**
     * random access, add n to index of this.
     */
    RowMajorIterator operator+(difference_type n) const
    {
        return {*m_matrix, m_i + n};
    }
    /**
     * random access, add n to index of iterator.
     */
    friend RowMajorIterator operator+(difference_type n, const RowMajorIterator& iter)
    {
        return {*iter.m_matrix, iter.m_i + n};
    }
    /**
     * add n to index of this.
     */
    RowMajorIterator& operator+=(difference_type n)
    {
        m_i += n;
        return *this;
    }
    /**
     * pre decrement operator.
     */
    RowMajorIterator& operator--()
    {
        --m_i;
        return *this;
    }
    /**
     * post decrement operator.
     */
    RowMajorIterator operator--(int)
    {
        auto cpy = *this;
        --m_i;
        return cpy;
    }
    /**
     * take n from the index of this.
     */
    RowMajorIterator operator-(difference_type n) const
    {
        return {*m_matrix, m_i - n};
    }
    /**
     * take n from the index of this.
     */
    RowMajorIterator& operator-=(difference_type n)
    {
        m_i -= n;
        return *this;
    }
    /**
     * calculate the distance between this and r.
     */
    difference_type operator-(RowMajorIterator r) const
    {
        return m_i - r.m_i;
    }

    /**
     * dereference this to get the element pointed to.
     * may return a temporary value instead of a reference to the element if the matrix is not evaluated 
     * in memory, e.g. in the case of Eigen::MatrixXd::Constant(r, c, v).
     */
    decltype(auto) operator*() const
    {
        return (*m_matrix)(m_i / m_matrix->cols(), m_i % m_matrix->cols());
    }

    /**
     * get a pointer to the element this iterator points to.
     * may return a proxy instead of a reference to the element if the matrix is not evaluated 
     * in memory, e.g. in the case of Eigen::MatrixXd::Constant(r, c, v).
     * The proxy stores a copy of the element and forwards the address of this copy.
     * @{
     */
    template <class Dummy = MatrixRef, std::enable_if_t<details::IsElementReference<Dummy>::value, void*> = nullptr>
    pointer operator->() const
    {
        return &(**this);
    }
    template <class Dummy = MatrixRef, std::enable_if_t<!details::IsElementReference<Dummy>::value, void*> = nullptr>
    pointer operator->() const
    {
        return Proxy{**this};
    }
    /**@}*/

    /**
     * lesser comparison operator.
     */
    bool operator<(const RowMajorIterator& other) const
    {
        return m_i < other.m_i;
    }
    /**
     * greater than comparison operator.
     */
    bool operator>(const RowMajorIterator& other) const
    {
        return other < *this;
    }
    /**
     * lesser equal comparison operator.
     */
    bool operator<=(const RowMajorIterator& other) const
    {
        return !(other < *this);
    }
    /**
     * greater equal comparison operator.
     */
    bool operator>=(const RowMajorIterator& other) const
    {
        return !(*this < other);
    }
    /**
     * equality comparison operator.
     */
    bool operator==(const RowMajorIterator& other) const
    {
        return !(*this < other) && !(other < *this);
    }
    /**
     * inequality comparison operator.
     */
    bool operator!=(const RowMajorIterator& other) const
    {
        return !(*this == other);
    }

private:
    MatrixPtr m_matrix;
    Eigen::Index m_i;
};

/**
 * create a non-const iterator to first element of the matrix m.
 * only enabled if the matrix is evaluated in memory, i.e. elements can be modified. 
 */
template <class M>
std::enable_if_t<conjunction_v<std::is_base_of<Eigen::EigenBase<M>, M>, details::IsElementReference<M>>,
                 RowMajorIterator<M, false>>
begin(M& m)
{
    return {m, 0};
}

/**
 * create a const iterator to first element of the matrix m.
 */
template <class M>
std::enable_if_t<std::is_base_of<Eigen::EigenBase<M>, M>::value, RowMajorIterator<M, true>> begin(const M& m)
{
    return {m, 0};
}

/**
 * create a const iterator to first element of the matrix m.
 */
template <class M>
std::enable_if_t<std::is_base_of<Eigen::EigenBase<M>, M>::value, RowMajorIterator<M, true>> cbegin(const M& m)
{
    return {m, 0};
}

/**
 * create a non-const end iterator for the matrix m.
 */
template <class M>
std::enable_if_t<conjunction_v<std::is_base_of<Eigen::EigenBase<M>, M>, details::IsElementReference<M>>,
                 RowMajorIterator<M, false>>
end(M& m)
{
    return {m, m.size()};
}

/**
 * create a const end iterator for the matrix m.
 */
template <class M>
std::enable_if_t<std::is_base_of<Eigen::EigenBase<M>, M>::value, RowMajorIterator<M, true>> end(const M& m)
{
    return {m, m.size()};
}

/**
 * create a non-const end iterator for the matrix m.
 */
template <class M>
std::enable_if_t<std::is_base_of<Eigen::EigenBase<M>, M>::value, RowMajorIterator<M, true>> cend(const M& m)
{
    return {m, m.size()};
}

} // namespace mio
