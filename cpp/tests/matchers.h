#ifndef EPI_TESTS_MATCHERS_H
#define EPI_TESTS_MATCHERS_H

#include "epidemiology/utils/compiler_diagnostics.h"
#include "gmock/gmock.h"

/**
 * @brief overload gtest printer function for eigen matrices.
 * @note see https://stackoverflow.com/questions/25146997/teach-google-test-how-to-print-eigen-matrix
 */
template <class M>
struct MatrixPrintWrap : public M {
    friend void PrintTo(const MatrixPrintWrap& m, std::ostream* os)
    {
        if (m.rows() == 1)
        {            
            //print row vector inline
            (*os) << m;
        }
        else if (m.cols() == 1)
        {
            //print col vector inline transposed
            (*os) << m.transpose() << " T";
        }
        else
        {
            //print matrix on its own
            (*os) << '\n' << m;
        }
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

/**
 * gmock matcher that checks where the elements of a container are linearly spaced.
 * @param b minimum value
 * @param e maximum value
 * @param num_points number of linearly spaced points in [b, e]
 * @return matcher that accepts a stl container
 */
template <class T>
auto ElementsAreLinspace(T b, T e, size_t num_points)
{
    assert(num_points >= 2);

    std::vector<T> values(num_points);
    for (size_t i = 0; i < num_points; i++) {
        values[i] = b + i * (e - b) / (num_points - 1);
    }
    return testing::ElementsAreArray(values);
}

/**
 * gmock matcher, checks if each element of two eigen matrices are within tolerance.
 * @param other matrix to compare
 * @param rtol relative tolerance
 * @param atol absolute tolerance
 * @return matcher that accepts eigen matrix types
 */
MATCHER_P3(MatrixNear, other, rtol, atol,
           "approx. equal to " + testing::PrintToString(print_wrap(other)) + " (rtol = " + testing::PrintToString(rtol) +
               ", atol = " + testing::PrintToString(atol) + ")")
{
    epi::unused(result_listener);
    return ((arg - other).array().abs() <= (atol + rtol * other.array().abs())).all();
}

/**
 * gmock matcher, checks if each element of two eigen matrices are close.
 * @param other matrix to compare
 * @return matcher that accepts eigen matrix types
 */
MATCHER_P(MatrixNear, other,
           "approx. equal to " + testing::PrintToString(print_wrap(other)) + " (rtol = 1e-15, atol = 1e-15)")
{
    epi::unused(result_listener);
    return ((arg - other).array().abs() <= (1e-15 + 1e-15 * other.array().abs())).all();
}

#endif //EPI_TESTS_MATCHERS_H