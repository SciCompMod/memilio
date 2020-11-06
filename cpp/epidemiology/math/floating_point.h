#ifndef EPI_MATH_FLOATING_POINT_H
#define EPI_MATH_FLOATING_POINT_H

namespace epi
{
    /**
     * maximum absolute value of two numbers.
     * @param v1 first number
     * @param v2 second number
     * @return maximum absolute value between v1 and v2
     */
    template<class T>
    T abs_max(T v1, T v2)
    {
        return std::max(std::abs(v1), std::abs(v2));
    }

    /**
     * compare two floating points for equality with tolerances.
     * @param v1 first floating point value
     * @param v2 second floating point value
     * @param abs_tol maximum allowed absolute difference, default 0. The default will not work for comparison with 0.
     * @param rel_tol maximum allowed relative difference, default numeric_limits::min.
     * @return true if v1 and v2 are within tolerance of each other.
     */
    template<class T>
    bool floating_point_equal(T v1, T v2, T abs_tol = 0, T rel_tol = std::numeric_limits<T>::min())
    {
        auto diff = std::abs(v1 - v2);
        return diff <= abs_tol || diff <= abs_max(v1, v2) * rel_tol;
    }
}

#endif