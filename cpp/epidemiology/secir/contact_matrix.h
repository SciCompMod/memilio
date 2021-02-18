#ifndef EPI_SECIR_CONTACT_FREQUENCY_MATRIX_H
#define EPI_SECIR_CONTACT_FREQUENCY_MATRIX_H

#include "epidemiology/utils/eigen.h"
#include "epidemiology/utils/stl_util.h"
#include "epidemiology/secir/damping.h"

#include <vector>
#include <numeric>
#include <ostream>

namespace epi
{

/**
 * represents the coefficient wise matrix (or vector) expression B - D * M
 * where B is a baseline, M is a minimum and D is some time dependent complex damping factor.
 * Base class for e.g. time dependent contact matrices.
 * Coefficient wise expression, so B, D, M matrices must have the same shape.
 * @see Damping
 * @see ContactMatrix
 * @tparam D instance of Dampings or compatible type
 */
template<class D>
class DampingMatrixExpression
{    
public:
    using DampingsType = D;
    using Shape        = typename DampingsType::Shape;
    using Matrix       = typename Shape::Matrix;

    /**
     * construct with baseline and minimum contacts and zero dampings.
     * @param baseline matrix expression 
     * @param minimum matrix expression, must be same size as baseline
     * @tparam M, K matrix expressions compatible with Matrix type
     */
    template <class M, class K,
              class = std::enable_if_t<
                  (is_matrix_expression<std::decay_t<M>>::value && is_matrix_expression<std::decay_t<K>>::value), void>>
    DampingMatrixExpression(M&& baseline, K&& minimum)
        : m_baseline(std::forward<M>(baseline))
        , m_minimum(std::forward<K>(minimum))
        , m_dampings(Shape::get_shape_of(m_baseline))
    {
        assert(Shape::get_shape_of(m_minimum) == Shape::get_shape_of(m_baseline));
        m_dampings.finalize();
    }

    /**
     * construct with only baseline, minimum and dampings is zero.
     * @param baseline matrix expression 
     * @tparam M matrix expressions compatible with Matrix type
     */
    template <class M, class = std::enable_if_t<is_matrix_expression<std::decay_t<M>>::value, void>>
    explicit DampingMatrixExpression(M&& baseline)
        : DampingMatrixExpression(std::forward<M>(baseline), Matrix::Zero(Shape::get_shape_of(baseline).rows(),
                                                                          Shape::get_shape_of(baseline).cols()))
    {
    }

    /**
     * construct with shape.
     * baseline, minimum and dampings all zero.
     * @param shape_args shape arguments.
     * @tparam T shape arguments.
     */
    template <class... T, class = std::enable_if_t<std::is_constructible<Shape, T...>::value, void>>
    explicit DampingMatrixExpression(T... shape_args)
        : DampingMatrixExpression(Matrix::Zero(Shape(shape_args...).rows(), Shape(shape_args...).cols()))
    {
    }

    /**
     * adds a damping.
     * @see Dampings::add
     */
    template <class... T>
    void add_damping(T&&... t)
    {
        m_dampings.add(std::forward<T>(t)...);
    }

    /**
     * list of dampings.
     */
    auto get_dampings() const
    {
        return make_range(m_dampings.begin(), m_dampings.end());
    }
    auto get_dampings()
    {
        return make_range(m_dampings.begin(), m_dampings.end());
    }

    /**
     * get the baseline matrix.
     */
    const Matrix& get_baseline() const
    {
        return m_baseline;
    }
    Matrix& get_baseline()
    {
        return m_baseline;
    }

    /**
     * get the minimum matrix.
     */
    const Matrix& get_minimum() const
    {
        return m_minimum;
    }
    Matrix& get_minimum()
    {
        return m_minimum;
    }

    /**
     * dimensions of the matrix.
     */
    Shape get_shape() const
    {
        return Shape::get_shape_of(m_baseline);
    }

    /**
     * equality operators.
     */
    bool operator==(const DampingMatrixExpression& other) const
    {
        return m_baseline == other.m_baseline && m_minimum == other.m_minimum && m_dampings == other.m_dampings;
    }
    bool operator!=(const DampingMatrixExpression& other) const
    {
        return !(*this == other);
    }

    /**
     * Applies dampings to compute the real contact frequency at a point in time.
     * Uses lazy evaluation, coefficients are calculated on indexed access.
     * @param t time in the simulation
     * @return matrix expression (num_groups x num_groups)
     */
    auto get_matrix_at(SimulationTime t) const
    {
        return m_baseline - (m_dampings.get_matrix_at(t).array() * (m_baseline - m_minimum).array()).matrix();
    }
    auto get_matrix_at(double t) const
    {
        return get_matrix_at(SimulationTime(t));
    }

    /**
     * gtest printer.
     */
    friend void PrintTo(const DampingMatrixExpression& self, std::ostream* os)
    {
        *os << '\n' << self.m_baseline;
        *os << '\n' << self.m_minimum;
        PrintTo(self.m_dampings, os);
    }

private:
    Matrix m_baseline;
    Matrix m_minimum;
    DampingsType m_dampings;
};

/**
 * represents a collection of DampingMatrixExpressions that are summed up.
 * @tparam E some instance of DampingMatrixExpression or compatible type.
 */
template<class E>
class DampingMatrixExpressionGroup
{
public:
    using value_type      = E;
    using reference       = value_type&;
    using const_reference = const value_type&;
    using iterator        = typename std::vector<value_type>::iterator;
    using const_iterator  = typename std::vector<value_type>::const_iterator;
    
    using Shape           = typename value_type::Shape;
    using Matrix          = typename value_type::Matrix;
    using DampingsType    = typename value_type::DampingsType;

    /**
     * create a collection.
     * @param num_groups number of groups.
     * @param num_matrices number of matrices.
     */
    template<class... T, class = std::enable_if_t<std::is_constructible<Shape, T...>::value, int>>
    explicit DampingMatrixExpressionGroup(size_t num_matrices, T... shape_args)
        : m_matrices(num_matrices, value_type{shape_args...})
    {
        assert(num_matrices > 0);
    }

    /**
     * create a collection that contains these matrices.
     * @param il initializer list of matrices, must all be the same size.
     */
    DampingMatrixExpressionGroup(std::initializer_list<value_type> il)
        : m_matrices(il)
    {
        assert(il.size() > 0);
    }

    /**
     * access one matrix.
     */
    reference operator[](size_t i)
    {
        return m_matrices[i];
    }
    const_reference operator[](size_t i) const
    {
        return m_matrices[i];
    }

    /**
     * get the number of matrices.
     */
    size_t get_num_matrices() const
    {
        return m_matrices.size();
    }

    /**
     * get the number of groups.
     */
    Shape get_shape() const
    {
        return m_matrices[0].get_shape();
    }

    /** 
     * equality operators.
     */
    bool operator==(const DampingMatrixExpressionGroup& other) const
    {
        return m_matrices == other.m_matrices;
    }
    bool operator!=(const DampingMatrixExpressionGroup& other) const
    {
        return !(*this == other);
    }

    /**
     * add the same damping to all matrices.
     * @see ContactMatrix::add_damping
     */
    template <class... T>
    void add_damping(T&&... t)
    {
        for (auto& m : *this) {
            m.add_damping(std::forward<T>(t)...);
        }
    }

    /**
     * get the real contact frequency at a point in time.
     * sum of all contained matrices.
     * @param t point in time
     * @return matrix expression of size num_groups x num_groups
     */
    template <class T>
    auto get_matrix_at(T t) const
    {
        return Eigen::MatrixXd::NullaryExpr(
            get_shape().rows(), get_shape().cols(), [t, this](Eigen::Index i, Eigen::Index j) {
                return std::accumulate(m_matrices.begin(), m_matrices.end(), 0.0, [t, i, j](double s, auto& m) {
                    return s + m.get_matrix_at(t)(i, j);
                });
            });
    }

    /**
     * STL iterators over matrices.
     */
    iterator begin()
    {
        return m_matrices.begin();
    }
    iterator end()
    {
        return m_matrices.end();
    }
    const_iterator begin() const
    {
        return m_matrices.begin();
    }
    const_iterator end() const
    {
        return m_matrices.end();
    }

    /**
     * gtest printer.
     */
    friend void PrintTo(const DampingMatrixExpressionGroup& self, std::ostream* os)
    {
        for (auto& m : self.m_matrices) {
            PrintTo(m, os);
            *os << '\n';
        }
    }

private:
    std::vector<value_type> m_matrices;
};

/**
 * represents time dependent contact frequencies between groups.
 * consists of constant baseline and irreducible minimum contacts.
 * actual contacts are adjusted over time by dampings.
 * The effective contacts are B - D * (B - M), where B is the baseline, D are
 * combined dampings and M is the minimum.
 * The minimum is not necessarily smaller than the baseline in every entry,
 * but it reflects the state at maximum lockdown. In places where the 
 * minimum is greater than the baseline, a positive damping increases
 * instead of reduces contacts.
 * All these members are matrix valued, e.g. B_ij are the normal contacts
 * that one person in group i has with persons in group j.
 */
class ContactMatrix : public DampingMatrixExpression<SquareDampings>
{
public:
    using Base = DampingMatrixExpression<SquareDampings>;
    using Base::Base;

    /**
     * get the number of groups.
     */
    Eigen::Index get_num_groups() const
    {
        return Base::get_shape().rows();
    }
};

/**
 * represents a collection of contact frequency matrices that whose sum is the total
 * number of contacts.
 * can separate matrices of contacts in different contexts, e.g. work, leisure, etc. 
 */
class ContactMatrixGroup : public DampingMatrixExpressionGroup<ContactMatrix>
{
public:
    using Base = DampingMatrixExpressionGroup<ContactMatrix>;
    using Base::Base;
    
    /**
     * get the number of groups.
     */
    Eigen::Index get_num_groups() const
    {
        return (*this)[0].get_num_groups();
    }
};

} // namespace epi

#endif //EPI_SECIR_CONTACT_FREQUENCY_MATRIX_H