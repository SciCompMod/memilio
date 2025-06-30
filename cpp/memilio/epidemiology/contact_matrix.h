/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Martin J. Kuehn, Daniel Abele
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
#ifndef EPI_ODE_CONTACT_FREQUENCY_MATRIX_H
#define EPI_ODE_CONTACT_FREQUENCY_MATRIX_H

#include "memilio/epidemiology/damping.h"
#include "memilio/math/matrix_shape.h"
#include "memilio/utils/stl_util.h"
#include "memilio/utils/logging.h"

#include <numeric>
#include <ostream>
#include <vector>

namespace mio
{

/**
 * represents the coefficient wise matrix (or vector) expression B - D * M
 * where B is a baseline, M is a minimum and D is some time dependent complex damping factor.
 * Base class for e.g. time dependent contact matrices.
 * Coefficient wise expression, so B, D, M matrices must have the same shape.
 * @see Damping
 * @see ContactMatrix
 * @tparam FP floating point, e.g. double
 * @tparam D instance of Dampings or compatible type
 */
template <typename FP, class D>
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
    template <class M, class K>
    DampingMatrixExpression(const Eigen::MatrixBase<M>& baseline, const Eigen::MatrixBase<K>& minimum)
        : m_baseline(baseline)
        , m_minimum(minimum)
        , m_dampings(Shape::get_shape_of(m_baseline))
    {
        assert(Shape::get_shape_of(m_minimum) == Shape::get_shape_of(m_baseline));
    }

    /**
     * construct with only baseline, minimum and dampings is zero.
     * @param baseline matrix expression
     * @tparam M matrix expressions compatible with Matrix type
     */
    template <class M>
    explicit DampingMatrixExpression(const Eigen::MatrixBase<M>& baseline)
        : DampingMatrixExpression(
              baseline, Matrix::Zero(Shape::get_shape_of(baseline).rows(), Shape::get_shape_of(baseline).cols()))
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
     * remove a damping.
     * @param i index to remove.
     */
    void remove_damping(size_t i)
    {
        m_dampings.remove(i);
    }

    /**
     * remove all dampings.
     */
    void clear_dampings()
    {
        m_dampings.clear();
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
     * Enable/disable automatic cache update of the dampings.
     * @see Dampings::set_automatic_cache_update
     */
    void set_automatic_cache_update(bool b)
    {
        m_dampings.set_automatic_cache_update(b);
    }

    /**
     * Applies dampings to compute the real contact frequency at a point in time.
     * Uses lazy evaluation, coefficients are calculated on indexed access.
     * @param t time in the simulation
     * @return matrix of size num_groups x num_groups
     */
    Matrix get_matrix_at(SimulationTime<FP> t) const
    {
        Matrix damping_at_t = m_dampings.get_matrix_at(t);
        assert(Shape::get_shape_of(damping_at_t) == Shape::get_shape_of(m_baseline));

        if (damping_at_t.rows() != m_baseline.rows() || damping_at_t.cols() != m_baseline.cols()) {
            mio::log_error("DampingMatrixExpression::get_matrix_at: Damping matrix at time {} has shape ({}, {}), "
                           "expected ({}, {}). ",
                           t.get(), damping_at_t.rows(), damping_at_t.cols(), m_baseline.rows(), m_baseline.cols());
        }

        return (m_baseline - damping_at_t.cwiseProduct(m_baseline - m_minimum)).eval();
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

    /**
     * serialize this.
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("DampingMatrixExpression");
        obj.add_element("Baseline", get_baseline());
        obj.add_element("Minimum", get_minimum());
        obj.add_list("Dampings", get_dampings().begin(), get_dampings().end());
    }

protected:
    /**
     * deserialize an object of a class derived from this.
     */
    template <class IOContext, class Derived>
    static IOResult<Derived> deserialize(IOContext& io, Tag<Derived>)
    {
        auto obj = io.expect_object("DampingMatrixExpression");
        auto b   = obj.expect_element("Baseline", Tag<Matrix>{});
        auto m   = obj.expect_element("Minimum", Tag<Matrix>{});
        auto d   = obj.expect_list("Dampings", Tag<typename DampingsType::value_type>{});
        return apply(
            io,
            [](auto&& b_, auto&& m_, auto&& d_) -> IOResult<Derived> {
                if (Shape::get_shape_of(b_) != Shape::get_shape_of(m_)) {
                    return failure(StatusCode::InvalidValue, "Baseline and Minimum must have the same shape.");
                }
                auto r = Derived(b_, m_);
                for (auto&& e : d_) {
                    if (e.get_shape() != Shape::get_shape_of(b_)) {
                        return failure(StatusCode::InvalidValue, "Dampings must have the same shape as the Baseline.");
                    }
                    r.add_damping(e);
                }
                return success(r);
            },
            b, m, d);
    }

public:
    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<DampingMatrixExpression> deserialize(IOContext& io)
    {
        return deserialize(io, Tag<DampingMatrixExpression>{});
    }

private:
    Matrix m_baseline;
    Matrix m_minimum;
    DampingsType m_dampings;
};

/**
 * represents a collection of DampingMatrixExpressions that are summed up.
 * @param FP floating point, e.g. double
 * @tparam E some instance of DampingMatrixExpression or compatible type.
 */
template <typename FP, class E>
class DampingMatrixExpressionGroup
{
public:
    using value_type      = E;
    using reference       = value_type&;
    using const_reference = const value_type&;
    using iterator        = typename std::vector<value_type>::iterator;
    using const_iterator  = typename std::vector<value_type>::const_iterator;

    using Shape        = typename value_type::Shape;
    using Matrix       = typename value_type::Matrix;
    using DampingsType = typename value_type::DampingsType;

    /**
     * create a collection.
     * @param num_groups number of groups.
     * @param num_matrices number of matrices.
     */
    template <class... T, class = std::enable_if_t<std::is_constructible<Shape, T...>::value, int>>
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

    DampingMatrixExpressionGroup(const std::vector<value_type>& v)
        : m_matrices(v)
    {
        assert(v.size() > 0);
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
     * dimensions of the matrices.
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
     * remove all dampings from all matrices.
     */
    void clear_dampings()
    {
        for (auto& m : *this) {
            m.clear_dampings();
        }
    }

    /**
     * Enable/disable automatic cache update for each contained matrix.
     * @see Dampings::set_automatic_cache_update
     */
    void set_automatic_cache_update(bool b)
    {
        for (auto& m : *this) {
            m.set_automatic_cache_update(b);
        }
    }

    /**
     * get the real contact frequency at a point in time.
     * sum of all contained matrices.
     * @param t point in time
     * @return matrix of size num_groups x num_groups
     */
    Matrix get_matrix_at(SimulationTime<FP> t) const
    {
        const auto shape = get_shape();
        Matrix result    = Matrix::Zero(shape.rows(), shape.cols());
        for (const auto& m : m_matrices) {
            result += m.get_matrix_at(t);
        }
        return result;
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

    /**
     * serialize this.
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("DampingMatrixExpressionGroup");
        obj.add_list("Matrices", begin(), end());
    }

protected:
    /**
     * deserialize an object of a class derived from this.
     * @see mio::deserialize
     */
    template <class IOContext, class Derived>
    static IOResult<Derived> deserialize(IOContext& io, Tag<Derived>)
    {
        auto obj = io.expect_object("DampingMatrixExpressionGroup");
        auto m   = obj.expect_list("Matrices", Tag<value_type>{});
        return apply(
            io,
            [](auto&& m_) -> IOResult<Derived> {
                //validation
                if (m_.empty()) {
                    return failure(StatusCode::InvalidValue,
                                   "DampingMatrixExpressionGroup must have at least one matrix.");
                }
                auto shape = m_[0].get_shape();
                for (size_t i = 1; i < m_.size(); ++i) {
                    if (m_[i].get_shape() != shape) {
                        return failure(StatusCode::InvalidValue,
                                       "Elements of DampingMatrixExpressionGroup must all have the same shape.");
                    }
                }

                return success(Derived(m_));
            },
            m);
    }

public:
    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<DampingMatrixExpressionGroup> deserialize(IOContext& io)
    {
        return deserialize(io, Tag<DampingMatrixExpressionGroup>{});
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
template <typename FP>
class ContactMatrix : public DampingMatrixExpression<FP, SquareDampings<FP>>
{
public:
    using Base = DampingMatrixExpression<FP, SquareDampings<FP>>;
    using Base::Base;

    /**
     * get the number of groups.
     */
    Eigen::Index get_num_groups() const
    {
        return Base::get_shape().rows();
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<ContactMatrix> deserialize(IOContext& io)
    {
        return Base::deserialize(io, Tag<ContactMatrix>{});
    }
};

/**
 * represents a collection of contact frequency matrices that whose sum is the total
 * number of contacts.
 * can separate matrices of contacts in different contexts, e.g. work, leisure, etc.
 */
template <typename FP>
class ContactMatrixGroup : public DampingMatrixExpressionGroup<FP, ContactMatrix<FP>>
{
public:
    using Base = DampingMatrixExpressionGroup<FP, ContactMatrix<FP>>;
    using Base::Base;

    /**
     * get the number of groups.
     */
    Eigen::Index get_num_groups() const
    {
        return (*this)[0].get_num_groups();
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<ContactMatrixGroup> deserialize(IOContext& io)
    {
        return Base::deserialize(io, Tag<ContactMatrixGroup>{});
    }
};

} // namespace mio

#endif //EPI_ODE_CONTACT_FREQUENCY_MATRIX_H
