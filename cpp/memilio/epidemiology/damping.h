/* 
* Copyright (C) 2020-2024 MEmilio
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
#ifndef DAMPING_H
#define DAMPING_H

#include "memilio/math/eigen.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/type_safe.h"
#include "memilio/utils/stl_util.h"
#include "memilio/math/matrix_shape.h"
#include "memilio/math/smoother.h"
#include "memilio/math/floating_point.h"

#include <tuple>
#include <vector>
#include <algorithm>
#include <ostream>

namespace mio
{

/**
 * integer damping level.
 */
DECL_TYPESAFE(int, DampingLevel);

/**
 * integer damping type.
 */
DECL_TYPESAFE(int, DampingType);

/**
 * double simulation time.
 */
class MEMILIO_ENABLE_EBO SimulationTime : public TypeSafe<double, SimulationTime>,
                                          public OperatorAdditionSubtraction<SimulationTime>,
                                          public OperatorScalarMultiplicationDivision<SimulationTime, double>,
                                          public OperatorComparison<SimulationTime>
{
public:
    using TypeSafe<double, SimulationTime>::TypeSafe;
};

/**
 * represent interventions or effects that affect contact frequencies between multiple groups.
 * Dampings have a level and a type and are active from a certain point in time forward.
 * Dampings are square matrix valued, coefficient d_ij affects the contacts from group i to group j.
 * @tparam S Matrix shape type
 */
template <class S>
class Damping : public std::tuple<typename S::Matrix, DampingLevel, DampingType, SimulationTime>
{
public:
    using Shape  = S;
    using Matrix = typename Shape::Matrix;
    using Base   = std::tuple<Matrix, DampingLevel, DampingType, SimulationTime>;

    /**
     * create a default Damping.
     * @param shape_args arguments to construct the shape of the damping matrix (can be Shape itself, copy ctor)
     * @tparam T constructor arguments of Damping::Shape.
     */
    template <class... T, class = std::enable_if_t<std::is_constructible<Shape, T...>::value, int>>
    explicit Damping(T... shape_args)
        : Base(Matrix::Zero(Shape(shape_args...).rows(), Shape(shape_args...).cols()), {}, {}, {})
    {
    }

    /**
     * create a Damping.
     * @param m matrix of damping coefficients
     * @param level damping level
     * @param type damping type
     * @param t time at which the damping becomes active
     * @tparam ME matrix expression, must be compatible with Shape
     */
    template <class ME>
    Damping(const Eigen::MatrixBase<ME>& m, DampingLevel level, DampingType type, SimulationTime t)
        : Base(m, level, type, t)
    {
        assert((get_coeffs().array() <= 1.).all() && "damping coefficient out of range");
    }

    /**
     * create a Damping with constant coefficients.
     * @param shape_args arguments to construct the shape of the damping matrix (can be Shape itself, copy ctor)
     * @param d damping coefficient for all groups.
     * @param level damping level
     * @param type damping type
     * @param t time at which the damping becomes active
     * @tparam T Shape constructor arguments.
     */
    template <class... T, class = std::enable_if_t<std::is_constructible<Shape, T...>::value, void>>
    Damping(double d, DampingLevel level, DampingType type, SimulationTime t, T... shape_args)
        : Damping(Matrix::Constant(Shape(shape_args...).rows(), Shape(shape_args...).cols(), d), level, type, t)
    {
    }

    /**
     * create a Damping at level and type zero
     * @param m damping coefficients
     * @param t time at which the damping becomes active
     * @tparam ME matrix expression, must be compatible with Damping::Matrix
     */
    template <class ME>
    Damping(const Eigen::MatrixBase<ME>& m, SimulationTime t)
        : Damping(m, DampingLevel(0), DampingType(0), t)
    {
    }

    /**
     * create a Damping with constant coefficients and zero level and type.
     * @param shape_args arguments to construct the shape of the damping matrix (can be Shape itself, copy ctor)
     * @param d damping coefficient for all groups.
     * @param t time at which the damping becomes active
     * @tparam T Shape constructor arguments.
     */
    template <class... T, class = std::enable_if_t<std::is_constructible<Shape, T...>::value, void>>
    Damping(double d, SimulationTime t, T... shape_args)
        : Damping(d, DampingLevel(0), DampingType(0), t, shape_args...)
    {
    }

    /**
     * the time this damping becomes active.
     */
    SimulationTime& get_time()
    {
        return std::get<SimulationTime>(*this);
    }
    const SimulationTime& get_time() const
    {
        return std::get<SimulationTime>(*this);
    }

    /**
     * the level of this damping.
     */
    DampingLevel& get_level()
    {
        return std::get<DampingLevel>(*this);
    }
    const DampingLevel& get_level() const
    {
        return std::get<DampingLevel>(*this);
    }

    /**
     * the type of this damping.
     */
    DampingType& get_type()
    {
        return std::get<DampingType>(*this);
    }
    const DampingType& get_type() const
    {
        return std::get<DampingType>(*this);
    }

    /**
     * the coefficients of this damping.
     */
    const Matrix& get_coeffs() const
    {
        return std::get<Matrix>(*this);
    }
    Matrix& get_coeffs()
    {
        return std::get<Matrix>(*this);
    }

    /**
     * shape of the damping matrix.
     */
    Shape get_shape() const
    {
        return Shape::get_shape_of(get_coeffs());
    }

    /**
     * GTest printer.
     */
    friend void PrintTo(const Damping& self, std::ostream* os)
    {
        *os << '[' << std::get<SimulationTime>(self) << ',' << std::get<DampingType>(self) << ','
            << std::get<DampingLevel>(self) << ']';
        *os << '\n' << std::get<Matrix>(self);
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Damping");
        obj.add_element("Time", get_time());
        obj.add_element("Type", get_type());
        obj.add_element("Level", get_level());
        obj.add_element("Coeffs", get_coeffs());
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Damping> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Damping");
        auto ti  = obj.expect_element("Time", Tag<SimulationTime>{});
        auto ty  = obj.expect_element("Type", Tag<DampingType>{});
        auto l   = obj.expect_element("Level", Tag<DampingLevel>{});
        auto c   = obj.expect_element("Coeffs", Tag<Matrix>{});
        return apply(
            io,
            [](auto&& ti_, auto&& ty_, auto&& l_, auto&& c_) {
                return Damping(c_, l_, ty_, ti_);
            },
            ti, ty, l, c);
    }
};

/**
 * collection of dampings at different time points.
 * combination of dampings is computed in the way described at get_matrix_at.
 * @see get_matrix_at 
 * @tparam D an instance of Damping template or compatible type. 
 */
template <class D>
class Dampings
{
public:
    using Shape           = typename D::Shape;
    using Matrix          = typename Shape::Matrix;
    using value_type      = D;
    using reference       = value_type&;
    using const_reference = const value_type&;
    using iterator        = typename std::vector<value_type>::iterator;
    using const_iterator  = typename std::vector<value_type>::const_iterator;

    /**
     * create damping collection.
     * @param shape_args shape constructor arguments.
     * @param num_dampings number of initial elements in the collection
     * @tparam T Shape constructor arguments.
     */
    template <class... T, class = std::enable_if_t<std::is_constructible<Shape, T...>::value, void>>
    explicit Dampings(T... shape_args)
        : m_dampings()
        , m_shape(shape_args...)
    {
        update_cache();
    }

    /**
     * create damping collection.
     * @param il initializer list of Dampings
     */
    Dampings(std::initializer_list<value_type> il)
        : m_dampings(il)
    {
        assert(il.size() > 0);
        m_shape = il.begin()->get_shape();
        update_cache();
    }

    /**
     * add a damping.
     * @param damping a Damping
     */
    void add(const value_type& damping)
    {
        add_(damping);
    }
    template <class... T>
    void add(T&&... t)
    {
        add_(value_type(std::forward<T>(t)...));
    }
    template <class... T>
    void add(double d, T... t)
    {
        add_(value_type(d, std::forward<T>(t)..., m_shape));
    }

    /**
     * remove the damping at index i.
     */
    void remove(size_t i)
    {
        assert(m_dampings.size() > i);
        m_dampings.erase(m_dampings.begin() + i);
        automatic_cache_update();
    }

    /**
     * remove all dampings.
     */
    void clear()
    {
        m_dampings.clear();
        automatic_cache_update();
    }

    /**
    * Disable the internal cache to speed up multiple modifications in a row. 
    * This class has an internal cache where all dampings are combined into a single time series of matrices for faster lookup.
    * By default, the cache is automatically updated when dampings are added or removed, but this can be expensive.
    * This function can be used to disable the cache so that add() or remove() can be called multiple times
    * in a row more efficiently. Afterwards, the cache needs to be reenabled.
    * @param b True if cache updates are enabled.
    */
    void set_automatic_cache_update(bool b)
    {
        m_automatic_cache_update = b;
        automatic_cache_update();
    }

    /**
     * Computes the real contact frequency at a point in time.
     * Combines the dampings that are active at a point in time according to
     * type and level rules. Dampings on different levels apply "multiplicatively".
     * Dampings on the same level apply additively.
     * e.g.
     * Two dampings a and b on different levels combine to (1 - (1 - a)(1 - b)) or (a + b - ab).
     * Two dampings a and b on the same level combine to (a + b). 
     * 
     * Transitions between different contact frequencies are smoothed out over one day to avoid discontinuities.
     * Uses lazy evaluation, coefficients are calculated on indexed access.
     * @param t time in the simulation
     * @return matrix expression 
     */
    auto get_matrix_at(SimulationTime t) const
    {
        assert(!m_accumulated_dampings_cached.empty() &&
               "Cache is not current. Did you disable the automatic cache update?");
        auto ub =
            std::upper_bound(m_accumulated_dampings_cached.begin(), m_accumulated_dampings_cached.end(),
                             std::make_tuple(t), [](auto&& tup1, auto&& tup2) {
                                 return double(std::get<SimulationTime>(tup1)) < double(std::get<SimulationTime>(tup2));
                             });
        auto damping =
            smoother_cosine(double(t), double(std::get<SimulationTime>(*ub)) - 1, double(std::get<SimulationTime>(*ub)),
                            std::get<Matrix>(*(ub - 1)), std::get<Matrix>(*ub));
        return damping;
    }
    auto get_matrix_at(double t) const
    {
        return get_matrix_at(SimulationTime(t));
    }

    /**
     * access one damping in this collection.
     */
    reference operator[](size_t i)
    {
        return m_dampings[i];
    }
    const_reference operator[](size_t i) const
    {
        return m_dampings[i];
    }

    /**
     * equality operators.
     */
    bool operator==(const Dampings& other) const
    {
        return m_dampings == other.m_dampings;
    }
    bool operator!=(const Dampings& other) const
    {
        return !(*this == other);
    }

    /**
     * get the number of matrices.
     */
    size_t get_num_dampings() const
    {
        return m_dampings.size();
    }

    /**
     * dimensions of the damping matrix.
     */
    Shape get_shape() const
    {
        return m_shape;
    }

    /**
     * STL iterators over matrices.
     */
    iterator begin()
    {
        return m_dampings.begin();
    }
    iterator end()
    {
        return m_dampings.end();
    }
    const_iterator begin() const
    {
        return m_dampings.begin();
    }
    const_iterator end() const
    {
        return m_dampings.end();
    }

    /**
     * GTest printer.
     */
    friend void PrintTo(const Dampings& self, std::ostream* os)
    {
        for (auto& d : self.m_dampings) {
            *os << '\n'
                << '[' << std::get<SimulationTime>(d) << ',' << std::get<DampingType>(d) << ','
                << std::get<DampingLevel>(d) << ']';
            *os << '\n' << std::get<Matrix>(d);
        }
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Dampings");
        obj.add_element("Shape", get_shape());
        obj.add_list("Dampings", begin(), end());
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Dampings> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Dampings");
        auto s   = obj.expect_element("Shape", Tag<Shape>{});
        auto d   = obj.expect_list("Dampings", Tag<value_type>{});
        return apply(
            io,
            [](auto&& s_, auto&& d_) -> IOResult<Dampings> {
                Dampings dampings(s_);
                for (auto&& i : d_) {
                    if (i.get_shape() != s_) {
                        return failure(StatusCode::InvalidValue, "Dampings must all have the same shape.");
                    }
                    dampings.add(i);
                }
                return success(dampings);
            },
            s, d);
    }

private:
    /**
     * internal add.
     */
    void add_(const value_type& damping);

    /**
     * compute the cache of accumulated dampings.
     * if this is used after adding dampings, all subsequent calls to get_matrix_at()
     * are quick and threadsafe. Otherwise the cache is updated automatically on the first call.
     */
    void update_cache();

    /**
     * updates the cache if automatic cache updates are enabled.
     * @see update_cache(), set_automatic_cache_update()
     */
    void automatic_cache_update()
    {
        if (m_automatic_cache_update) {
            update_cache();
        }
    }

    /**
     * replace matrices of the same type, sum up matrices on the same level.
     * add new types/levels if necessary.
     */
    static void update_active_dampings(
        const value_type& damping,
        std::vector<std::tuple<std::reference_wrapper<const Matrix>, DampingLevel, DampingType>>& active_by_type,
        std::vector<std::tuple<Matrix, DampingLevel>>& sum_by_level);

    /**
     * e.g. inclusive_exclusive_sum({A, B, C}) = A + B + C - AB - BC - AC + ABC
     * equal to but more efficient than 1 - (1 - A)(1 - B)(1 - C))
     */
    template <class Iter>
    static void inclusive_exclusive_sum_rec(Iter b, Iter e, Matrix& sum)
    {
        if (b != e) {
            auto& mat_b   = std::get<Matrix>(*b);
            auto mat_prod = (sum.array() * mat_b.array()).matrix();
            sum           = sum + mat_b - mat_prod;
            inclusive_exclusive_sum_rec(++b, e, sum);
        }
    }
    template <class Tuple>
    static Matrix inclusive_exclusive_sum(const std::vector<Tuple>& v)
    {
        assert(!v.empty());
        Matrix sum = std::get<Matrix>(v.front());
        inclusive_exclusive_sum_rec(v.begin() + 1, v.end(), sum);
        return sum;
    }

private:
    std::vector<value_type> m_dampings;
    Shape m_shape;
    std::vector<std::tuple<Matrix, SimulationTime>> m_accumulated_dampings_cached;
    bool m_automatic_cache_update = true;
};

template <class D>
void Dampings<D>::update_cache()
{
    using std::get;

    if (m_accumulated_dampings_cached.empty()) {
        m_accumulated_dampings_cached.emplace_back(Matrix::Zero(m_shape.rows(), m_shape.cols()),
                                                   SimulationTime(std::numeric_limits<double>::lowest()));

        std::vector<std::tuple<std::reference_wrapper<const Matrix>, DampingLevel, DampingType>> active_by_type;
        std::vector<std::tuple<Matrix, DampingLevel>> sum_by_level;
        for (auto& damping : m_dampings) {
            //update active damping
            update_active_dampings(damping, active_by_type, sum_by_level);
            auto combined_damping = inclusive_exclusive_sum(sum_by_level);
            assert((combined_damping.array() <= 1).all() && "unexpected error, accumulated damping out of range.");
            if (floating_point_equal(double(get<SimulationTime>(damping)),
                                     double(get<SimulationTime>(m_accumulated_dampings_cached.back())), 1e-15, 1e-15)) {
                std::get<Matrix>(m_accumulated_dampings_cached.back()) = combined_damping;
            }
            else {
                m_accumulated_dampings_cached.emplace_back(combined_damping, get<SimulationTime>(damping));
            }
        }

        m_accumulated_dampings_cached.emplace_back(get<Matrix>(m_accumulated_dampings_cached.back()),
                                                   SimulationTime(std::numeric_limits<double>::max()));
    }
}

template <class D>
void Dampings<D>::add_(const value_type& damping)
{
    assert(damping.get_shape() == m_shape && "Inconsistent matrix shape.");
    insert_sorted_replace(m_dampings, damping, [](auto& tup1, auto& tup2) {
        return std::make_tuple(tup1.get_time(), int(tup1.get_type()), int(tup1.get_level())) <
               std::make_tuple(tup2.get_time(), int(tup2.get_type()), int(tup2.get_level()));
    });
    m_accumulated_dampings_cached.clear();
    if (m_automatic_cache_update) {
        update_cache();
    }
}

template <class S>
void Dampings<S>::update_active_dampings(
    const value_type& damping,
    std::vector<std::tuple<std::reference_wrapper<const Matrix>, DampingLevel, DampingType>>& active_by_type,
    std::vector<std::tuple<Matrix, DampingLevel>>& sum_by_level)
{
    using std::get;

    const int MatrixIdx = 0;

    //find active with same type and level if existent
    auto iter_active_same_type = std::find_if(active_by_type.begin(), active_by_type.end(), [&damping](auto& active) {
        return get<DampingLevel>(active) == get<DampingLevel>(damping) &&
               get<DampingType>(active) == get<DampingType>(damping);
    });
    if (iter_active_same_type != active_by_type.end()) {
        //replace active of the same type and level
        auto& active_same_type = *iter_active_same_type;
        //find active with the same level
        auto& sum_same_level = *std::find_if(sum_by_level.begin(), sum_by_level.end(), [&damping](auto& sum) {
            return get<DampingLevel>(sum) == get<DampingLevel>(damping);
        });
        //remove active with the same type and level and add new one
        get<MatrixIdx>(sum_same_level) += get<MatrixIdx>(damping) - get<MatrixIdx>(active_same_type).get();
        get<MatrixIdx>(active_same_type) = get<MatrixIdx>(damping);
    }
    else {
        //add new type to active
        active_by_type.emplace_back(get<MatrixIdx>(damping), get<DampingLevel>(damping), get<DampingType>(damping));
        //find damping with same level if existent
        auto iter_sum_same_level = std::find_if(sum_by_level.begin(), sum_by_level.end(), [&damping](auto& sum) {
            return get<DampingLevel>(sum) == get<DampingLevel>(damping);
        });
        if (iter_sum_same_level != sum_by_level.end()) {
            //add to existing level
            get<MatrixIdx>(*iter_sum_same_level) += get<MatrixIdx>(damping);
        }
        else {
            //add new level
            sum_by_level.emplace_back(get<MatrixIdx>(damping), get<DampingLevel>(damping));
        }
    }
}

/**
 * aliases for common damping specializations.
 */
using SquareDamping  = Damping<SquareMatrixShape>;
using SquareDampings = Dampings<SquareDamping>;
using VectorDamping  = Damping<ColumnVectorShape>;
using VectorDampings = Dampings<VectorDamping>;

} // namespace mio

#endif // DAMPING_H
