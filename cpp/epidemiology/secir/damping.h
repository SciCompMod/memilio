#ifndef DAMPING_H
#define DAMPING_H

#include "epidemiology/utils/eigen.h"
#include "epidemiology/utils/type_safe.h"
#include "epidemiology/utils/stl_util.h"
#include "epidemiology/math/smoother.h"

#include <tuple>
#include <vector>
#include <algorithm>
#include <ostream>

namespace epi
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
class SimulationTime : public TypeSafe<double, SimulationTime>,
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
 * @tparam M Matrix type
 */
template<class M>
class Damping : public std::tuple<M, DampingLevel, DampingType, SimulationTime>
{
public:
    using Matrix = M;
    using Base = std::tuple<Matrix, DampingLevel, DampingType, SimulationTime>;

    /**
     * create a default Damping.
     * @param rows number of rows of the damping matrix
     * @param cols number of cols of the damping matrix
     */
    Damping(Eigen::Index rows = 1, Eigen::Index cols = 1)
        : Base(Matrix::Zero(rows, cols), DampingLevel{}, DampingType{}, SimulationTime{})
    {
    }

    /**
     * create a Damping.
     * @param m matrix of damping coefficients
     * @param level damping level
     * @param type damping type
     * @param t time at which the damping becomes active
     * @tparam ME matrix expression, must be compatible with Damping::Matrix
     */
    template <class ME>
    Damping(ME&& m, DampingLevel level, DampingType type, SimulationTime t)
        : Base(std::forward<ME>(m), level, type, t)
    {
        assert((get_coeffs().array() >= 0.).all() && (get_coeffs().array() <= 1.).all() && "damping coefficient out of range");
    }

    /**
     * create a Damping with constant coefficients.
     * @param rows number of rows of the damping matrix.
     * @param cols number of cols of the damping matrix.
     * @param d damping coefficient for all groups.
     * @param level damping level
     * @param type damping type
     * @param t time at which the damping becomes active
     */
    Damping(Eigen::Index rows, Eigen::Index cols, double d, DampingLevel level, DampingType type, SimulationTime t)
        : Damping(Matrix::Constant(rows, cols, d), level, type, t)
    {
    }

    /**
     * create a Damping at level and type zero
     * @param m damping coefficients
     * @param t time at which the damping becomes active
     * @tparam ME matrix expression, must be compatible with Damping::Matrix
     */
    template <class ME>
    Damping(ME&& m, SimulationTime t)
        : Damping(std::forward<ME>(m), DampingLevel(0), DampingType(0), t)
    {
    }

    /**
     * create a Damping with constant coefficients and zero level and type.
     * @param rows number of rows of the damping matrix.
     * @param cols number of cols of the damping matrix.
     * @param d damping coefficient for all groups.
     * @param t time at which the damping becomes active
     */
    Damping(Eigen::Index rows, Eigen::Index cols, double d, SimulationTime t)
        : Damping(rows, cols, d, DampingLevel(0), DampingType(0), t)
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
     * dimensions of the damping matrix.
     */
    Eigen::Index get_num_rows() const
    {
        return get_coeffs().rows();
    }
    Eigen::Index get_num_cols() const
    {
        return get_coeffs().cols();
    }
};

/**
 * collection of dampings at different time points.
 */
template<class M>
class Dampings
{
public:
    using Matrix          = M;
    using value_type      = Damping<M>;
    using reference       = value_type&;
    using const_reference = const value_type&;
    using iterator        = typename std::vector<value_type>::iterator;
    using const_iterator  = typename std::vector<value_type>::const_iterator;

    /**
     * create damping collection.
     * @param rows number of rows of the damping matrices.
     * @param cols number of cols of the damping matrices.
     * @param num_dampings number of initial elements in the collection
     */
    Dampings(Eigen::Index rows = 1, Eigen::Index cols = 1, size_t num_dampings = 0)
        : m_dampings(num_dampings, value_type(rows, cols))
        , m_num_rows(rows)
        , m_num_cols(cols)
    {
    }

    /**
     * create damping collection.
     * @param il initializer list of Dampings
     */
    Dampings(std::initializer_list<value_type> il)
        : m_dampings(il)
    {
        assert(il.size() > 0);
        m_num_rows = m_dampings.front().get_coeffs().rows();
        m_num_cols = m_dampings.front().get_coeffs().rows();
    }

    /**
     * add a damping.
     * @param damping a Damping
     */
    void add(const Damping<M>& damping)
    {
        add_(damping);
    }
    template<class... T>
    void add(T&&... t)
    {
        add_(Damping<M>{std::forward<T>(t)...});
    }
    template<class...T>
    void add(double d, T... t)
    {
        add_(Damping<M>{m_num_rows, m_num_cols, d, std::forward<T>(t)...});
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
        finalize();
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
     * compute the cache of accumulated dampings.
     * if this is used after adding dampings, all subsequent calls to get_matrix_at()
     * are quick and threadsafe. Otherwise the cache is updated automatically on the first call.
     */
    void finalize() const;

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
    Eigen::Index get_num_rows() const
    {
        return m_num_rows;
    }
    Eigen::Index get_num_cols() const
    {
        return m_num_cols;
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

    friend void PrintTo(const Dampings& self, std::ostream* os)
    {
        for (auto& d : self.m_dampings) {
            *os << '\n'
                << '[' << std::get<SimulationTime>(d) << ',' << std::get<DampingType>(d) << ','
                << std::get<DampingLevel>(d) << ']';
            *os << '\n' << std::get<Matrix>(d);
        }
    }

private:
    /**
     * internal add.
     */
    void add_(const Damping<M>& damping);

    /**
     * replace matrices of the same type, sum up matrices on the same level.
     * add new types/levels if necessary.
     */
    static void update_active_dampings(
        const Damping<M>& damping,
        std::vector<std::tuple<std::reference_wrapper<const Matrix>, DampingLevel, DampingType>>&
            active_by_type,
        std::vector<std::tuple<Matrix, DampingLevel>>& sum_by_level);

    /**
     * e.g. inclusive_exclusive_sum({A, B, C}) = A + B + C - AB - BC - AC + ABC
     * equal to but more efficient than 1 - (1 - A)(1 - B)(1 - C))
     */
    template <class Iter>
    static void inclusive_exclusive_sum_rec(Iter b, Iter e, Matrix& sum)
    {
        if (b != e) {
            sum = sum + std::get<Matrix>(*b) - (sum.array() * std::get<Matrix>(*b).array()).matrix();
            inclusive_exclusive_sum_rec(++b, e, sum);
        }
    }
    template <class Tuple>
    static Matrix inclusive_exclusive_sum(const std::vector<Tuple>& v)
    {
        assert(!v.empty());
        auto& m  = std::get<Matrix>(v.front());
        auto sum = m.eval();
        inclusive_exclusive_sum_rec(v.begin() + 1, v.end(), sum);
        return sum;
    }

private:
    std::vector<Damping<M>> m_dampings;
    Eigen::Index m_num_rows;
    Eigen::Index m_num_cols;
    mutable std::vector<std::tuple<Matrix, SimulationTime>> m_accumulated_dampings_cached;
};

template<class M>
void Dampings<M>::finalize() const
{
    using std::get;

    if (m_accumulated_dampings_cached.empty()) {
        m_accumulated_dampings_cached.emplace_back(Matrix::Zero(m_num_rows, m_num_cols),
                                                   SimulationTime(std::numeric_limits<double>::lowest()));

        std::vector<std::tuple<std::reference_wrapper<const Matrix>, DampingLevel, DampingType>>
            active_by_type;
        std::vector<std::tuple<Matrix, DampingLevel>> sum_by_level;
        for (auto& damping : m_dampings) {
            update_active_dampings(damping, active_by_type, sum_by_level);
            m_accumulated_dampings_cached.emplace_back(inclusive_exclusive_sum(sum_by_level),
                                                       get<SimulationTime>(damping));
            assert((get<Matrix>(m_accumulated_dampings_cached.back()).array() <= 1).all() &&
                   (get<Matrix>(m_accumulated_dampings_cached.back()).array() >= 0).all() &&
                   "unexpected error, accumulated damping out of range.");
        }

        m_accumulated_dampings_cached.emplace_back(get<Matrix>(m_accumulated_dampings_cached.back()),
                                                   SimulationTime(std::numeric_limits<double>::max()));
    }
}

template<class M>
void Dampings<M>::add_(const Damping<M>& damping)
{
    assert(damping.get_coeffs().rows() == m_num_rows && damping.get_coeffs().cols() == m_num_cols);
    insert_sorted_replace(m_dampings, damping, [](auto& tup1, auto& tup2) {
        return double(std::get<SimulationTime>(tup1)) < double(std::get<SimulationTime>(tup2));
    });
    m_accumulated_dampings_cached.clear();
}

template<class M>
void Dampings<M>::update_active_dampings(
    const Damping<M>& damping,
    std::vector<std::tuple<std::reference_wrapper<const Matrix>, DampingLevel, DampingType>>& active_by_type,
    std::vector<std::tuple<Matrix, DampingLevel>>& sum_by_level)
{
    using std::get;

    const int MatrixIdx = 0;

    auto iter_active_same_type = std::find_if(active_by_type.begin(), active_by_type.end(), [&damping](auto& active) {
        return get<DampingLevel>(active) == get<DampingLevel>(damping) &&
               get<DampingType>(active) == get<DampingType>(damping);
    });
    if (iter_active_same_type != active_by_type.end()) {
        //replace active of the same type and level
        auto& active_same_type = *iter_active_same_type;
        auto& sum_same_level   = *std::find_if(sum_by_level.begin(), sum_by_level.end(), [&damping](auto& sum) {
            return get<DampingLevel>(sum) == get<DampingLevel>(damping);
        });
        get<MatrixIdx>(sum_same_level) += get<MatrixIdx>(damping) - get<MatrixIdx>(active_same_type).get();
        get<MatrixIdx>(active_same_type) = get<MatrixIdx>(damping);
    }
    else {
        //add new type
        active_by_type.emplace_back(get<MatrixIdx>(damping), get<DampingLevel>(damping), get<DampingType>(damping));

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
 * specialization of Damping with square matrix damping coefficients.
 */
class SquareDamping : public Damping<Eigen::MatrixXd>
{
public:
    using Base = Damping<Eigen::MatrixXd>;

    /**
     * create a default Damping.
     * @param num_groups number of groups (i.e. size of matrix)
     */
    SquareDamping(Eigen::Index num_groups = 1)
        : Base(num_groups, num_groups)
    {
    }

    /**
     * create a Damping.
     * @param m square matrix of damping coefficients
     * @param level damping level
     * @param type damping type
     * @param t time at which the damping becomes active
     * @tparam ME matrix expression
     */
    template <class ME>
    SquareDamping(ME&& m, DampingLevel level, DampingType type, SimulationTime t)
        : Base(std::forward<ME>(m), level, type, t)
    {
    }

    /**
     * create a Damping with constant coefficients.
     * @param num_groups number of groups.
     * @param d damping coefficient for all groups.
     * @param level damping level
     * @param type damping type
     * @param t time at which the damping becomes active
     */
    SquareDamping(Eigen::Index num_groups, double d, DampingLevel level, DampingType type, SimulationTime t)
        : Base(num_groups, num_groups, d, level, type, t)
    {
    }

    /**
     * create a Damping at level and type zero
     * @param m square damping coefficients
     * @param t time at which the damping becomes active
     * @tparam ME matrix expression
     */
    template <class ME>
    SquareDamping(ME&& m, SimulationTime t)
        : Base(std::forward<ME>(m), t)
    {
    }

    /**
     * create a Damping with constant coefficients and zero level and type.
     * @param num_groups number of groups.
     * @param d damping coefficient for all groups.
     * @param t time at which the damping becomes active
     */
    SquareDamping(Eigen::Index num_groups, double d, SimulationTime t)
        : Base(num_groups, num_groups, d, t)
    {
    }

    /**
     * dimensions of the damping matrix.
     */
    Eigen::Index get_num_groups() const 
    {
        return get_num_rows();
    }
};

/**
 * specialization of Dampings with square matrix damping coefficients.
 */
class SquareDampings : public Dampings<Eigen::MatrixXd>
{
public:
    using Matrix = Eigen::MatrixXd;
    using Base = Dampings<Matrix>;

    /**
     * create damping collection.
     * @param num_groups number of groups (i.e. size of damping matrices)
     * @param num_dampings number of initial elements in the collection
     */ 
    SquareDampings(Eigen::Index num_groups = 1, size_t num_dampings = 0)
        : Base(num_groups, num_groups, num_dampings)
    {
    }

    /**
     * create damping collection.
     * @param il initializer list of Dampings
     */
    SquareDampings(std::initializer_list<Damping<Matrix>> il)
        : Base(il)
    {
    }

    /**
     * dimensions of the damping matrix.
     */
    Eigen::Index get_num_groups() const
    {
        return get_num_rows();
    }
};

/**
 * specialization of Damping with vector damping coefficients.
 */
class VectorDamping : public Damping<Eigen::VectorXd>
{
public:
    using Base = Damping<Eigen::VectorXd>;

    /**
     * create a default Damping.
     * @param num_coeffs number of groups (i.e. size of vector)
     */
    VectorDamping(Eigen::Index num_coeffs = 1)
        : Base(num_coeffs, 1)
    {
    }

    /**
     * create a Damping.
     * @param v vector of damping coefficients
     * @param level damping level
     * @param type damping type
     * @param t time at which the damping becomes active
     * @tparam VE vector expression
     */
    template <class VE>
    VectorDamping(VE&& v, DampingLevel level, DampingType type, SimulationTime t)
        : Base(std::forward<VE>(), level, type, t)
    {
    }

    /**
     * create a Damping with constant coefficients.
     * @param num_coeffs number of coefficients.
     * @param d damping coefficient for all groups.
     * @param level damping level
     * @param type damping type
     * @param t time at which the damping becomes active
     */
    VectorDamping(Eigen::Index num_coeffs, double d, DampingLevel level, DampingType type, SimulationTime t)
        : Base(num_coeffs, 1, d, level, type, t)
    {
    }

    /**
     * create a Damping at level and type zero.
     * @param v vector of damping coefficients.
     * @param t time at which the damping becomes active.
     * @tparam VE matrix expression.
     */
    template <class VE>
    VectorDamping(VE&& v, SimulationTime t)
        : Base(std::forward<VE>(v), t)
    {
    }

    /**
     * create a Damping with constant coefficients and zero level and type.
     * @param num_coeffs number of coefficients.
     * @param d damping coefficient for all groups.
     * @param t time at which the damping becomes active
     */
    VectorDamping(Eigen::Index num_coeffs, double d, SimulationTime t)
        : Base(num_coeffs, 1, d, t)
    {
    }
};

/**
 * specialization of Damping with vector damping coefficients.
 */
class VectorDampings : public Dampings<Eigen::VectorXd>
{
public:
    using Matrix = Eigen::VectorXd;
    using Base = Dampings<Matrix>;

    /**
     * create damping collection.
     * @param num_coeffs number of coefficients (i.e. size of damping vector)
     * @param num_dampings number of initial elements in the collection
     */ 
    VectorDampings(Eigen::Index num_coeffs = 1, size_t num_dampings = 0)
        : Base(num_coeffs, num_coeffs, num_dampings)
    {
    }

    /**
     * create damping collection.
     * @param il initializer list of Dampings
     */
    VectorDampings(std::initializer_list<Damping<Matrix>> il)
        : Base(il)
    {
    }
};

} // namespace epi

#endif // DAMPING_H
