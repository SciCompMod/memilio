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
DECL_TYPESAFE(DampingLevel, int);

/**
 * integer damping type.
 */
DECL_TYPESAFE(DampingType, int);

/**
 * double simulation time.
 */
DECL_TYPESAFE(SimulationTime, double);

/**
 * represent interventions or effects that affect contact frequencies between multiple groups.
 * Dampings have a level and a type and are active from a certain point in time forward.
 * Dampings are square matrix valued, coefficient d_ij affects the contacts from group i to group j.
 */
class Damping : public std::tuple<Eigen::MatrixXd, DampingLevel, DampingType, SimulationTime>
{
public:
    using Base = std::tuple<Eigen::MatrixXd, DampingLevel, DampingType, SimulationTime>;

    /**
     * create a default Damping.
     * @param num_groups number of groups (i.e. size of matrix)
     */
    Damping(Eigen::Index num_groups = 1)
        : Base(Eigen::MatrixXd::Zero(num_groups, num_groups), DampingLevel{}, DampingType{}, SimulationTime{})
    {
    }

    /**
     * create a Damping.
     * @param m matrix of damping coefficients
     * @param level damping level
     * @param type damping type
     * @param t time at which the damping becomes active
     */
    template <class M>
    Damping(M&& m, DampingLevel level, DampingType type, SimulationTime t)
        : Base(std::forward<M>(m), level, type, t)
    {
        assert((get_coeffs().array() >= 0.).all() && (get_coeffs().array() <= 1.).all() && "damping coefficient out of range");
    }

    /**
     * create a Damping with constant coefficients.
     * @param num_groups number of groups.
     * @param d damping coefficient for all groups.
     * @param level damping level
     * @param type damping type
     * @param t time at which the damping becomes active
     */
    Damping(Eigen::Index num_groups, double d, DampingLevel level, DampingType type, SimulationTime t)
        : Damping(Eigen::MatrixXd::Constant(num_groups, num_groups, d), level, type, t)
    {
    }

    /**
     * create a Damping at level and type zero
     * @param m damping coefficients
     * @param t time at which the damping becomes active
     */
    template <class M>
    Damping(M&& m, SimulationTime t)
        : Damping(std::forward<M>(m), DampingLevel(0), DampingType(0), t)
    {
    }

    /**
     * create a Damping with constant coefficients and zero level and type.
     * @param num_groups number of groups.
     * @param d damping coefficient for all groups.
     * @param t time at which the damping becomes active
     */
    Damping(Eigen::Index num_groups, double d, SimulationTime t)
        : Damping(num_groups, d, DampingLevel(0), DampingType(0), t)
    {
    }

    SimulationTime& get_time()
    {
        return std::get<SimulationTime>(*this);
    }
    const SimulationTime& get_time() const
    {
        return std::get<SimulationTime>(*this);
    }

    DampingLevel& get_level()
    {
        return std::get<DampingLevel>(*this);
    }
    const DampingLevel& get_level() const
    {
        return std::get<DampingLevel>(*this);
    }

    DampingType& get_type()
    {
        return std::get<DampingType>(*this);
    }
    const DampingType& get_type() const
    {
        return std::get<DampingType>(*this);
    }

    const Eigen::MatrixXd& get_coeffs() const
    {
        return std::get<Eigen::MatrixXd>(*this);
    }
    Eigen::MatrixXd& get_coeffs()
    {
        return std::get<Eigen::MatrixXd>(*this);
    }
};

/**
 * collection of dampings at different time points.
 */
class Dampings
{
public:
    using value_type      = Damping;
    using reference       = value_type&;
    using const_reference = const value_type&;
    using iterator        = std::vector<value_type>::iterator;
    using const_iterator  = std::vector<value_type>::const_iterator;

    /**
     * create damping collection.
     * @param num_groups number of groups (i.e. size of damping matrices)
     * @param num_dampings number of initial elements in the collection
     */ 
    Dampings(Eigen::Index num_groups = 1, size_t num_dampings = 0)
        : m_dampings(num_dampings, Damping(num_groups))
        , m_num_groups(num_groups)
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
        m_num_groups = m_dampings.front().get_coeffs().rows();
    }

    /**
     * add a damping.
     * @param damping a Damping
     */
    void add(const Damping& damping)
    {
        add_(damping);
    }
    template<class... T>
    void add(T&&... t)
    {
        add_(Damping{std::forward<T>(t)...});
    }
    template<class...T>
    void add(double d, T... t)
    {
        add_(Damping{m_num_groups, d, std::forward<T>(t)...});
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
     * @return matrix expression (num_groups x num_groups)
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
                            std::get<Eigen::MatrixXd>(*(ub - 1)), std::get<Eigen::MatrixXd>(*ub));
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
     * access one matrix.
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

    Eigen::Index get_num_groups() const
    {
        return m_num_groups;
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
            *os << '\n' << std::get<Eigen::MatrixXd>(d);
        }
    }

private:
    /**
     * internal add.
     */
    void add_(const Damping& damping);

    /**
     * replace matrices of the same type, sum up matrices on the same level.
     * add new types/levels if necessary.
     */
    static void update_active_dampings(
        const Damping& damping,
        std::vector<std::tuple<std::reference_wrapper<const Eigen::MatrixXd>, DampingLevel, DampingType>>&
            active_by_type,
        std::vector<std::tuple<Eigen::MatrixXd, DampingLevel>>& sum_by_level);

    /**
     * e.g. inclusive_exclusive_sum({A, B, C}) = A + B + C - AB - BC - AC + ABC
     * equal to but more efficient than 1 - (1 - A)(1 - B)(1 - C))
     */
    template <class Iter>
    static void inclusive_exclusive_sum_rec(Iter b, Iter e, Eigen::MatrixXd& sum)
    {
        if (b != e) {
            sum = sum + std::get<Eigen::MatrixXd>(*b) - (sum.array() * std::get<Eigen::MatrixXd>(*b).array()).matrix();
            inclusive_exclusive_sum_rec(++b, e, sum);
        }
    }
    template <class Tuple>
    static Eigen::MatrixXd inclusive_exclusive_sum(const std::vector<Tuple>& v)
    {
        assert(!v.empty());
        auto& m  = std::get<Eigen::MatrixXd>(v.front());
        auto sum = m.eval();
        inclusive_exclusive_sum_rec(v.begin() + 1, v.end(), sum);
        return sum;
    }

private:
    std::vector<Damping> m_dampings;
    Eigen::Index m_num_groups;
    mutable std::vector<std::tuple<Eigen::MatrixXd, SimulationTime>> m_accumulated_dampings_cached;
};

} // namespace epi

#endif // DAMPING_H
