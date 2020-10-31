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
 * represents time dependent contact frequencies between groups.
 * consists of constant baseline and irreducible minimum contacts.
 * actual contacts are adjusted over time by dampings.
 * The effective contacts are B - D * (B - M), where B is the baseline, D are
 * combined dampings and M is the minimum.
 * All these members are matrix valued, e.g. B_ij are the normal contacts
 * that one person in group i has with persons in group j.
 */
class ContactMatrix
{
public:
    /**
     * construct with baseline and minimum contacts.
     * @param baseline square matrix expression 
     * @param minimum square matrix expression, must be same size as baseline
     */
    template <class M, class K>
    ContactMatrix(M&& baseline, K&& minimum)
        : m_baseline(baseline)
        , m_minimum(minimum)
        , m_dampings(get_num_groups())
    {
        assert(m_baseline.rows() == m_baseline.cols() && m_baseline.rows() > 0);
        assert(m_minimum.rows() == m_baseline.rows() && m_minimum.cols() == m_baseline.cols());
        m_dampings.finalize();
    }
    /**
     * construct contacts with baseline, minimum is set to zero.
     * @param baseline Matrix expression.
     */
    template <class M,
              std::enable_if_t<std::is_base_of<Eigen::MatrixBase<std::decay_t<M>>, std::decay_t<M>>::value, int> = 0>
    explicit ContactMatrix(M&& baseline)
        : ContactMatrix(baseline, Eigen::MatrixXd::Zero(baseline.rows(), baseline.cols()))
    {
    }
    /**
     * constructor with zero baseline and minimum.
     * @param num_groups number of groups that have contacts with each other.
     */
    explicit ContactMatrix(Eigen::Index num_groups = 1)
        : ContactMatrix(Eigen::MatrixXd::Zero(num_groups, num_groups))
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
    const Eigen::MatrixXd& get_baseline() const
    {
        return m_baseline;
    }
    Eigen::MatrixXd& get_baseline()
    {
        return m_baseline;
    }

    /**
     * get the minimum matrix.
     */
    const Eigen::MatrixXd& get_minimum() const
    {
        return m_minimum;
    }
    Eigen::MatrixXd& get_minimum()
    {
        return m_minimum;
    }

    /**
     * get the number of groups.
     */
    Eigen::Index get_num_groups() const
    {
        return m_baseline.rows();
    }

    /**
     * equality operators.
     */
    bool operator==(const ContactMatrix& other) const
    {
        return m_baseline == other.m_baseline && m_minimum == other.m_minimum && m_dampings == other.m_dampings;
    }
    bool operator!=(const ContactMatrix& other) const
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
    friend void PrintTo(const ContactMatrix& self, std::ostream* os)
    {
        *os << '\n' << self.m_baseline;
        *os << '\n' << self.m_minimum;
        PrintTo(self.m_dampings, os);
    }

private:
    Eigen::MatrixXd m_baseline;
    Eigen::MatrixXd m_minimum;
    Dampings m_dampings;
};

/**
 * represents a collection of contact frequency matrices that whose sum is the total
 * number of contacts.
 * can separate matrices of contacts in different contexts, e.g. work, leisure, etc. 
 */
class ContactMatrixGroup
{
public:
    using value_type      = ContactMatrix;
    using reference       = value_type&;
    using const_reference = const value_type&;
    using iterator        = std::vector<value_type>::iterator;
    using const_iterator  = std::vector<value_type>::const_iterator;

    /**
     * create a collection.
     * @param num_groups number of groups.
     * @param num_matrices number of matrices.
     */
    ContactMatrixGroup(Eigen::Index num_groups = 1, size_t num_matrices = 1)
        : m_matrices(num_matrices, ContactMatrix(num_groups))
    {
        assert(num_matrices > 0);
    }

    /**
     * create a collection that contains these matrices.
     * @param il initializer list of matrices, must all be the same size.
     */
    ContactMatrixGroup(std::initializer_list<ContactMatrix> il)
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
    Eigen::Index get_num_groups() const
    {
        return m_matrices[0].get_num_groups();
    }

    /** 
     * equality operators.
     */
    bool operator==(const ContactMatrixGroup& other) const
    {
        return m_matrices == other.m_matrices;
    }
    bool operator!=(const ContactMatrixGroup& other) const
    {
        return !(*this == other);
    }

    /**
     * add the damping to all matrices.
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
            get_num_groups(), get_num_groups(), [t, this](Eigen::Index i, Eigen::Index j) {
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
    friend void PrintTo(const ContactMatrixGroup& self, std::ostream* os)
    {
        for (auto& m : self.m_matrices) {
            PrintTo(m, os);
            *os << '\n';
        }
    }

private:
    std::vector<ContactMatrix> m_matrices;
};

} // namespace epi

#endif //EPI_SECIR_CONTACT_FREQUENCY_MATRIX_H