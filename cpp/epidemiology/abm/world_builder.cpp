#include "epidemiology/abm/world_builder.h"

#include "epidemiology/utils/eigen.h"

#include <unsupported/Eigen/NonLinearOptimization>
#include <cmath>

namespace epi
{
struct MyFunc {
    enum
    {
        InputsAtCompileTime = Eigen::Dynamic,
        ValuesAtCompileTime = Eigen::Dynamic
    };

    using Scalar       = double;
    using InputType    = Eigen::VectorXd;
    using ValueType    = Eigen::VectorXd;
    using JacobianType = Eigen::MatrixXd;

    Eigen::VectorXd m_y;
    int m_n;
    int m_l;

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& f) const
    {
        // reshape x to get a matrix for easier notation
        Eigen::Map<const Eigen::MatrixXd> M(x.data(), m_n, m_l);
        // at every location contains a certain number of people
        // at the moment: every location has the same number of people
        f.segment(0, m_l) = M.colwise().sum().transpose() - m_y.segment(0, m_l);

        // the number of people per age group is given
        f.segment(m_l, m_n) = M.rowwise().sum() - m_y.segment(m_l, m_n);

        //the average of the contact matrices per location is given
        // go through every entry of the contact matrix
        for (int i = 0; i < m_n; i++) {
            for (int j = 0; j < m_n; j++) {
                //sum over enty (i,j) of the contact matrices of every location
                double sum = 0;
                for (int k = 0; k < m_l; k++) {
                    if (i != j) {
                        sum += M(j, k) * M(i, k);
                    }
                    else {
                        sum += (M(i, k) - 1) * M(i, k);
                    }
                }

                //compare average of local contact matrices to original contact matrix
                f(m_n + m_l + i + m_n * j) = sum - m_y(m_n + m_l + i + m_n * j) * M.row(i).sum();
            }
        }

        f.segment(m_y.size(), inputs()) = 1e2 * x.array().min(0);

        return 0;
    }

    int inputs() const
    {
        return m_n * m_l;
    }
    int values() const
    {
        return m_y.size() + inputs();
    }
};

Eigen::VectorXd find_optimal_locations(Eigen::VectorXd& num_people_sorted, int num_locs,
                                       Eigen::MatrixXd& contact_matrix)
{
    int n = contact_matrix.rows();
    int l = num_locs;

    Eigen::VectorXd y((n + 1) * n + l);

    // size of the new locations
    y.segment(0, l) = Eigen::VectorXd::Constant(l, num_people_sorted.sum() / l);
    // number of people per age group
    y.segment(l, n) = num_people_sorted;
    // entries of contact matrix (entry (i,j) goes to (j*n + i))
    Eigen::Map<Eigen::VectorXd> v(contact_matrix.data(), n * n);
    y.segment(l + n, n * n) = v;

    MyFunc func{y, n, l};

    Eigen::NumericalDiff<MyFunc> func_with_num_diff(func);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MyFunc>> lm(func_with_num_diff);
    Eigen::VectorXd x;

    //the following starting values provide a rough fit.
    x.setConstant(n * l, num_people_sorted.sum() / (n * l));

    lm.minimize(x);

    return x;
}

void create_locations(uint32_t num_locs, LocationType type, World& world, Eigen::MatrixXd& contact_matrix)
{
    uint32_t num_ages = contact_matrix.rows();
    // sort people according to their age groups
    std::vector<std::vector<uint32_t>> people_sorted(num_ages);
    uint32_t index = 0;
    for (auto& person : world.get_persons()) {
        people_sorted[size_t(person.get_age())].push_back(index);
        ++index;
    }

    // get vector with number of people per age group
    Eigen::VectorXd num_people_sorted(num_ages);
    for (size_t i = 0; i < num_ages; i++) {
        num_people_sorted(i) = double(people_sorted[i].size());
    }

    // find optimal number of people per age group for every new location by minimizing nonlinear optimization problem
    Eigen::VectorXd sol = find_optimal_locations(num_people_sorted, num_locs, contact_matrix);

    // create locations and assign them to people
    Eigen::VectorXi current_index = Eigen::VectorXi::Zero(num_ages);
    for (size_t i = 0; i < num_locs; i++) {
        auto loc = world.add_location(type);
        //add enough people of each age group
        for (size_t j = 0; j < num_ages; j++) {
            int num = (int)round(sol(j + num_ages * i));
            for (int counter = 0; counter < num && current_index(j) < num_people_sorted(j); counter++) {
                auto person = world.get_persons().begin() + people_sorted[j][current_index(j)];
                person->set_assigned_location(loc);
                current_index(j)++;
            }
        }
    }
}

} // namespace epi
