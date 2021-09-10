#include "epidemiology/abm/world_builder.h"

#include "epidemiology/utils/eigen.h"

#include <unsupported/Eigen/NonLinearOptimization>

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
        // at every location live a certain number of people
        for (auto i = 0; i < m_l; i++) {
            f(i) = x.segment(m_n * m_l, (m_n + 1) * m_l - 1).sum() - m_y(i);
        }
        return 0;
    }

    int inputs() const
    {
        return m_n * m_l;
    }
    int values() const
    {
        return m_y.size();
    }
};

void minimize()
{
    Eigen::VectorXd y(4);
    y << 1, 2, 3, 4;
    int n = 2, l = 2, num_persons = 30;
    //size of the locations as a constraint
    y.segment(0, l - 1) = Eigen::VectorXd::Constant(l, num_persons / l);

    MyFunc func{y, n, l};

    Eigen::NumericalDiff<MyFunc> func_with_num_diff(func);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MyFunc>> lm(func_with_num_diff);
    Eigen::VectorXd x(n * l);
    lm.minimize(x);

    std::cout << x;
}

/*void create_locations(int num_locs, LocationType type, World& world)
{
    // sort people according to their age groups
}
*/

} // namespace epi
