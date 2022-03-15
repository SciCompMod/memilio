#include "IDE/IDE.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"

#include <vector>
#include <iostream>

int main()
{
    using Vec = mio::TimeSeries<double>::Vector;
    // example to demonstrate how the IDE-Model can be used
    int tmax  = 15;
    int N     = 810000;
    double dt = 1.0 / 10.0;
    mio::TimeSeries<double> result(1);

    result.add_time_point<Eigen::VectorXd>(-16.5, Vec::Constant(1, (double)800000));

    while (result.get_last_time() < 0) {
        result.add_time_point(round((result.get_last_time() + dt) * 100.0) / 100.0,
                              Vec::Constant(1, (double)result.get_last_value()[0] + result.get_last_time() / 10.0));
    }

    mio::IdeModel model(std::move(result), dt, N);

    mio::ContactMatrix& contact_matrix = model.get_contact_matrix();
    contact_matrix =
        mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 1 / 8.2), Eigen::MatrixXd::Constant(1, 1, 0.5 / 8.2));
    contact_matrix.add_damping(0.7, mio::SimulationTime(3.0));

    model.simulate(tmax);
    model.calculate_EIR();
    model.print_result(true);
}