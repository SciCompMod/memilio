#ifndef IDE_H
#define IDE_H

#include "memilio/math/eigen.h"

#include <vector>
#include <array>
#include <numeric>

namespace mio
{

class IdeModel{
    public:
        IdeModel();
    private:
        double Beta(double tau, int p=3, int q=10) const;
        double S_derivative(int idx) const;
        double num_integration_inner_integral(double time) const;

        double timelatency=3.3;
        double timestepinfectious=8.2;
        double dt;
        std::vector<Eigen::VectorXd> result;
        
        };
}// namespace mio
#endif