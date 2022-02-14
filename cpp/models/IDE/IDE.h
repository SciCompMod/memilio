#ifndef IDE_H
#define IDE_H

#include "memilio/math/eigen.h"

#include <vector>

namespace mio
{

class IdeModel{
    public:
        IdeModel(std::vector<Eigen::Vector2d> init, int length_init, double dt_init, int N_init);
        void set_latencytime(double latency);
        void set_infectioustime(double infectious);
        std::vector<Eigen::Vector2d> simulate(int t_max);
        void print_result() const;
        void add_damping(double time, double R0t_time);

    private:
        double Beta(double tau, double p=3.0, double q=10.0) const;
        double S_derivative(int idx) const;
        double num_integration_inner_integral(int idx) const;

        double timelatency=3.3;
        double timeinfectious=8.2;

        // vector containing one time Step per entry stored in an Eigen Vector (time, number of Susceptible at time t, R0t)
        std::vector<Eigen::Vector2d> result; 
        int length_result;

        std::vector<Eigen::Vector2d> R0t; 
        int length_R0t=0;

        double dt;
        int k;
        int l;
        int N;
        };
}// namespace mio
#endif