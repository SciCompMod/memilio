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
        IdeModel(std::vector<Eigen::Vector2d> init, int length_init, double dt_init, int N_init);
        void set_latencytime(double latency);
        void set_infectioustime(double infectious);
        std::vector<Eigen::Vector2d> simulate(int duration);
        void print_result() const;

    private:
        double Beta(double tau, double p=3.0, double q=10.0) const;
        double S_derivative(int idx) const;
        double num_integration_inner_integral(int idx) const;

        double timelatency=3.3;
        double timeinfectious=8.2;

        std::vector<Eigen::Vector2d> result; // vector mit Eigen 2 (t, S) und davon dann im vec bel viele einträge
        int length;

        double dt;
        int k;
        int l;
        int N;
        };
}// namespace mio
#endif