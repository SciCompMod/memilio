#ifndef IDE_H
#define IDE_H

#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/contact_matrix.h"

#include <vector>

namespace mio
{
class IdeModel{
    public:
        IdeModel(TimeSeries<double> init, double dt_init, int N_init);

        ContactMatrix& get_contact_matrix();
        void set_latencytime(double latency);
        void set_infectioustime(double infectious);

        const TimeSeries<double>& simulate(int t_max);
        const TimeSeries<double>& calculate_EIR();
        void print_result(bool calculated_SEIR=false) const;
        
    private:
        double Beta(double tau, double p=3.0, double q=10.0) const;
        double S_derivative(int idx) const;
        double num_integration_inner_integral(int idx) const;

        // vector containing one time Step per entry stored in an Eigen Vector (time, number of Susceptible at time t, R0t)
        TimeSeries<double> result; 
        TimeSeries<double> result_SEIR=TimeSeries<double>(4); 

        ContactMatrix contact_matrix=ContactMatrix(Eigen::MatrixXd::Constant(1,1,10),Eigen::MatrixXd::Constant(1,1,1));

        double latencyTime=3.3;
        double infectiousTime=8.2;

        double dt=0;
        int k=0;
        int l=0;
        int N=0;
        };
}// namespace mio
#endif