#include "IDE/IDE.h"
#include "memilio/epidemiology/contact_matrix.h"
#include <iostream>

namespace mio{

    using Vec = TimeSeries<double>::Vector;

    /**
    * Constructor of IDEModel
    * @param init TimeSeries with 3d Vectors in its entries containing the initial (time, quantity of susceptible, R0t) values
    *   the time values should have a distance of dt_init; we need time values from -(k-1)*dt until 0 -> condition: length_init>= k
    *   (with k=std::ceil((infectiousTime+latencyTime)/dt)) with default values: take 11.6 and you are safe
    * @param dt_init size of the time step used for numerical integration
    * @param N population of the considered Region 
    */
    IdeModel::IdeModel(TimeSeries<double> init, double dt_init, int N_init)
    :result(init), dt(dt_init),N(N_init){
        l=(int) std::floor(latencyTime/dt);
        k=(int) std::ceil((infectiousTime+latencyTime)/dt);
    }

    /**
    * Change latencytime, default value: 3.3
    * @param latency The value of the new latency time
    */
    void IdeModel::set_latencytime(double latency){
        latencyTime=latency;
        //initialize parameters k and  l used for numerical calculations in the simulation
        l=(int) std::floor(latencyTime/dt);
        k=(int) std::ceil((infectiousTime+latencyTime)/dt);
    }

    /**
    * Change infectioustime, default value: 8.2
    * @param infectious The value of the new infectious time
    */
    void IdeModel::set_infectioustime(double infectious){
        infectiousTime=infectious;
        k=(int) std::ceil((infectiousTime+latencyTime)/dt);
    }

    /**
    * Get Reference to used ContactMatrix
    */
    ContactMatrix& IdeModel::get_contact_matrix(){
        return contact_matrix;
    }

    /**
    * Density of the generalized beta distribution used for the function f_{beta} of the IDGL model.
    * 
    * @param tau Functionparameter
    * @param p parameter p of the generalized Beta distribution
    * @param q parameter q of the generalized Beta distribution
    */
    double IdeModel::Beta(double tau, double p, double q) const{
        if((latencyTime<tau) && (infectiousTime+latencyTime>tau)){
            return tgamma(p+q)*pow(tau-latencyTime,p-1)*pow(infectiousTime+latencyTime-tau,q-1)/
                (tgamma(p)*tgamma(q)*pow(infectiousTime,p+q-1));
        }
        return 0.0;
    }

    /**
    * Numerical differentiation of S using a central difference quotient.
    * @param idx Time index at which the numerical differentiation is to be performed
    * @return S'(t[idx]) numerically approximated
    */
    double IdeModel::S_derivative(int idx) const{
        return (result[idx+1][0]-result[idx-1][0])/(2*dt);
    }

    /**
    * Numerical integration of the inner integral of the integro-differential equation for the group S using a trapezoidal sum.
    * @param idx Index of the time in the inner integral 
    * @return Result of the numerical integration
    */
    double IdeModel::num_integration_inner_integral(int idx) const{
        double res=0.5 * (Beta(result.get_time(idx) - result.get_time(idx-k)) * S_derivative(k) + 
                Beta(result.get_time(idx) - result.get_time(idx-l)) * S_derivative(idx-l));
        int i=idx-k+1;
        while (i<=idx-l-2){
            res+=(Beta(result.get_time(idx) - result.get_time(i))*S_derivative(i));
            ++i;
        }
        return res;
    }

    /**
    * Simulation of the course of infection
    * @param t_max last simulation day
    * @return result of the simulation, stored in a TimeSeries with simulation time and associated number of susceptibles (S).
    */
    const TimeSeries<double>& IdeModel::simulate(int t_max){
        while (result.get_last_time()<t_max){
            result.add_time_point(round((result.get_last_time()+dt)*10000.0)/10000.0);
            Eigen::Index idx = result.get_num_time_points();

            //R0t = effective contactfrequency at time t * timeinfectious
            auto R0t1 = contact_matrix.get_matrix_at(result.get_time(idx-2))(0,0) * infectiousTime;
            auto R0t2 = contact_matrix.get_matrix_at(result.get_last_time())(0,0) * infectiousTime;

            result.get_last_value()=Vec::Constant(1,
                    result[idx-2][0] * exp( dt * (0.5 * 1/N) * ( 
                        R0t1 * num_integration_inner_integral(idx-2) + 
                        R0t2 * num_integration_inner_integral(idx-1) )));
        }
        return result;
    }

    /**
    * Calculate the distribution of the population in E,I and R based on the calculated values for S.
    * Here average values are calculated using the average latency and infectious time and not the distributions 
    * actually used by the model. (because they are not explicitly calculable)
    * @return result of the calculation stored in an TimeSeries. The TimeSeries contains the simulation time and an associated Vector 
    * with values for S, E, I and R 
    */
    const TimeSeries<double>& IdeModel::calculate_EIR(){
        Eigen::Index num_points=result.get_num_time_points();
        double S,E,I,R;
        for (int i = k; i < num_points; ++i) {
            S=result[i][0];
            E=result[i-l][0]-S;
            I=result[i-k][0]-result[i-l][0];
            R=N-S-E-I;
            result_SEIR.add_time_point(result.get_time(i),(Vec(4) << S,E,I,R).finished());
        }
        return result_SEIR;
    }

    /**
    * print simulation result
    * @param calculated_SEIR if the values for E,I and R are calculated before with function "calculate_EIR",
    * one may want to print out these values too.
    */
    void IdeModel::print_result(bool calculated_SEIR) const{
        if (calculated_SEIR){
            std::cout<<"time  |  S  |  E  |  I  |  R"<<std::endl;
            Eigen::Index num_points=result_SEIR.get_num_time_points();
            for (int i = 0; i < num_points; ++i) {
                std::cout<<result_SEIR.get_time(i)<<"  |  "<<result_SEIR[i][0]<<"  |  "<<result_SEIR[i][1]
                    <<"  |  "<<result_SEIR[i][2]<<"  |  "<<result_SEIR[i][3]<<std::endl;
            }
        }else{
            std::cout<<"time  |  number of susceptibles"<<std::endl;
            Eigen::Index num_points=result.get_num_time_points();
            for (int i = 0; i < num_points; ++i) {
                std::cout<<result.get_time(i)<<"  |  "<<result[i][0]<<std::endl;
            }
        }
    }

}// namespace mio