#include "IDE/IDE.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace mio{
    //TODO: R0 einbauen, austesten ob funktioniert( in examples), check if all includes necessary(array?)
    /**
    * Constructor of IDEModel
    * @param init Vector with 2d Vectors in its entries containing the initial (time, quantity of susceptible) values
    *   the time values should have a distance of dt_init; we need time values from -(k-1)*dt until 0 -> condition: length_init>= k
    *   (with k=std::ceil((timeinfectious+timelatency)/dt)) with default values: take 11.6 and you are safe
    * @param length_init length of the vector containing the inital values
    * @param dt_init size of the time step used for numerical integration
    * @param N population of the considered Region 
    */
    IdeModel::IdeModel(std::vector<Eigen::Vector2d> init, int length_init, double dt_init, int N_init)
    :result(init),length(length_init), dt(dt_init),N(N_init){
        l=(int) std::floor(timelatency/dt);
        k=(int) std::ceil((timeinfectious+timelatency)/dt);
    }

    /**
    * Change latencytime, default value: 3.3
    * @param latency The value of the new latency time
    */
    void IdeModel::set_latencytime(double latency){
        timelatency=latency;
        l=(int) std::floor(timelatency/dt);
        k=(int) std::ceil((timeinfectious+timelatency)/dt);
    }

    /**
    * Change infectioustime, default value: 8.2
    * @param infectious The value of the new infectious time
    */
    void IdeModel::set_infectioustime(double infectious){
        timeinfectious=infectious;
        k=(int) std::ceil((timeinfectious+timelatency)/dt);
    }

    /**
    * Density of the generalized beta distribution used for the function f_{beta} of the IDGL model.
    * 
    * @param tau: Function  of 
    * @param p parameter p of the generalized Beta distribution
    * @param q parameter q of the generalized Beta distribution
    */
    double IdeModel::Beta(double tau, double p, double q) const{
        if((timelatency<tau) && (timeinfectious+timelatency>tau)){
            return tgamma(p+q)*pow(tau-timelatency,p-1)*pow(timeinfectious+timelatency-tau,q-1)/(tgamma(p)*tgamma(q)*pow(timeinfectious,p+q-1));
        }
        return 0.0;
    }

    /**
    * Numerical differentiation of S using a central difference quotient.
    * @param idx Time index at which the numerical differentiation is to be performed
    * @return S'(t[idx]) numerically approximated
    */
    double IdeModel::S_derivative(int idx) const{
        return (result[idx+1][1]-result[idx-1][1])/(2*dt);
    }

    /**
    * Numerical integration of the inner integral of the integro-differential equation for the group S using a trapezoidal sum.
    * @param idx Index of the time in the inner integral 
    * @return Result of the numerical integration
    */
    double IdeModel::num_integration_inner_integral(int idx) const{
        double res=0.5 * (Beta(result[idx][0] -result[idx-k][0]) * S_derivative(k) + 
                Beta(result[idx][0]-result[idx-l][0]) * S_derivative(idx-l));
        int i=idx-k+1;
        while (i<=idx-l-2){
            res+=(Beta(result[idx][0]-result[i][0])*S_derivative(i));
            i++;
        }
        return res;
    }

    /**
    * Simulation of the course of infection
    * @param duration duration of the simulation in days
    * @return result of the simulation, stored in an vector with 2d vectors in each entry. 
    *   These 2d vectors contains the time values in the first entry and the S (susceptible) values in the second one.
    */
    std::vector<Eigen::Vector2d> IdeModel::simulate(int duration){
        while (result[length-1][0]<duration){
            result.push_back(Eigen::Vector2d (round((result[length-1][0]+dt)*10000.0)/10000.0,0));
            length++;
            result[length-1][1]=result[length-2][1] * exp(dt/(2*N)*(num_integration_inner_integral(length-1)+num_integration_inner_integral(length-2)));
        }
        return result;
    }

    /**
    * print simulation result
    */
    void IdeModel::print_result() const{
        std::cout<<"time  |  number of susceptibles"<<std::endl;
        for (int i=0;i<length;i++){
            std::cout<<result[i][0]<<"  |  "<<result[i][1]<<std::endl;
        }
        
    }
}