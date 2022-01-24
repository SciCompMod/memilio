#include "IDE/IDE.h"
#include <math.h>
#include <cmath>

namespace mio{

    IdeModel::IdeModel(std::vector<Eigen::VectorXd> init, int length_init, double dt_init)
    :result(init),length(length_init), dt(dt_init){
        l=(int) std::floor(timelatency/dt);
        k=(int) std::ceil((timeinfectious+timelatency)/dt);
    }

    void IdeModel::set_latencytime(double latency){
        timelatency=latency;
        l=(int) std::floor(timelatency/dt);
        k=(int) std::ceil((timeinfectious+timelatency)/dt);
    }

    void IdeModel::set_infectioustime(double infectious){
        timeinfectious=infectious;
        k=(int) std::ceil((timeinfectious+timelatency)/dt);
    }

    double IdeModel::Beta(double tau, double p, double q) const{
        if((timelatency<tau) && (timeinfectious+timelatency>tau)){
            return tgamma(p+q)*pow(tau-timelatency,p-1)*pow(timeinfectious+timelatency-tau,q-1)/(tgamma(p)*tgamma(q)*pow(timeinfectious,p+q-1));
        }
        return 0.0;
    }

    double IdeModel::S_derivative(int idx) const{
        return (result[idx+1][1]-result[idx-1][1])/(2*dt);
    }

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
}