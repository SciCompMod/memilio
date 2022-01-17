#include "IDE/IDE.h"
#include <math.h>

namespace mio{

    IdeModel::Beta(double tau, double p=3, double q=10) const{
        if((timelatency<tau) && (timeinfectious+timelatency>tau)){
            return tgamma(p+q)*pow(tau-timelatency,p-1)*pow(timeinfectious+timelatency-tau,q-1)/(tgamma(p)*tgamma(q)*pow(timeinfectious,p+q-1));
        }
        return 0;
    }

    IdeModel::S_derivative(int idx) const{
        return (s[idx+1]-S[idx-1])/(2*dt);
    }

    IdeModel::num_integration_inner_integral(double time) const{
        double a=time-(timelatency+timeinfectious);
        double b=time-timelatency;

        int k=0,l=0,i=0;
        while (t[i+1]<=tx0l){
            if (t[i]<=a){
                k=i;
                l=i;
            } else if(t[i]<=b){
                l=i;
            }
            i++;
        }
        double res=0.5 * (IdeModel::Beta(time - t[k], latent_period, x1) * S_derivative(k,dt_inverse) + 
                IdeModel::Beta(time - t[l], latent_period, x1) * S_derivative(l,dt_inverse));
        k++;
        while (k<=l-1){
            res+=(IdeModel::Beta(time-t[k],latent_period,x1)*S_derivative(k,dt_inverse));
            k++;
        }
    
    return res;

    }
}