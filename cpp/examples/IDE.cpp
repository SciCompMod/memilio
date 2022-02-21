#include "IDE/IDE.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
//#include "secir/secir_result_io.h"

#include <vector>
#include <iostream>

int main(){
    using Vec = mio::TimeSeries<double>::Vector;
    // example to demonstrate how the IDE-Model can be used 
    int tmax=15;
    int N=810000;
    double dt=1.0/10.0;
    mio::TimeSeries<double> result(1);
    
    result.add_time_point<Eigen::VectorXd>(-16.5,Vec::Constant(1, (double)800000));

    while(result.get_last_time()<0){
        result.add_time_point(round((result.get_last_time()+dt)*100.0)/100.0,
            Vec::Constant(1,(double)result.get_last_value()[0]+result.get_last_time()/10.0));
    }

    mio::IdeModel model(result, dt, N);
    model.add_damping(5.0,0.5);
    model.add_damping(7,3.5);
    model.simulate(tmax);
    model.calculate_EIR(); //uto res =
    model.print_result(true);
   // mio::save_result(res, "test_result.h5");
}