#include "IDE/IDE.h"
#include "memilio/math/eigen.h"

#include <vector>

int main(){
    // example to demonstrate how the IDE-Model can be used 
    int tmax=15;
    int N=800000;
    double dt=1.0/10.0;
    std::vector<Eigen::Vector2d> result;
    int length_res=0;

    result.push_back(Eigen::Vector2d (-16.5,(double)N));
    length_res++;

    while(result[length_res-1][0]<0){
        result.push_back(Eigen::Vector2d (round((result[length_res-1][0]+dt)*100.0)/100.0,result[length_res-1][1]+result[length_res-1][0]/10.0));
        length_res++;
    }
    mio::IdeModel model(result, length_res, dt, N);
    model.add_damping(5.0,0.5);
    model.add_damping(7,3.5);
    model.simulate(tmax);
    model.print_result();
}