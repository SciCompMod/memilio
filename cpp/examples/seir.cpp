#include <epidemiology/seir.h>

int main()
{
    double t0 = 0;
    double tmax = 200;
    double dt = 0.1;

    seirParam<double> params = seirParam<double>();

    printSeirParams(params);

    std::vector<std::vector<double> > seir(0);

    simulate_seir(t0, tmax, dt, params, seir);

}
