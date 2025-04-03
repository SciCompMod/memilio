#define R123_NO_CUDA_DEVICE_RANDOM 1
// #define CUDA_TEST_SIZE 1024
#include "boost/predef/compiler/nvcc.h"
#include <cuda_runtime.h>

#include "memilio/math/integrator.h"
#include "memilio/math/adapt_rk.h"

#include <iostream>
#include <memory>
#include <cmath>

Eigen::VectorX<double> t_offset, t_scale, amplitude;

void set_params(Eigen::Index size) {
    t_offset = Eigen::VectorX<double>::Random(size);
    t_scale = Eigen::VectorX<double>::Random(size);
    amplitude = Eigen::VectorX<double>::Random(size);
}


void rhs(Eigen::Ref<const Eigen::VectorX<double>> x, double t, Eigen::Ref<Eigen::VectorX<double>> dxdt) {
    dxdt = amplitude.array() * (t * t_scale + t_offset).array().sin();
}

int main()
{   
    using namespace mio;
    
    const int size = 10;

    // Guard the CUDA test with proper CUDA error handling
    cudaError_t cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        std::cerr << "CUDA initialization failed: " << cudaGetErrorString(cudaStatus) << std::endl;
        std::cout << "CUDA test failed! Continuing without CUDA." << std::endl;
    }
    else {
        std::cout << "CUDA initialization succeeded." <<  std::endl;
        
        cudaDeviceReset();
    } 

    // TODO: nvidia-x-markers??

    set_params(size);

    auto core = std::make_shared<RKIntegratorCore<double>>(1e-3, 1e-8, 0.001, 1.);
    OdeIntegrator<double> integrator(core);

    Eigen::VectorX<double> init(size);
    init.fill(1.0); 
    TimeSeries<double> results(0, init);
    results.print_table();

    const double tmax = 6.29;
    double dt = 0.1;

    integrator.advance(rhs, tmax, dt, results);

    results.print_table();

    return 0;
}
