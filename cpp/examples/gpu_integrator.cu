#define R123_NO_CUDA_DEVICE_RANDOM 1
// #define CUDA_TEST_SIZE 1024
#include "boost/predef/compiler/nvcc.h"
#include <cuda_runtime.h>

#include "memilio/math/integrator.h"
#include "memilio/math/adapt_rk.h"
#include "memilio/math/stepper_wrapper.h"

#include <iostream>
#include <memory>
#include <cmath>

uint64_t seed                  = 147634;
// Numerical Recipes, ranqd1
const uint64_t rand_modulus    = (uint64_t(1) << 32);
const uint64_t rand_multiplier = 1664525;
const uint64_t rand_increment  = 1013904223;

uint64_t randc() {
    seed = (rand_multiplier * seed + rand_increment) & (rand_modulus - 1);
    return seed; 
}

double uniform_rand(double min, double max) {
    auto val = randc() / ((double)rand_modulus) * (max - min) + min;
    mio::log_info("rand {}", val);
    return val;
}

Eigen::VectorXd t_offset, t_scale;
// Eigen::VectorXd amplitude;
Eigen::MatrixXd amplitude_lincomb;

void set_params(Eigen::Index size, Eigen::Index band_width, double min, double max) {
    t_offset = Eigen::VectorXd::Zero(size);
    t_scale =  Eigen::VectorXd::Zero(size);
    // amplitude = Eigen::VectorXd(size);
    amplitude_lincomb = Eigen::MatrixXd::Zero(size, size);


    for (Eigen::Index i = 0; i < size; i++) {
        t_offset[i] = uniform_rand(min, max);
        t_scale[i] = uniform_rand(min, max);
        // amplitude[i] = uniform_rand(min, max);
        
        for (Eigen::Index j = i - (band_width / 2); j < i + ((band_width + 1) / 2) && j < size; j++) {
            if (j >= 0)
                amplitude_lincomb(i,j) = uniform_rand(min, max);
        } 
    }
}


void rhs(Eigen::Ref<const Eigen::VectorXd> x, double t, Eigen::Ref<Eigen::VectorXd> dxdt) {
    // dxdt = amplitude.array() * (t * t_scale + t_offset).array().sin();
    dxdt = amplitude_lincomb.matrix() * (t * t_scale + t_offset).array().sin().matrix();
}

int main()
{   
    using namespace mio;
    set_log_level(LogLevel::off);
    
    mio::log_debug("Enter Main");
    
    const int size = 3000;

    // Guard the CUDA test with proper CUDA error handling
    // cudaError_t cudaStatus = cudaSetDevice(0);
    // if (cudaStatus != cudaSuccess) {
    //     std::cerr << "CUDA initialization failed: " << cudaGetErrorString(cudaStatus) << std::endl;
    //     std::cout << "CUDA test failed! Continuing without CUDA." << std::endl;
    // }
    // else {
    //     std::cout << "CUDA initialization succeeded." <<  std::endl;
        
    //     cudaDeviceReset();
    // } 

    // TODO: nvidia-x-markers??

    set_params(size, 7, -3.0, 3.0);
    mio::log_debug("Params Set");

    // std::cout << amplitude_lincomb << "\n";

    const double abs_tol = 1e-3, rel_tol = 1e-8, min_dt = 1e-2, max_dt = 1e+2;
    // auto core = std::make_shared<mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_cash_karp54>>(abs_tol, rel_tol, min_dt, max_dt);
    auto core = std::make_shared<mio::RKIntegratorCore<double>>(abs_tol, rel_tol, min_dt, max_dt);

    mio::log_debug("Core Set");

    OdeIntegrator<double> integrator(core);

    mio::log_debug("Integrator Set");

    TimeSeries<double> results(0, Eigen::VectorXd::Zero(size));

    const double tmax = 100 * M_PI;
    double dt = 0.1;

    mio::log_debug("Results Set");
    mio::log_debug("Integrating...");

    integrator.advance(rhs, tmax, dt, results);
    
    mio::log_debug("Integration Finished");

    if (size < 5)
        results.print_table();
    else
        std::cout << "Num time steps: " << results.get_num_time_points() << "\n"; 

    mio::log_debug("Exit Main");
    return 0;
}
