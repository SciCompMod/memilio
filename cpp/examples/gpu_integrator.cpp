#define R123_NO_CUDA_DEVICE_RANDOM 1
// #define CUDA_TEST_SIZE 1024

#include <cuda_runtime.h>

#include "gpu_integrator.h"

// #include "memilio/math/integrator.h"
// #include "memilio/math/adapt_rk.h"

#include <numeric>
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
    return val;
}

std::vector<double> t_offset, t_scale;
// std::vector<double> amplitude;
std::vector<std::vector<double>> amplitude_lincomb;

void set_params(size_t size, size_t band_width, double min, double max) {
    t_offset = std::vector<double>(size, 0.);
    t_scale =  std::vector<double>(size, 0.);
    // amplitude = std::vector<double>(size, 0.); (void)band_width;
    amplitude_lincomb = std::vector<std::vector<double>>(size, std::vector<double>(size, 0.));


    for (size_t i = 0; i < size; i++) {
        t_offset[i] = uniform_rand(min, max);
        t_scale[i] = uniform_rand(min, max);
        // amplitude[i] = uniform_rand(min, max);
        
        for (int j = - ((int)band_width / 2); j < (((int)band_width + 1) / 2); j++) {
            if ((int)i + j >= 0 && i + j < size)
                amplitude_lincomb[i][i + j] = uniform_rand(min, max);
        } 
    }
}


void rhs(const std::vector<double>& x, double t, std::vector<double>& dxdt) {
    // for (size_t i = 0; i < x.size(); i++) {
    //     dxdt[i] = 0;
    // }
    const size_t size = x.size();

    // dxdt[i] = amplitude[i] * std::sin(t * t_scale[i] + t_offset[i]);
    #pragma acc parallel loop present(amplitude_lincomb, t_scale, t_offset) default(present) // copy(x, t, dxdt)
    for (size_t i = 0; i < size; i++) {
        
        double sum = 0;
        #pragma acc loop reduction(+:sum)
        for (size_t j = 0; j < size; j++) {
            sum += amplitude_lincomb[i][j] * std::sin(t * t_scale[j] + t_offset[j]);
        }
        dxdt[i] = sum;
    }
}

// void rhs2(Eigen::Ref<const Eigen::VectorXd> x, double t, Eigen::Ref<Eigen::VectorXd> dxdt) {
//     for (size_t i = 0; i < x.size(); i++) {
//         // dxdt[i] = amplitude[i] * std::sin(t * t_scale[i] + t_offset[i]);
        
//         dxdt[i] = 0;
//         for (size_t j = 0; j < x.size(); j++) {
//             dxdt[i] += amplitude_lincomb[i][j] * std::sin(t * t_scale[j] + t_offset[j]);
//         }
//     }
// }



namespace mio {
    void log_debug(std::string_view s) {
        std::cout << s << "\n";
    }
}

int main()
{   

    // using namespace mio;
    // set_log_level(LogLevel::off);
    
    mio::log_debug("Enter Main");

    const int size = 10000, band_width = 2;

    // // Guard the CUDA test with proper CUDA error handling
    // // cudaError_t cudaStatus = cudaSetDevice(0);
    // // if (cudaStatus != cudaSuccess) {
    // //     std::cerr << "CUDA initialization failed: " << cudaGetErrorString(cudaStatus) << std::endl;
    // //     std::cout << "CUDA test failed! Continuing without CUDA." << std::endl;
    // // }
    // // else {
    // //     std::cout << "CUDA initialization succeeded." <<  std::endl;
        
    // //     cudaDeviceReset();
    // // } 

    set_params(size, band_width, -3.0, 3.0);
    mio::log_debug("Params Set");

    #pragma acc declare copyin(readonly: amplitude_lincomb, t_offset, t_scale)
    {}

    // print(t_offset);
    // print(t_scale);
    // print(amplitude_lincomb);
    // std::cout << "\n";

    // // std::cout << amplitude_lincomb << "\n";

    const double abs_tol = 1e-3, rel_tol = 1e-8, min_dt = 1e-2, max_dt = 1e+2;
    // // auto core = std::make_shared<mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_cash_karp54>>(abs_tol, rel_tol, min_dt, max_dt);
    // auto core = std::make_shared<mio::RKIntegratorCore<double>>(abs_tol, rel_tol, min_dt, max_dt);

    Monstrosity stepper{abs_tol, rel_tol, min_dt, max_dt, std::vector<double>(size), std::vector<double>(size), std::vector<std::vector<double>>(tableau().entries_low.size(), std::vector<double>(size)), tableau()};

    mio::log_debug("Core Set");

    // OdeIntegrator<double> integrator(core);

    mio::log_debug("Integrator Set");

    // TimeSeries<double> results(0, Eigen::VectorXd::Zero(size));

    const double tmax = 100 * M_PI;
    double dt = 0.1;

    double t = 0.0;
    std::vector<double> x(size, 0.0);
    std::vector<double> x2(size, 0.0);

    std::cout << "\n";

    mio::log_debug("Results Set");    
    mio::log_debug("Integrating...");


    while (t < tmax) {
        stepper.step(&rhs, x, t, dt, x2);
        // print(x);
        // print(x2);
        for (size_t i = 0; i < size; i++) {
            x[i] = x2[i];
            x2[i] = 0;
        }
        // std::cin.ignore();
    }

    std::cout << t << " " << dt << " " << std::accumulate(x.begin(), x.end(), 0.0) << "\n";
    // integrator.advance(rhs, tmax, dt, results);
    
    mio::log_debug("Integration Finished");

    // if (size < 5)
    //     results.print_table();
    // else
    //     std::cout << "Num time steps: " << results.get_num_time_points() << "\n"; 

    mio::log_debug("Exit Main");
    // return 0;
}
