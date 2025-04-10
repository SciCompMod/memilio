/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Julia Bicker, Sascha Heinz Korf, Carlotta von Gerstein
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#include <vector>
#include <cuda_runtime.h>
#include "abm/interface_cuda.h"
#include "abm/location_type.h"
#include "abm/infection_state.h"
#include "abm/time.h"
#include <cmath>
#include <curand_kernel.h>
#include <stdio.h>
#include <chrono>  // Add this for timing measurements

namespace mio {
namespace abm {

__constant__ int weekend_cutoff = 5 ;
__constant__ int event_gotimeweekend = 10 ;
__constant__ int event_gotime_weekday = 19 ;
__constant__ int event_comebacktime = 20 ;


__device__ LocationType get_buried(const GPurson& person, int t){
    auto current_loc = person.current_loc;
    if (person.infection_state == InfectionState::Dead) {
        return LocationType::Cemetery;
    }
    return current_loc;
}

__device__ LocationType return_home_when_recovered(const GPurson& person, int t){
    auto current_loc = person.current_loc;
    if ((current_loc == LocationType::Hospital || current_loc == LocationType::ICU) &&
        person.infection_state == InfectionState::Recovered) {
        return LocationType::Home;
    }
    return current_loc;
}

__device__ LocationType go_to_hospital(const GPurson& person, int t){
    auto current_loc = person.current_loc;
    if (person.infection_state == InfectionState::InfectedSevere) {
        return LocationType::Hospital;
    }
    return current_loc;
}

__device__ LocationType go_to_icu(const GPurson& person, int t){
    auto current_loc = person.current_loc;
    if (person.infection_state == InfectionState::InfectedCritical) {
        return LocationType::ICU;
    }
    return current_loc;
}

__device__ bool random_transition(curandState_t* rng_state, double dt_days, double rate){
    float u_exp = curand_uniform(rng_state);
    double v = -std::logf(u_exp) / rate;
    return v < dt_days;
}

__device__ LocationType go_to_event(const GPurson& person, curandState_t* rng_state, int t, double dt_days, double rate){
    auto current_loc = person.current_loc;
    if(current_loc == LocationType::Home){
        if(random_transition(rng_state, dt_days, rate)){
            return LocationType::SocialEvent;
        }
        // if(t%24 >= event_gotime_weekday){
        //     return LocationType::SocialEvent;
        // }
    }
    else if(current_loc == LocationType::SocialEvent){
        if(t%24 >= event_comebacktime && person.time_at_location_hours >= 2.0){
            return LocationType::Home;
        }
    }
    return current_loc;
}
    
__device__ LocationType try_mobility_rule(const GPurson& person, curandState_t* rng_state, int t, double dt_days, double rate){
    auto loc_type = get_buried(person,t);
    if(loc_type != person.current_loc){
        return loc_type;
    }
    loc_type = return_home_when_recovered(person,t);
    if(loc_type != person.current_loc){
        return loc_type;
    }
    loc_type = go_to_hospital(person,t);
    if(loc_type != person.current_loc){
        return loc_type;
    }
    loc_type = go_to_icu(person,t);
    if(loc_type != person.current_loc){
        return loc_type;
    }
    loc_type = go_to_event(person, rng_state, t, dt_days, rate);
    if(loc_type != person.current_loc){
        return loc_type;
    }
    return person.current_loc;
}

// CUDA kernel for mobility rule get_buried
__global__ void next_loc(const GPurson* persons, LocationType* results, int num_persons, int t, double dt_days, unsigned long long seed, curandState_t* states, double rate) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num_persons) {
        curand_init(seed, idx, 0 &states[idx]);
        results[persons[idx].id] = try_mobility_rule(persons[idx], states[idx], t, dt_days, rate);
    }
}


// Helper function to measure elapsed time
double elapsedMilliseconds(const std::chrono::high_resolution_clock::time_point& start) {
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

// CUDA implementation for LogTimeAtLocationForEachPerson
std::vector<LocationType> mobility_rules(const std::vector<GPurson>& gPursons, int num_persons, int t, double dt_days, double rate, unsigned long long seed) 
{
    // Start timing
    auto start_total = std::chrono::high_resolution_clock::now();
    
    // Create results vector
    std::vector<LocationType> next_locs(num_persons, LocationType(0));
    
    if (num_persons == 0 || gPursons.empty()) {
        return next_locs;
    }
    
    printf("CUDA Performance Report:\n");
    printf("Number of persons: %d\n", num_persons);
    
    // Get device properties to verify CUDA is working
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    printf("Using GPU: %s with %d multiprocessors\n", prop.name, prop.multiProcessorCount);
    
    // Start memory allocation timing
    auto start_alloc = std::chrono::high_resolution_clock::now();
    
    // Allocate device memory
    GPurson* d_persons = nullptr;
    LocationType* d_results = nullptr;
    curandState_t* d_rng_states = nullptr;
    
    cudaError_t err = cudaMalloc(&d_persons, num_persons * sizeof(GPurson));
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (cudaMalloc persons): %s\n", cudaGetErrorString(err));
        return next_locs;
    }
    
    err = cudaMalloc(&d_results, num_persons * sizeof(LocationType));
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (cudaMalloc results): %s\n", cudaGetErrorString(err));
        cudaFree(d_persons);
        return next_locs;
    }

    err = cudaMalloc(&d_rng_states, num_persons * sizeof(curandState_t));
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (cudaMalloc rng_states): %s\n", cudaGetErrorString(err));
        cudaFree(d_rng_states);
        return next_locs;
    }
    
    // Initialize results to zero
    err = cudaMemset(d_results, 0, num_persons * sizeof(LocationType));
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (cudaMemset): %s\n", cudaGetErrorString(err));
        cudaFree(d_persons);
        cudaFree(d_results);
        return next_locs;
    }
    
    printf("Memory allocation time: %.3f ms\n", elapsedMilliseconds(start_alloc));
    
    // Start memory copy timing
    auto start_copy = std::chrono::high_resolution_clock::now();
    
    // Copy data to device
    err = cudaMemcpy(d_persons, gPursons.data(), gPursons.size() * sizeof(GPurson), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (cudaMemcpy to device): %s\n", cudaGetErrorString(err));
        cudaFree(d_persons);
        cudaFree(d_results);
        return next_locs;
    }
    
    printf("H2D memory copy time: %.3f ms\n", elapsedMilliseconds(start_copy));
    
    // Start kernel timing
    auto start_kernel = std::chrono::high_resolution_clock::now();
    
    // Launch kernel
    int blockSize = 256;
    int numBlocks = (num_persons + blockSize - 1) / blockSize;
    printf("Launching kernel with %d blocks of %d threads\n", numBlocks, blockSize);
    
    next_loc<<<numBlocks, blockSize>>>(d_persons, d_results, num_persons, t, dt_days, seed, d_rng_states, rate);
    
    // Check for kernel launch errors
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (kernel launch): %s\n", cudaGetErrorString(err));
        cudaFree(d_persons);
        cudaFree(d_results);
        return next_locs;
    }
    
    // Wait for kernel to finish
    err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (device sync): %s\n", cudaGetErrorString(err));
        cudaFree(d_persons);
        cudaFree(d_results);
        return next_locs;
    }
    
    printf("Kernel execution time: %.3f ms\n", elapsedMilliseconds(start_kernel));
    
    // Start copy back timing
    auto start_copyback = std::chrono::high_resolution_clock::now();
    
    // Copy results back to host
    err = cudaMemcpy(next_locs.data(), d_results, num_persons * sizeof(LocationType), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (cudaMemcpy to host): %s\n", cudaGetErrorString(err));
    }
    
    printf("D2H memory copy time: %.3f ms\n", elapsedMilliseconds(start_copyback));
    
    // Start cleanup timing
    auto start_cleanup = std::chrono::high_resolution_clock::now();
    
    // Free device memory
    cudaFree(d_persons);
    cudaFree(d_results);
    
    printf("Cleanup time: %.3f ms\n", elapsedMilliseconds(start_cleanup));
    printf("Total CUDA time: %.3f ms\n", elapsedMilliseconds(start_total));
    
    // // For comparison, let's also measure the CPU version
    // auto start_cpu = std::chrono::high_resolution_clock::now();
    // std::vector<int> cpu_results(num_persons, 0.0);
    
    // // Sequential CPU version
    // for (const auto& person : gPursons) {
    //     cpu_results[person.id] = person.time_at_location_hours;
    // }
    
    // printf("CPU sequential time: %.3f ms\n", elapsedMilliseconds(start_cpu));
    // printf("Speedup: %.2fx\n", elapsedMilliseconds(start_cpu) / elapsedMilliseconds(start_total));
    
    // // Verify results match between CPU and GPU
    // bool results_match = true;
    // for (int i = 0; i < num_persons; i++) {
    //     if (std::abs(tal[i] - cpu_results[i]) > 1e-6) {
    //         results_match = false;
    //         break;
    //     }
    // }
    // printf("Results match between CPU and GPU: %s\n", results_match ? "YES" : "NO");
    
    return next_locs;
}

}} // namespace mio::abm


