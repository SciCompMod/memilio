/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Sascha Korf
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
#include <stdio.h>
#include <chrono>  // Add this for timing measurements

namespace mio {
namespace abm {

// CUDA kernel to compute time at location for each person
__global__ void computeTimeAtLocationKernel(const CudaPerson* persons, double* results, int num_persons) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num_persons) {
        results[persons[idx].id] = persons[idx].time_at_location_hours;
    }
}

// Helper function to measure elapsed time
double elapsedMilliseconds(const std::chrono::high_resolution_clock::time_point& start) {
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

// CUDA implementation for LogTimeAtLocationForEachPerson
std::vector<double> logTimeAtLocationCuda(const std::vector<CudaPerson>& cuda_persons, int num_persons) 
{
    // Start timing
    auto start_total = std::chrono::high_resolution_clock::now();
    
    // Create results vector
    std::vector<double> tal(num_persons, 0.0);
    
    if (num_persons == 0 || cuda_persons.empty()) {
        return tal;
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
    CudaPerson* d_persons = nullptr;
    double* d_results = nullptr;
    
    cudaError_t err = cudaMalloc(&d_persons, num_persons * sizeof(CudaPerson));
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (cudaMalloc persons): %s\n", cudaGetErrorString(err));
        return tal;
    }
    
    err = cudaMalloc(&d_results, num_persons * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (cudaMalloc results): %s\n", cudaGetErrorString(err));
        cudaFree(d_persons);
        return tal;
    }
    
    // Initialize results to zero
    err = cudaMemset(d_results, 0, num_persons * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (cudaMemset): %s\n", cudaGetErrorString(err));
        cudaFree(d_persons);
        cudaFree(d_results);
        return tal;
    }
    
    printf("Memory allocation time: %.3f ms\n", elapsedMilliseconds(start_alloc));
    
    // Start memory copy timing
    auto start_copy = std::chrono::high_resolution_clock::now();
    
    // Copy data to device
    err = cudaMemcpy(d_persons, cuda_persons.data(), cuda_persons.size() * sizeof(CudaPerson), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (cudaMemcpy to device): %s\n", cudaGetErrorString(err));
        cudaFree(d_persons);
        cudaFree(d_results);
        return tal;
    }
    
    printf("H2D memory copy time: %.3f ms\n", elapsedMilliseconds(start_copy));
    
    // Start kernel timing
    auto start_kernel = std::chrono::high_resolution_clock::now();
    
    // Launch kernel
    int blockSize = 256;
    int numBlocks = (num_persons + blockSize - 1) / blockSize;
    printf("Launching kernel with %d blocks of %d threads\n", numBlocks, blockSize);
    
    computeTimeAtLocationKernel<<<numBlocks, blockSize>>>(d_persons, d_results, num_persons);
    
    // Check for kernel launch errors
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (kernel launch): %s\n", cudaGetErrorString(err));
        cudaFree(d_persons);
        cudaFree(d_results);
        return tal;
    }
    
    // Wait for kernel to finish
    err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error (device sync): %s\n", cudaGetErrorString(err));
        cudaFree(d_persons);
        cudaFree(d_results);
        return tal;
    }
    
    printf("Kernel execution time: %.3f ms\n", elapsedMilliseconds(start_kernel));
    
    // Start copy back timing
    auto start_copyback = std::chrono::high_resolution_clock::now();
    
    // Copy results back to host
    err = cudaMemcpy(tal.data(), d_results, num_persons * sizeof(double), cudaMemcpyDeviceToHost);
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
    
    // For comparison, let's also measure the CPU version
    auto start_cpu = std::chrono::high_resolution_clock::now();
    std::vector<double> cpu_results(num_persons, 0.0);
    
    // Sequential CPU version
    for (const auto& person : cuda_persons) {
        cpu_results[person.id] = person.time_at_location_hours;
    }
    
    printf("CPU sequential time: %.3f ms\n", elapsedMilliseconds(start_cpu));
    printf("Speedup: %.2fx\n", elapsedMilliseconds(start_cpu) / elapsedMilliseconds(start_total));
    
    // Verify results match between CPU and GPU
    bool results_match = true;
    for (int i = 0; i < num_persons; i++) {
        if (std::abs(tal[i] - cpu_results[i]) > 1e-6) {
            results_match = false;
            break;
        }
    }
    printf("Results match between CPU and GPU: %s\n", results_match ? "YES" : "NO");
    
    return tal;
}

}} // namespace mio::abm

