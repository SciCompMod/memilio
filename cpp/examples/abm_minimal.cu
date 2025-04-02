/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Khoa Nguyen
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
#define MEMILIO_WITH_CUDA TRUE
#ifdef MEMILIO_WITH_CUDA
#define R123_NO_CUDA_DEVICE_RANDOM 1
#endif

// Note: Place all includes after the defines
// #include "abm/household.h"
// #include "abm/lockdown_rules.h"
// #include "abm/model.h"
// #include "abm/common_abm_loggers.h"

#include <iostream>


#ifdef MEMILIO_WITH_CUDA
#include <cuda_runtime.h>
#include <iostream>
#include <vector>

// CUDA test functions remain unchanged...


// Define a larger test size for multi-core testing
#define CUDA_TEST_SIZE 1024

// Simple CUDA kernel to verify CUDA functionality
__global__ void testCudaKernel(int* result) 
{
    *result = 42;
}

// More complex kernel to test parallel computing on multiple cores
__global__ void testParallelKernel(float* input, float* output, int size) 
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        // Simple computation to verify parallel processing
        float value = input[idx];
        // Do some work to ensure the GPU is actually utilized
        for(int i = 0; i < 1000; i++) {
            value = sinf(value) * cosf(value) + sqrtf(fabs(value));
        }
        output[idx] = value;
    }
}

// Helper function to test if CUDA is working properly
bool testCuda(int size)
{
    int deviceCount = 0;
    cudaError_t error = cudaGetDeviceCount(&deviceCount);
    
    if (error != cudaSuccess || deviceCount == 0) {
        std::cout << "CUDA test: No CUDA devices found!" << std::endl;
        return false;
    }
    
    std::cout << "CUDA test: Found " << deviceCount << " CUDA device(s)" << std::endl;
    
    // Display information about the GPU
    cudaDeviceProp deviceProp;
    for (int device = 0; device < deviceCount; device++) {
        cudaGetDeviceProperties(&deviceProp, device);
        std::cout << "  Device " << device << ": " << deviceProp.name << std::endl;
        std::cout << "   Compute capability: " << deviceProp.major << "." << deviceProp.minor << std::endl;
        std::cout << "    Multiprocessors: " << deviceProp.multiProcessorCount << std::endl;
        std::cout << "    Max threads per block: " << deviceProp.maxThreadsPerBlock << std::endl;
    }
    
    // Basic test - single value
    int* d_result;
    int h_result = 0;
    
    // Allocate device memory
    cudaMalloc((void**)&d_result, sizeof(int));
    
    // Launch kernel
    testCudaKernel<<<1, size>>>(d_result);
    
    // Copy result back
    cudaMemcpy(&h_result, d_result, sizeof(int), cudaMemcpyDeviceToHost);
    
    // Free device memory
    cudaFree(d_result);
    
    if (h_result != 42) {
        std::cout << "CUDA test: Failed! Basic kernel did not produce expected result." << std::endl;
        return false;
    }
    
    std::cout << "CUDA test: Basic single-thread test passed." << std::endl;
    
    // Extended test - multiple cores
    std::cout << "CUDA test: Running multi-core performance test..." << std::endl;
    
    // Create input data
    std::vector<float> h_input(CUDA_TEST_SIZE);
    std::vector<float> h_output(CUDA_TEST_SIZE);
    
    // Initialize input data
    for (int i = 0; i < CUDA_TEST_SIZE; i++) {
        h_input[i] = static_cast<float>(i) * 0.01f;
    }
    
    // Allocate device memory
    float* d_input;
    float* d_output;
    cudaMalloc((void**)&d_input, CUDA_TEST_SIZE * sizeof(float));
    cudaMalloc((void**)&d_output, CUDA_TEST_SIZE * sizeof(float));
    
    // Copy input data to device
    cudaMemcpy(d_input, h_input.data(), CUDA_TEST_SIZE * sizeof(float), cudaMemcpyHostToDevice);
    
    // Create CUDA events for timing
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    
    // Record start time
    cudaEventRecord(start);
    
    // Launch parallel kernel (use 256 threads per block)
    int blockSize = 256;
    int numBlocks = (CUDA_TEST_SIZE + blockSize - 1) / blockSize;
    testParallelKernel<<<numBlocks, blockSize>>>(d_input, d_output, CUDA_TEST_SIZE);
    
    // Record end time
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    
    // Calculate elapsed time
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    
    // Copy result back
    cudaMemcpy(h_output.data(), d_output, CUDA_TEST_SIZE * sizeof(float), cudaMemcpyDeviceToHost);
    
    // Clean up
    cudaFree(d_input);
    cudaFree(d_output);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    // Verify results (just check that they're not all zeros)
    bool hasNonZeroResults = false;
    for (int i = 0; i < CUDA_TEST_SIZE; i++) {
        if (h_output[i] != 0.0f) {
            hasNonZeroResults = true;
            break;
        }
    }
    
    if (!hasNonZeroResults) {
        std::cout << "CUDA test: Failed! Parallel kernel did not produce valid results." << std::endl;
        return false;
    }
    
    std::cout << "CUDA test: Multi-core test passed." << std::endl;
    std::cout << "CUDA test: Processing time: " << milliseconds << " ms" << std::endl;
    std::cout << "CUDA test: Success! CUDA is working properly with multiple cores." << std::endl;
    
    return true;
}
#endif


int abm_minimal_main(int size)
{
    // Test CUDA if enabled
    #ifdef MEMILIO_WITH_CUDA
    std::cout << "Testing CUDA capabilities..." << std::endl;
    
    // Guard the CUDA test with proper CUDA error handling
    cudaError_t cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        std::cerr << "CUDA initialization failed: " << cudaGetErrorString(cudaStatus) << std::endl;
        std::cout << "CUDA test failed! Continuing without CUDA." << std::endl;
    }
    else {
        bool cudaWorking = testCuda(size);
        std::cout << "CUDA test " << (cudaWorking ? "passed!" : "failed!") << std::endl;
        
        // Reset the device before running non-CUDA code
        cudaDeviceReset();
    }
    #else
    std::cout << "CUDA support is not enabled." << std::endl;
    #endif

    // Run ABM simulation (this doesn't use CUDA and shouldn't cause conflicts)
    // runABMSimulation();

    return 0;
}
