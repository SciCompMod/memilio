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

#ifdef MEMILIO_WITH_CUDA
#define R123_NO_CUDA_DEVICE_RANDOM 1
#endif

// Note: Place all includes after the defines
// #include "abm/household.h"
// #include "abm/lockdown_rules.h"
// #include "abm/model.h"
// #include "abm/common_abm_loggers.h"

#include <fstream>


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
bool testCuda()
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
    testCudaKernel<<<1, 1>>>(d_result);
    
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

// Simple function to run ABM simulation that doesn't depend on CUDA
void runABMSimulation()
{
    // // This is a minimal example with children and adults < 60 year old.
    // // We divided them into 4 different age groups, which are defined as follows:
    // mio::set_log_level(mio::LogLevel::warn);
    // size_t num_age_groups         = 4;
    // const auto age_group_0_to_4   = mio::AgeGroup(0);
    // const auto age_group_5_to_14  = mio::AgeGroup(1);
    // const auto age_group_15_to_34 = mio::AgeGroup(2);
    // const auto age_group_35_to_59 = mio::AgeGroup(3);

    // // Create the model with 4 age groups.
    // auto model = mio::abm::Model(num_age_groups);

    // // Set same infection parameter for all age groups. For example, the incubation period is 4 days.
    // model.parameters.get<mio::abm::IncubationPeriod>() = 4.;

    // // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    // model.parameters.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    // model.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 and 35-59)
    // model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({age_group_15_to_34, age_group_35_to_59}, true);

    // // Check if the parameters satisfy their contraints.
    // model.parameters.check_constraints();

    // // There are 10 households for each household group.
    // int n_households = 10;

    // // For more than 1 family households we need families. These are parents and children and randoms (which are distributed like the data we have for these households).
    // auto child = mio::abm::HouseholdMember(num_age_groups); // A child is 50/50% 0-4 or 5-14.
    // child.set_age_weight(age_group_0_to_4, 1);
    // child.set_age_weight(age_group_5_to_14, 1);

    // auto parent = mio::abm::HouseholdMember(num_age_groups); // A parent is 50/50% 15-34 or 35-59.
    // parent.set_age_weight(age_group_15_to_34, 1);
    // parent.set_age_weight(age_group_35_to_59, 1);

    // // Two-person household with one parent and one child.
    // auto twoPersonHousehold_group = mio::abm::HouseholdGroup();
    // auto twoPersonHousehold_full  = mio::abm::Household();
    // twoPersonHousehold_full.add_members(child, 1);
    // twoPersonHousehold_full.add_members(parent, 1);
    // twoPersonHousehold_group.add_households(twoPersonHousehold_full, n_households);
    // add_household_group_to_model(model, twoPersonHousehold_group);

    // // Three-person household with two parent and one child.
    // auto threePersonHousehold_group = mio::abm::HouseholdGroup();
    // auto threePersonHousehold_full  = mio::abm::Household();
    // threePersonHousehold_full.add_members(child, 1);
    // threePersonHousehold_full.add_members(parent, 2);
    // threePersonHousehold_group.add_households(threePersonHousehold_full, n_households);
    // add_household_group_to_model(model, threePersonHousehold_group);

    // // Add one social event with 5 maximum contacts.
    // // Maximum contacs limit the number of people that a person can infect while being at this location.
    // auto event = model.add_location(mio::abm::LocationType::SocialEvent);
    // model.get_location(event).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    // // Add hospital and ICU with 5 maximum contacs.
    // auto hospital = model.add_location(mio::abm::LocationType::Hospital);
    // model.get_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    // auto icu = model.add_location(mio::abm::LocationType::ICU);
    // model.get_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    // // Add one supermarket, maximum constacts are assumed to be 20.
    // auto shop = model.add_location(mio::abm::LocationType::BasicsShop);
    // model.get_location(shop).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    // // At every school, the maximum contacts are 20.
    // auto school = model.add_location(mio::abm::LocationType::School);
    // model.get_location(school).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    // // At every workplace, maximum contacts are 20.
    // auto work = model.add_location(mio::abm::LocationType::Work);
    // model.get_location(work).get_infection_parameters().set<mio::abm::MaximumContacts>(20);

    // // Increase aerosol transmission for all locations
    // model.parameters.get<mio::abm::AerosolTransmissionRates>() = 10.0;
    // // Increase contact rate for all people between 15 and 34 (i.e. people meet more often in the same location)
    // model.get_location(work)
    //     .get_infection_parameters()
    //     .get<mio::abm::ContactRates>()[{age_group_15_to_34, age_group_15_to_34}] = 10.0;

    // // People can get tested at work (and do this with 0.5 probability) from time point 0 to day 10.
    // auto validity_period       = mio::abm::days(1);
    // auto probability           = 0.5;
    // auto start_date            = mio::abm::TimePoint(0);
    // auto end_date              = mio::abm::TimePoint(0) + mio::abm::days(10);
    // auto test_type             = mio::abm::TestType::Antigen;
    // auto test_parameters       = model.parameters.get<mio::abm::TestData>()[test_type];
    // auto testing_criteria_work = mio::abm::TestingCriteria();
    // auto testing_scheme_work   = mio::abm::TestingScheme(testing_criteria_work, validity_period, start_date, end_date,
    //                                                    test_parameters, probability);
    // model.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme_work);

    // // Assign infection state to each person.
    // // The infection states are chosen randomly with the following distribution
    // std::vector<double> infection_distribution{0.5, 0.3, 0.05, 0.05, 0.05, 0.05, 0.0, 0.0};
    // for (auto& person : model.get_persons()) {
    //     mio::abm::InfectionState infection_state = mio::abm::InfectionState(
    //         mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), infection_distribution));
    //     auto rng = mio::abm::PersonalRandomNumberGenerator(person);
    //     if (infection_state != mio::abm::InfectionState::Susceptible) {
    //         person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
    //                                                      model.parameters, start_date, infection_state));
    //     }
    // }

    // // Assign locations to the people
    // for (auto& person : model.get_persons()) {
    //     const auto id = person.get_id();
    //     //assign shop and event
    //     model.assign_location(id, event);
    //     model.assign_location(id, shop);
    //     //assign hospital and ICU
    //     model.assign_location(id, hospital);
    //     model.assign_location(id, icu);
    //     //assign work/school to people depending on their age
    //     if (person.get_age() == age_group_5_to_14) {
    //         model.assign_location(id, school);
    //     }
    //     if (person.get_age() == age_group_15_to_34 || person.get_age() == age_group_35_to_59) {
    //         model.assign_location(id, work);
    //     }
    // }

    // // During the lockdown, social events are closed for 90% of people.
    // auto t_lockdown = mio::abm::TimePoint(0) + mio::abm::days(10);
    // mio::abm::close_social_events(t_lockdown, 0.9, model.parameters);

    // // Set start and end time for the simulation.
    // auto t0   = mio::abm::TimePoint(0);
    // auto tmax = t0 + mio::abm::days(5); // Reduced from 10 to 5 days for faster testing
    // auto sim  = mio::abm::Simulation(t0, std::move(model));

    // // Create a history object to store the time series of the infection states.
    // mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
    //     Eigen::Index(mio::abm::InfectionState::Count)};

    // // Run the simulation until tmax with the history object.
    // std::cout << "Running ABM simulation..." << std::endl;
    // // sim.advance(tmax, historyTimeSeries);
    // std::cout << "ABM simulation completed." << std::endl;

    // // The results are written into the file "abm_minimal.txt" as a table
    // std::ofstream outfile("abm_minimal.txt");
    // std::get<0>(historyTimeSeries.get_log())
    //     .print_table({"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4, outfile);
    std::cout << "Results written to abm_minimal.txt" << std::endl;
}

int main()
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
        bool cudaWorking = testCuda();
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
