#include <cstddef>
#include <iostream>
#include <chrono>
#include <utility>
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/parameter_distribution_wrapper.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/random_number_generator.h"

int main()
{
    mio::ParameterDistributionExponential exp1(1.);
    mio::ParameterDistributionWrapper distWrapper = mio::ParameterDistributionWrapper(exp1);
    mio::ParameterDistributionAbstraction distAbstraction1 =
        mio::ParameterDistributionAbstraction(mio::ParameterDistributionExponential(1.));
    mio::ParameterDistributionAbstraction distAbstraction(std::move(distAbstraction1));
    // std::chrono::steady_clock::time_point begin_wrapper = std::chrono::steady_clock::now();

    // for (size_t i = 0; i < 10000000; ++i) {
    //     auto v = distWrapper.get(mio::thread_local_rng());
    //     mio::unused(v);
    // }
    // std::chrono::steady_clock::time_point end_wrapper = std::chrono::steady_clock::now();
    // std::cout << "Time difference wrapper = "
    //           << std::chrono::duration_cast<std::chrono::seconds>(end_wrapper - begin_wrapper).count() << "[s]"
    //           << std::endl;
    // std::cout << "Time difference wrapper = "
    //           << std::chrono::duration_cast<std::chrono::milliseconds>(end_wrapper - begin_wrapper).count() << "[ms]"
    //           << std::endl;

    std::chrono::steady_clock::time_point begin_abstraction = std::chrono::steady_clock::now();
    for (size_t i = 0; i < 10000000; ++i) {
        auto v = distAbstraction.get(mio::thread_local_rng());
        mio::unused(v);
    }
    std::chrono::steady_clock::time_point end_abstraction = std::chrono::steady_clock::now();
    std::cout << "Time difference abstraction = "
              << std::chrono::duration_cast<std::chrono::seconds>(end_abstraction - begin_abstraction).count() << "[s]"
              << std::endl;
    std::cout << "Time difference abstraction = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_abstraction - begin_abstraction).count()
              << "[ms]" << std::endl;
    return 0;
}
