#include "memilio/compartments/simulation.h"
#include "models/ode_seir_metapop/model.h"
#include "models/ode_seir_metapop/parameters.h"
#include "models/ode_seir_metapop/infection_state.h"
#include "memilio/geography/regions.h"
#include "memilio/math/stepper_wrapper.h"
#include "memilio/math/euler.h"
#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

/*
 * Benchmark for SEIR metapopulation model runtime scaling.
 * Varies number of regions, keeps one age group, fixed parameters, simulates to tmax=500.
 * Writes CSV to compare_results/cpp_benchmark.csv with columns:
 * Regions,TimeSteps,RuntimeSeconds
 */

double median(std::vector<double> values)
{
    if (values.empty()) {
        return 0.0;
    }
    auto mid = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    double med = values[mid];
    if (values.size() % 2 == 0) {
        auto lower_max = *std::max_element(values.begin(), values.begin() + mid);
        med = 0.5 * (med + lower_max);
    }
    return med;
}

struct BenchmarkConfig
{
    std::string name;
    std::string output_path;
    std::function<std::shared_ptr<mio::IntegratorCore<double>>()> make_integrator;
};

int main(int argc, char** argv)
{
    using FP                     = double;
    double dt                    = 0.1; // fixed step (overridable via --dt)
    double t0                    = 0.0;
    double tmax                  = 500.0; // overridable via --tmax
    std::vector<int> region_grid = {1, 2, 4, 8, 16, 32, 64, 128, 256};
    int runs                     = 40;

    // CLI parsing: flags --dt <val>, --tmax <val>, positional integers = regions
    if (argc > 1) {
        region_grid.clear();
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--dt" && i + 1 < argc) {
                dt = std::stod(argv[++i]);
            }
            else if (arg.rfind("--dt=", 0) == 0) {
                dt = std::stod(arg.substr(5));
            }
            else if (arg == "--tmax" && i + 1 < argc) {
                tmax = std::stod(argv[++i]);
            }
            else if (arg.rfind("--tmax=", 0) == 0) {
                tmax = std::stod(arg.substr(7));
            }
            else if (arg == "--runs" && i + 1 < argc) {
                runs = std::max(1, std::stoi(argv[++i]));
            }
            else if (arg.rfind("--runs=", 0) == 0) {
                runs = std::max(1, std::stoi(arg.substr(7)));
            }
            else {
                try {
                    region_grid.push_back(std::stoi(arg));
                }
                catch (...) {
                    std::cerr << "Ignoring unrecognized arg: " << arg << "\n";
                }
            }
        }
        if (region_grid.empty()) {
            region_grid = {1, 2, 4, 8, 16, 32, 64, 128, 256};
        }
    }

    std::vector<BenchmarkConfig> configs;
    configs.push_back({"RK4", "/localdata1/code_2025/memilio/compare_results/cpp_benchmark_rk4.csv",
                       []() {
                           return std::make_shared<
                               mio::ExplicitStepperWrapper<double, boost::numeric::odeint::runge_kutta4>>();
                       }});
    configs.push_back({"Euler", "/localdata1/code_2025/memilio/compare_results/cpp_benchmark_euler.csv",
                       []() { return std::make_shared<mio::EulerIntegratorCore<double>>(); }});

    runs = std::max(1, runs);

    for (const auto& config : configs) {
        std::ofstream ofs(config.output_path);
        if (!ofs.is_open()) {
            std::cerr << "Failed to open output file: " << config.output_path << "\n";
            return 1;
        }
        ofs << "Regions,TimeSteps,RuntimeSeconds,TotalNoIOSeconds\n";

        for (auto n_regions : region_grid) {
            std::vector<double> sim_durations;
            std::vector<double> total_durations;
            sim_durations.reserve(runs);
            total_durations.reserve(runs);
            int steps = 0;

            for (int run_idx = 0; run_idx < runs; ++run_idx) {
                auto total_start = std::chrono::steady_clock::now();

                std::shared_ptr<mio::IntegratorCore<FP>> integrator = config.make_integrator();

                mio::oseirmetapop::Model<FP> model(n_regions, 1);
                for (int r = 0; r < n_regions; ++r) {
                    model.populations[{mio::regions::Region(r), mio::AgeGroup(0),
                                       mio::oseirmetapop::InfectionState::Susceptible}] = 10000;
                }
                model.populations[{mio::regions::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Exposed}] +=
                    100;
                model.populations[{mio::regions::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Susceptible}] -=
                    100;

                Eigen::MatrixXd commuting = Eigen::MatrixXd::Identity(n_regions, n_regions);
                model.set_commuting_strengths(commuting);

                model.parameters.template get<mio::oseirmetapop::ContactPatterns<>>()
                    .get_cont_freq_mat()[0]
                    .get_baseline()
                    .setConstant(2.7);
                model.parameters.set<mio::oseirmetapop::TimeExposed<>>(3.335);
                model.parameters.set<mio::oseirmetapop::TimeInfected<>>(8.097612257);
                model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<>>(0.07333);

                auto sim_start                      = std::chrono::steady_clock::now();
                auto result                         = mio::simulate(t0, tmax, dt, model, integrator);
                auto sim_end                        = std::chrono::steady_clock::now();
                auto total_end                      = std::chrono::steady_clock::now();
                std::chrono::duration<double> sim_d = sim_end - sim_start;
                std::chrono::duration<double> tot_d = total_end - total_start;

                sim_durations.push_back(sim_d.count());
                total_durations.push_back(tot_d.count());
                steps = result.get_num_time_points();
            }

            double sim_median = median(sim_durations);
            double tot_median = median(total_durations);
            ofs << n_regions << "," << steps << "," << sim_median << "," << tot_median << "\n";
            std::cout << "[" << config.name << "] Regions=" << n_regions << " steps=" << steps
                      << " median_sim=" << sim_median << "s median_total_no_io=" << tot_median << "s (dt=" << dt
                      << ", tmax=" << tmax << ", runs=" << runs << ")\n";
        }
        std::cout << "Benchmark CSV written to: " << config.output_path << "\n";
    }
    return 0;
}
