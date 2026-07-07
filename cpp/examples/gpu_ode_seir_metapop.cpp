/***********************************************************************************************************************
########################################################################################################################
README
########################################################################################################################

Installation on sc-030247l (team server):

from "memilio/" run:

# clean build env
rm -rf cpp/build/CMake*
# load nvc
module load PrgEnv/nvhpc24-openmpi
# generate, build, run
cmake -S cpp -B cpp/build/ -DMEMILIO_ENABLE_HDF5=OFF
cmake --build cpp/build -t gpu_ode_seir_metapop_example

########################################################################################################################
***********************************************************************************************************************/

#include <iomanip>
#include <iostream>
#include <cmath>
#include <vector>

// #define CUDA_TEST_SIZE 1024
#ifdef MEMILIO_WITH_CUDA
#include <cuda_runtime.h>
#endif

#include "gpu_integrator.h"

// #include "memilio/math/integrator.h"
// #include "memilio/math/adapt_rk.h"

enum InfectionState
{
    Susceptible,
    Exposed,
    Infected,
    Recovered,
    Count
};

enum Flows
{
    StoE,
    EtoI,
    ItoR
};

size_t n_age_groups;
size_t n_regions;
std::vector<std::vector<double>>
    commuting_strengths; // param CommutingStrengths, lost time dependency, dim=(n_regions, n_regions)
std::vector<std::vector<double>>
    contact_patterns; // param ContactPatterns, lost time dependency, dim=(n_age_groups, n_age_groups)
std::vector<std::vector<double>>
    transmission_probability_on_contact; // param TransmissionProbabilityOnContact, dim=(n_regions, n_age_groups)
std::vector<std::vector<double>>
    population_after_commuting; // param PopulationAfterCommuting, element-wise inverted, dim=(n_regions, n_age_groups)
std::vector<std::vector<double>>
    time_exposed_inv; // param TimeExposed, element-wise inverted, dim=(n_regions, n_age_groups)
std::vector<std::vector<double>>
    time_infected_inv; // param TimeInfected, element-wise inverted, dim=(n_regions, n_age_groups)

struct Index {
    size_t region, age_group, state;
};
size_t product(Index i)
{
    return i.region * i.age_group * i.state;
}
size_t flatten(Index i, Index dims)
{
    return (i.region * dims.age_group + i.age_group) * dims.state + i.state;
}

void rhs(const std::vector<double>& x, double t, std::vector<double>& dxdt)
{
    (void)(t); // unused

    using enum InfectionState;
    const Index dims{n_regions, n_age_groups, Count};
    const auto flat_index = [dims](const Index& i) -> size_t {
        return flatten(i, dims);
    };
    for (size_t region_n = 0; region_n < n_regions; region_n++) {
        std::vector<double> infections_due_commuting(n_age_groups, 0.0); // for region_n
        for (size_t age_i = 0; age_i < n_age_groups; age_i++) {
            // based on matrix calculations from cpp/models/ode_seir_metapop/model.h
            // note: for this comment, '*' is a matrix-matrix product, ':' is element-wise division
            // infections_due_commuting
            //  = commuting_strengths * infectious_share_per_region
            //  = commuting_strengths * ((commuting_strengths.transpose() * infected_pop) : PopulationAfterCommuting)
            //                        ^ outer product                     ^ inner product
            // product indezes are (region_n, region_outer) * ( ( (region_outer, region_inner).T * (region_inner, age_i) ) : (region_outer, age_i) )
            for (size_t r_outer = 0; r_outer < n_regions; r_outer++) {
                double infectious_share_per_region = 0; // per "outer" region (and age_i)
                for (size_t r_inner = 0; r_inner < n_regions; r_inner++) {
                    const double& infected_pop                   = x[flat_index({r_inner, age_i, Infected})];
                    const double& commuting_strengths_transposed = commuting_strengths[r_inner][r_outer];
                    infectious_share_per_region += commuting_strengths_transposed * infected_pop;
                }
                infections_due_commuting[age_i] += commuting_strengths[region_n][r_outer] *
                                                   infectious_share_per_region /
                                                   population_after_commuting[r_outer][age_i];
            }
        }
        for (size_t age_i = 0; age_i < n_age_groups; age_i++) {
            double& flow_S_to_E = dxdt[flat_index({region_n, age_i, Flows::StoE})];

            flow_S_to_E = 0;
            for (size_t age_j = 0; age_j < n_age_groups; age_j++) {
                const size_t Ejn = flat_index({region_n, age_j, Exposed});
                const size_t Ijn = flat_index({region_n, age_j, Infected});
                const size_t Rjn = flat_index({region_n, age_j, Recovered});
                const size_t Sjn = flat_index({region_n, age_j, Susceptible});

                const double N     = x[Sjn] + x[Ejn] + x[Ijn] + x[Rjn];
                const double div_N = (N < 1e-12) ? double(0.0) : double(1.0 / N);
                const double coeffStoE =
                    0.5 * contact_patterns[age_i][age_j] * transmission_probability_on_contact[region_n][age_i];
                flow_S_to_E += (x[Ijn] * div_N + infections_due_commuting[age_j]) * coeffStoE *
                               x[flat_index({region_n, age_i, Susceptible})];
            }
            dxdt[flat_index({region_n, age_i, Flows::EtoI})] =
                x[flat_index({region_n, age_i, Exposed})] * time_exposed_inv[region_n][age_i];
            dxdt[flat_index({region_n, age_i, Flows::ItoR})] =
                x[flat_index({region_n, age_i, Exposed})] * time_infected_inv[region_n][age_i];
        }
    }
}

// originally "set_commuting_strengths"
auto set_population_after_commuting(const std::vector<double>& population, const Index& dims)
{
    // assumes population_after_commuting = 0
    for (size_t region_n = 0; region_n < n_regions; ++region_n) {
        for (size_t age = 0; age < n_age_groups; ++age) {
            double population_n = 0;
            for (size_t state = 0; state < (size_t)InfectionState::Count; state++) {
                population_n += population[flatten({region_n, age, state}, dims)];
            }
            population_after_commuting[region_n][age] += population_n;
            for (size_t region_m = 0; region_m < n_regions; ++region_m) {
                const double commuters = commuting_strengths[region_n][region_m] * population_n;
                population_after_commuting[region_n][age] -= commuters;
                population_after_commuting[region_m][age] += commuters;
            }
        }
    }
}

namespace mio
{
void log_debug(std::string_view s)
{
    std::cout << s << "\n";
}
} // namespace mio

inline std::ostream& set_ostream_format(std::ostream& out, size_t width, size_t precision, char fill = ' ')
{
    // Note: ostream& operator<< returns a reference to itself
    return out << std::setw(width) << std::fixed << std::setprecision(precision) << std::setfill(fill);
}
void print_table(std::vector<std::vector<double>> results, std::ostream& out,
                 const std::vector<std::string>& column_labels = {}, size_t width = 16, size_t precision = 5,
                 char separator = ' ', const std::string header_prefix = "\n")
{
    // Note: input manipulators (like std::setw, std::left) are consumed by the first argument written to the stream
    // print column labels
    const auto w = width, p = precision;
    out << header_prefix;
    set_ostream_format(out, w, p) << std::left << "Time";
    for (size_t k = 0; k < static_cast<size_t>(results[0].size() - 1); k++) {
        out << separator;
        if (k < column_labels.size()) {
            set_ostream_format(out, w, p) << std::left << column_labels[k];
        }
        else {
            set_ostream_format(out, w, p) << std::left << "C" + std::to_string(k + 1);
        }
    }
    // print values as table
    auto num_points = static_cast<size_t>(results.size());
    for (size_t i = 0; i < num_points; i++) {
        out << "\n";
        set_ostream_format(out, w, p) << std::right << results[i][0];
        for (size_t j = 0; j < static_cast<size_t>(results[i].size() - 1); j++) {
            out << separator;
            set_ostream_format(out, w, p) << std::right << results[i][j + 1];
        }
    }
    out << "\n";
}

int main()
{
    // using namespace mio;
    // set_log_level(LogLevel::off);

    mio::log_debug("Enter Main");

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

    // dimensions and sizes
    n_age_groups           = 1;
    n_regions              = 3;
    const Index dimensions = {n_regions, n_age_groups, InfectionState::Count};
    const size_t size      = product(dimensions);
    std::vector<double> initial_values(size); // initial condition / "population"
    // initialization scope
    {
        const auto vec_2d = [](size_t x, size_t y) {
            return std::vector<std::vector<double>>(x, std::vector<double>(y));
        };
        // apply dimensions to parameters / rhs constants
        commuting_strengths                 = vec_2d(n_regions, n_regions);
        contact_patterns                    = vec_2d(n_age_groups, n_age_groups);
        transmission_probability_on_contact = vec_2d(n_regions, n_age_groups);
        population_after_commuting          = vec_2d(n_regions, n_age_groups);
        time_exposed_inv                    = vec_2d(n_regions, n_age_groups);
        time_infected_inv                   = vec_2d(n_regions, n_age_groups);
        // initialize values
        const std::vector<double> pop{10000 - 100, 100, 0, 0};
        int T_E = 3, T_I = 7;

        for (size_t a = 0; a < n_age_groups; a++) {
            for (size_t b = 0; b < n_age_groups; b++) {
                contact_patterns[a][b] = 2.7;
            }
        }

        commuting_strengths = {{0.4, 0.3, 0.3}, {0.2, 0.7, 0.1}, {0.4, 0.1, 0.5}};

        for (size_t r = 0; r < n_regions; r++) {
            for (size_t a = 0; a < n_age_groups; a++) {
                for (size_t c = 0; c < InfectionState::Count; c++) {
                    initial_values[flatten({r, a, c}, dimensions)] = pop[c];
                }
                transmission_probability_on_contact[r][a] = 0.07333;
                population_after_commuting[r][a]          = 0; // 0-init, then use set_population_after_commuting
                time_exposed_inv[r][a]                    = 1. / T_E++;
                time_infected_inv[r][a]                   = 1. / T_I++;
            }
        }

        set_population_after_commuting(initial_values, dimensions); // originally "set_commuting_strengths"
    }
    mio::log_debug("Params Set");

#pragma acc declare copyin(readonly : n_age_groups, n_regions, commuting_strengths, contact_patterns,                  \
                           transmission_probability_on_contact, population_after_commuting, time_exposed_inv,          \
                           time_infected_inv)
    {
    }

    // print(t_offset);
    // print(t_scale);
    // print(amplitude_lincomb);
    // std::cout << "\n";

    // // std::cout << amplitude_lincomb << "\n";

    const double abs_tol = 1e-3, rel_tol = 1e-8, min_dt = 1e-2, max_dt = 1e+2;
    // // auto core = std::make_shared<mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_cash_karp54>>(abs_tol, rel_tol, min_dt, max_dt);
    // auto core = std::make_shared<mio::RKIntegratorCore<double>>(abs_tol, rel_tol, min_dt, max_dt);

    Monstrosity stepper{abs_tol,
                        rel_tol,
                        min_dt,
                        max_dt,
                        std::vector<double>(size),
                        std::vector<double>(size),
                        std::vector<std::vector<double>>(tableau().entries_low.size(), std::vector<double>(size)),
                        tableau()};

    mio::log_debug("Core Set");

    // OdeIntegrator<double> integrator(core);

    mio::log_debug("Integrator Set");

    // TimeSeries<double> results(0, Eigen::VectorXd::Zero(size));
    std::vector<std::vector<double>> results;

    const double tmax = 10;
    double dt         = 0.1;

    double t = 0.0;

    results.reserve(size_t(tmax / dt));
    results.push_back(initial_values);

    std::cout << "\n";

    mio::log_debug("Results Set");
    mio::log_debug("Integrating...");

    while (t < tmax) {
        results.emplace_back(size, 0.0);
        const size_t n = results.size() - 1;
        stepper.step(&rhs, results[n - 1], t, dt, results[n]);
        // std::cin.ignore();
    }

    // integrator.advance(rhs, tmax, dt, results);

    mio::log_debug("Integration Finished");

    // if (size < 5)
    //     results.print_table();
    // else
    //     std::cout << "Num time steps: " << results.get_num_time_points() << "\n";
    print_table(results, std::cout);

    mio::log_debug("Exit Main");
    // return 0;
}
