#include <epidemiology/seir.h>
#include <epidemiology/logging.h>

int main()
{
    epi::set_log_level(epi::LogLevel::debug);

    double t0   = 0;
    double tmax = 1;
    double dt   = 0.1;

    epi::log_info("Simulating SEIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    epi::SeirParams params;

    print_seir_params(params);

    std::vector<std::vector<double>> seir(0);

    simulate(t0, tmax, dt, params, seir);
}