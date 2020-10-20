#include <epidemiology/secir/seir.h>
#include <epidemiology/model/simulation.h>
#include <epidemiology/utils/logging.h>

int main()
{
    epi::set_log_level(epi::LogLevel::debug);

    double t0   = 0;
    double tmax = 1;
    double dt   = 0.001;

    epi::log_info("Simulating SEIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    epi::SeirModel model;

    double total_population = 10000;
    model.populations.set(100, epi::SeirInfType::E);
    model.populations.set(100, epi::SeirInfType::I);
    model.populations.set(100, epi::SeirInfType::R);
    model.populations.set(total_population - model.populations.get(epi::SeirInfType::E) -
                              model.populations.get(epi::SeirInfType::I) - model.populations.get(epi::SeirInfType::R),
                          epi::SeirInfType::S);
    // suscetible now set with every other update
    // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
    model.parameters.times.set_incubation(5.2);
    model.parameters.times.set_cont_freq(0.4);
    model.parameters.times.set_infectious(6);

    print_seir_params(model);

    auto seir = simulate(t0, tmax, dt, model);

    printf("\n number total: %f\n",
           seir.get_last_value()[0] + seir.get_last_value()[1] + seir.get_last_value()[2] + seir.get_last_value()[3]);
}
