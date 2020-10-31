#include <epidemiology/secir/seir.h>
#include <epidemiology/utils/logging.h>

int main()
{
    epi::set_log_level(epi::LogLevel::debug);

    double t0   = 0;
    double tmax = 1;
    double dt   = 0.001;

    epi::log_info("Simulating SEIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    epi::SeirParams params;

    double total_population = 10000;
    params.populations.set({epi::SeirCompartments::E}, 100);
    params.populations.set({epi::SeirCompartments::I}, 100);
    params.populations.set({epi::SeirCompartments::R}, 100);
    params.populations.set({epi::SeirCompartments::S}, total_population -
                                                           params.populations.get({epi::SeirCompartments::E}) -
                                                           params.populations.get({epi::SeirCompartments::I}) -
                                                           params.populations.get({epi::SeirCompartments::R}));
    // suscetible now set with every other update
    // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
    params.times.set_incubation(5.2);
    params.times.set_infectious(6);
    params.contact_frequency.get_baseline()(0, 0) = 0.4;

    print_seir_params(params);

    auto seir = simulate(t0, tmax, dt, params);

    printf("\n number total: %f\n",
           seir.get_last_value()[0] + seir.get_last_value()[1] + seir.get_last_value()[2] + seir.get_last_value()[3]);
}
