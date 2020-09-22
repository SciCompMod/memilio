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
    params.times.set_cont_freq(0.4);
    params.times.set_infectious(6);

    print_seir_params(params);

    std::vector<Eigen::VectorXd> seir(0);

    simulate(t0, tmax, dt, params, seir);

    printf("\n number total: %f\n",
           seir[seir.size() - 1][0] + seir[seir.size() - 1][1] + seir[seir.size() - 1][2] + seir[seir.size() - 1][3]);
}
