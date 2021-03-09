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
    model.populations[{epi::SeirInfType::E}] = 100;
    model.populations[{epi::SeirInfType::I}] = 100;
    model.populations[{epi::SeirInfType::R}] = 100;
    model.populations[{epi::SeirInfType::S}] = total_population - model.populations[{epi::SeirInfType::E}] - model.populations[{epi::SeirInfType::I}] - model.populations[{epi::SeirInfType::R}];
    // suscetible now set with every other update
    // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
    model.parameters.set<epi::StageTimeIncubationInv>(1./5.2);
    model.parameters.set<epi::StageTimeInfectiousInv>(1./6);
    model.parameters.set<epi::TransmissionRisk>(0.04);
    model.parameters.get<epi::ContactFrequency>().get_baseline()(0, 0) = 10;

    print_seir_params(model);

    auto seir = simulate(t0, tmax, dt, model);

    printf("\n number total: %f\n",
           seir.get_last_value()[0] + seir.get_last_value()[1] + seir.get_last_value()[2] + seir.get_last_value()[3]);
}
