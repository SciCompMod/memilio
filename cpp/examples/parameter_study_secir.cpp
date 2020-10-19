#include <epidemiology_io/secir_parameters_io.h>
#include <epidemiology/secir/parameter_space.h>
#include <epidemiology/secir/parameter_studies.h>
#include <epidemiology_io/secir_parameters_io.h>
#include <epidemiology/migration/migration.h>

#include <tixi.h>

int main(int argc, char* argv[])
{
    epi::set_log_level(epi::LogLevel::debug);

    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double tinc    = 5.2, // R_2^(-1)+R_3^(-1)
        tinfmild   = 6, // 4-14  (=R4^(-1))
        tserint    = 4.2, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        thosp2home = 12, // 7-16 (=R5^(-1))
        thome2hosp = 5, // 2.5-7 (=R6^(-1))
        thosp2icu  = 2, // 1-3.5 (=R7^(-1))
        ticu2home  = 8, // 5-16 (=R8^(-1))
        tinfasy    = 6.2, // (=R9^(-1)=R_3^(-1)+0.5*R_4^(-1))
        ticu2death = 5; // 3.5-7 (=R5^(-1))

    double tinfasy2 = 1.0 / (0.5 / (tinfmild - tserint) + 0.5 / tinfmild);
    if (fabs(tinfasy2 - tinfasy) > 0) {
        epi::log_warning("----> TODO / To consider: In the HZI paper, tinfasy (the asymptomatic infectious time) or "
                         "R9^(-1)=R_3^(-1)+0.5*R_4^(-1) is directly given by R_3 and R_4 and maybe should not be an "
                         "'additional parameter'");
    }

    double cont_freq = 10, // see Polymod study
        inf_prob = 0.05, carr_infec = 0.67,
           alpha = 0.09, // 0.01-0.16
        beta     = 0.25, // 0.05-0.5
        delta    = 0.3, // 0.15-0.77
        rho      = 0.2, // 0.1-0.35
        theta    = 0.25; // 0.15-0.4

    double num_total_t0 = 10000, num_exp_t0 = 100, num_inf_t0 = 50, num_car_t0 = 50, num_hosp_t0 = 20, num_icu_t0 = 10,
           num_rec_t0 = 10, num_dead_t0 = 0;

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    epi::SecirModel<epi::AgeGroup1> model;
    int num_groups = model.parameters.get_num_groups();
    double fact    = 1.0 / (double)num_groups;

    auto& params = model.parameters;

    params.set_icu_capacity(std::numeric_limits<double>::max());
    params.set_start_day(0);
    params.set_seasonality(0);

    for (size_t i = 0; i < num_groups; i++) {
        params.times[i].set_incubation(tinc);
        params.times[i].set_infectious_mild(tinfmild);
        params.times[i].set_serialinterval(tserint);
        params.times[i].set_hospitalized_to_home(thosp2home);
        params.times[i].set_home_to_hospitalized(thome2hosp);
        params.times[i].set_hospitalized_to_icu(thosp2icu);
        params.times[i].set_icu_to_home(ticu2home);
        params.times[i].set_infectious_asymp(tinfasy);
        params.times[i].set_icu_to_death(ticu2death);

        model.populations.set(fact * num_exp_t0, (epi::AgeGroup1)i, epi::InfectionType::E);
        model.populations.set(fact * num_car_t0, (epi::AgeGroup1)i, epi::InfectionType::C);
        model.populations.set(fact * num_inf_t0, (epi::AgeGroup1)i, epi::InfectionType::I);
        model.populations.set(fact * num_hosp_t0, (epi::AgeGroup1)i, epi::InfectionType::H);
        model.populations.set(fact * num_icu_t0, (epi::AgeGroup1)i, epi::InfectionType::U);
        model.populations.set(fact * num_rec_t0, (epi::AgeGroup1)i, epi::InfectionType::R);
        model.populations.set(fact * num_dead_t0, (epi::AgeGroup1)i, epi::InfectionType::D);
        model.populations.set_difference_from_group_total(fact * num_total_t0, (epi::AgeGroup1)i, (epi::AgeGroup1)i,
                                                          epi::InfectionType::S);

        params.probabilities[i].set_infection_from_contact(inf_prob);
        params.probabilities[i].set_carrier_infectability(carr_infec);
        params.probabilities[i].set_asymp_per_infectious(alpha);
        params.probabilities[i].set_risk_from_symptomatic(beta);
        params.probabilities[i].set_hospitalized_per_infectious(rho);
        params.probabilities[i].set_icu_per_hospitalized(theta);
        params.probabilities[i].set_dead_per_icu(delta);
    }

    epi::ContactFrequencyMatrix& cont_freq_matrix = params.get_contact_patterns();
    epi::Damping dummy(30., 0.3);
    for (int i = 0; i < num_groups; i++) {
        for (int j = i; j < num_groups; j++) {
            cont_freq_matrix.set_cont_freq(fact * cont_freq, i, j);
        }
    }

    epi::set_params_distributions_normal(model, t0, tmax, 0.2);

    // write parameter space without parameter study
    std::string path2 = "/Parameters2";
    TixiDocumentHandle handle2;
    tixiCreateDocument("Parameters2", &handle2);
    epi::write_parameter_space(handle2, path2, model, 0, 2);
    tixiSaveDocument(handle2, "Parameters2.xml");
    tixiCloseDocument(handle2);

    // draw sample and write same parameter space but with different current values
    epi::draw_sample(model);
    std::string path3 = "/Parameters3";
    TixiDocumentHandle handle3;
    tixiCreateDocument("Parameters3", &handle3);
    epi::write_parameter_space(handle3, path3, model, 0, 2);
    tixiSaveDocument(handle3, "Parameters3.xml");
    tixiCloseDocument(handle3);

    // create study
    epi::ParameterStudy<epi::SecirModel<epi::AgeGroup1>> parameter_study(model, t0, tmax, 0.2, 1);

    // write and run study
    std::string path = "/Parameters";
    TixiDocumentHandle handle;
    tixiCreateDocument("Parameters", &handle);

    int io_mode = 2;
    epi::write_parameter_study(handle, path, parameter_study, io_mode);
    tixiSaveDocument(handle, "Parameters.xml");
    tixiCloseDocument(handle);

    tixiOpenDocument("Parameters.xml", &handle);
    epi::ParameterStudy<epi::SecirModel<epi::AgeGroup1>> read_study =
        epi::read_parameter_study<epi::SecirModel<epi::AgeGroup1>>(handle, path);
    int run     = 0;
    auto lambda = [&run, t0, tmax](const auto& model, const auto& secir_result, int node) {
        epi::write_single_run_params(run++, model, t0, tmax, secir_result, node);
    };
    auto results = read_study.run(lambda);

#if 0
    if (argc > 1) {
        // If provided, the first argument is the input file
        input_filename = argv[1];
    }
    else {
        // If not provided, we use a sample input file
        input_filename = "parameter_studies_example_input.txt";
    }

    // Create parameter study
    ParameterStudy parameter_study(input_filename);
#endif
    return 0;
}
