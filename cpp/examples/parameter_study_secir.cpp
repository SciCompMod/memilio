#include <epidemiology_io/secir_parameters_io.h>
#include <epidemiology/secir/parameter_space.h>
#include <epidemiology/secir/parameter_studies.h>
#include <epidemiology/migration/migration.h>

int main()
{
    epi::set_log_level(epi::LogLevel::debug);

    double t0   = 0;
    double tmax = 50;

    double tinc    = 5.2, // R_2^(-1)+R_3^(-1)
        tinfmild   = 6, // 4-14  (=R4^(-1))
        tserint    = 4.2, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        thosp2home = 12, // 7-16 (=R5^(-1))
        thome2hosp = 5, // 2.5-7 (=R6^(-1))
        thosp2icu  = 2, // 1-3.5 (=R7^(-1))
        ticu2home  = 8, // 5-16 (=R8^(-1))
        // tinfasy    = 6.2, // (=R9^(-1)=R_3^(-1)+0.5*R_4^(-1))
        ticu2death = 5; // 3.5-7 (=R5^(-1))

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

    epi::SecirModel model(1);
    epi::AgeGroup num_groups = model.parameters.get_num_groups();
    double fact    = 1.0 / (double)(size_t)num_groups;

    auto& params = model.parameters;

    params.set<epi::ICUCapacity>(std::numeric_limits<double>::max());
    params.set<epi::StartDay>(0);
    params.set<epi::Seasonality>(0);

    for (auto i = epi::AgeGroup(0); i < num_groups; i++) {
        params.get<epi::IncubationTime>()[i] = tinc;
        params.get<epi::InfectiousTimeMild>()[i] = tinfmild;
        params.get<epi::SerialInterval>()[i] = tserint;
        params.get<epi::HospitalizedToHomeTime>()[i] = thosp2home;
        params.get<epi::HomeToHospitalizedTime>()[i] = thome2hosp;
        params.get<epi::HospitalizedToICUTime>()[i] = thosp2icu;
        params.get<epi::ICUToHomeTime>()[i] = ticu2home;
        params.get<epi::ICUToDeathTime>()[i] = ticu2death;

        model.populations[{i,epi::InfectionState::Exposed}] = fact * num_exp_t0;
        model.populations[{i,epi::InfectionState::Carrier}] = fact * num_car_t0;
        model.populations[{i,epi::InfectionState::Infected}] = fact * num_inf_t0;
        model.populations[{i,epi::InfectionState::Hospitalized}] = fact * num_hosp_t0;
        model.populations[{i,epi::InfectionState::ICU}] = fact * num_icu_t0;
        model.populations[{i,epi::InfectionState::Recovered}] = fact * num_rec_t0;
        model.populations[{i,epi::InfectionState::Dead}] = fact * num_dead_t0;
        model.populations.set_difference_from_group_total<epi::AgeGroup>({i, epi::InfectionState::Susceptible},
                                                                         fact * num_total_t0);

        params.get<epi::InfectionProbabilityFromContact>()[i] = inf_prob;
        params.get<epi::RelativeCarrierInfectability>()[i] = carr_infec;
        params.get<epi::AsymptoticCasesPerInfectious>()[i] = alpha;
        params.get<epi::RiskOfInfectionFromSympomatic>()[i] = beta;
        params.get<epi::HospitalizedCasesPerInfectious>()[i] = rho;
        params.get<epi::ICUCasesPerHospitalized>()[i] = theta;
        params.get<epi::DeathsPerHospitalized>()[i] = delta;
    }

    params.apply_constraints();

    epi::ContactMatrixGroup& contact_matrix = params.get<epi::ContactPatterns>();
    contact_matrix[0] = epi::ContactMatrix(Eigen::MatrixXd::Constant((size_t)num_groups,
                                                                     (size_t)num_groups, fact * cont_freq));

    epi::set_params_distributions_normal(model, t0, tmax, 0.2);

    auto write_parameters_status = epi::write_json("parameters.json", model);
    if (!write_parameters_status)
    {
        std::cout << "Error writing parameters: " << write_parameters_status.error().formatted_message(); 
        return -1;
    }

    //create study
    epi::ParameterStudy<epi::SecirSimulation<>> parameter_study(model, t0, tmax, 0.2, 1);

    //run study
    int run                        = 0;
    auto lambda                    = [&run](auto&& graph) {
        auto write_result_status = epi::write_single_run_result(run++, graph);
        if (!write_result_status) {            
            std::cout << "Error writing result: " << write_result_status.error().formatted_message(); 
        }
    };
    parameter_study.run(lambda);

    return 0;
}
