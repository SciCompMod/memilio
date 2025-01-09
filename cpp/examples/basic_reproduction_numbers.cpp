
#include "models/ode_seir/model.h"
#include "models/ode_seir_mobility/model.h"
// #include "models/ode_seir_mobility_improved/model.h"

#include "memilio/math/euler.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/custom_index_array.h"

Eigen::MatrixXd get_contact_matrix()
{
    Eigen::MatrixXd contact_matrix_eigen(6, 6);
    contact_matrix_eigen << 3.9547, 1.1002, 2.9472, 2.05, 0.3733, 0.0445, 0.3327, 3.5892, 1.236, 1.9208, 0.2681, 0.0161,
        0.246, 0.7124, 5.6518, 3.2939, 0.2043, 0.0109, 0.1742, 0.8897, 3.3124, 4.5406, 0.4262, 0.0214, 0.0458, 0.1939,
        0.5782, 1.3825, 1.473, 0.0704, 0.1083, 0.1448, 0.4728, 0.9767, 0.6266, 0.1724;

    return contact_matrix_eigen;
}

const ScalarType TimeExposed[]                      = {3.335, 3.335, 3.335, 3.335, 3.335, 3.335};
const ScalarType TimeInfected[]                     = {8.0096875, 8.0096875, 8.2182, 8.1158, 8.033, 7.985};
const ScalarType TransmissionProbabilityOnContact[] = {0.03, 0.06, 0.06, 0.06, 0.09, 0.175};

void seir(size_t number_regions, ScalarType tmax)
{
    mio::set_log_level(mio::LogLevel::off);
    ScalarType t0                = 0.;
    ScalarType dt                = 0.1;
    ScalarType number_age_groups = 6;

    mio::oseir::Model<ScalarType> model(number_age_groups);
    auto& population = model.populations;

    for (size_t j = 0; j < number_age_groups; j++) {

        population[{mio::AgeGroup(j), mio::oseir::InfectionState::Susceptible}] = number_regions * 10000;
    }
    population[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] += 100;
    population[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] -= 100;

    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.template get<mio::oseir::ContactPatterns<>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline() = get_contact_matrix();

    for (size_t j = 0; j < number_age_groups; j++) {
        model.parameters.template get<mio::oseir::TimeExposed<>>()[mio::AgeGroup(j)]  = TimeExposed[j];
        model.parameters.template get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(j)] = TimeInfected[j];
        model.parameters.template get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(j)] =
            TransmissionProbabilityOnContact[j];
    }

    std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator = std::make_shared<mio::EulerIntegratorCore<>>();

    auto result = simulate(t0, tmax, dt, model, integrator);

    auto basic_reproduction_number = model.get_reproduction_number(t0, result).value();
    std::cout << "\"SEIR\": " << basic_reproduction_number << ", " << std::endl;
}

void wang(size_t number_regions, ScalarType tmax)
{
    mio::set_log_level(mio::LogLevel::off);
    ScalarType t0                = 0.;
    ScalarType dt                = 0.1;
    ScalarType number_age_groups = 6;

    mio::oseirmobility::Model<ScalarType> model(number_regions, number_age_groups);
    auto& population = model.populations;

    for (size_t j = 0; j < number_age_groups; j++) {
        for (size_t i = 0; i < number_regions; i++) {
            population[{mio::oseirmobility::Region(i), mio::AgeGroup(j),
                        mio::oseirmobility::InfectionState::Susceptible}] = 10000;
        }
    }
    population[{mio::oseirmobility::Region(0), mio::AgeGroup(0), mio::oseirmobility::InfectionState::Exposed}] += 100;
    population[{mio::oseirmobility::Region(0), mio::AgeGroup(0), mio::oseirmobility::InfectionState::Susceptible}] -=
        100;

    double fraction_commuter = 1. / (2 * number_regions);
    Eigen::MatrixXd mobility_data_commuter =
        Eigen::MatrixXd::Constant(number_regions, number_regions, fraction_commuter) -
        fraction_commuter *
            Eigen::MatrixXd::Identity(number_regions, number_regions); // Ensure that the diagonal is zero
    for (size_t county_idx_i = 0; county_idx_i < number_regions; ++county_idx_i) {
        mobility_data_commuter(county_idx_i, county_idx_i) = 1 - mobility_data_commuter.rowwise().sum()(county_idx_i);
    }
    model.parameters.template get<mio::oseirmobility::CommutingStrengths<>>().get_cont_freq_mat()[0].get_baseline() =
        mobility_data_commuter;

    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.template get<mio::oseirmobility::ContactPatterns<>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline() = get_contact_matrix();

    for (size_t j = 0; j < number_age_groups; j++) {
        model.parameters.template get<mio::oseirmobility::TimeExposed<>>()[mio::AgeGroup(j)]  = TimeExposed[j];
        model.parameters.template get<mio::oseirmobility::TimeInfected<>>()[mio::AgeGroup(j)] = TimeInfected[j];
        model.parameters.template get<mio::oseirmobility::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(j)] =
            TransmissionProbabilityOnContact[j];
    }

    std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator = std::make_shared<mio::EulerIntegratorCore<>>();

    auto result = simulate(t0, tmax, dt, model, integrator);

    auto basic_reproduction_number = model.get_reproduction_number(t0, result).value();
    std::cout << "\"Wang\": " << basic_reproduction_number << "}" << std::endl;
}

// void metapopulation(size_t number_regions, ScalarType tmax)
// {
//     mio::set_log_level(mio::LogLevel::off);
//     ScalarType t0                = 0.;
//     ScalarType dt                = 0.1;
//     ScalarType number_age_groups = 6;

//     mio::oseirmobilityimproved::Model<ScalarType> model(number_regions, number_age_groups);
//     auto& population = model.populations;

//     for (size_t j = 0; j < number_age_groups; j++) {
//         for (size_t i = 0; i < number_regions; i++) {
//             population[{mio::oseirmobilityimproved::Region(i), mio::AgeGroup(j),
//                         mio::oseirmobilityimproved::InfectionState::Susceptible}] = 10000;
//         }
//     }
//     population[{mio::oseirmobilityimproved::Region(0), mio::AgeGroup(0),
//                 mio::oseirmobilityimproved::InfectionState::Exposed}] += 100;
//     population[{mio::oseirmobilityimproved::Region(0), mio::AgeGroup(0),
//                 mio::oseirmobilityimproved::InfectionState::Susceptible}] -= 100;

//     double fraction_commuter = 1. / (2 * number_regions);
//     Eigen::MatrixXd mobility_data_commuter =
//         Eigen::MatrixXd::Constant(number_regions, number_regions, fraction_commuter) -
//         fraction_commuter *
//             Eigen::MatrixXd::Identity(number_regions, number_regions); // Ensure that the diagonal is zero
//     for (size_t county_idx_i = 0; county_idx_i < number_regions; ++county_idx_i) {
//         mobility_data_commuter(county_idx_i, county_idx_i) = 1 - mobility_data_commuter.rowwise().sum()(county_idx_i);
//     }
//     model.parameters.template get<mio::oseirmobilityimproved::CommutingStrengths<>>()
//         .get_cont_freq_mat()[0]
//         .get_baseline() = mobility_data_commuter;

//     mio::ContactMatrixGroup& contact_matrix =
//         model.parameters.template get<mio::oseirmobilityimproved::ContactPatterns<>>().get_cont_freq_mat();
//     contact_matrix[0].get_baseline() = get_contact_matrix();

//     for (size_t j = 0; j < number_age_groups; j++) {
//         model.parameters.template get<mio::oseirmobilityimproved::TimeExposed<>>()[mio::AgeGroup(j)]  = TimeExposed[j];
//         model.parameters.template get<mio::oseirmobilityimproved::TimeInfected<>>()[mio::AgeGroup(j)] = TimeInfected[j];
//         model.parameters
//             .template get<mio::oseirmobilityimproved::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(j)] =
//             TransmissionProbabilityOnContact[j];
//     }

//     mio::ContactMatrixGroup& commuting_strengths =
//         model.parameters.template get<mio::oseirmobilityimproved::CommutingStrengths<>>().get_cont_freq_mat();

//     auto& population_after_commuting = model.m_population_after_commuting;
//     for (size_t region_n = 0; region_n < number_regions; ++region_n) {
//         for (size_t age = 0; age < number_age_groups; ++age) {
//             double population_n = 0;
//             for (size_t state = 0; state < (size_t)mio::oseirmobilityimproved::InfectionState::Count; state++) {
//                 population_n += population[{mio::oseirmobilityimproved::Region(region_n), mio::AgeGroup(age),
//                                             mio::oseirmobilityimproved::InfectionState(state)}];
//             }
//             population_after_commuting[{mio::oseirmobilityimproved::Region(region_n), mio::AgeGroup(age)}] +=
//                 population_n;
//             for (size_t region_m = 0; region_m < number_regions; ++region_m) {
//                 population_after_commuting[{mio::oseirmobilityimproved::Region(region_n), mio::AgeGroup(age)}] -=
//                     commuting_strengths[0].get_baseline()(region_n, region_m) * population_n;
//                 population_after_commuting[{mio::oseirmobilityimproved::Region(region_m), mio::AgeGroup(age)}] +=
//                     commuting_strengths[0].get_baseline()(region_n, region_m) * population_n;
//             }
//         }
//     }

//     std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator = std::make_shared<mio::EulerIntegratorCore<>>();

//     auto result = simulate(t0, tmax, dt, model, integrator);

//     auto basic_reproduction_number = model.get_reproduction_number(t0, result).value();
//     std::cout << "\"Metapopulation\": " << basic_reproduction_number << "}" << std::endl;
// }

int main()
{
    const ScalarType tmax = 1.;
    size_t num_regions    = 150;

    std::cout << "{ \"Regions\": " << num_regions << ", " << std::endl;

    seir(num_regions, tmax);
    wang(num_regions, tmax);
    // metapopulation(num_regions, tmax);
    return 0;
}