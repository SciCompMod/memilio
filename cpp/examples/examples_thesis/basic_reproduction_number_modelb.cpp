
#include "models/ode_metapop_wang/model.h"

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

void calculate_basic_reproduction_number(size_t number_regions, ScalarType tmax)
{
    mio::set_log_level(mio::LogLevel::off);
    ScalarType t0                = 0.;
    ScalarType dt                = 0.1;
    ScalarType number_age_groups = 6;

    mio::oseirmetapopwang::Model<ScalarType> model(number_regions, number_age_groups);
    auto& population = model.populations;

    for (size_t j = 0; j < number_age_groups; j++) {
        for (size_t i = 0; i < number_regions; i++) {
            population[{mio::oseirmetapopwang::Region(i), mio::AgeGroup(j),
                        mio::oseirmetapopwang::InfectionState::Susceptible}] = 10000;
        }
    }

    double fraction_commuter = 1. / (2 * number_regions);
    Eigen::MatrixXd mobility_data_commuter =
        Eigen::MatrixXd::Constant(number_regions, number_regions, fraction_commuter) -
        fraction_commuter *
            Eigen::MatrixXd::Identity(number_regions, number_regions); // Ensure that the diagonal is zero
    for (size_t county_idx_i = 0; county_idx_i < number_regions; ++county_idx_i) {
        mobility_data_commuter(county_idx_i, county_idx_i) = 1 - mobility_data_commuter.rowwise().sum()(county_idx_i);
    }
    model.parameters.template get<mio::oseirmetapopwang::CommutingStrengths<>>().get_cont_freq_mat()[0].get_baseline() =
        mobility_data_commuter;

    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.template get<mio::oseirmetapopwang::ContactPatterns<>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline() = get_contact_matrix();

    for (size_t j = 0; j < number_age_groups; j++) {
        model.parameters.template get<mio::oseirmetapopwang::TimeExposed<>>()[mio::AgeGroup(j)]  = TimeExposed[j];
        model.parameters.template get<mio::oseirmetapopwang::TimeInfected<>>()[mio::AgeGroup(j)] = TimeInfected[j];
        model.parameters.template get<mio::oseirmetapopwang::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(j)] =
            TransmissionProbabilityOnContact[j];
    }

    std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator = std::make_shared<mio::EulerIntegratorCore<>>();

    auto result = simulate(t0, tmax, dt, model, integrator);

    auto basic_reproduction_number = model.get_reproduction_number(t0, result).value();
    std::cout << "\"Model B\": " << basic_reproduction_number << "}" << std::endl;
}

int main()
{
    const ScalarType tmax = 1.;
    size_t num_regions    = 1;

    std::cout << "{ \"Regions\": " << num_regions << ", " << std::endl;

    calculate_basic_reproduction_number(num_regions, tmax);
    return 0;
}