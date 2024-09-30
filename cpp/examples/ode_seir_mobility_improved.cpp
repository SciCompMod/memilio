
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/io/mobility_io.h"
#include "models/ode_seir_mobility_improved/infection_state.h"
#include "models/ode_seir_mobility_improved/model.h"
#include "models/ode_seir_mobility_improved/parameters.h"
#include "models/ode_seir_mobility_improved/regions.h"
#include "memilio/io/io.h"
#include "memilio/io/result_io.h"
#include "Eigen/Sparse"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0.;
    ScalarType tmax = 15.;
    ScalarType dt   = 0.5;

    std::vector<int> region_ids  = {1001, 1002};
    ScalarType number_regions    = region_ids.size();
    ScalarType number_age_groups = 1;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    const std::string& mobility_data   = "";
    const std::string& trip_chain_data = "";

    mio::oseirmobilityimproved::Model<ScalarType> model(number_regions, number_age_groups);
    model.populations[{mio::oseirmobilityimproved::Region(0), mio::AgeGroup(0),
                       mio::oseirmobilityimproved::InfectionState::Exposed}]     = 10;
    model.populations[{mio::oseirmobilityimproved::Region(0), mio::AgeGroup(0),
                       mio::oseirmobilityimproved::InfectionState::Susceptible}] = 9990;
    model.populations[{mio::oseirmobilityimproved::Region(1), mio::AgeGroup(0),
                       mio::oseirmobilityimproved::InfectionState::Exposed}]     = 0;
    model.populations[{mio::oseirmobilityimproved::Region(1), mio::AgeGroup(0),
                       mio::oseirmobilityimproved::InfectionState::Susceptible}] = 10000;

    model.parameters.set<mio::oseirmobilityimproved::TransmissionProbabilityOnContact<>>(1.);

    model.parameters.set<mio::oseirmobilityimproved::TimeExposed<>>(3.);
    model.parameters.set<mio::oseirmobilityimproved::TimeInfected<>>(5.);

    model.parameters.set<mio::oseirmobilityimproved::ImpactTransmissionDuringCommuting<>>(0.);
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseirmobilityimproved::ContactPatterns<>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(2.7);
    // contact_matrix[0].add_damping(0.5, mio::SimulationTime(5));

    mio::ContactMatrixGroup& commuting_strengths =
        model.parameters.get<mio::oseirmobilityimproved::CommutingStrengths<>>().get_cont_freq_mat();
    Eigen::MatrixXd values(2, 2);
    values(0, 0)                          = 0.95;
    values(0, 1)                          = 0.05;
    values(1, 0)                          = 0.01;
    values(1, 1)                          = 0.99;
    commuting_strengths[0].get_baseline() = values;

    auto& population = model.parameters.get<mio::oseirmobilityimproved::PopulationSizes<>>();
    for (int n = 0; n < number_regions; ++n) {
        population[{mio::oseirmobilityimproved::Region(n)}] +=
            model.populations.get_group_total(mio::oseirmobilityimproved::Region(n));
        for (int m = 0; m < number_regions; ++m) {
            population[{mio::oseirmobilityimproved::Region(n)}] -=
                values(n, m) * model.populations.get_group_total(mio::oseirmobilityimproved::Region(n));
            population[{mio::oseirmobilityimproved::Region(m)}] +=
                values(n, m) * model.populations.get_group_total(mio::oseirmobilityimproved::Region(n));
        }
    }
    using DefaultIntegratorCore =
        mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>;

    std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator = std::make_shared<DefaultIntegratorCore>();

    model.check_constraints();

    auto result_from_sim = simulate(t0, tmax, dt, model, integrator);

    auto save_result_status =
        mio::save_result({result_from_sim}, region_ids, number_regions * number_age_groups, "ode_result_improved.h5");

    // bool print_to_terminal = true;

    // result_from_sim.print_table();

    // if (print_to_terminal) {

    //     std::vector<std::string> vars = {"S", "E", "I", "R"};
    //     printf("\n # t");
    //     for (size_t i = 0; i < (size_t)model.parameters.get_num_regions(); i++) {
    //         for (size_t k = 0; k < (size_t)mio::oseirmobilityimproved::InfectionState::Count; k++) {
    //             printf(" %s_%d", vars[k].c_str(), (int)i);
    //         }
    //     }

    //     auto num_points = static_cast<size_t>(result_from_sim.get_num_time_points());
    //     for (size_t i = 0; i < num_points; i++) {
    //         printf("\n%.14f ", result_from_sim.get_time(i));
    //         for (size_t k = 0; k < (size_t)model.parameters.get_num_regions(); k++) {
    //             for (size_t j = 0; j < (size_t)mio::oseirmobilityimproved::InfectionState::Count; j++) {
    //                 printf(" %.14f", result_from_sim.get_value(
    //                                      i)[j + (size_t)mio::oseirmobilityimproved::InfectionState::Count * (int)k]);
    //             }
    //         }
    //     }
    //     printf("\n");
    // }
}
