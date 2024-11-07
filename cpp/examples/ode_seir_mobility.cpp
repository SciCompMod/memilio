
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/io/mobility_io.h"
#include "models/ode_seir_mobility/infection_state.h"
#include "models/ode_seir_mobility/model.h"
#include "models/ode_seir_mobility/parameters.h"
#include "models/ode_seir_mobility/regions.h"
#include "memilio/io/io.h"
#include "memilio/io/result_io.h"

template <typename FP = ScalarType>
mio::IOResult<void> set_mobility_weights(const std::string& mobility_data, mio::oseirmobility::Model<FP>& model,
                                         size_t number_regions)
{
    // mobility between nodes
    BOOST_OUTCOME_TRY(auto&& mobility_data_commuter,
                      mio::read_mobility_plain(mobility_data + "/mobility" + "/commuter_migration_test.txt"));
    if (mobility_data_commuter.rows() != Eigen::Index(number_regions) ||
        mobility_data_commuter.cols() != Eigen::Index(number_regions)) {
        return mio::failure(mio::StatusCode::InvalidValue,
                            "Mobility matrices do not have the correct size. You may need to run "
                            "transformMobilitydata.py from pycode memilio epidata package.");
    }

    for (size_t county_idx_i = 0; county_idx_i < number_regions; ++county_idx_i) {
        auto population_i = model.populations.get_group_total(mio::oseirmobility::Region(county_idx_i));
        mobility_data_commuter.row(county_idx_i) /= population_i;
    }
    model.parameters.template get<mio::oseirmobility::CommutingStrengths<>>().get_cont_freq_mat()[0].get_baseline() =
        mobility_data_commuter;

    return mio::success();
}

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0.;
    ScalarType tmax = 15.;
    ScalarType dt   = 0.1;

    ScalarType number_regions = 2;
    std::vector<int> region_ids(number_regions);
    iota(region_ids.begin(), region_ids.end(), 1);
    ScalarType number_age_groups = 1;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    const std::string& mobility_data = "";

    mio::oseirmobility::Model<ScalarType> model(number_regions, number_age_groups);
    model.populations[{mio::oseirmobility::Region(0), mio::AgeGroup(0), mio::oseirmobility::InfectionState::Exposed}] =
        10;
    model.populations[{mio::oseirmobility::Region(0), mio::AgeGroup(0),
                       mio::oseirmobility::InfectionState::Susceptible}] = 9990;
    for (int i = 1; i < number_regions; i++) {
        model.populations[{mio::oseirmobility::Region(i), mio::AgeGroup(0),
                           mio::oseirmobility::InfectionState::Exposed}]     = 0;
        model.populations[{mio::oseirmobility::Region(i), mio::AgeGroup(0),
                           mio::oseirmobility::InfectionState::Susceptible}] = 10000;
    }

    model.parameters.set<mio::oseirmobility::TransmissionProbabilityOnContact<>>(1.);

    model.parameters.set<mio::oseirmobility::TimeExposed<>>(3.);
    model.parameters.set<mio::oseirmobility::TimeInfected<>>(5.);

    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseirmobility::ContactPatterns<>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(2.7);
    // contact_matrix[0].add_damping(0.5, mio::SimulationTime(5));

    auto result_preprocess = set_mobility_weights(mobility_data, model, number_regions);

    using DefaultIntegratorCore =
        mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>;

    std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator = std::make_shared<DefaultIntegratorCore>();

    model.check_constraints();

    auto result_from_sim = simulate(t0, tmax, dt, model, integrator);

    auto save_result_status =
        mio::save_result({result_from_sim}, region_ids, number_regions * number_age_groups, "ode_result_standard.h5");

    // bool print_to_terminal = true;

    // result_from_sim.print_table();

    // if (print_to_terminal) {

    //     std::vector<std::string> vars = {"S", "E", "I", "R"};
    //     printf("\n # t");
    //     for (size_t i = 0; i < (size_t)model.parameters.get_num_regions(); i++) {
    //         for (size_t k = 0; k < (size_t)mio::oseirmobility::InfectionState::Count; k++) {
    //             printf(" %s_%d", vars[k].c_str(), (int)i);
    //         }
    //     }

    //     auto num_points = static_cast<size_t>(sir.get_num_time_points());
    //     for (size_t i = 0; i < num_points; i++) {
    //         printf("\n%.14f ", sir.get_time(i));
    //         for (size_t k = 0; k < (size_t)model.parameters.get_num_regions(); k++) {
    //             for (size_t j = 0; j < (size_t)mio::oseirmobility::InfectionState::Count; j++) {
    //                 printf(" %.14f", sir.get_value(i)[j + (size_t)mio::oseirmobility::InfectionState::Count * (int)k]);
    //             }
    //         }
    //     }
    //     printf("\n");
    // }
}
