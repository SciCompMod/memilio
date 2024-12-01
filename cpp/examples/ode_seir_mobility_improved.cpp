
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
#include "memilio/io/epi_data.h"

#include <chrono>

template <typename FP = ScalarType>
mio::IOResult<void> set_population_data(mio::oseirmobilityimproved::Model<FP>& model,
                                        const std::string& population_data_path)
{
    BOOST_OUTCOME_TRY(auto&& node_ids, mio::get_node_ids(population_data_path, true, true));

    BOOST_OUTCOME_TRY(const auto&& population_data, mio::read_population_data(population_data_path, true));

    for (auto&& entry : population_data) {
        auto it = std::find_if(node_ids.begin(), node_ids.end(), [&entry](auto r) {
            return r == 0 ||
                   (entry.county_id && mio::regions::StateId(r) == mio::regions::get_state_id(int(*entry.county_id))) ||
                   (entry.county_id && mio::regions::CountyId(r) == *entry.county_id) ||
                   (entry.district_id && mio::regions::DistrictId(r) == *entry.district_id);
        });
        if (it != node_ids.end()) {
            auto region_idx = size_t(it - node_ids.begin());
            for (size_t age = 0; age < (size_t)model.parameters.get_num_agegroups(); age++) {
                model.populations[{mio::oseirmobilityimproved::Region(region_idx), mio::AgeGroup(age),
                                   mio::oseirmobilityimproved::InfectionState::Susceptible}] =
                    entry.population[mio::AgeGroup(age)];
            }
        }
    }

    return mio::success();
}

template <typename FP = ScalarType>
mio::IOResult<void> set_mobility_weights(mio::oseirmobilityimproved::Model<FP>& model, const std::string& mobility_data)
{
    size_t number_regions = (size_t)model.parameters.get_num_regions();
    if (number_regions == 1) {
        model.parameters.template get<mio::oseirmobilityimproved::CommutingStrengths<>>()
            .get_cont_freq_mat()[0]
            .get_baseline()
            .setConstant(1.0);

        return mio::success();
    }
    else {
        // mobility between nodes
        BOOST_OUTCOME_TRY(auto&& mobility_data_commuter, mio::read_mobility_plain(mobility_data));
        if (mobility_data_commuter.rows() != Eigen::Index(number_regions) ||
            mobility_data_commuter.cols() != Eigen::Index(number_regions)) {
            return mio::failure(mio::StatusCode::InvalidValue,
                                "Mobility matrices do not have the correct size. You may need to run "
                                "transformMobilitydata.py from pycode memilio epidata package.");
        }

        for (size_t county_idx_i = 0; county_idx_i < number_regions; ++county_idx_i) {
            auto population_i = model.populations.get_group_total(mio::oseirmobilityimproved::Region(county_idx_i));
            auto test         = population_i;
            mobility_data_commuter.row(county_idx_i) /= population_i;
            mio::unused(test);
            mobility_data_commuter(county_idx_i, county_idx_i) =
                1 - mobility_data_commuter.rowwise().sum()(county_idx_i);
        }
        model.parameters.template get<mio::oseirmobilityimproved::CommutingStrengths<>>()
            .get_cont_freq_mat()[0]
            .get_baseline() = mobility_data_commuter;

        return mio::success();
    }
}

template <typename FP = ScalarType>
mio::IOResult<void> set_parameters_and_population(mio::oseirmobilityimproved::Model<FP>& model,
                                                  const std::string& population_data_path,
                                                  const std::string& mobility_data)
{
    BOOST_OUTCOME_TRY(set_population_data(model, population_data_path));

    auto& populations = model.populations;
    auto& parameters  = model.parameters;

    size_t number_regions    = (size_t)parameters.get_num_regions();
    size_t number_age_groups = (size_t)parameters.get_num_agegroups();

    BOOST_OUTCOME_TRY(set_mobility_weights(model, mobility_data));

    populations[{mio::oseirmobilityimproved::Region(0), mio::AgeGroup(3),
                 mio::oseirmobilityimproved::InfectionState::Susceptible}] -= 100;
    populations[{mio::oseirmobilityimproved::Region(0), mio::AgeGroup(3),
                 mio::oseirmobilityimproved::InfectionState::Exposed}] += 100;

    mio::ContactMatrixGroup& contact_matrix =
        parameters.template get<mio::oseirmobilityimproved::ContactPatterns<>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(7.95);

    parameters.template set<mio::oseirmobilityimproved::TransmissionProbabilityOnContact<>>(1.);

    parameters.template set<mio::oseirmobilityimproved::TimeExposed<>>(3.);
    parameters.template set<mio::oseirmobilityimproved::TimeInfected<>>(5.);

    mio::ContactMatrixGroup& commuting_strengths =
        parameters.template get<mio::oseirmobilityimproved::CommutingStrengths<>>().get_cont_freq_mat();

    auto& population_after_commuting = model.m_population_after_commuting;
    for (size_t region_n = 0; region_n < number_regions; ++region_n) {
        for (size_t age = 0; age < number_age_groups; ++age) {
            double population_n = 0;
            for (size_t state = 0; state < (size_t)mio::oseirmobilityimproved::InfectionState::Count; state++) {
                population_n += populations[{mio::oseirmobilityimproved::Region(region_n), mio::AgeGroup(age),
                                             mio::oseirmobilityimproved::InfectionState(state)}];
            }
            population_after_commuting[{mio::oseirmobilityimproved::Region(region_n), mio::AgeGroup(age)}] +=
                population_n;
            for (size_t region_m = 0; region_m < number_regions; ++region_m) {
                population_after_commuting[{mio::oseirmobilityimproved::Region(region_n), mio::AgeGroup(age)}] -=
                    commuting_strengths[0].get_baseline()(region_n, region_m) * population_n;
                population_after_commuting[{mio::oseirmobilityimproved::Region(region_m), mio::AgeGroup(age)}] +=
                    commuting_strengths[0].get_baseline()(region_n, region_m) * population_n;
            }
        }
    }

    return mio::success();
}

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0.;
    ScalarType tmax = 0.2;
    ScalarType dt   = 0.1;

    ScalarType number_regions = 53;
    std::vector<int> region_ids(number_regions);
    iota(region_ids.begin(), region_ids.end(), 1);
    ScalarType number_age_groups = 6;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    const std::string& mobility_data   = "";
    const std::string& population_data = "";

    mio::oseirmobilityimproved::Model<ScalarType> model(number_regions, number_age_groups);
    auto result_prepare_simulation = set_parameters_and_population(model, population_data, mobility_data);

    // using DefaultIntegratorCore =
    //     mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>;

    std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator = std::make_shared<mio::EulerIntegratorCore<>>();

    model.check_constraints();

    printf("Start Simulation\n");
    auto t1              = std::chrono::high_resolution_clock::now();
    auto result_from_sim = simulate(t0, tmax, dt, model, integrator);
    auto t2              = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> ms_double = t2 - t1;

    printf("Runtime: %f\n", ms_double.count());
    result_from_sim.print_table({"S", "E", "I", "R"});

    auto save_result_status =
        mio::save_result({result_from_sim}, region_ids, number_regions * number_age_groups, "ode_result_test.h5");

    auto reproduction_numbers = model.get_reproduction_numbers(result_from_sim);
    std::cout << "\nbasis reproduction number: " << reproduction_numbers[0] << "\n";
}
