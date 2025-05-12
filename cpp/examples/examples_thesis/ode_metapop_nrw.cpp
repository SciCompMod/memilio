
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/io/mobility_io.h"
#include "models/ode_metapop/infection_state.h"
#include "models/ode_metapop/model.h"
#include "models/ode_metapop/parameters.h"
#include "models/ode_metapop/regions.h"
#include "memilio/io/io.h"
#include "memilio/io/result_io.h"
#include "memilio/io/epi_data.h"

/**
 * Set epidemiological parameters of Sars-CoV-2 for a immune-naive
 * population and wild type variant.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_covid_parameters(mio::oseirmetapop::Parameters<double>& params)
{
    params.template set<mio::oseirmetapop::TimeExposed<>>(3.335);
    params.get<mio::oseirmetapop::TimeInfected<>>()[mio::AgeGroup(0)] = 8.0096875;
    params.get<mio::oseirmetapop::TimeInfected<>>()[mio::AgeGroup(1)] = 8.0096875;
    params.get<mio::oseirmetapop::TimeInfected<>>()[mio::AgeGroup(2)] = 8.2182;
    params.get<mio::oseirmetapop::TimeInfected<>>()[mio::AgeGroup(3)] = 8.1158;
    params.get<mio::oseirmetapop::TimeInfected<>>()[mio::AgeGroup(4)] = 8.033;
    params.get<mio::oseirmetapop::TimeInfected<>>()[mio::AgeGroup(5)] = 7.985;

    params.get<mio::oseirmetapop::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(0)] = 0.03;
    params.get<mio::oseirmetapop::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(1)] = 0.06;
    params.get<mio::oseirmetapop::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(2)] = 0.06;
    params.get<mio::oseirmetapop::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(3)] = 0.06;
    params.get<mio::oseirmetapop::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(4)] = 0.09;
    params.get<mio::oseirmetapop::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(5)] = 0.175;
    printf("Setting epidemiological parameters successful.\n");
    return mio::success();
}

/**
 * indices of contact matrix corresponding to locations where contacts occur.
 */
enum class ContactLocation
{
    Home = 0,
    School,
    Work,
    Other,
    Count,
};

static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}};

/**
 * Set contact matrices.
 * Reads contact matrices from files in the data directory.
 * @param data_dir data directory.
 * @param params Object that the contact matrices will be added to.
 * @returns any io errors that happen during reading of the files.
 */
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::oseirmetapop::Parameters<double>& params)
{
    //TODO: io error handling
    auto contact_matrices = mio::ContactMatrixGroup(contact_locations.size(), size_t(params.get_num_agegroups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(
            auto&& baseline,
            mio::read_mobility_plain(
                (data_dir / "Germany" / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));
        BOOST_OUTCOME_TRY(
            auto&& minimum,
            mio::read_mobility_plain(
                (data_dir / "Germany" / "contacts" / ("minimum_" + contact_location.second + ".txt")).string()));
        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = minimum;
    }
    params.get<mio::oseirmetapop::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);

    printf("Setting contact matrices successful.\n");
    return mio::success();
}

template <typename FP = ScalarType>
mio::IOResult<void> set_population_data(mio::oseirmetapop::Model<FP>& model, const fs::path& data_dir)
{
    BOOST_OUTCOME_TRY(
        auto&& node_ids,
        mio::get_node_ids((data_dir / "Germany" / "pydata" / "county_current_population_nrw.json").string(), true,
                          true));

    BOOST_OUTCOME_TRY(const auto&& population_data,
                      mio::read_population_data(
                          (data_dir / "Germany" / "pydata" / "county_current_population_nrw.json").string(), true));

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
                model.populations[{mio::oseirmetapop::Region(region_idx), mio::AgeGroup(age),
                                   mio::oseirmetapop::InfectionState::Susceptible}] =
                    entry.population[mio::AgeGroup(age)];
            }
        }
    }

    printf("Setting population data successful.\n");
    return mio::success();
}

template <typename FP = ScalarType>
mio::IOResult<void> set_mobility_weights(mio::oseirmetapop::Model<FP>& model, const fs::path& data_dir)
{
    size_t number_regions = (size_t)model.parameters.get_num_regions();
    if (number_regions == 1) {
        model.set_commuting_strengths();

        return mio::success();
    }
    else {
        // mobility between nodes
        BOOST_OUTCOME_TRY(
            auto&& mobility_data_commuter,
            mio::read_mobility_plain((data_dir / "Germany" / "mobility" / "commuter_mobility_2022_nrw.txt").string()));
        if (mobility_data_commuter.rows() != Eigen::Index(number_regions) ||
            mobility_data_commuter.cols() != Eigen::Index(number_regions)) {
            return mio::failure(mio::StatusCode::InvalidValue,
                                "Mobility matrices do not have the correct size. You may need to run "
                                "transformMobilitydata.py from pycode memilio epidata package.");
        }

        for (size_t county_idx_i = 0; county_idx_i < number_regions; ++county_idx_i) {
            auto population_i = model.populations.get_group_total(mio::oseirmetapop::Region(county_idx_i));
            mobility_data_commuter.row(county_idx_i) /= population_i;
            mobility_data_commuter(county_idx_i, county_idx_i) =
                1 - mobility_data_commuter.rowwise().sum()(county_idx_i);
        }

        model.set_commuting_strengths(mobility_data_commuter);

        printf("Setting mobility weights successful.\n");
        return mio::success();
    }
}

template <typename FP = ScalarType>
mio::IOResult<void> set_parameters_and_population(mio::oseirmetapop::Model<FP>& model, const fs::path& data_dir)
{
    auto& populations = model.populations;
    auto& parameters  = model.parameters;

    BOOST_OUTCOME_TRY(set_population_data(model, data_dir));
    populations[{mio::oseirmetapop::Region(27), mio::AgeGroup(4), mio::oseirmetapop::InfectionState::Susceptible}] -=
        100;
    populations[{mio::oseirmetapop::Region(27), mio::AgeGroup(4), mio::oseirmetapop::InfectionState::Exposed}] += 100;

    BOOST_OUTCOME_TRY(set_mobility_weights(model, data_dir));

    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, parameters))

    BOOST_OUTCOME_TRY(set_covid_parameters(parameters));

    return mio::success();
}

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0.;
    ScalarType tmax = 100.;
    ScalarType dt   = 0.1;

    ScalarType number_regions    = 53;
    ScalarType number_age_groups = 6;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    const std::string& data_dir = "";

    mio::oseirmetapop::Model<ScalarType> model(number_regions, number_age_groups);
    auto result_prepare_simulation = set_parameters_and_population(model, data_dir);

    model.check_constraints();

    printf("Start Simulation\n");
    auto result_from_sim = simulate(t0, tmax, dt, model);

    auto result = mio::interpolate_simulation_result(result_from_sim);

    auto save_result_status = mio::save_result({result}, {1}, number_regions * number_age_groups, "ode_result_nrw.h5");
}
