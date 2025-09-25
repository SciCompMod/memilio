#include "ad/ad.hpp"
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/parameter_space.h"
#include "ode_secirvvs/parameters_io.h"
#include "memilio/compartments/simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/io/mobility_io.h"
#include "memilio/io/result_io.h"
#include "memilio/math/eigen.h"
#include "memilio/data/analyze_result.h"
#include "ode_secirvvs/analyze_result.h"

#include "boost/outcome/try.hpp"
#include "boost/outcome/result.hpp"
#include "boost/optional.hpp"
#include "boost/filesystem.hpp"

namespace fs = boost::filesystem;


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

/**
 * different types of NPI, used as DampingType.
 */
enum class Intervention
{
    Home,
    SchoolClosure,
    HomeOffice,
    GatheringBanFacilitiesClosure,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Count,
};

/**
 * different level of NPI, used as DampingLevel.
 */
enum class InterventionLevel
{
    Main,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Holidays,
    Count,
};

static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}
                                                                        };

using FP = double;
// using FP = typename ad::ga1s<double>::type; // algorithmic differentiation data type: scalar tangent-linear mode

mio::TimeSeries<FP> aggregate_result(mio::TimeSeries<FP> results, int num_groups)
{
    // Aggregate results for each compartment
    size_t num_days = results.get_num_time_points();
    mio::TimeSeries<FP> results_aggregated = mio::TimeSeries<FP>::zero(num_days, (size_t)mio::osecirvvs::InfectionState::Count);

    size_t day = 0;
    for (auto t : results.get_times()) {
        results_aggregated.get_time(day) = t;
        for (size_t infection_state = 0; infection_state < (size_t)mio::osecirvvs::InfectionState::Count; infection_state++){
            for (size_t age = 0; age < (size_t)num_groups; age++) {

                auto age_group_offset = age * (size_t)mio::osecirvvs::InfectionState::Count;
                results_aggregated[day](infection_state) += results[day](infection_state + age_group_offset);
            }
        }
        day ++;
    }
    return results_aggregated;
}


/**
 * Set contact matrices.
 * Reads contact matrices from files in the data directory.
 * @param data_dir data directory.
 * @param params Object that the contact matrices will be added to.
 * @returns any io errors that happen during reading of the files.
 */
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::osecirvvs::Parameters<FP>& params)
{
    //TODO: io error handling
    auto contact_matrices = mio::ContactMatrixGroup<FP>(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));

        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline.cast<FP>();
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic>::Zero(6, 6);
    }
    params.get<mio::osecirvvs::ContactPatterns<FP>>() = mio::UncertainContactMatrix<FP>(contact_matrices);

    return mio::success();
}

mio::IOResult<void> set_synthetic_contact_matrices(mio::osecirvvs::Parameters<FP>& params)
{  
    auto contact_matrix_group = mio::ContactMatrixGroup<FP>(contact_locations.size(), size_t(1));
    contact_matrix_group[0].get_baseline().setConstant(1);
    contact_matrix_group[1].get_baseline().setConstant(1);
    contact_matrix_group[2].get_baseline().setConstant(1);
    contact_matrix_group[3].get_baseline().setConstant(1);
    params.get<mio::osecirvvs::ContactPatterns<FP>>() = mio::UncertainContactMatrix<FP>(contact_matrix_group);

    return mio::success();
}

mio::IOResult<void> set_population_data(const fs::path& filename, mio::osecirvvs::Model<FP>& model)
{  
    // BOOST_OUTCOME_TRY(auto&& num_lines, mio::count_lines((filename).string()));

    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return mio::failure(mio::StatusCode::FileNotFound, (filename).string());
    }

    try {
        std::string tp;
        size_t linenumber = 0;
        while (getline(file, tp)) {
            auto line = mio::split(tp, '\t');

            for (size_t j = 0; j < size_t(line.size()); j++) {
                model.populations[{(mio::AgeGroup)linenumber, (mio::osecirvvs::InfectionState)j}] = std::stod(line[j]);
            }
            linenumber++;
        }
    }
    catch (std::runtime_error& ex) {
        return mio::failure(mio::StatusCode::InvalidFileFormat, (filename).string() + ": " + ex.what());
    }
    // mio::osecirvvs::draw_sample_demographics(model);

    // std::cout << ad::value(model.populations[{mio::AgeGroup(5), mio::osecirvvs::InfectionState::SusceptibleNaive}].value()) << "\n";
    return mio::success();
}

mio::IOResult<void> set_synthetic_population_data(mio::osecirvvs::Model<FP>& model)
{
    for (mio::AgeGroup i = 0; i < model.parameters.get_num_groups(); i++) {
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedNaive}]                                = 10;
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedImprovedImmunity}]                     = 11;
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedPartialImmunity}]                      = 12;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}]                     = 13;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}]            = 13;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}]           = 14;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]  = 14;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}]          = 15;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 15;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}]                       = 5;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}]              = 5;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}]             = 6;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]    = 6;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}]            = 7;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]   = 7;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereNaive}]                         = 8;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereImprovedImmunity}]              = 1;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSeverePartialImmunity}]               = 2;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalNaive}]                       = 3;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalPartialImmunity}]             = 4;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalImprovedImmunity}]            = 5;
        model.populations[{i, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}]                 = 6;
        model.populations[{i, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}]                  = 7;
        model.populations[{i, mio::osecirvvs::InfectionState::DeadNaive}]                    = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::DeadPartialImmunity}]          = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::DeadImprovedImmunity}]         = 0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>(
            {i, mio::osecirvvs::InfectionState::SusceptibleNaive}, FP(100000));
    }

    return mio::success();
}

mio::IOResult<void> set_covid_parameters(mio::osecirvvs::Parameters<FP>& params, FP tmax)
{
    params.get<mio::osecirvvs::ICUCapacity<FP>>()          = 100;
    params.get<mio::osecirvvs::TestAndTraceCapacity<FP>>() = 0.0143;
    const size_t daily_vaccinations                                      = 10;
    params.get<mio::osecirvvs::DailyPartialVaccinations<FP>>().resize(
        mio::SimulationDay((size_t)ad::value(tmax) + 1));
    params.get<mio::osecirvvs::DailyFullVaccinations<FP>>().resize(mio::SimulationDay((size_t)ad::value(tmax) + 1));

    for (mio::AgeGroup i = 0; i < params.get_num_groups(); i++) {
        for (size_t j = 0; j < tmax + 1; ++j) {
            auto num_vaccinations = static_cast<FP>(j * daily_vaccinations);
            params
                .get<mio::osecirvvs::DailyPartialVaccinations<FP>>()[{i, mio::SimulationDay(j)}] =
                num_vaccinations;
            params
                .get<mio::osecirvvs::DailyFullVaccinations<FP>>()[{i, mio::SimulationDay(j)}] =
                num_vaccinations;
        }
        params.get<mio::osecirvvs::DynamicNPIsImplementationDelay<FP>>() = 7;

        //times
        params.get<mio::osecirvvs::TimeExposed<FP>>()[i]            = 3.33;
        params.get<mio::osecirvvs::TimeInfectedNoSymptoms<FP>>()[i] = 1.87;
        params.get<mio::osecirvvs::TimeInfectedSymptoms<FP>>()[i]   = 7;
        params.get<mio::osecirvvs::TimeInfectedSevere<FP>>()[i]     = 6;
        params.get<mio::osecirvvs::TimeInfectedCritical<FP>>()[i]   = 7;

        //probabilities
        params.get<mio::osecirvvs::TransmissionProbabilityOnContact<FP>>()[i] = 0.15;
        params.get<mio::osecirvvs::RelativeTransmissionNoSymptoms<FP>>()[i]   = 0.5;
        // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
        // depends on incidence and test and trace capacity
        params.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<FP>>()[i]    = 0.0;
        params.get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<FP>>()[i] = 0.4;
        params.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<FP>>()[i]    = 0.2;
        params.get<mio::osecirvvs::SeverePerInfectedSymptoms<FP>>()[i]         = 0.1;
        params.get<mio::osecirvvs::CriticalPerSevere<FP>>()[i]                 = 0.1;
        params.get<mio::osecirvvs::DeathsPerCritical<FP>>()[i]                 = 0.1;

        params.get<mio::osecirvvs::ReducExposedPartialImmunity<FP>>()[i]           = 0.8;
        params.get<mio::osecirvvs::ReducExposedImprovedImmunity<FP>>()[i]          = 0.331;
        params.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<FP>>()[i]  = 0.65;
        params.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<FP>>()[i] = 0.243;
        params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[i] =
            0.1;
        params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[i] =
            0.091;
        params.get<mio::osecirvvs::ReducTimeInfectedMild<FP>>()[i] = 0.9;
    }

    params.get<mio::osecirvvs::Seasonality<FP>>() = 0.2;

    return mio::success();
}

mio::IOResult<void> set_npis(mio::osecirvvs::Parameters<FP>& params, FP tmax)
{
    auto damping_helper = [=](mio::SimulationTime<FP> t, FP min, FP max, mio::DampingLevel damping_level, mio::DampingType damping_type, const std::vector<size_t> location,
                Eigen::VectorX<FP> group_weights) {
        auto p = mio::UncertainValue<FP>(0.5 * (max + min));
        p.set_distribution(mio::ParameterDistributionUniform(ad::value(min), ad::value(max)));
        return mio::DampingSampling<FP>(p, damping_level, damping_type, t, location, group_weights);
    };

    auto group_weights_all     = Eigen::VectorX<FP>::Constant(size_t(params.get_num_groups()), 1.0);

    auto school_closure = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(int(InterventionLevel::Main)),
                                            mio::DampingType(int(Intervention::SchoolClosure)),
                                            {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto home_office = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(int(InterventionLevel::Main)),
                                            mio::DampingType(int(Intervention::HomeOffice)),
                                            {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto physical_distancing_school = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto physical_distancing_work = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto physical_distancing_other = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::Other)}, group_weights_all);
    };

    auto& contacts         = params.get<mio::osecirvvs::ContactPatterns<FP>>();
    auto& contact_dampings = contacts.get_dampings();

    FP controllValue = 1.;
    for(int t = 0; t < ad::value(tmax); t++)
    {
        contact_dampings.push_back(school_closure(mio::SimulationTime<FP>(t), FP(1. * controllValue))); 
        contact_dampings.push_back(home_office(mio::SimulationTime<FP>(t), FP(0.25 * controllValue)));
        contact_dampings.push_back(physical_distancing_school(mio::SimulationTime<FP>(t), FP(0.25 * controllValue)));
        contact_dampings.push_back(physical_distancing_work(mio::SimulationTime<FP>(t), FP(0.25 * controllValue)));
        contact_dampings.push_back(physical_distancing_other(mio::SimulationTime<FP>(t), FP(0.35 * controllValue)));
    }
    contacts.make_matrix();
    // std::cout << ad::value(contact_matrix.get_matrix_at(FP(4.0))) << "\n";
    // std::cout << ad::value(contact_matrix.get_matrix_at(FP(5.0))) << "\n";
    // std::cout << ad::value(contact_matrix.get_matrix_at(FP(6.0))) << "\n";

    return mio::success();
}

mio::IOResult<void> run(const fs::path& data_dir)
{
    FP t0   = 0;
    FP tmax = 56;
    FP dt   = 0.1;
    size_t num_age_groups = 6;

    mio::log_info("Simulating SECIRVVS; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::osecirvvs::Parameters<FP> params(num_age_groups);
    
    BOOST_OUTCOME_TRY(set_covid_parameters(params, tmax));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));
    // BOOST_OUTCOME_TRY(set_npis(params, tmax));


    mio::osecirvvs::Model<FP> model(num_age_groups);
    model.parameters = params;

    // BOOST_OUTCOME_TRY(set_synthetic_population_data(model));
    BOOST_OUTCOME_TRY(set_population_data((data_dir / ".." / "compartment_initialization_2025-05-26" / "initialization_sum.txt"), model));
    model.apply_constraints();

    // use default Cash-Karp adaptive integrator
    mio::osecirvvs::Simulation<FP> sim(model, t0, dt);
    BOOST_OUTCOME_TRY(set_npis(sim.get_model().parameters, tmax));

    sim.advance(tmax);
    mio::TimeSeries<FP> result = sim.get_result();

    // mio::TimeSeries<FP> result = mio::osecirvvs::simulate<FP>(t0, tmax, dt, model);

    // bool print_to_terminal = true;

    // auto contact_matrix         = model.parameters.get<mio::osecirvvs::ContactPatterns<FP>>().get_cont_freq_mat();
    // std::cout << ad::value(contact_matrix.get_matrix_at(FP(4.0))) << "\n";
    // if (print_to_terminal) {
    //     printf("%.14f\n", ad::value(result.get_last_time()));
    //     for (size_t j = 0; j < (size_t)mio::osecirvvs::InfectionState::Count; j++) {
    //         printf("compartment %d: %.14f\n", (int)j, ad::value(result.get_last_value()[j]));
    //     }
    //     printf("%.14f\n", ad::value(result.get_last_time()));
    //     for (size_t j = 0; j < (size_t)mio::osecirvvs::InfectionState::Count; j++) {
    //         printf("compartment %d: %.14f\n", (int)j, ad::value(result.get_last_value()[j]));
    //     }
    // }

    // auto result_interpolated = mio::interpolate_simulation_result(result);

    std::ofstream outFileResult("Result.csv");
    std::vector<std::string> vars = {"S_n", "S_p", "E_n", "E_p", "E_i", "C_n", "C_p", "C_i", "C_confirmed_n", "C_confirmed_p", "C_confirmed_i",  "I_n", "I_p", "I_i", "I_confirmed_n", "I_confirmed_p", "I_confirmed_i", "H_n", "H_p", "H_i", "U_n", "U_p", "U_i", "S_i", "D_n", "D_p", "D_i"};
    mio::TimeSeries<FP> result_interpolated = mio::interpolate_simulation_result(result);
    mio::TimeSeries<FP> result_aggregated = aggregate_result(result_interpolated, num_age_groups);
    result_aggregated.print_table(outFileResult, vars, 21, 10);
    outFileResult.close();

    FP constraint = 0;
    for (int t = 0; t < tmax; t++) {
        // current_constraint += result_aggregated.get_value(t)[(int)mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed] + result_aggregated.get_value(t)[(int)mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed] + result_aggregated.get_value(t)[(int)mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed] + 
        //     result_aggregated.get_value(t)[(int)mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed] + result_aggregated.get_value(t)[(int)mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed] + result_aggregated.get_value(t)[(int)mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed];
        FP current_constraint = result_aggregated.get_value(t)[(int)mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive] + result_aggregated.get_value(t)[(int)mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity] + result_aggregated.get_value(t)[(int)mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity] + 
            result_aggregated.get_value(t)[(int)mio::osecirvvs::InfectionState::InfectedSymptomsNaive] + result_aggregated.get_value(t)[(int)mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity] + result_aggregated.get_value(t)[(int)mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity];
        
        std::cout << "Value at " << t << " is " << result_aggregated.get_value(t)[(int)mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive] << "\n";
        constraint = std::max(constraint, current_constraint);
    }
    std::cout << "Max: " << constraint << "\n";
    return mio::success();
}

int main(int argc, char** argv)
{
    mio::set_log_level(mio::LogLevel::info);
    std::string data_dir;

    if (argc == 2){
        data_dir   = argv[1];
    }
    else{
        printf("Input Arguments: {data_dir}\n");
        return 0;
    }

    auto result = run(data_dir);
    return 0;
}