#include "epidemiology/secir/parameter_studies.h"
#include "epidemiology/utils/regions.h"
#include "epidemiology_io/secir_parameters_io.h"
#include "epidemiology_io/secir_result_io.h"
#include "epidemiology_io/mobility_io.h"
#include "boost/filesystem.hpp"
#include <cstdio>
#include <iomanip>

namespace fs = boost::filesystem;

//contact matrix indices
enum class ContactLocation
{
    Home = 0,
    School,
    Work,
    Other,
    Count,
};

//types of damping
enum class Intervention
{
    Home = 0,
    SchoolClosure,
    HomeOffice,
    GatheringBanFacilitiesClosure,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
};

void assign_uniform_distribution(epi::UncertainValue& p, double min, double max)
{
    p   = epi::UncertainValue(0.5 * (max + min));
    min = max = double(p);
    p.set_distribution(epi::ParameterDistributionUniform(min, max));
}

template <size_t N>
void array_assign_uniform_distribution(epi::CustomIndexArray<epi::UncertainValue, epi::AgeGroup>& array,
                                       const double (&min)[N], const double (&max)[N])
{
    assert(N == array.numel());
    for (auto i = epi::AgeGroup(0); i < epi::AgeGroup(N); ++i) {
        assign_uniform_distribution(array[i], min[size_t(i)], max[size_t(i)]);
    }
}

void array_assign_uniform_distribution(epi::CustomIndexArray<epi::UncertainValue, epi::AgeGroup>& array, double min,
                                       double max)
{
    for (auto i = epi::AgeGroup(0); i < array.size<epi::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min, max);
    }
}

epi::IOResult<void> set_covid_parameters(epi::SecirParams& params)
{
    //times
    const double tinc             = 5.2; // R_2^(-1)+R_3^(-1)
    const double tserint_min      = 0.5 * 2.67 + 0.5 * 5.2; // R_2^(-1)+0.5*R_3^(-1)
    const double tserint_max      = 0.5 * 4.00 + 0.5 * 5.2;
    const double t_inf_rec_min    = 5.6; // R4^(-1) = T_I^R
    const double t_inf_rec_max    = 8.4;
    const double t_inf_hosp_min[] = {9, 9, 9, 5, 5, 5}; // R6^(-1) = T_I^H
    const double t_inf_hosp_max[] = {12, 12, 12, 7, 7, 7};
    const double t_hosp_rec_min[] = {4, 4, 5, 7, 9, 13}; // R5^(-1) = T_H^R
    const double t_hosp_rec_max[] = {6, 6, 7, 9, 11, 17};
    const double t_hosp_icu_min   = 3; // R7^(-1) = T_H^U
    const double t_hosp_icu_max   = 7;
    const double t_icu_rec_min[]  = {5, 5, 5, 14, 14, 10}; // R8^(-1) = T_U^R
    const double t_icu_rec_max[]  = {9, 9, 9, 21, 21, 15};
    const double t_icu_dead_min[] = {4, 4, 4, 15, 15, 10}; // 5-16 (=R8^(-1) = T_U^R)
    const double t_icu_dead_max[] = {8, 8, 8, 18, 18, 12};

    array_assign_uniform_distribution(params.get<epi::IncubationTime>(), tinc, tinc);
    array_assign_uniform_distribution(params.get<epi::SerialInterval>(), tserint_min, tserint_max);
    array_assign_uniform_distribution(params.get<epi::InfectiousTimeMild>(), t_inf_rec_min, t_inf_rec_max);
    array_assign_uniform_distribution(params.get<epi::HomeToHospitalizedTime>(), t_inf_hosp_min, t_inf_hosp_max);
    array_assign_uniform_distribution(params.get<epi::HospitalizedToHomeTime>(), t_hosp_rec_min, t_hosp_rec_max);
    array_assign_uniform_distribution(params.get<epi::HospitalizedToICUTime>(), t_hosp_icu_min, t_hosp_icu_max);
    array_assign_uniform_distribution(params.get<epi::ICUToHomeTime>(), t_icu_rec_min, t_icu_rec_max);
    array_assign_uniform_distribution(params.get<epi::ICUToDeathTime>(), t_icu_dead_min, t_icu_dead_max);

    //probabilities
    const double transmission_risk_min[] = {0.02, 0.05, 0.05, 0.05, 0.08, 0.15};
    const double transmission_risk_max[] = {0.04, 0.07, 0.07, 0.07, 0.10, 0.20};
    const double carr_infec_min          = 1;
    const double carr_infec_max          = 1;
    const double beta_low_incidenc_min   = 0.1; // beta (depends on incidence and test and trace capacity)
    const double beta_low_incidenc_max   = 0.3;
    const double beta_high_incidence_min = 0.3;
    const double beta_high_incidence_max = 0.5;
    const double prob_car_rec_min[]      = {0.2, 0.2, 0.15, 0.15, 0.15, 0.15}; // alpha
    const double prob_car_rec_max[]      = {0.3, 0.3, 0.25, 0.25, 0.25, 0.25};
    const double prob_inf_hosp_min[]     = {0.006, 0.006, 0.015, 0.049, 0.15, 0.20}; // rho
    const double prob_inf_hosp_max[]     = {0.009, 0.009, 0.023, 0.074, 0.18, 0.25};
    const double prob_hosp_icu_min[]     = {0.05, 0.05, 0.05, 0.10, 0.25, 0.35}; // theta
    const double prob_hosp_icu_max[]     = {0.10, 0.10, 0.10, 0.20, 0.35, 0.45};
    const double prob_icu_dead_min[]     = {0.00, 0.00, 0.10, 0.10, 0.30, 0.5}; // delta
    const double prob_icu_dead_max[]     = {0.10, 0.10, 0.18, 0.18, 0.50, 0.7};

    array_assign_uniform_distribution(params.get<epi::InfectionProbabilityFromContact>(), transmission_risk_min,
                                      transmission_risk_max);
    array_assign_uniform_distribution(params.get<epi::RelativeCarrierInfectability>(), carr_infec_min, carr_infec_max);
    array_assign_uniform_distribution(params.get<epi::RiskOfInfectionFromSympomatic>(), beta_low_incidenc_min,
                                      beta_low_incidenc_max);
    array_assign_uniform_distribution(params.get<epi::MaxRiskOfInfectionFromSympomatic>(), beta_high_incidence_min,
                                      beta_high_incidence_max);
    array_assign_uniform_distribution(params.get<epi::AsymptoticCasesPerInfectious>(), prob_car_rec_min,
                                      prob_car_rec_max);
    array_assign_uniform_distribution(params.get<epi::HospitalizedCasesPerInfectious>(), prob_inf_hosp_min,
                                      prob_inf_hosp_max);
    array_assign_uniform_distribution(params.get<epi::ICUCasesPerHospitalized>(), prob_hosp_icu_min, prob_hosp_icu_max);
    array_assign_uniform_distribution(params.get<epi::DeathsPerHospitalized>(), prob_icu_dead_min, prob_icu_dead_max);

    //sasonality
    const double seasonality_min = 0.1;
    const double seasonality_max = 0.3;

    assign_uniform_distribution(params.get<epi::Seasonality>(), seasonality_min, seasonality_max);

    return epi::success();
}

static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}};

epi::IOResult<void> set_contact_matrices(const fs::path& data_dir, epi::SecirParams& params)
{
    //TODO: io error handling
    auto contact_matrices = epi::ContactMatrixGroup(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(baseline,
                          epi::read_mobility_plain(
                              (data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));
        BOOST_OUTCOME_TRY(minimum,
                          epi::read_mobility_plain(
                              (data_dir / "contacts" / ("minimum_" + contact_location.second + ".txt")).string()));
        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = minimum;
    }
    params.get<epi::ContactPatterns>() = epi::UncertainContactMatrix(contact_matrices);

    return epi::success();
}

epi::IOResult<void> set_npis(const fs::path& data_dir, epi::Date start_date, epi::Date end_date,
                             epi::SecirParams& params)
{
    auto& contacts         = params.get<epi::ContactPatterns>();
    auto& contact_dampings = contacts.get_dampings();

    //weights for age groups affected by an NPI
    auto group_weights_all     = Eigen::VectorXd::Constant(size_t(params.get_num_groups()), 1.0);
    auto group_weights_seniors = Eigen::VectorXd::NullaryExpr(size_t(params.get_num_groups()), [](auto&& i) {
        return i == 5 ? 1.0 : i == 4 ? 0.5 : 0.0; //65-80 only partially
    });

    //helper functions that create dampings for specific NPIs
    auto contacts_at_home = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(1), epi::DampingType(int(Intervention::Home)), t,
                                    {size_t(ContactLocation::Home)}, group_weights_all);
    };
    auto school_closure = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(1), epi::DampingType(int(Intervention::SchoolClosure)), t,
                                    {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto home_office = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(1), epi::DampingType(int(Intervention::HomeOffice)), t,
                                    {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto social_events = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(1),
                                    epi::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), t,
                                    {size_t(ContactLocation::Other)}, group_weights_all);
    };
    auto social_events_work = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(1),
                                    epi::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), t,
                                    {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto physical_distancing_home_school = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(
            v, epi::DampingLevel(2), epi::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
            {size_t(ContactLocation::Home), size_t(ContactLocation::School)}, group_weights_all);
    };
    auto physical_distancing_work_other = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(2),
                                    epi::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
                                    {size_t(ContactLocation::Work), size_t(ContactLocation::Other)}, group_weights_all);
    };
    auto senior_awareness = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(3), epi::DampingType(int(Intervention::SeniorAwareness)), t,
                                    {size_t(ContactLocation::Home), size_t(ContactLocation::Other)},
                                    group_weights_seniors);
    };

    //SPRING 2020 LOCKDOWN SCENARIO
    auto start_spring_date = epi::Date(2020, 3, 18);
    if (start_spring_date < end_date) {
        auto start_spring = epi::SimulationTime(epi::get_offset_in_days(start_spring_date, start_date));
        contact_dampings.push_back(contacts_at_home(start_spring, 0.6, 0.8));
        contact_dampings.push_back(school_closure(start_spring, 1.0, 1.0));
        contact_dampings.push_back(home_office(start_spring, 0.2, 0.3));
        contact_dampings.push_back(social_events(start_spring, 0.6, 0.8));
        contact_dampings.push_back(social_events_work(start_spring, 0.1, 0.2));
        contact_dampings.push_back(physical_distancing_home_school(start_spring, 0.4, 0.6));
        contact_dampings.push_back(physical_distancing_work_other(start_spring, 0.4, 0.6));
        contact_dampings.push_back(senior_awareness(start_spring, 0.0, 0.0));
    }

    // SUMMER 2020 SCENARIO
    auto start_summer_date = epi::Date(2020, 5, 15);
    if (start_summer_date < end_date) {
        auto start_summer       = epi::SimulationTime(epi::get_offset_in_days(start_summer_date, start_date));
        auto school_reopen_time = epi::SimulationTime(epi::get_offset_in_days(epi::Date(2020, 6, 15), start_date));
        contact_dampings.push_back(contacts_at_home(start_summer, 0.0, 0.2));
        contact_dampings.push_back(school_closure(start_summer, 0.5, 0.5)); //schools partially reopened
        contact_dampings.push_back(school_closure(school_reopen_time, 0.0, 0.0)); //school fully reopened
        contact_dampings.push_back(home_office(start_summer, 0.2, 0.3));
        contact_dampings.push_back(social_events(start_summer, 0.0, 0.2));
        contact_dampings.push_back(social_events_work(start_summer, 0.0, 0.05));
        contact_dampings.push_back(physical_distancing_home_school(start_summer, 0.0, 0.2));
        contact_dampings.push_back(physical_distancing_work_other(start_summer, 0.0, 0.2));
        contact_dampings.push_back(senior_awareness(start_summer, 0.0, 0.0));
    }

    //autumn enforced attention
    auto start_autumn_date = epi::Date(2020, 10, 1);
    if (start_autumn_date < end_date) {
        auto start_autumn = epi::SimulationTime(epi::get_offset_in_days(start_autumn_date, start_date));
        contact_dampings.push_back(contacts_at_home(start_autumn, 0.2, 0.4));
        contact_dampings.push_back(physical_distancing_home_school(start_autumn, 0.2, 0.4));
        contact_dampings.push_back(physical_distancing_work_other(start_autumn, 0.2, 0.4));
    }

    //autumn lockdown light
    auto start_autumn_lockdown_date = epi::Date(2020, 11, 1);
    if (start_autumn_lockdown_date < end_date) {
        auto start_autumn_lockdown =
            epi::SimulationTime(epi::get_offset_in_days(start_autumn_lockdown_date, start_date));
        contact_dampings.push_back(contacts_at_home(start_autumn_lockdown, 0.4, 0.6));
        contact_dampings.push_back(school_closure(start_autumn_lockdown, 0.0, 0.0));
        contact_dampings.push_back(home_office(start_autumn_lockdown, 0.2, 0.3));
        contact_dampings.push_back(social_events(start_autumn_lockdown, 0.6, 0.8));
        contact_dampings.push_back(social_events_work(start_autumn_lockdown, 0.0, 0.1));
        contact_dampings.push_back(physical_distancing_home_school(start_autumn_lockdown, 0.2, 0.4));
        contact_dampings.push_back(physical_distancing_work_other(start_autumn_lockdown, 0.4, 0.6));
        contact_dampings.push_back(senior_awareness(start_autumn_lockdown, 0.0, 0.0));
    }

    //winter lockdown
    auto start_winter_lockdown_date = epi::Date(2020, 12, 16);
    if (start_winter_lockdown_date < end_date) {
        double min = 0.6, max = 0.8; //for strictest scenario: 0.8 - 1.0
        auto start_winter_lockdown =
            epi::SimulationTime(epi::get_offset_in_days(start_winter_lockdown_date, start_date));
        contact_dampings.push_back(contacts_at_home(start_winter_lockdown, min, max));
        contact_dampings.push_back(school_closure(start_winter_lockdown, 1.0, 1.0));
        contact_dampings.push_back(home_office(start_winter_lockdown, 0.2, 0.3));
        contact_dampings.push_back(social_events(start_winter_lockdown, min, max));
        contact_dampings.push_back(social_events_work(start_winter_lockdown, 0.1, 0.2));
        contact_dampings.push_back(physical_distancing_home_school(start_winter_lockdown, 0.2, 0.4));
        contact_dampings.push_back(physical_distancing_work_other(start_winter_lockdown, min, max));
        contact_dampings.push_back(senior_awareness(start_winter_lockdown, 0.0, 0.0));

        //relaxing of restrictions over christmas days
        auto xmas_date = epi::Date(2020, 12, 24);
        auto xmas      = epi::SimulationTime(epi::get_offset_in_days(xmas_date, start_date));
        contact_dampings.push_back(contacts_at_home(xmas, 0.0, 0.0));
        contact_dampings.push_back(home_office(xmas, 0.4, 0.5));
        contact_dampings.push_back(social_events(xmas, 0.4, 0.6));
        contact_dampings.push_back(physical_distancing_home_school(xmas, 0.0, 0.0));
        contact_dampings.push_back(physical_distancing_work_other(xmas, 0.4, 0.6));

        // after christmas
        auto after_xmas_date = epi::Date(2020, 12, 27);
        auto after_xmas      = epi::SimulationTime(epi::get_offset_in_days(after_xmas_date, start_date));
        contact_dampings.push_back(contacts_at_home(after_xmas, min, max));
        contact_dampings.push_back(home_office(after_xmas, 0.2, 0.3));
        contact_dampings.push_back(social_events(after_xmas, 0.6, 0.8));
        contact_dampings.push_back(physical_distancing_home_school(after_xmas, 0.2, 0.4));
        contact_dampings.push_back(physical_distancing_work_other(after_xmas, min, max));
    }

    //local dynamic NPIs
    auto& dynamic_npis        = params.get<epi::DynamicNPIsInfected>();
    auto dynamic_npi_dampings = std::vector<epi::DampingSampling>();
    dynamic_npi_dampings.push_back(contacts_at_home(epi::SimulationTime(0), 0.2, 0.4));
    dynamic_npi_dampings.push_back(school_closure(epi::SimulationTime(0), 1.0, 1.0)); //0.25 - 0.25 in autumn
    dynamic_npi_dampings.push_back(home_office(epi::SimulationTime(0), 0.2, 0.3));
    dynamic_npi_dampings.push_back(social_events(epi::SimulationTime(0), 0.2, 0.4));
    dynamic_npi_dampings.push_back(social_events_work(epi::SimulationTime(0), 0.0, 0.0));
    dynamic_npi_dampings.push_back(physical_distancing_home_school(epi::SimulationTime(0), 0.2, 0.4));
    dynamic_npi_dampings.push_back(physical_distancing_work_other(epi::SimulationTime(0), 0.2, 0.4));
    dynamic_npi_dampings.push_back(senior_awareness(epi::SimulationTime(0), 0.0, 0.0));
    dynamic_npis.set_interval(epi::SimulationTime(3.0));
    dynamic_npis.set_duration(epi::SimulationTime(14.0));
    dynamic_npis.set_base_value(100'000);
    dynamic_npis.set_threshold(10.0, dynamic_npi_dampings);

    //school holidays (holiday periods are set per node, see set_nodes)
    auto school_holiday_value = epi::UncertainValue();
    assign_uniform_distribution(school_holiday_value, 1.0, 1.0);
    contacts.get_school_holiday_damping() = epi::DampingSampling(
        school_holiday_value, epi::DampingLevel(5), epi::DampingType(int(Intervention::SchoolClosure)),
        epi::SimulationTime(0.0), {size_t(ContactLocation::School)}, group_weights_all);

    return epi::success();
}

epi::IOResult<void> set_nodes(const epi::SecirParams& params, epi::Date start_date, epi::Date end_date,
                              const fs::path& data_dir,
                              epi::Graph<epi::SecirModel, epi::MigrationParameters>& params_graph)
{
    namespace de = epi::regions::de;

    BOOST_OUTCOME_TRY(county_ids, epi::get_county_ids((data_dir / "pydata" / "Germany").string()));
    std::vector<epi::SecirModel> counties(county_ids.size(), epi::SecirModel(size_t(params.get_num_groups())));
    for (auto& county : counties) {
        county.parameters = params;
    }
    // auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 2.5);
    // auto scaling_factor_icu      = 1.0;
    // BOOST_OUTCOME_TRY(epi::read_population_data_county(counties, start_date, county_ids, scaling_factor_infected,
    //                                                    scaling_factor_icu, (data_dir / "pydata" / "Germany").string()));

    for (size_t county_idx = 0; county_idx < counties.size(); ++county_idx) {
        double nb_total_t0 = 10000, nb_exp_t0 = 2, nb_inf_t0 = 0, nb_car_t0 = 0, nb_hosp_t0 = 0, nb_icu_t0 = 0,
               nb_rec_t0 = 0, nb_dead_t0 = 0;

        nb_exp_t0 = (county_idx % 10 + 1) * 3;

        for (epi::AgeGroup i = 0; i < params.get_num_groups(); i++) {
            counties[county_idx].populations[{i, epi::InfectionState::Exposed}]      = nb_exp_t0;
            counties[county_idx].populations[{i, epi::InfectionState::Carrier}]      = nb_car_t0;
            counties[county_idx].populations[{i, epi::InfectionState::Infected}]     = nb_inf_t0;
            counties[county_idx].populations[{i, epi::InfectionState::Hospitalized}] = nb_hosp_t0;
            counties[county_idx].populations[{i, epi::InfectionState::ICU}]          = nb_icu_t0;
            counties[county_idx].populations[{i, epi::InfectionState::Recovered}]    = nb_rec_t0;
            counties[county_idx].populations[{i, epi::InfectionState::Dead}]         = nb_dead_t0;
            counties[county_idx].populations.set_difference_from_group_total<epi::AgeGroup>(
                {i, epi::InfectionState::Susceptible}, nb_total_t0);
        }

        //local parameters
        auto tnt_capacity = counties[county_idx].populations.get_total() * 7.5 / 100000.;
        assign_uniform_distribution(counties[county_idx].parameters.get<epi::TestAndTraceCapacity>(),
                                    0.8 * tnt_capacity, 1.2 * tnt_capacity);

        //holiday periods (damping set globally, see set_npis)
        auto holiday_periods =
            de::get_holidays(de::get_state_id(de::CountyId(county_ids[county_idx])), start_date, end_date);
        auto& contacts = counties[county_idx].parameters.get<epi::ContactPatterns>();
        contacts.get_school_holidays() =
            std::vector<std::pair<epi::SimulationTime, epi::SimulationTime>>(holiday_periods.size());
        std::transform(
            holiday_periods.begin(), holiday_periods.end(), contacts.get_school_holidays().begin(), [=](auto& period) {
                return std::make_pair(epi::SimulationTime(epi::get_offset_in_days(period.first, start_date)),
                                      epi::SimulationTime(epi::get_offset_in_days(period.second, start_date)));
            });

        //uncertainty in populations
        //TODO: do we need uncertainty in age groups as well?
        for (auto i = epi::AgeGroup(0); i < params.get_num_groups(); i++) {
            for (auto j = epi::Index<epi::InfectionState>(0); j < epi::InfectionState::Count; ++j) {
                auto& compartment_value = counties[county_idx].populations[{i, j}];
                assign_uniform_distribution(compartment_value, 0.9 * double(compartment_value),
                                            1.1 * double(compartment_value));
            }
        }

        params_graph.add_node(county_ids[county_idx], counties[county_idx]);
    }
    return epi::success();
}

epi::IOResult<void> set_edges(const fs::path& data_dir,
                              epi::Graph<epi::SecirModel, epi::MigrationParameters>& params_graph)
{
    //migration between nodes
    BOOST_OUTCOME_TRY(migration_data_commuter,
                      epi::read_mobility_plain((data_dir / "migration" / "commuter_migration_scaled.txt").string()));
    BOOST_OUTCOME_TRY(migration_data_twitter,
                      epi::read_mobility_plain((data_dir / "migration" / "twitter_scaled_1252.txt").string()));
    if (migration_data_commuter.rows() != params_graph.nodes().size() ||
        migration_data_commuter.cols() != params_graph.nodes().size() ||
        migration_data_twitter.rows() != params_graph.nodes().size() ||
        migration_data_twitter.cols() != params_graph.nodes().size()) {
        return epi::failure(epi::StatusCode::InvalidValue, "Migration matrices not the correct size.");
    }

    auto migrating_compartments = {epi::InfectionState::Susceptible, epi::InfectionState::Exposed,
                                   epi::InfectionState::Carrier, epi::InfectionState::Infected,
                                   epi::InfectionState::Recovered};
    for (size_t county_idx_i = 0; county_idx_i < params_graph.nodes().size(); ++county_idx_i) {
        for (size_t county_idx_j = 0; county_idx_j < params_graph.nodes().size(); ++county_idx_j) {
            auto& populations = params_graph.nodes()[county_idx_i].property.populations;
            //migration coefficients have the same number of components as the contact matrices.
            //so that the same NPIs/dampings can be used for both (e.g. more home office => fewer commuters)
            auto migration_coeffs = epi::MigrationCoefficientGroup(contact_locations.size(), populations.numel());

            //commuters
            auto working_population = 0.0;
            auto min_commuter_age   = epi::AgeGroup(2);
            auto max_commuter_age   = epi::AgeGroup(4); //this group is partially retired, only partially commutes
            for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                working_population += populations.get_group_total(age) * (age == max_commuter_age ? 0.33 : 1.0);
            }
            auto commuter_coeff_ij = migration_data_commuter(county_idx_i, county_idx_j) /
                                     working_population; //data is absolute numbers, we need relative
            for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                for (auto compartment : migrating_compartments) {
                    auto coeff_index = populations.get_flat_index({age, compartment});
                    migration_coeffs[size_t(ContactLocation::Work)].get_baseline()[coeff_index] =
                        commuter_coeff_ij * (age == max_commuter_age ? 0.33 : 1.0);
                }
            }
            //others
            auto total_population = populations.get_total();
            auto twitter_coeff    = migration_data_twitter(county_idx_i, county_idx_j) /
                                 total_population; //data is absolute numbers, we need relative
            for (auto age = epi::AgeGroup(0); age < populations.size<epi::AgeGroup>(); ++age) {
                for (auto compartment : migrating_compartments) {
                    auto coeff_idx = populations.get_flat_index({age, compartment});
                    migration_coeffs[size_t(ContactLocation::Other)].get_baseline()[coeff_idx] = twitter_coeff;
                }
            }

            //only add edges with migration above thresholds for performance
            //thresholds are chosen empirically so that more than 99% of migration is covered, approx. 1/3 of the edges
            if (commuter_coeff_ij > 4e-5 || twitter_coeff > 1e-5) {
                params_graph.add_edge(county_idx_i, county_idx_j, std::move(migration_coeffs));
            }
        }
    }

    return epi::success();
}

epi::IOResult<void> run(const fs::path& data_dir)
{
    const auto start_date   = epi::Date(2020, 12, 12);
    const auto start_day    = epi::get_day_in_year(start_date);
    const auto num_days_sim = 20.0;
    const auto end_date     = epi::offset_date_by_days(start_date, int(std::ceil(num_days_sim)));
    const auto num_runs     = 1;

    //global parameters
    const int num_age_groups = 6;
    epi::SecirParams params(num_age_groups);
    params.get<epi::StartDay>() = start_day;
    BOOST_OUTCOME_TRY(set_covid_parameters(params));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));
    BOOST_OUTCOME_TRY(set_npis(data_dir, start_date, end_date, params));

    //graph of counties with populations and local parameters
    //and migration between counties
    epi::Graph<epi::SecirModel, epi::MigrationParameters> params_graph;
    BOOST_OUTCOME_TRY(set_nodes(params, start_date, end_date, data_dir, params_graph));
    BOOST_OUTCOME_TRY(set_edges(data_dir, params_graph));

    std::cout << std::setprecision(16);
    //parameter study
    auto parameter_study  = epi::ParameterStudy<epi::SecirSimulation<>>{params_graph, 0.0, num_days_sim, 0.5, num_runs};
    auto ensemble_results = std::vector<std::vector<epi::TimeSeries<double>>>{};
    auto ensemble_params  = std::vector<std::vector<epi::SecirModel>>{};
    ensemble_results.reserve(size_t(num_runs));
    ensemble_params.reserve(size_t(num_runs));
    parameter_study.run([&ensemble_results, &ensemble_params](auto results_graph) {
        ensemble_results.push_back(epi::interpolate_simulation_result(results_graph));

        ensemble_params.emplace_back();
        ensemble_params.back().reserve(results_graph.nodes().size());
        std::transform(results_graph.nodes().begin(), results_graph.nodes().end(),
                       std::back_inserter(ensemble_params.back()), [](auto&& node) {
                           return node.property.get_simulation().get_model();
                       });
    });

    auto& params_out = ensemble_params[0][0].parameters;
    std::cout << "matrix at 0\n";
    std::cout << params_out.get<epi::ContactPatterns>().get_cont_freq_mat().get_matrix_at(0.0) << '\n';
    std::cout << "matrix at 20\n";
    std::cout << params_out.get<epi::ContactPatterns>().get_cont_freq_mat().get_matrix_at(20.0) << '\n';
    // for (size_t i = 0; i < contact_locations.size(); ++i) {
    //     std::cout << "matrix " << i << " at 0\n";
    //     std::cout << params_out.get<epi::ContactPatterns>().get_cont_freq_mat()[i].get_matrix_at(0.0) << '\n';
    // }
    // int i = 0;
    // for (auto&& d : params_out.get<epi::ContactPatterns>().get_cont_freq_mat()[1].get_dampings()) {
    //     std::cout << "\nDamping " << ++i << ":\n";
    //     PrintTo(d, &std::cout);
    // }

    for (auto& n : ensemble_results[0]) {
        for (auto i : {Eigen::Index(0), n.get_num_time_points() - 1}) {
            const auto& v = n.get_value(i);
            std::cout << v.transpose() << '\n';
        }
    }

    return epi::success();
}

int main(int argc, char** argv)
{
    epi::set_log_level(epi::LogLevel::warn);

    epi::thread_local_rng().seed({1993982109, 2219896052, 2121723912, 560069401, 4132934777, 1129445404});
    printf("Seeds: ");
    for (auto s : epi::thread_local_rng().get_seeds()) {
        printf("%u, ", s);
    }
    printf("\n");

    auto result = run(argv[1]);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
    }
}