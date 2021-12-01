/**
* Simulation application that was used to produce results of the following publication(s):
* M. J. Kühn et al, 2021: Assessment of effective mitigation and prediction of the spread of SARS-CoV-2 in Germany 
* using demographic information and spatial resolution
* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#include "epidemiology/secir/parameter_studies_vaccinated.h"
#include "epidemiology/utils/regions.h"
#include "epidemiology_io/secir_parameters_io.h"
#include "epidemiology_io/secir_result_io.h"
#include "epidemiology_io/mobility_io.h"
#include "boost/filesystem.hpp"
#include <cstdio>
#include <iomanip>

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
};

/**
 * Set a value and distribution of an UncertainValue.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution.
 * @param p uncertain value to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void assign_uniform_distribution(epi::UncertainValue& p, double min, double max)
{
    p = epi::UncertainValue(0.5 * (max + min));
    p.set_distribution(epi::ParameterDistributionUniform(min, max));
}

/**
 * Set a value and distribution of an array of UncertainValues.
 * Assigns average of min[i] and max[i] as a value and UNIFORM(min[i], max[i]) as a distribution for
 * each element i of the array.
 * @param array array of UncertainValues to set.
 * @param min minimum of distribution for each element of array.
 * @param max minimum of distribution for each element of array.
 */
template <size_t N>
void array_assign_uniform_distribution(epi::CustomIndexArray<epi::UncertainValue, epi::AgeGroup>& array,
                                       const double (&min)[N], const double (&max)[N])
{
    assert(N == array.numel());
    for (auto i = epi::AgeGroup(0); i < epi::AgeGroup(N); ++i) {
        assign_uniform_distribution(array[i], min[size_t(i)], max[size_t(i)]);
    }
}

/**
 * Set a value and distribution of an array of UncertainValues.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution to every element of the array.
 * @param array array of UncertainValues to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void array_assign_uniform_distribution(epi::CustomIndexArray<epi::UncertainValue, epi::AgeGroup>& array, double min,
                                       double max)
{
    for (auto i = epi::AgeGroup(0); i < array.size<epi::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min, max);
    }
}

/**
 * Set epidemiological parameters of Covid19.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
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
    array_assign_uniform_distribution(params.get<epi::SerialInterval>(), tserint_max, tserint_max);
    array_assign_uniform_distribution(params.get<epi::InfectiousTimeMild>(), t_inf_rec_max, t_inf_rec_max);
    array_assign_uniform_distribution(params.get<epi::HomeToHospitalizedTime>(), t_inf_hosp_max, t_inf_hosp_max);
    array_assign_uniform_distribution(params.get<epi::HospitalizedToHomeTime>(), t_hosp_rec_max, t_hosp_rec_max);
    array_assign_uniform_distribution(params.get<epi::HospitalizedToICUTime>(), t_hosp_icu_max, t_hosp_icu_max);
    array_assign_uniform_distribution(params.get<epi::ICUToHomeTime>(), t_icu_rec_max, t_icu_rec_max);
    array_assign_uniform_distribution(params.get<epi::ICUToDeathTime>(), t_icu_dead_max, t_icu_dead_max);

    //probabilities
    double fac_variant                   = 1.4;
    const double transmission_risk_min[] = {0.02 * fac_variant, 0.05 * fac_variant, 0.05 * fac_variant,
                                            0.05 * fac_variant, 0.08 * fac_variant, 0.10 * fac_variant};

    const double transmission_risk_max[] = {0.04 * fac_variant, 0.07 * fac_variant, 0.07 * fac_variant,
                                            0.07 * fac_variant, 0.10 * fac_variant, 0.15 * fac_variant};
    const double carr_infec_min          = 0.5;
    const double carr_infec_max          = 0.5;
    const double beta_low_incidenc_min   = 0.0; // beta (depends on incidence and test and trace capacity)
    const double beta_low_incidenc_max   = 0.2;
    const double beta_high_incidence_min = 0.4;
    const double beta_high_incidence_max = 0.5;
    const double prob_car_rec_min[]      = {0.2, 0.2, 0.15, 0.15, 0.15, 0.15}; // alpha
    const double prob_car_rec_max[]      = {0.3, 0.3, 0.25, 0.25, 0.25, 0.25};
    const double prob_inf_hosp_min[]     = {0.006, 0.006, 0.015, 0.049, 0.15, 0.20}; // rho
    const double prob_inf_hosp_max[]     = {0.009, 0.009, 0.023, 0.074, 0.18, 0.25};
    const double prob_hosp_icu_min[]     = {0.05, 0.05, 0.05, 0.10, 0.25, 0.35}; // theta
    const double prob_hosp_icu_max[]     = {0.10, 0.10, 0.10, 0.20, 0.35, 0.45};
    const double prob_icu_dead_min[]     = {0.00, 0.00, 0.10, 0.10, 0.30, 0.5}; // delta
    const double prob_icu_dead_max[]     = {0.10, 0.10, 0.18, 0.18, 0.50, 0.7};

    const double reduc_vacc_exp_min        = 0.75; //partial vacc
    const double reduc_vacc_exp_max        = 0.85;
    const double reduc_immune_exp_min      = 0.281; //full vacc or recov
    const double reduc_immune_exp_max      = 0.381;
    const double reduc_exp_inf_min         = 0.6;
    const double reduc_exp_inf_max         = 0.7;
    const double reduc_immune_exp_inf_min  = 0.193;
    const double reduc_immune_exp_inf_max  = 0.293;
    const double reduc_inf_hosp_min        = 0.05;
    const double reduc_inf_hosp_max        = 0.15;
    const double reduc_immune_inf_hosp_min = 0.041;
    const double reduc_immune_inf_hosp_max = 0.141;

    const double reduc_time = 0.5;

    array_assign_uniform_distribution(params.get<epi::InfectionProbabilityFromContact>(), transmission_risk_min,
                                      transmission_risk_min);
    array_assign_uniform_distribution(params.get<epi::RelativeCarrierInfectability>(), carr_infec_max, carr_infec_max);
    array_assign_uniform_distribution(params.get<epi::RiskOfInfectionFromSympomatic>(), beta_low_incidenc_max,
                                      beta_low_incidenc_max);
    array_assign_uniform_distribution(params.get<epi::MaxRiskOfInfectionFromSympomatic>(), beta_high_incidence_max,
                                      beta_high_incidence_max);
    array_assign_uniform_distribution(params.get<epi::AsymptoticCasesPerInfectious>(), prob_car_rec_max,
                                      prob_car_rec_max);
    array_assign_uniform_distribution(params.get<epi::HospitalizedCasesPerInfectious>(), prob_inf_hosp_max,
                                      prob_inf_hosp_max);
    array_assign_uniform_distribution(params.get<epi::ICUCasesPerHospitalized>(), prob_hosp_icu_max, prob_hosp_icu_max);
    array_assign_uniform_distribution(params.get<epi::DeathsPerHospitalized>(), prob_icu_dead_max, prob_icu_dead_max);

    array_assign_uniform_distribution(params.get<epi::ReducVaccExp>(), reduc_vacc_exp_min, reduc_vacc_exp_max);
    array_assign_uniform_distribution(params.get<epi::ReducImmuneExp>(), reduc_immune_exp_min, reduc_immune_exp_max);
    array_assign_uniform_distribution(params.get<epi::ReducExpInf>(), reduc_exp_inf_min, reduc_exp_inf_max);
    array_assign_uniform_distribution(params.get<epi::ReducImmuneExpInf>(), reduc_immune_exp_inf_min,
                                      reduc_immune_exp_inf_max);
    array_assign_uniform_distribution(params.get<epi::ReducInfHosp>(), reduc_inf_hosp_min, reduc_inf_hosp_max);
    array_assign_uniform_distribution(params.get<epi::ReducImmuneInfHosp>(), reduc_immune_inf_hosp_min,
                                      reduc_immune_inf_hosp_max);
    array_assign_uniform_distribution(params.get<epi::ReducTime>(), reduc_time,
                                      reduc_time);

    //sasonality
    const double seasonality_min = 0.1;
    const double seasonality_max = 0.1;

    assign_uniform_distribution(params.get<epi::Seasonality>(), seasonality_max, seasonality_max);

    return epi::success();
}

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
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::MatrixXd::Zero(6,6);;
    }
    params.get<epi::ContactPatterns>() = epi::UncertainContactMatrix(contact_matrices);

    return epi::success();
}

/**
 * Set NPIs.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param params Object that the NPIs will be added to.
 * @returns Currently generates no errors.
 */
epi::IOResult<void> set_npis(epi::Date start_date, epi::Date end_date, epi::SecirParams& params, bool late, bool masks)
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
        return epi::DampingSampling(v, epi::DampingLevel(int(InterventionLevel::Main)),
                                    epi::DampingType(int(Intervention::Home)), t, {size_t(ContactLocation::Home)},
                                    group_weights_all);
    };
    auto school_closure = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(int(InterventionLevel::Main)),
                                    epi::DampingType(int(Intervention::SchoolClosure)), t,
                                    {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto home_office = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(int(InterventionLevel::Main)),
                                    epi::DampingType(int(Intervention::HomeOffice)), t, {size_t(ContactLocation::Work)},
                                    group_weights_all);
    };
    auto social_events = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(int(InterventionLevel::Main)),
                                    epi::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), t,
                                    {size_t(ContactLocation::Other)}, group_weights_all);
    };
    auto social_events_work = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(int(InterventionLevel::Main)),
                                    epi::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), t,
                                    {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto physical_distancing_home = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                    epi::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
                                    {size_t(ContactLocation::Home)},
                                    group_weights_all);
    };
    auto physical_distancing_work_other_school = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                    epi::DampingType(int(Intervention::PhysicalDistanceAndMasks)), t,
                                    {size_t(ContactLocation::Work), size_t(ContactLocation::Other), size_t(ContactLocation::School)}, group_weights_all);
    };
    auto senior_awareness = [=](auto t, auto min, auto max) {
        auto v = epi::UncertainValue();
        assign_uniform_distribution(v, min, max);
        return epi::DampingSampling(v, epi::DampingLevel(int(InterventionLevel::SeniorAwareness)),
                                    epi::DampingType(int(Intervention::SeniorAwareness)), t,
                                    {size_t(ContactLocation::Home), size_t(ContactLocation::Other)},
                                    group_weights_seniors);
    };

    //INITIAL
    auto t0 = epi::SimulationTime(0);
    contact_dampings.push_back(contacts_at_home(t0, 0.0, 0.1));
    contact_dampings.push_back(school_closure(t0, 0.0, 0.1));
    contact_dampings.push_back(home_office(t0, 0.2, 0.4));
    contact_dampings.push_back(social_events(t0, 0.2, 0.4));
    contact_dampings.push_back(social_events_work(t0, 0.0, 0.0));
    contact_dampings.push_back(physical_distancing_home(t0, 0.0, 0.0));
    contact_dampings.push_back(physical_distancing_work_other_school(t0, 0.1, 0.3));
    contact_dampings.push_back(senior_awareness(t0, 0.2, 0.4));

    //5.11. Vorwarnstufe
    //max 10 Pers.
    auto date_vorwarnstufe = epi::Date(2021, 11, 5);
    auto t_vorwarnstufe = epi::SimulationTime(epi::get_offset_in_days(date_vorwarnstufe, start_date));
    contact_dampings.push_back(contacts_at_home(t_vorwarnstufe, 0.2, 0.4));
    // contact_dampings.push_back(school_closure(t_vorwarnstufe, 0.0, 0.1));
    contact_dampings.push_back(home_office(t_vorwarnstufe, 0.3, 0.5));
    contact_dampings.push_back(social_events(t_vorwarnstufe, 0.4, 0.6));
    // contact_dampings.push_back(social_events_work(t_vorwarnstufe, 0.0, 0.0));
    // contact_dampings.push_back(physical_distancing_home(t_vorwarnstufe, 0.0, 0.0));
    // contact_dampings.push_back(physical_distancing_work_other_school(t_vorwarnstufe, 0.1, 0.3));
    // contact_dampings.push_back(senior_awareness(t_vorwarnstufe, 0.2, 0.4));

    //8.11. Neue Verordnung
    //max 10 Pers.
    //2G Veranstaltungen
    //FFP2 ÖPNV
    //Testpflicht für Pflegepersonal
    auto date_verordnung = epi::Date(2021, 11, 8);
    auto t_verordnung = epi::SimulationTime(epi::get_offset_in_days(date_verordnung, start_date));
    // contact_dampings.push_back(contacts_at_home(t_verordnung, 0.2, 0.4));
    // contact_dampings.push_back(school_closure(t_verordnung, 0.0, 0.1));
    // contact_dampings.push_back(home_office(t_verordnung, 0.3, 0.5));
    contact_dampings.push_back(social_events(t_verordnung, 0.5, 0.7));
    // contact_dampings.push_back(social_events_work(t_verordnung, 0.0, 0.0));
    // contact_dampings.push_back(physical_distancing_home(t_verordnung, 0.0, 0.0));
    contact_dampings.push_back(physical_distancing_work_other_school(t_verordnung, 0.3, 0.5));
    contact_dampings.push_back(senior_awareness(t_verordnung, 0.3, 0.5));

    //19.11. überlastungsstufe
    //1+1
    //3G->2G
    //teilweise Schulschließungen?
    auto date_uberlastung = epi::Date(2021, 11, 19);
    auto t_uberlastung = epi::SimulationTime(epi::get_offset_in_days(date_uberlastung, start_date));
    contact_dampings.push_back(contacts_at_home(t_uberlastung, 0.3, 0.5));
    contact_dampings.push_back(school_closure(t_uberlastung, 0.1, 0.3));
    // contact_dampings.push_back(home_office(t_uberlastung, 0.3, 0.5));
    // contact_dampings.push_back(social_events(t_uberlastung, 0.5, 0.7));
    // contact_dampings.push_back(social_events_work(t_uberlastung, 0.0, 0.0));
    // contact_dampings.push_back(physical_distancing_home(t_uberlastung, 0.0, 0.0));
    // contact_dampings.push_back(physical_distancing_work_other_school(t_uberlastung, 0.3, 0.5));
    // contact_dampings.push_back(senior_awareness(t_uberlastung, 0.3, 0.5));

    //22.11. Notfallveordnung (PDF)
    //1+1
    //Clubs, Sport, Großveranstaltungen, u.a. Social Events geschlossen, max 10 Pers.
    //Gastronomie, Handle (außer Grundvers.) 2G, sonst FFP2 und qm-Vorgabe
    //teilweise Schulschließungen
    auto date_notfallvo = epi::Date(2021, 11, 22);
    auto t_notfallvo = epi::SimulationTime(epi::get_offset_in_days(date_notfallvo, start_date));
    contact_dampings.push_back(contacts_at_home(t_notfallvo, 0.4, 0.6));
    contact_dampings.push_back(school_closure(t_notfallvo, 0.3, 0.5));
    contact_dampings.push_back(home_office(t_notfallvo, 0.5, 0.7));
    contact_dampings.push_back(social_events(t_notfallvo, 0.6, 0.8));
    // contact_dampings.push_back(social_events_work(t_notfallvo, 0.0, 0.0));
    // contact_dampings.push_back(physical_distancing_home(t_notfallvo, 0.0, 0.0));
    // contact_dampings.push_back(physical_distancing_work_other_school(t_notfallvo, 0.3, 0.5));
    // contact_dampings.push_back(senior_awareness(t_notfallvo, 0.3, 0.5));

    //23.11. InfSG
    //Homeofficepflicht
    //3G Arbeit und ÖPNV
    auto date_infsg = epi::Date(2021, 11, 23);
    auto t_infsg = epi::SimulationTime(epi::get_offset_in_days(date_infsg, start_date));
    contact_dampings.push_back(contacts_at_home(t_infsg, 0.5, 0.7));
    // contact_dampings.push_back(school_closure(t_infsg, 0.3, 0.5));
    contact_dampings.push_back(home_office(t_infsg, 0.6, 0.8));
    // contact_dampings.push_back(social_events(t_infsg, 0.6, 0.8));
    // contact_dampings.push_back(social_events_work(t_infsg, 0.0, 0.0));
    // contact_dampings.push_back(physical_distancing_home(t_infsg, 0.0, 0.0));
    // contact_dampings.push_back(physical_distancing_work_other_school(t_infsg, 0.3, 0.5));
    // contact_dampings.push_back(senior_awareness(t_infsg, 0.3, 0.5));

    //26.11. vermehrt Schuschließungen   
    auto date_schools = epi::Date(2021, 11, 26);
    auto t_schools = epi::SimulationTime(epi::get_offset_in_days(date_schools, start_date));
    // contact_dampings.push_back(contacts_at_home(t_schools, 0.5, 0.7));
    contact_dampings.push_back(school_closure(t_schools, 0.4, 0.6));
    // contact_dampings.push_back(home_office(t_schools, 0.6, 0.8));
    // contact_dampings.push_back(social_events(t_schools, 0.6, 0.8));
    // contact_dampings.push_back(social_events_work(t_schools, 0.0, 0.0));
    // contact_dampings.push_back(physical_distancing_home(t_schools, 0.0, 0.0));
    // contact_dampings.push_back(physical_distancing_work_other_school(t_schools, 0.3, 0.5));
    // contact_dampings.push_back(senior_awareness(t_schools, 0.3, 0.5));

    //29.11. 
    //Eingeschränkter Schulbetrieb (PDF)
    auto date_schools2 = epi::Date(2021, 11, 29);
    auto t_schools2 = epi::SimulationTime(epi::get_offset_in_days(date_schools2, start_date));
    // contact_dampings.push_back(contacts_at_home(t_schools2, 0.5, 0.7));
    contact_dampings.push_back(school_closure(t_schools2, 0.5, 0.7));
    // contact_dampings.push_back(home_office(t_schools2, 0.6, 0.8));
    // contact_dampings.push_back(social_events(t_schools2, 0.6, 0.8));
    // contact_dampings.push_back(social_events_work(t_schools2, 0.0, 0.0));
    // contact_dampings.push_back(physical_distancing_home(t_schools2, 0.0, 0.0));
    // contact_dampings.push_back(physical_distancing_work_other_school(t_schools2, 0.3, 0.5));
    // contact_dampings.push_back(senior_awareness(t_schools2, 0.3, 0.5));

    //Dynamic NPIS
    // double reduced_effect = 0.1;
    // narrow                = 0.0;
    //local dynamic NPIs
    // auto& dynamic_npis        = params.get<epi::DynamicNPIsInfected>();
    // auto dynamic_npi_dampings = std::vector<epi::DampingSampling>();
    // dynamic_npi_dampings.push_back(
    //     contacts_at_home(epi::SimulationTime(0), 0.3 + narrow - reduced_effect, 0.3 - narrow - reduced_effect));
    // dynamic_npi_dampings.push_back(school_closure(epi::SimulationTime(0), 0.4 + narrow - reduced_effect,
    //                                               0.4 - narrow - reduced_effect)); //0.25 - 0.25 in autumn
    // dynamic_npi_dampings.push_back(
    //     home_office(epi::SimulationTime(0), 0.2 + narrow - reduced_effect, 0.2 - narrow - reduced_effect));
    // dynamic_npi_dampings.push_back(
    //     social_events(epi::SimulationTime(0), 0.6 + narrow - reduced_effect, 0.6 - narrow - reduced_effect));
    // dynamic_npi_dampings.push_back(social_events_work(epi::SimulationTime(0), 0.0, 0.0));
    // dynamic_npi_dampings.push_back(physical_distancing_home_school(
    //     epi::SimulationTime(0), 0.1 + narrow - reduced_effect, 0.1 - narrow - reduced_effect));
    // dynamic_npi_dampings.push_back(physical_distancing_work_other(epi::SimulationTime(0), 0.4 + narrow - reduced_effect,
    //                                                               0.4 - narrow - reduced_effect));
    // dynamic_npi_dampings.push_back(senior_awareness(epi::SimulationTime(0), 0.0, 0.0));

    // auto dynamic_npi_dampings2 = std::vector<epi::DampingSampling>();
    // dynamic_npi_dampings2.push_back(
    //     contacts_at_home(epi::SimulationTime(0), 0.5 + narrow - reduced_effect, 0.5 - narrow - reduced_effect));
    // dynamic_npi_dampings2.push_back(school_closure(epi::SimulationTime(0), 0.4 + narrow - reduced_effect,
    //                                                0.4 - narrow - reduced_effect)); //0.25 - 0.25 in autumn
    // dynamic_npi_dampings2.push_back(
    //     home_office(epi::SimulationTime(0), 0.3 + narrow - reduced_effect, 0.3 - narrow - reduced_effect));
    // dynamic_npi_dampings2.push_back(
    //     social_events(epi::SimulationTime(0), 0.6 + narrow - reduced_effect, 0.6 - narrow - reduced_effect));
    // dynamic_npi_dampings2.push_back(social_events_work(epi::SimulationTime(0), 0.0, 0.0));
    // dynamic_npi_dampings2.push_back(physical_distancing_home_school(
    //     epi::SimulationTime(0), 0.1 + narrow - reduced_effect, 0.1 - narrow - reduced_effect));
    // dynamic_npi_dampings2.push_back(physical_distancing_work_other(
    //     epi::SimulationTime(0), 0.4 + narrow - reduced_effect, 0.4 - narrow - reduced_effect));
    // dynamic_npi_dampings2.push_back(senior_awareness(epi::SimulationTime(0), 0.0, 0.0));

    // dynamic_npis.set_interval(epi::SimulationTime(1.0));
    // dynamic_npis.set_duration(epi::SimulationTime(14.0));
    // dynamic_npis.set_base_value(100'000);
    // dynamic_npis.set_threshold(35.0, dynamic_npi_dampings);
    // dynamic_npis.set_threshold(100.0, dynamic_npi_dampings2);

    //school holidays (holiday periods are set per node, see set_nodes)
    auto school_holiday_value = epi::UncertainValue();
    assign_uniform_distribution(school_holiday_value, 1.0, 1.0);
    contacts.get_school_holiday_damping() =
        epi::DampingSampling(school_holiday_value, epi::DampingLevel(int(InterventionLevel::Holidays)),
                             epi::DampingType(int(Intervention::SchoolClosure)), epi::SimulationTime(0.0),
                             {size_t(ContactLocation::School)}, group_weights_all);

    return epi::success();
}

/**
 * Set synthetic population data for testing.
 * Same total populaton but different spread of infection in each county.
 * @param counties parameters for each county.
 */
void set_synthetic_population_data(std::vector<epi::SecirModelV>& counties)
{
    for (size_t county_idx = 0; county_idx < counties.size(); ++county_idx) {
        double nb_total_t0 = 10000, nb_exp_t0 = 10, nb_inf_t0 = 10, nb_car_t0 = 10, nb_hosp_t0 = 10, nb_icu_t0 = 10,
               nb_rec_t0 = 100, nb_dead_t0 = 10;

        //nb_exp_t0 = (county_idx % 10 + 1) * 3;

        for (epi::AgeGroup i = 0; i < counties[county_idx].parameters.get_num_groups(); i++) {
            counties[county_idx].populations[{i, epi::InfectionStateV::Exposed}]        = nb_exp_t0;
            counties[county_idx].populations[{i, epi::InfectionStateV::Carrier}]        = nb_car_t0;
            counties[county_idx].populations[{i, epi::InfectionStateV::Infected}]       = nb_inf_t0;
            counties[county_idx].populations[{i, epi::InfectionStateV::Hospitalized}]   = nb_hosp_t0;
            counties[county_idx].populations[{i, epi::InfectionStateV::ICU}]            = nb_icu_t0;
            counties[county_idx].populations[{i, epi::InfectionStateV::Recovered}]      = nb_rec_t0;
            counties[county_idx].populations[{i, epi::InfectionStateV::Dead}]           = nb_dead_t0;
            counties[county_idx].populations[{i, epi::InfectionStateV::ExposedV1}]      = nb_exp_t0;
            counties[county_idx].populations[{i, epi::InfectionStateV::CarrierV1}]      = nb_car_t0;
            counties[county_idx].populations[{i, epi::InfectionStateV::InfectedV1}]     = nb_inf_t0;
            counties[county_idx].populations[{i, epi::InfectionStateV::HospitalizedV1}] = nb_hosp_t0;
            counties[county_idx].populations[{i, epi::InfectionStateV::ICUV1}]          = nb_icu_t0;
            counties[county_idx].populations[{i, epi::InfectionStateV::SusceptibleV1}]  = 500;
            counties[county_idx].populations.set_difference_from_group_total<epi::AgeGroup>(
                {i, epi::InfectionStateV::Susceptible}, nb_total_t0);
        }
    }
}

enum class RegionLevel
{
    County,
    State
};

struct RegionSpec
{
    RegionLevel level;
    int id;
};

/**
 * Adds county nodes to graph.
 * Reads list counties and populations from files in the data directory. 
 * @param params Parameters that are shared between all nodes.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param data_dir data directory.
 * @param params_graph graph object that the nodes will be added to.
 * @returns any io errors that happen during reading of the files.
 */
epi::IOResult<void> set_nodes(const epi::SecirParams& params, RegionSpec region_spec, epi::Date start_date,
                              epi::Date end_date, const fs::path& data_dir,
                              epi::Graph<epi::SecirModelV, epi::MigrationParameters>& params_graph)
{
    namespace de = epi::regions::de;

    BOOST_OUTCOME_TRY(county_ids, epi::get_county_ids((data_dir / "pydata" / "Germany").string()));

    auto county_ids_new_end = std::remove_if(county_ids.begin(), county_ids.end(), [region_spec](auto&& county_id) {
        return (region_spec.level != RegionLevel::County || county_id != region_spec.id) &&
               (region_spec.level != RegionLevel::State || de::get_state_id(de::CountyId(county_id)) != de::StateId(region_spec.id));
    });
    county_ids.erase(county_ids_new_end, county_ids.end());
    
    std::vector<epi::SecirModelV> counties(county_ids.size(), epi::SecirModelV(int(size_t(params.get_num_groups()))));
    for (auto& county : counties) {
        county.parameters = params;
    }

    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 1.0);
    auto scaling_factor_icu      = 1.0;
    BOOST_OUTCOME_TRY((epi::read_population_data_county<epi::SecirModelV, epi::InfectionStateV>(
        counties, start_date, county_ids, scaling_factor_infected, scaling_factor_icu,
        (data_dir / "pydata" / "Germany").string())));
    BOOST_OUTCOME_TRY(epi::read_vaccine_data(counties, start_date, county_ids,
                                             (data_dir / "pydata" / "Germany").string(),
                                             epi::get_offset_in_days(end_date, start_date)));
    //set_synthetic_population_data(counties);

    epi::SecirModelV region(int(size_t(params.get_num_groups())));
    region.parameters = params;
    region.populations.array() = 0.0;
    for (auto& county : counties) {
        region.populations.array() += county.populations.array();
    }
    for (auto age = epi::AgeGroup(0); age < region.parameters.get_num_groups(); ++age) {
        region.parameters.get<epi::DailyFirstVaccination>()[age] = std::vector<double>(counties[0].parameters.get<epi::DailyFirstVaccination>()[age].size(), 0.0);
        region.parameters.get<epi::DailyFullVaccination>()[age] = std::vector<double>(counties[0].parameters.get<epi::DailyFullVaccination>()[age].size(), 0.0);
        for (auto& county : counties) {
            for (auto t = size_t(0); t < region.parameters.get<epi::DailyFirstVaccination>()[age].size(); ++t) {
                region.parameters.get<epi::DailyFirstVaccination>()[age][t] += county.parameters.get<epi::DailyFirstVaccination>()[age][t];
            }
            for (auto t = size_t(0); t < region.parameters.get<epi::DailyFullVaccination>()[age].size(); ++t) {
                region.parameters.get<epi::DailyFullVaccination>()[age][t] += county.parameters.get<epi::DailyFullVaccination>()[age][t];
            }
        }
    }

    //local parameters
    auto tnt_capacity = region.populations.get_total() * 1.43 / 100000.;
    assign_uniform_distribution(region.parameters.get<epi::TestAndTraceCapacity>(),
                                1.0 * tnt_capacity, 1.0 * tnt_capacity);

    //holiday periods (damping set globally, see set_npis)
    auto holiday_periods = de::get_holidays(region_spec.level == RegionLevel::County ? de::get_state_id(de::CountyId(region_spec.id)) : de::StateId(region_spec.id), start_date, end_date);
    auto& contacts = region.parameters.get<epi::ContactPatterns>();
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
        for (auto j = epi::Index<epi::InfectionStateV>(0); j < epi::InfectionStateV::Count; ++j) {
            auto& compartment_value = region.populations[{i, j}];
            assign_uniform_distribution(compartment_value, 1.0 * double(compartment_value),
                                        1.0 * double(compartment_value));
        }
    }

    params_graph.add_node(region_spec.id, region);
    return epi::success();
}

/**
 * Adds edges to graph.
 * Edges represent commuting and other migration between counties.
 * Reads migration from files in the data directory.
 * @param data_dir data directory.
 * @param params_graph graph object that the nodes will be added to.
 * @returns any io errors that happen during reading of the files.
 */
epi::IOResult<void> set_edges(const fs::path& data_dir,
                              epi::Graph<epi::SecirModelV, epi::MigrationParameters>& params_graph)
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

    auto migrating_compartments = {epi::InfectionStateV::Susceptible, epi::InfectionStateV::Exposed,
                                   epi::InfectionStateV::Carrier,     epi::InfectionStateV::Infected,
                                   epi::InfectionStateV::Recovered,   epi::InfectionStateV::SusceptibleV1,
                                   epi::InfectionStateV::ExposedV1,   epi::InfectionStateV::CarrierV1,
                                   epi::InfectionStateV::InfectedV1,  epi::InfectionStateV::Recovered};
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

/**
 * Create the input graph for the parameter study.
 * Reads files from the data directory.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param data_dir data directory.
 * @returns created graph or any io errors that happen during reading of the files.
 */
epi::IOResult<epi::Graph<epi::SecirModelV, epi::MigrationParameters>>
create_graph(RegionSpec region_spec, epi::Date start_date, epi::Date end_date, const fs::path& data_dir, bool late, bool masks)
{
    const auto start_day = epi::get_day_in_year(start_date);
    int start_summer;
    if (late) {
        start_summer = epi::get_day_in_year(epi::Date(2021, 8, 1));
    }
    else {
        start_summer = epi::get_day_in_year(epi::Date(2021, 7, 1));
    }

    //global parameters
    const int num_age_groups = 6;
    epi::SecirParams params(num_age_groups);
    params.get<epi::StartDay>()    = start_day;
    params.get<epi::StartSummer>() = start_summer;
    params.get<epi::AvgSecondVaccDay>().array().setConstant(epi::get_offset_in_days(epi::Date(2021, 7, 1), start_date));
    // params.get<epi::AvgSecondVaccDay>()[epi::AgeGroup(0)] = epi::get_offset_in_days(epi::Date(2021, 9, 1), start_date);
    // params.get<epi::AvgSecondVaccDay>()[epi::AgeGroup(1)] = epi::get_offset_in_days(epi::Date(2021, 9, 1), start_date);
    // params.get<epi::AvgSecondVaccDay>()[epi::AgeGroup(2)] = epi::get_offset_in_days(epi::Date(2021, 7, 1), start_date);
    // params.get<epi::AvgSecondVaccDay>()[epi::AgeGroup(3)] = epi::get_offset_in_days(epi::Date(2021, 7, 1), start_date);
    // params.get<epi::AvgSecondVaccDay>()[epi::AgeGroup(4)] = epi::get_offset_in_days(epi::Date(2021, 7, 1), start_date);
    // params.get<epi::AvgSecondVaccDay>()[epi::AgeGroup(5)] = epi::get_offset_in_days(epi::Date(2021, 10, 1), start_date);
    BOOST_OUTCOME_TRY(set_covid_parameters(params));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));
    BOOST_OUTCOME_TRY(set_npis(start_date, end_date, params, late, masks));

    //graph of counties with populations and local parameters
    //and migration between counties
    epi::Graph<epi::SecirModelV, epi::MigrationParameters> params_graph;
    BOOST_OUTCOME_TRY(set_nodes(params, region_spec, start_date, end_date, data_dir, params_graph));
    //BOOST_OUTCOME_TRY(set_edges(data_dir, params_graph));

    return epi::success(params_graph);
}

/**
 * Load the input graph for the parameter study that was previously saved.
 * @param save_dir directory where the graph was saved.
 * @returns created graph or any io errors that happen during reading of the files.
 */
epi::IOResult<epi::Graph<epi::SecirModelV, epi::MigrationParameters>> load_graph(const fs::path& save_dir)
{
    return epi::read_graph<epi::SecirModelV>(save_dir.string());
}

/**
 * Save the input graph for the parameter study.
 * @param save_dir directory where the graph will be saved.
 * @returns any io errors that happen during writing of the files.
 */
epi::IOResult<void> save_graph(const epi::Graph<epi::SecirModelV, epi::MigrationParameters>& params_graph,
                               const fs::path& save_dir)
{
    return epi::write_graph(params_graph, save_dir.string());
}

/**
 * Create an unconnected graph.
 * Can be used to save space on disk when writing parameters if the edges are not required.
 * @param params parameters for each county node.
 * @param county_ids id of each county node.
 * @return graph with county nodes but no edges.
 */
auto make_graph_no_edges(const std::vector<epi::SecirModelV>& params, const std::vector<int>& county_ids)
{
    //make a graph without edges for writing to file
    auto graph = epi::Graph<epi::SecirModelV, epi::MigrationParameters>();
    for (auto i = size_t(0); i < county_ids.size(); ++i) {
        graph.add_node(county_ids[i], params[i]);
    }
    return graph;
}

/**
 * Save the result of a single parameter study run.
 * Creates a new subdirectory for this run.
 * @param result result of the simulation run.
 * @param params parameters used for the simulation run.
 * @param county_ids ids of the county nodes.
 * @param result_dir top level directory for all results of the parameter study.
 * @param run_idx index of the run.
 * @return any io errors that occur during writing of the files.
 */
epi::IOResult<void> save_result(const std::vector<epi::TimeSeries<double>>& result,
                                const std::vector<epi::SecirModelV>& params, const std::vector<int>& county_ids,
                                const fs::path& result_dir, size_t run_idx)
{
    auto result_dir_run = result_dir / ("run" + std::to_string(run_idx));
    BOOST_OUTCOME_TRY(epi::create_directory(result_dir_run.string()));
    BOOST_OUTCOME_TRY(
        (epi::save_result<epi::InfectionStateV>(result, county_ids, (result_dir_run / "Result.h5").string())));
    BOOST_OUTCOME_TRY(
        epi::write_graph(make_graph_no_edges(params, county_ids), result_dir_run.string(), epi::IOF_OmitDistributions));
    return epi::success();
}

/**
 * Save the results of a parameter study.
 * Stores different percentiles and sums of the results and parameters. 
 * @param ensemble_results result of each simulation run.
 * @param ensemble_params parameters used for each simulation run.
 * @param county_ids ids of the county nodes.
 * @param result_dir top level directory for all results of the parameter study.
 * @return any io errors that occur during writing of the files.
 */
epi::IOResult<void> save_results(const std::vector<std::vector<epi::TimeSeries<double>>>& ensemble_results,
                                 const std::vector<std::vector<epi::SecirModelV>>& ensemble_params,
                                 const std::vector<int>& county_ids, const fs::path& result_dir)
{
    //save results and sum of results over nodes
    auto ensemble_result_sum = epi::sum_nodes(ensemble_results);
    for (size_t i = 0; i < ensemble_result_sum.size(); ++i) {
        BOOST_OUTCOME_TRY((epi::save_result<epi::InfectionStateV>(
            ensemble_result_sum[i], {0}, (result_dir / ("results_run" + std::to_string(i) + "_sum.h5")).string())));
        BOOST_OUTCOME_TRY((epi::save_result<epi::InfectionStateV>(
            ensemble_results[i], county_ids, (result_dir / ("results_run" + std::to_string(i) + ".h5")).string())));
    }

    //make directories for percentiles
    auto result_dir_p05 = result_dir / "p05";
    auto result_dir_p25 = result_dir / "p25";
    auto result_dir_p50 = result_dir / "p50";
    auto result_dir_p75 = result_dir / "p75";
    auto result_dir_p95 = result_dir / "p95";
    BOOST_OUTCOME_TRY(epi::create_directory(result_dir_p05.string()));
    BOOST_OUTCOME_TRY(epi::create_directory(result_dir_p25.string()));
    BOOST_OUTCOME_TRY(epi::create_directory(result_dir_p50.string()));
    BOOST_OUTCOME_TRY(epi::create_directory(result_dir_p75.string()));
    BOOST_OUTCOME_TRY(epi::create_directory(result_dir_p95.string()));

    //save percentiles of results, summed over nodes
    {
        auto ensemble_results_sum_p05 = epi::ensemble_percentile(ensemble_result_sum, 0.05);
        auto ensemble_results_sum_p25 = epi::ensemble_percentile(ensemble_result_sum, 0.25);
        auto ensemble_results_sum_p50 = epi::ensemble_percentile(ensemble_result_sum, 0.50);
        auto ensemble_results_sum_p75 = epi::ensemble_percentile(ensemble_result_sum, 0.75);
        auto ensemble_results_sum_p95 = epi::ensemble_percentile(ensemble_result_sum, 0.95);

        BOOST_OUTCOME_TRY((epi::save_result<epi::InfectionStateV>(ensemble_results_sum_p05, {0},
                                                                  (result_dir_p05 / "Results_sum.h5").string())));
        BOOST_OUTCOME_TRY((epi::save_result<epi::InfectionStateV>(ensemble_results_sum_p25, {0},
                                                                  (result_dir_p25 / "Results_sum.h5").string())));
        BOOST_OUTCOME_TRY((epi::save_result<epi::InfectionStateV>(ensemble_results_sum_p50, {0},
                                                                  (result_dir_p50 / "Results_sum.h5").string())));
        BOOST_OUTCOME_TRY((epi::save_result<epi::InfectionStateV>(ensemble_results_sum_p75, {0},
                                                                  (result_dir_p75 / "Results_sum.h5").string())));
        BOOST_OUTCOME_TRY((epi::save_result<epi::InfectionStateV>(ensemble_results_sum_p95, {0},
                                                                  (result_dir_p95 / "Results_sum.h5").string())));
    }

    //save percentiles of results
    {
        auto ensemble_results_p05 = epi::ensemble_percentile(ensemble_results, 0.05);
        auto ensemble_results_p25 = epi::ensemble_percentile(ensemble_results, 0.25);
        auto ensemble_results_p50 = epi::ensemble_percentile(ensemble_results, 0.50);
        auto ensemble_results_p75 = epi::ensemble_percentile(ensemble_results, 0.75);
        auto ensemble_results_p95 = epi::ensemble_percentile(ensemble_results, 0.95);

        BOOST_OUTCOME_TRY((epi::save_result<epi::InfectionStateV>(ensemble_results_p05, county_ids,
                                                                  (result_dir_p05 / "Results.h5").string())));
        BOOST_OUTCOME_TRY((epi::save_result<epi::InfectionStateV>(ensemble_results_p25, county_ids,
                                                                  (result_dir_p25 / "Results.h5").string())));
        BOOST_OUTCOME_TRY((epi::save_result<epi::InfectionStateV>(ensemble_results_p50, county_ids,
                                                                  (result_dir_p50 / "Results.h5").string())));
        BOOST_OUTCOME_TRY((epi::save_result<epi::InfectionStateV>(ensemble_results_p75, county_ids,
                                                                  (result_dir_p75 / "Results.h5").string())));
        BOOST_OUTCOME_TRY((epi::save_result<epi::InfectionStateV>(ensemble_results_p95, county_ids,
                                                                  (result_dir_p95 / "Results.h5").string())));
    }

    //save percentiles of parameters
    {
        auto ensemble_params_p05 =
            epi::ensemble_params_percentile<epi::SecirModelV, epi::InfectionStateV>(ensemble_params, 0.05);
        auto ensemble_params_p25 =
            epi::ensemble_params_percentile<epi::SecirModelV, epi::InfectionStateV>(ensemble_params, 0.25);
        auto ensemble_params_p50 =
            epi::ensemble_params_percentile<epi::SecirModelV, epi::InfectionStateV>(ensemble_params, 0.50);
        auto ensemble_params_p75 =
            epi::ensemble_params_percentile<epi::SecirModelV, epi::InfectionStateV>(ensemble_params, 0.75);
        auto ensemble_params_p95 =
            epi::ensemble_params_percentile<epi::SecirModelV, epi::InfectionStateV>(ensemble_params, 0.95);

        auto make_graph = [&county_ids](auto&& params) {
            return make_graph_no_edges(params, county_ids);
        };
        BOOST_OUTCOME_TRY(
            epi::write_graph(make_graph(ensemble_params_p05), result_dir_p05.string(), epi::IOF_OmitDistributions));
        BOOST_OUTCOME_TRY(
            epi::write_graph(make_graph(ensemble_params_p25), result_dir_p25.string(), epi::IOF_OmitDistributions));
        BOOST_OUTCOME_TRY(
            epi::write_graph(make_graph(ensemble_params_p50), result_dir_p50.string(), epi::IOF_OmitDistributions));
        BOOST_OUTCOME_TRY(
            epi::write_graph(make_graph(ensemble_params_p75), result_dir_p75.string(), epi::IOF_OmitDistributions));
        BOOST_OUTCOME_TRY(
            epi::write_graph(make_graph(ensemble_params_p95), result_dir_p95.string(), epi::IOF_OmitDistributions));
    }
    return epi::success();
}

/**
 * Different modes for running the parameter study.
 */
enum class RunMode
{
    Load,
    Save,
};

epi::IOResult<void> generate_rki_series(RegionSpec region_spec, epi::Date start_date,
    epi::Date end_date, bool late, bool masks, const fs::path& data_dir, const fs::path& result_dir)
{
    namespace de = epi::regions::de;

    const auto start_day = epi::get_day_in_year(start_date);
    int start_summer;
    if (late) {
        start_summer = epi::get_day_in_year(epi::Date(2021, 8, 1));
    }
    else {
        start_summer = epi::get_day_in_year(epi::Date(2021, 7, 1));
    }

    //global parameters
    const int num_age_groups = 6;
    epi::SecirParams params(num_age_groups);
    params.get<epi::StartDay>()    = start_day;
    params.get<epi::StartSummer>() = start_summer;
    BOOST_OUTCOME_TRY(set_covid_parameters(params));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));
    BOOST_OUTCOME_TRY(set_npis(start_date, end_date, params, late, masks));

    BOOST_OUTCOME_TRY(county_ids, epi::get_county_ids((data_dir / "pydata" / "Germany").string()));

    auto county_ids_new_end = std::remove_if(county_ids.begin(), county_ids.end(), [region_spec](auto&& county_id) {
        return (region_spec.level != RegionLevel::County || county_id != region_spec.id) &&
               (region_spec.level != RegionLevel::State || de::get_state_id(de::CountyId(county_id)) != de::StateId(region_spec.id));
    });
    county_ids.erase(county_ids_new_end, county_ids.end());
    
    std::vector<epi::SecirModelV> counties(county_ids.size(), epi::SecirModelV(int(size_t(params.get_num_groups()))));
    for (auto& county : counties) {
        county.parameters = params;
    }

    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 1.0);
    auto scaling_factor_icu      = 1.0;

    auto r = epi::extrapolate_rki_results<epi::SecirModelV, epi::InfectionStateV>(
        counties, (data_dir / "pydata" / "Germany").string(), result_dir.string(), county_ids, start_date,
        scaling_factor_infected, scaling_factor_icu, epi::get_offset_in_days(end_date, start_date));
    BOOST_OUTCOME_TRY(r);

    return epi::success();
}

/**
 * Run the parameter study.
 * Load a previously stored graph or create a new one from data. 
 * The graph is the input for the parameter study.
 * A newly created graph is saved and can be reused.
 * @param mode Mode for running the parameter study.
 * @param data_dir data directory. Not used if mode is RunMode::Load.
 * @param save_dir directory where the graph is loaded from if mode is RunMOde::Load or save to if mode is RunMode::Save.
 * @param result_dir directory where all results of the parameter study will be stored.
 * @returns any io error that occurs during reading or writing of files.
 */
epi::IOResult<void> run(RunMode mode, const fs::path& data_dir, const fs::path& save_dir, const fs::path& result_dir,
                        bool late, bool masks)
{
    const auto start_date   = epi::Date(2021, 10, 15);
    const auto num_days_sim = 70;
    const auto end_date     = epi::offset_date_by_days(start_date, int(std::ceil(num_days_sim)));
    const auto num_runs     = 500;
    const auto region       = RegionSpec{RegionLevel::State, 14};

    //create or load graph
    epi::Graph<epi::SecirModelV, epi::MigrationParameters> params_graph;
    if (mode == RunMode::Save) {
        BOOST_OUTCOME_TRY(created, create_graph(region, start_date, end_date, data_dir, late, masks));
        BOOST_OUTCOME_TRY(save_graph(created, save_dir));
        params_graph = created;
    }
    else {
        BOOST_OUTCOME_TRY(loaded, load_graph(save_dir));
        params_graph = loaded;
    }

    std::vector<int> county_ids(params_graph.nodes().size());
    std::transform(params_graph.nodes().begin(), params_graph.nodes().end(), county_ids.begin(), [](auto& n) {
        return n.id;
    });

    //run parameter study
    auto parameter_study =
        epi::ParameterStudyV<epi::SecirSimulationV<>>{params_graph, 0.0, num_days_sim, 0.5, num_runs};
    auto ensemble_results = std::vector<std::vector<epi::TimeSeries<double>>>{};
    ensemble_results.reserve(size_t(num_runs));
    auto ensemble_params = std::vector<std::vector<epi::SecirModelV>>{};
    ensemble_params.reserve(size_t(num_runs));
    auto save_result_result = epi::IOResult<void>(epi::success());
    auto run_idx            = size_t(0);
    parameter_study.run([&](auto results_graph) {
        ensemble_results.push_back(epi::interpolate_simulation_result(results_graph));

        ensemble_params.emplace_back();
        ensemble_params.back().reserve(results_graph.nodes().size());
        std::transform(results_graph.nodes().begin(), results_graph.nodes().end(),
                       std::back_inserter(ensemble_params.back()), [](auto&& node) {
                           return node.property.get_simulation().get_model();
                       });

        if (save_result_result) {
            save_result_result =
                save_result(ensemble_results.back(), ensemble_params.back(), county_ids, result_dir, run_idx);
        }
        std::cout << "run " << run_idx << " complete.\n";
        ++run_idx;
    }, true);
    BOOST_OUTCOME_TRY(save_result_result);

    //post process to get detected cases per day
    auto ensemble_pp_result = std::vector<std::vector<epi::TimeSeries<double>>>{};
    for (auto run_idx = size_t(0); run_idx < ensemble_params.size(); ++run_idx)
    {
        ensemble_pp_result.emplace_back();
        for (auto node_idx = size_t(0); node_idx < ensemble_params[run_idx].size(); ++node_idx)
        {
            ensemble_pp_result.back().emplace_back(ensemble_results[run_idx][node_idx].get_num_elements());
            for (auto t_idx = Eigen::Index(0); t_idx < ensemble_results[run_idx][node_idx].get_num_time_points(); ++t_idx)
            {
                auto y = ensemble_results[run_idx][node_idx][t_idx];
                auto t = ensemble_results[run_idx][node_idx].get_time(t_idx);
                auto dydt = ensemble_pp_result.back().back().add_time_point(t);
                ensemble_params[run_idx][node_idx].get_derivatives_postprocess(y, t, dydt);
            }
        }
    }

    std::cout << "Saving Results\n";
    BOOST_OUTCOME_TRY(save_results(ensemble_pp_result, ensemble_params, county_ids, result_dir));

    // std::cout << "Generate rki series.\n";
    // BOOST_OUTCOME_TRY(generate_rki_series(region, start_date, end_date, late, masks, data_dir, result_dir));

    return epi::success();
}

int main(int argc, char** argv)
{
    //TODO: proper command line interface to set:
    //- number of runs
    //- start and end date (may be incompatible with runmode::load)
    //- seeds
    //- log level
    //- ...

    epi::set_log_level(epi::LogLevel::warn);

    bool late  = false;
    bool masks = false;

    RunMode mode;
    std::string save_dir;
    std::string data_dir;
    std::string result_dir;
    if (argc == 4) {
        mode       = RunMode::Save;
        data_dir   = argv[1];
        save_dir   = argv[2];
        result_dir = argv[3];
        printf("Reading data from \"%s\", saving graph to \"%s\".\n", data_dir.c_str(), save_dir.c_str());
    }
    else if (argc == 3) {
        mode       = RunMode::Load;
        save_dir   = argv[1];
        result_dir = argv[2];
        data_dir   = "";
        printf("Loading graph from \"%s\".\n", save_dir.c_str());
    }
    else {
        printf("Usage:\n");
        printf("paper1 <data_dir> <save_dir> <result_dir>\n");
        printf("\tMake graph with data from <data_dir> and save at <save_dir>, then run the simulation.\n");
        printf("\tStore the results in <result_dir>\n");
        printf("paper1 <load_dir> <result_dir>\n");
        printf("\tLoad graph from <load_dir>, then run the simulation.\n");
        return 0;
    }
    
    boost::filesystem::path dir(result_dir);
    bool created = boost::filesystem::create_directories(dir);

    if (created) {
        epi::log_info("Directory '{:s}' was created.", dir.string());
    }
    printf("Saving results to \"%s\".\n", result_dir.c_str());

    // epi::thread_local_rng().seed({}); //set seeds, e.g., for debugging
    printf("Seeds: ");
    for (auto s : epi::thread_local_rng().get_seeds()) {
        printf("%u, ", s);
    }
    printf("\n");

    auto result = run(mode, data_dir, save_dir, result_dir, late, masks);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    } else {
        printf("Simulation succesful.\n");
    }
    return 0;
}
