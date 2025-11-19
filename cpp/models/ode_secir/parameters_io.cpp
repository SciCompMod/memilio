/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kuehn, Lena Ploetzke
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

#include "memilio/utils/compiler_diagnostics.h"

//see below for line that causes this warning
GCC_CLANG_DIAGNOSTIC(push)
GCC_CLANG_DIAGNOSTIC(ignored "-Wmaybe-uninitialized")

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "ode_secir/parameters_io.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/parameters_io.h"
#include "memilio/io/io.h"
#include "memilio/utils/stl_util.h"
#include "memilio/utils/date.h"

namespace mio
{

namespace osecir
{

namespace details
{

IOResult<void>
compute_confirmed_cases_data(std::vector<ConfirmedCasesDataEntry>& case_data, const int region, Date date,
                             std::vector<ScalarType>& num_Exposed, std::vector<ScalarType>& num_InfectedNoSymptoms,
                             std::vector<ScalarType>& num_InfectedSymptoms, std::vector<ScalarType>& num_InfectedSevere,
                             std::vector<ScalarType>& num_icu, std::vector<ScalarType>& num_death,
                             std::vector<ScalarType>& num_rec, const std::vector<int>& t_Exposed,
                             const std::vector<int>& t_InfectedNoSymptoms, const std::vector<int>& t_InfectedSymptoms,
                             const std::vector<int>& t_InfectedSevere, const std::vector<int>& t_InfectedCritical,
                             const std::vector<ScalarType>& mu_C_R, const std::vector<ScalarType>& mu_I_H,
                             const std::vector<ScalarType>& mu_H_U, const std::vector<ScalarType>& scaling_factor_inf)
{
    auto max_date_entry = std::max_element(case_data.begin(), case_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == case_data.end()) {
        log_error("RKI data file is empty.");
        return failure(StatusCode::InvalidFileFormat, "RKI file is empty.");
    }
    auto max_date = max_date_entry->date;
    if (max_date < date) {
        log_error("Specified date does not exist in RKI data");
        return failure(StatusCode::OutOfRange, "Specified date does not exist in RKI data.");
    }
    auto days_surplus = std::min(get_offset_in_days(max_date, date) - 6, 0);

    for (auto&& entry : case_data) {

        auto date_df = entry.date;
        auto age     = size_t(entry.age_group);

        if (date_df == offset_date_by_days(date, 0)) {
            num_InfectedSymptoms[age] += scaling_factor_inf[age] * entry.num_confirmed;
            // We intentionally do NOT multiply recovered with the scaling_factor_inf here.
            // If we apply the scaling factor to recovered as well, we would implicitly
            // assume that this factor holds for the entire historical period up to t0, which is not valid.
            num_rec[age] += entry.num_confirmed;
        }
        if (date_df == offset_date_by_days(date, days_surplus)) {
            num_InfectedNoSymptoms[age] -= 1 / (1 - mu_C_R[age]) * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (date_df == offset_date_by_days(date, t_InfectedNoSymptoms[age] + days_surplus)) {
            num_InfectedNoSymptoms[age] += 1 / (1 - mu_C_R[age]) * scaling_factor_inf[age] * entry.num_confirmed;
            num_Exposed[age] -= 1 / (1 - mu_C_R[age]) * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (date_df == offset_date_by_days(date, t_Exposed[age] + t_InfectedNoSymptoms[age] + days_surplus)) {
            num_Exposed[age] += 1 / (1 - mu_C_R[age]) * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (date_df == offset_date_by_days(date, -t_InfectedSymptoms[age])) {
            num_InfectedSymptoms[age] -= scaling_factor_inf[age] * entry.num_confirmed;
            num_InfectedSevere[age] += mu_I_H[age] * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (date_df == offset_date_by_days(date, -t_InfectedSymptoms[age] - t_InfectedSevere[age])) {
            num_InfectedSevere[age] -= mu_I_H[age] * scaling_factor_inf[age] * entry.num_confirmed;
            num_icu[age] += mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * entry.num_confirmed;
        }
        if (date_df ==
            offset_date_by_days(date, -t_InfectedSymptoms[age] - t_InfectedSevere[age] - t_InfectedCritical[age])) {
            num_death[age] += entry.num_deaths;
            num_icu[age] -= mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * entry.num_confirmed;
        }
    }

    for (size_t i = 0; i < ConfirmedCasesDataEntry::age_group_names.size(); i++) {
        // R(t0) = ΣC(t0) − I(t0) − H(t0) − U(t0) − D(t0)
        // subtract currently infectious/hospitalized/ICU/dead
        // Note: We divide I/H/U/D by scaling_factor_inf to "unscale" these contributions back to the
        // reported level before subtracting from recovered. If we also applied the scaling factor to
        // recovered, we would implicitly assume that the same underreporting applies to the entire
        // history up to t0, which would be wrong. The scaling factor should reflect underreporting
        // around t0 only.
        num_rec[i] -= (num_InfectedSymptoms[i] / scaling_factor_inf[i] + num_InfectedSevere[i] / scaling_factor_inf[i] +
                       num_icu[i] / scaling_factor_inf[i] + num_death[i] / scaling_factor_inf[i]);

        auto try_fix_constraints = [region, i](ScalarType& value, ScalarType error, auto str) {
            if (value < error) {
                //this should probably return a failure
                //but the algorithm is not robust enough to avoid large negative values and there are tests that rely on it
                log_error("{:s} for age group {:s} is {:.4f} for region {:d}, exceeds expected negative value.", str,
                          ConfirmedCasesDataEntry::age_group_names[i], value, region);
                value = 0.0;
            }
            else if (value < 0) {
                log_info("{:s} for age group {:s} is {:.4f} for region {:d}, automatically corrected", str,
                         ConfirmedCasesDataEntry::age_group_names[i], value, region);
                value = 0.0;
            }
        };

        try_fix_constraints(num_InfectedSymptoms[i], -5, "InfectedSymptoms");
        try_fix_constraints(num_InfectedNoSymptoms[i], -5, "InfectedNoSymptoms");
        try_fix_constraints(num_Exposed[i], -5, "Exposed");
        try_fix_constraints(num_InfectedSevere[i], -5, "InfectedSevere");
        try_fix_constraints(num_death[i], -5, "Dead");
        try_fix_constraints(num_icu[i], -5, "InfectedCritical");
        try_fix_constraints(num_rec[i], -20, "Recovered");
    }

    return success();
}

IOResult<void> set_confirmed_cases_data(Model<ScalarType>& model, std::vector<ConfirmedCasesDataEntry>& case_data,
                                        const int region, Date date, const std::vector<ScalarType>& scaling_factor_inf)
{
    const size_t num_age_groups = ConfirmedCasesDataEntry::age_group_names.size();
    // allow single scalar scaling that is broadcast to all age groups
    assert(scaling_factor_inf.size() == 1 || scaling_factor_inf.size() == num_age_groups);

    // Set scaling factors to match num age groups
    std::vector<ScalarType> scaling_factor_inf_full;
    if (scaling_factor_inf.size() == 1) {
        scaling_factor_inf_full.assign(num_age_groups, scaling_factor_inf[0]);
    }
    else {
        scaling_factor_inf_full = scaling_factor_inf;
    }

    std::vector<int> t_InfectedNoSymptoms;
    std::vector<int> t_Exposed;
    std::vector<int> t_InfectedSymptoms;
    std::vector<int> t_InfectedSevere;
    std::vector<int> t_InfectedCritical;

    std::vector<ScalarType> mu_C_R;
    std::vector<ScalarType> mu_I_H;
    std::vector<ScalarType> mu_H_U;
    std::vector<ScalarType> mu_U_D;

    std::vector<ScalarType> num_InfectedSymptoms(num_age_groups, 0.0);
    std::vector<ScalarType> num_death(num_age_groups, 0.0);
    std::vector<ScalarType> num_rec(num_age_groups, 0.0);
    std::vector<ScalarType> num_Exposed(num_age_groups, 0.0);
    std::vector<ScalarType> num_InfectedNoSymptoms(num_age_groups, 0.0);
    std::vector<ScalarType> num_InfectedSevere(num_age_groups, 0.0);
    std::vector<ScalarType> num_icu(num_age_groups, 0.0);

    const size_t model_groups = (size_t)model.parameters.get_num_groups();
    assert(model_groups == 1 || model_groups == num_age_groups);
    for (size_t ag = 0; ag < num_age_groups; ag++) {
        // If the model has fewer groups than casedata entries available,
        // reuse group 0 parameters for all RKI age groups
        const size_t group = (model_groups == num_age_groups) ? ag : 0;

        t_Exposed.push_back(
            static_cast<int>(std::round(model.parameters.template get<TimeExposed<ScalarType>>()[(AgeGroup)group])));
        t_InfectedNoSymptoms.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedNoSymptoms<ScalarType>>()[(AgeGroup)group])));
        t_InfectedSymptoms.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSymptoms<ScalarType>>()[(AgeGroup)group])));
        t_InfectedSevere.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedSevere<ScalarType>>()[(AgeGroup)group])));
        t_InfectedCritical.push_back(static_cast<int>(
            std::round(model.parameters.template get<TimeInfectedCritical<ScalarType>>()[(AgeGroup)group])));

        mu_C_R.push_back(model.parameters.template get<RecoveredPerInfectedNoSymptoms<ScalarType>>()[(AgeGroup)group]);
        mu_I_H.push_back(model.parameters.template get<SeverePerInfectedSymptoms<ScalarType>>()[(AgeGroup)group]);
        mu_H_U.push_back(model.parameters.template get<CriticalPerSevere<ScalarType>>()[(AgeGroup)group]);
    }

    BOOST_OUTCOME_TRY(compute_confirmed_cases_data(
        case_data, region, date, num_Exposed, num_InfectedNoSymptoms, num_InfectedSymptoms, num_InfectedSevere, num_icu,
        num_death, num_rec, t_Exposed, t_InfectedNoSymptoms, t_InfectedSymptoms, t_InfectedSevere, t_InfectedCritical,
        mu_C_R, mu_I_H, mu_H_U, scaling_factor_inf_full));

    if (model_groups == num_age_groups) {
        for (size_t i = 0; i < model_groups; i++) {
            model.populations[{AgeGroup(i), InfectionState::Exposed}]                     = num_Exposed[i];
            model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptoms}]          = num_InfectedNoSymptoms[i];
            model.populations[{AgeGroup(i), InfectionState::InfectedNoSymptomsConfirmed}] = 0;
            model.populations[{AgeGroup(i), InfectionState::InfectedSymptoms}]            = num_InfectedSymptoms[i];
            model.populations[{AgeGroup(i), InfectionState::InfectedSymptomsConfirmed}]   = 0;
            model.populations[{AgeGroup(i), InfectionState::InfectedSevere}]              = num_InfectedSevere[i];
            // Only set the number of ICU patients here, if the date is not available in the data.
            if (!is_divi_data_available(date)) {
                model.populations[{AgeGroup(i), InfectionState::InfectedCritical}] = num_icu[i];
            }
            model.populations[{AgeGroup(i), InfectionState::Recovered}] = num_rec[i];
            model.populations[{AgeGroup(i), InfectionState::Dead}]      = num_death[i];
        }
    }
    else {
        const auto sum_vec = [](const std::vector<ScalarType>& v) {
            return std::accumulate(v.begin(), v.end(), ScalarType(0.0), [](const ScalarType& a, const ScalarType& b) {
                return evaluate_intermediate<ScalarType>(a + b);
            });
        };
        model.populations[{AgeGroup(0), InfectionState::Exposed}]                     = sum_vec(num_Exposed);
        model.populations[{AgeGroup(0), InfectionState::InfectedNoSymptoms}]          = sum_vec(num_InfectedNoSymptoms);
        model.populations[{AgeGroup(0), InfectionState::InfectedNoSymptomsConfirmed}] = 0;
        model.populations[{AgeGroup(0), InfectionState::InfectedSymptoms}]            = sum_vec(num_InfectedSymptoms);
        model.populations[{AgeGroup(0), InfectionState::InfectedSymptomsConfirmed}]   = 0;
        model.populations[{AgeGroup(0), InfectionState::InfectedSevere}]              = sum_vec(num_InfectedSevere);
        if (!is_divi_data_available(date)) {
            model.populations[{AgeGroup(0), InfectionState::InfectedCritical}] = sum_vec(num_icu);
        }
        model.populations[{AgeGroup(0), InfectionState::Dead}]      = sum_vec(num_death);
        model.populations[{AgeGroup(0), InfectionState::Recovered}] = sum_vec(num_rec);
    }
    if (std::accumulate(num_InfectedSymptoms.begin(), num_InfectedSymptoms.end(), ScalarType(0.0),
                        [](const ScalarType& a, const ScalarType& b) {
                            return evaluate_intermediate<ScalarType>(a + b);
                        }) == 0.0) {
        log_warning(
            "No infections for unvaccinated reported on date {} for region {}. Population data has not been set.", date,
            region);
    }
    return success();
}

IOResult<void> set_confirmed_cases_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path,
                                        Date date, const std::vector<ScalarType>& scaling_factor_inf)
{
    BOOST_OUTCOME_TRY(auto&& case_data, mio::read_confirmed_cases_data(path));

    // sort case_data into regions and ignore once with no region associated
    std::vector<std::vector<ConfirmedCasesDataEntry>> vcase_data{model.size()};
    for (auto&& entry : case_data) {
        auto it = std::find_if(model.begin(), model.end(), [&entry](auto m) {
            return m.id == 0 || get_region_id(entry) == m.id;
        });
        if (it != model.end()) {
            auto region_idx = size_t(it - model.begin());
            vcase_data[region_idx].push_back(entry);
        }
    }

    for (size_t region_idx = 0; region_idx < model.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(set_confirmed_cases_data(model[region_idx].property, vcase_data[region_idx],
                                                   model[region_idx].id, date, scaling_factor_inf));
    }
    return success();
}

IOResult<void> set_population_data(Model<ScalarType>& model, const std::vector<ScalarType>& num_population,
                                   const int region)
{
    if (std::accumulate(num_population.begin(), num_population.end(), ScalarType(0.0),
                        [](const ScalarType& a, const ScalarType& b) {
                            return evaluate_intermediate<ScalarType>(a + b);
                        }) <= 0) {
        log_warning("No population data available for region " + std::to_string(region) +
                    ". Population data has not been set.");
        return success();
    }

    auto num_groups  = model.parameters.get_num_groups();
    auto data_groups = num_population.size();

    if (data_groups == (size_t)num_groups) {
        for (auto i = AgeGroup(0); i < num_groups; i++) {
            model.populations.template set_difference_from_group_total<AgeGroup>({i, InfectionState::Susceptible},
                                                                                 num_population[(size_t)i]);
        }
    }
    else if ((size_t)num_groups == 1 && data_groups >= 1) {
        const ScalarType total = std::accumulate(num_population.begin(), num_population.end(), ScalarType(0.0),
                                                 [](const ScalarType& a, const ScalarType& b) {
                                                     return evaluate_intermediate<ScalarType>(a + b);
                                                 });
        model.populations.template set_difference_from_group_total<AgeGroup>({AgeGroup(0), InfectionState::Susceptible},
                                                                             total);
    }

    for (auto i = AgeGroup(0); i < AgeGroup(6); i++) {
        for (auto j = Index<InfectionState>(0); j < InfectionState::Count; ++j) {
            if (model.populations[{i, j}] < 0) {
                log_warning("Compartment at age group {}, infection state {}, is negative: {}", size_t(i), size_t(j),
                            model.populations[{i, j}]);
            }
        }
    }
    return success();
}

IOResult<void> set_population_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path)
{
    std::vector<int> vregion;
    std::transform(model.begin(), model.end(), std::back_inserter(vregion), [](const auto& m) {
        return m.id;
    });
    BOOST_OUTCOME_TRY(const auto&& num_population, read_population_data(path, vregion));

    for (size_t region_idx = 0; region_idx < model.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(
            set_population_data(model[region_idx].property, num_population[region_idx], model[region_idx].id));
    }
    return success();
}

IOResult<void> set_divi_data(Model<ScalarType>& model, const ScalarType num_icu, ScalarType scaling_factor_icu)
{
    ScalarType sum_mu_I_U = 0;
    std::vector<ScalarType> mu_I_U;
    auto num_groups = model.parameters.get_num_groups();
    for (auto i = AgeGroup(0); i < num_groups; i++) {
        sum_mu_I_U += model.parameters.template get<CriticalPerSevere<ScalarType>>()[i] *
                      model.parameters.template get<SeverePerInfectedSymptoms<ScalarType>>()[i];
        mu_I_U.push_back(model.parameters.template get<CriticalPerSevere<ScalarType>>()[i] *
                         model.parameters.template get<SeverePerInfectedSymptoms<ScalarType>>()[i]);
    }

    for (auto i = AgeGroup(0); i < num_groups; i++) {
        model.populations[{i, InfectionState::InfectedCritical}] =
            scaling_factor_icu * num_icu * mu_I_U[(size_t)i] / sum_mu_I_U;
    }

    return success();
}

IOResult<void> set_divi_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, const std::string& path, Date date,
                             ScalarType scaling_factor_icu)
{
    // DIVI dataset will no longer be updated from CW29 2024 on.
    if (!is_divi_data_available(date)) {
        log_warning("No DIVI data available for date: {}. "
                    "ICU compartment will be set based on Case data.",
                    date);
        return success();
    }

    std::vector<int> vregion;
    std::transform(model.begin(), model.end(), std::back_inserter(vregion), [](const auto& m) {
        return m.id;
    });
    BOOST_OUTCOME_TRY(auto&& num_icu, read_divi_data(path, vregion, date));

    for (size_t region_idx = 0; region_idx < model.size(); ++region_idx) {
        BOOST_OUTCOME_TRY(set_divi_data(model[region_idx].property, num_icu[region_idx], scaling_factor_icu));
    }

    return success();
}

} //namespace details

IOResult<void> read_input_data(const mio::VectorRange<Node<Model<ScalarType>>>& model, Date date,
                               const std::vector<ScalarType>& scaling_factor_inf, ScalarType scaling_factor_icu,
                               const mio::regions::de::EpidataFilenames& epidata_filenames)
{
    BOOST_OUTCOME_TRY(details::set_divi_data(model, epidata_filenames.divi_data_path, date, scaling_factor_icu));

    BOOST_OUTCOME_TRY(
        details::set_confirmed_cases_data(model, epidata_filenames.case_data_path, date, scaling_factor_inf));
    BOOST_OUTCOME_TRY(details::set_population_data(model, epidata_filenames.population_data_path));
    return success();
}

#ifdef MEMILIO_HAS_HDF5

IOResult<void> export_input_data_timeseries(const mio::VectorRange<Node<Model<ScalarType>>> model,
                                            const std::string& results_dir, Date date,
                                            const std::vector<ScalarType>& scaling_factor_inf,
                                            ScalarType scaling_factor_icu, int num_days,
                                            const mio::regions::de::EpidataFilenames& epidata_filenames)
{
    const auto num_age_groups = (size_t)model[0].property.parameters.get_num_groups();
    // allow scalar scaling factor as convenience for 1-group models
    assert(scaling_factor_inf.size() == 1 || scaling_factor_inf.size() == num_age_groups);

    std::vector<TimeSeries<ScalarType>> extrapolated_data(
        model.size(), TimeSeries<ScalarType>::zero(num_days + 1, (size_t)InfectionState::Count * num_age_groups));

    for (int t = 0; t <= num_days; ++t) {
        auto offset_day = offset_date_by_days(date, t);

        BOOST_OUTCOME_TRY(
            read_input_data(model, offset_day, scaling_factor_inf, scaling_factor_icu, epidata_filenames));

        for (size_t r = 0; r < model.size(); r++) {
            extrapolated_data[r][t] = model[r].property.get_initial_values();
        }
    }

    std::vector<int> vregion;
    std::transform(model.begin(), model.end(), std::back_inserter(vregion), [](const auto& m) {
        return m.id;
    });
    BOOST_OUTCOME_TRY(save_result(extrapolated_data, vregion, static_cast<int>(num_age_groups),
                                  path_join(results_dir, "Results_rki.h5")));

    auto extrapolated_rki_data_sum = sum_nodes(std::vector<std::vector<TimeSeries<ScalarType>>>{extrapolated_data});
    BOOST_OUTCOME_TRY(save_result({extrapolated_rki_data_sum[0][0]}, {0}, static_cast<int>(num_age_groups),
                                  path_join(results_dir, "Results_rki_sum.h5")));

    return success();
}
#else
IOResult<void> export_input_data_county_timeseries(const mio::VectorRange<Node<Model<ScalarType>>>, const std::string&,
                                                   Date, const std::vector<int>&, const std::vector<ScalarType>&,
                                                   const ScalarType, const int,
                                                   const mio::regions::de::EpidataFilenames&)
{
    mio::log_warning("HDF5 not available. Cannot export time series of extrapolated real data.");
    return success();
}
#endif // MEMILIO_HAS_HDF5

} // namespace osecir

} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

GCC_CLANG_DIAGNOSTIC(pop)
