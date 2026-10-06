/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Sascha Korf
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
#include "halle_model.h"

#include "abm/infection.h"
#include "abm/location_type.h"
#include "abm/personal_rng.h"
#include "memilio/io/mobility_io.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"

#include <algorithm>
#include <fstream>
#include <map>
#include <numeric>
#include <sstream>

namespace mio
{
namespace halle
{

namespace
{

/// @brief Convert a {mean, standard deviation} pair into the {mu, sigma} of the matching log-normal distribution.
std::pair<double, double> mu_and_sigma(double mean, double stddev)
{
    const double mu    = std::log(mean * mean / std::sqrt(mean * mean + stddev * stddev));
    const double sigma = std::sqrt(std::log(1 + stddev * stddev / (mean * mean)));
    return {mu, sigma};
}

/// @brief Set a log-normally distributed transition time from its mean and standard deviation in days.
template <class Parameter>
void set_transition_time(abm::Parameters& params, double mean, double stddev)
{
    const auto p            = mu_and_sigma(mean, stddev);
    params.get<Parameter>() = ParameterDistributionLogNormal(p.first, p.second);
}

/**
 * @brief The age group index of the input data, mapped onto the six age groups used here.
 *
 * The population file stores an index into an eleven-group scheme (0: 0-4, 1: 5-6, 2: 7-15, 3: 16-18,
 * 4: 19-21, 5: 22-35, 6: 36-60, 7: 61-65, 8: 66-75, 9: 76-80, 10: 81+), which aggregates onto the six
 * RKI age groups used by all parameter values in this simulation.
 */
IOResult<AgeGroup> age_group_from_input(int input_age_index)
{
    switch (input_age_index) {
    case 0:
        return success(AgeGroup(0)); // 0-4
    case 1:
    case 2:
        return success(AgeGroup(1)); // 5-14
    case 3:
    case 4:
    case 5:
        return success(AgeGroup(2)); // 15-34
    case 6:
        return success(AgeGroup(3)); // 35-59
    case 7:
    case 8:
    case 9:
        return success(AgeGroup(4)); // 60-79
    case 10:
        return success(AgeGroup(5)); // 80+
    default:
        return failure(StatusCode::InvalidValue,
                       "Age group index " + std::to_string(input_age_index) + " is outside the expected range 0-10.");
    }
}

/// @brief Split a line of comma separated values.
std::vector<std::string> split_csv_line(const std::string& line)
{
    std::vector<std::string> values;
    std::stringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) {
        values.push_back(value);
    }
    return values;
}

/// @brief One row of the population file, after parsing.
struct PersonRow {
    AgeGroup age{0};
    int home_id     = 0;
    int school_id   = 0;
    int event_id    = 0;
    int shopping_id = 0;
    int work_id     = 0;
};

/**
 * @brief Read the Halle population file.
 *
 * Expects the header columns age, home_id, school_id, event_id, shopping_id and work_id; any further
 * columns are ignored. Reports an error rather than returning a partial population, so that a run
 * against unreadable data cannot look successful.
 */
IOResult<std::vector<PersonRow>> read_population_file(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, "Could not open population file " + filename);
    }

    std::string line;
    if (!std::getline(file, line)) {
        return failure(StatusCode::InvalidFileFormat, "Population file " + filename + " is empty.");
    }
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

    std::map<std::string, size_t> column;
    const auto header = split_csv_line(line);
    for (size_t i = 0; i < header.size(); ++i) {
        column[header[i]] = i;
    }
    for (const auto& required : {"age", "home_id", "school_id", "event_id", "shopping_id", "work_id"}) {
        if (column.count(required) == 0) {
            return failure(StatusCode::InvalidFileFormat,
                           "Population file " + filename + " has no column \"" + required + "\".");
        }
    }

    std::vector<PersonRow> rows;
    size_t line_number = 1;
    while (std::getline(file, line)) {
        ++line_number;
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty()) {
            continue;
        }
        const auto values = split_csv_line(line);
        if (values.size() < header.size()) {
            return failure(StatusCode::InvalidFileFormat, "Population file " + filename +
                                                              " has too few values in line " +
                                                              std::to_string(line_number) + ".");
        }
        PersonRow row;
        BOOST_OUTCOME_TRY(auto&& age, age_group_from_input(std::stoi(values[column["age"]])));
        row.age         = age;
        row.home_id     = std::stoi(values[column["home_id"]]);
        row.school_id   = std::stoi(values[column["school_id"]]);
        row.event_id    = std::stoi(values[column["event_id"]]);
        row.shopping_id = std::stoi(values[column["shopping_id"]]);
        row.work_id     = std::stoi(values[column["work_id"]]);
        rows.push_back(row);
    }
    if (rows.empty()) {
        return failure(StatusCode::InvalidFileFormat, "Population file " + filename + " contains no persons.");
    }
    return success(rows);
}

/**
 * @brief Add the locations of the population file to the model and assign every Person to them.
 *
 * Location ids of the input data are mapped onto model locations lazily, so only locations that are
 * actually used are created. The input topology is taken exactly as given: one model Location per id,
 * never split. The sizes it implies are plausible as institutions (schools median 104 and at most 700
 * Person%s, workplaces median 11 with a tail up to 4464, which is the scale of Halle's large employers),
 * so there is nothing to correct here.
 */
void add_persons_and_locations(abm::Model& model, const std::vector<PersonRow>& rows)
{
    // Location id of the input data -> model location, one per id.
    std::map<int, abm::LocationId> homes, events, shops, schools, works;

    auto get_or_add = [&model](std::map<int, abm::LocationId>& known, int input_id, abm::LocationType type) {
        auto it = known.find(input_id);
        if (it == known.end()) {
            it = known.emplace(input_id, model.add_location(type)).first;
        }
        return it->second;
    };

    // The population file has no sex column, but the PAIS parameters are stratified by sex, so it is drawn
    // per Person from the age dependent female share of the Halle population.
    const std::array<double, num_age_groups> female_share{0.481, 0.468, 0.490, 0.490, 0.568, 0.622};
    auto& uniform = UniformDistribution<double>::get_instance();

    for (const auto& row : rows) {
        const auto sex =
            uniform(model.get_rng(), 0.0, 1.0) < female_share[row.age.get()] ? abm::Sex::Female : abm::Sex::Male;
        const auto home = get_or_add(homes, row.home_id, abm::LocationType::Home);
        const auto pid  = model.add_person(home, row.age, sex);
        model.assign_location(pid, home);
        model.assign_location(pid, get_or_add(events, row.event_id, abm::LocationType::SocialEvent));
        model.assign_location(pid, get_or_add(shops, row.shopping_id, abm::LocationType::BasicsShop));
        if (model.parameters.get<abm::AgeGroupGotoSchool>()[row.age]) {
            model.assign_location(pid, get_or_add(schools, row.school_id, abm::LocationType::School));
        }
        if (model.parameters.get<abm::AgeGroupGotoWork>()[row.age]) {
            model.assign_location(pid, get_or_add(works, row.work_id, abm::LocationType::Work));
        }
    }
}

/**
 * @brief Set the age resolved contact rates of every Location from the German baseline contact matrices.
 *
 * ContactRates defaults to an all-zero baseline, so without this there is no contact transmission at all
 * and the epidemic is driven purely by the seeded history.
 *
 * The matrices are read from data/Germany/contacts rather than hardcoded, which is the same source the
 * setups on the abmXpanvadere and AIMS branches copied their tables from. Those copies both carry the
 * same transcription error in the home matrix (0.0504 where the file has 0.4504), which reading the file
 * avoids. The per-location factors are those of set_local_parameters_ger() on abmXpanvadere: they scale a
 * daily contact count to the time a Person spends at that location type and to the contact intensity
 * there.
 *
 * @param[in,out] model The model whose locations are given contact rates.
 * @param[in] contact_dir Directory holding baseline_home.txt and the other baseline matrices.
 */
IOResult<void> set_contact_rates(abm::Model& model, const std::string& contact_dir)
{
    const auto size = static_cast<Eigen::Index>(num_age_groups);
    BOOST_OUTCOME_TRY(auto&& home, read_mobility_plain(path_join(contact_dir, "baseline_home.txt")));
    BOOST_OUTCOME_TRY(auto&& school, read_mobility_plain(path_join(contact_dir, "baseline_school_pf_aver.txt")));
    BOOST_OUTCOME_TRY(auto&& work, read_mobility_plain(path_join(contact_dir, "baseline_work.txt")));
    BOOST_OUTCOME_TRY(auto&& other, read_mobility_plain(path_join(contact_dir, "baseline_other.txt")));

    for (const auto& [name, matrix] :
         {std::make_pair("baseline_home.txt", &home), std::make_pair("baseline_school_pf_aver.txt", &school),
          std::make_pair("baseline_work.txt", &work), std::make_pair("baseline_other.txt", &other)}) {
        if (matrix->rows() != size || matrix->cols() != size) {
            return failure(StatusCode::InvalidValue, std::string(name) + " is not a " + std::to_string(num_age_groups) +
                                                         "x" + std::to_string(num_age_groups) + " matrix.");
        }
    }

    auto as_contact_matrix = [](const Eigen::MatrixXd& baseline, double factor) {
        auto matrix           = ContactMatrix<ScalarType>(baseline.rows());
        matrix.get_baseline() = baseline * factor;
        return matrix;
    };

    // Factors from set_local_parameters_ger(): hours spent at the location type times contact intensity.
    const auto rates_home   = as_contact_matrix(home, 1.4 * 15.0);
    const auto rates_school = as_contact_matrix(school, 4.8);
    const auto rates_work   = as_contact_matrix(work, 3.0 * 0.5);
    const auto rates_event  = as_contact_matrix(other, 1.2 * 2.0 * 6.0);
    const auto rates_shop   = as_contact_matrix(other, 0.8 * 0.33 * 12.0);
    // Hospital, ICU and the transport types keep one contact per age group pair per day.
    auto rates_default           = ContactMatrix<ScalarType>(size);
    rates_default.get_baseline() = Eigen::MatrixXd::Constant(size, size, 1.0);

    for (auto& location : model.get_locations()) {
        switch (location.get_type()) {
        case abm::LocationType::Home:
            location.get_infection_parameters().get<abm::ContactRates>() = rates_home;
            break;
        case abm::LocationType::School:
            location.get_infection_parameters().get<abm::ContactRates>() = rates_school;
            break;
        case abm::LocationType::Work:
            location.get_infection_parameters().get<abm::ContactRates>() = rates_work;
            break;
        case abm::LocationType::SocialEvent:
            location.get_infection_parameters().get<abm::ContactRates>() = rates_event;
            break;
        case abm::LocationType::BasicsShop:
            location.get_infection_parameters().get<abm::ContactRates>() = rates_shop;
            break;
        default:
            location.get_infection_parameters().get<abm::ContactRates>() = rates_default;
            break;
        }
    }
    return success();
}

/**
 * @brief Add a Hospital and an ICU and assign every Person to both.
 *
 * Must be called after all Person%s are added. The mobility rules send severe and critical cases to the
 * Hospital and the ICU, but Model::perform_mobility only moves a Person to a location type they are
 * assigned to. Without the assignment they stay where they are, mostly at home, and keep infecting their
 * household. The infection course, and with it the deaths, is drawn in advance and does not depend on this.
 *
 * No MaximumContacts is set on any location: the contact matrices of set_contact_rates() are used as they
 * are, which is also what set_local_parameters_ger() on the abmXpanvadere branch does.
 */
void add_infrastructure(abm::Model& model)
{
    const auto hospital = model.add_location(abm::LocationType::Hospital);
    const auto icu      = model.add_location(abm::LocationType::ICU);
    for (auto& person : model.get_persons()) {
        model.assign_location(person.get_id(), hospital);
        model.assign_location(person.get_id(), icu);
    }
}

/// @brief A daily time series per age group, read from a CSV with columns date, age_group and a value column.
struct DailySeries {
    Date first_date{2020, 1, 1}; ///< Date of index 0 of every per-age series.
    std::vector<std::vector<double>> values{}; ///< Indexed [age group][day offset from first_date].
};

/**
 * @brief Read a CSV with the columns date, age_group and @p value_column into a daily series per age group.
 *
 * The age_group column holds the index of one of the six age groups used here. Missing days are zero.
 */
IOResult<DailySeries> read_daily_series(const std::string& filename, const std::string& value_column)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, "Could not open data file " + filename);
    }
    std::string line;
    if (!std::getline(file, line)) {
        return failure(StatusCode::InvalidFileFormat, "Data file " + filename + " is empty.");
    }
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

    std::map<std::string, size_t> column;
    const auto header = split_csv_line(line);
    for (size_t i = 0; i < header.size(); ++i) {
        column[header[i]] = i;
    }
    for (const auto& required : {std::string("date"), std::string("age_group"), value_column}) {
        if (column.count(required) == 0) {
            return failure(StatusCode::InvalidFileFormat,
                           "Data file " + filename + " has no column \"" + required + "\".");
        }
    }

    // Collect the rows first, since the file is not required to be sorted by date.
    std::vector<std::tuple<Date, size_t, double>> entries;
    Date min_date{9999, 12, 31}, max_date{1, 1, 1};
    while (std::getline(file, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty()) {
            continue;
        }
        const auto values = split_csv_line(line);
        BOOST_OUTCOME_TRY(auto&& date, parse_date(values[column["date"]]));
        const auto age = static_cast<size_t>(std::stoi(values[column["age_group"]]));
        if (age >= num_age_groups) {
            return failure(StatusCode::InvalidValue,
                           "Data file " + filename + " has age_group " + std::to_string(age) + " outside 0-5.");
        }
        entries.emplace_back(date, age, std::stod(values[column[value_column]]));
        if (get_offset_in_days(date, min_date) < 0) {
            min_date = date;
        }
        if (get_offset_in_days(date, max_date) > 0) {
            max_date = date;
        }
    }
    if (entries.empty()) {
        return failure(StatusCode::InvalidFileFormat, "Data file " + filename + " contains no rows.");
    }

    DailySeries series;
    series.first_date   = min_date;
    const auto num_days = static_cast<size_t>(get_offset_in_days(max_date, min_date)) + 1;
    series.values.assign(num_age_groups, std::vector<double>(num_days, 0.0));
    for (const auto& [date, age, value] : entries) {
        series.values[age][static_cast<size_t>(get_offset_in_days(date, min_date))] += value;
    }
    return success(series);
}

/**
 * @brief Draw @p count distinct indices from [0, @p size) without replacement.
 *
 * Uses a partial Fisher-Yates shuffle, so the cost is proportional to @p count rather than to @p size.
 */
std::vector<size_t> sample_without_replacement(RandomNumberGenerator& rng, size_t size, size_t count)
{
    count = std::min(count, size);
    std::map<size_t, size_t> swapped; // Sparse representation of the shuffled prefix.
    std::vector<size_t> result(count);
    auto& uniform_int = UniformIntDistribution<size_t>::get_instance();
    for (size_t i = 0; i < count; ++i) {
        const size_t j = uniform_int(rng, i, size - 1);
        auto value_at  = [&swapped](size_t k) {
            const auto it = swapped.find(k);
            return it == swapped.end() ? k : it->second;
        };
        result[i]  = value_at(j);
        swapped[j] = value_at(i);
    }
    return result;
}

/**
 * @brief Seed the infection and vaccination history of the population.
 *
 * For every age group, the number of Person%s that were infected in the lookback window before @p t0 is the
 * number of reported cases scaled by the fitted dark figure; the date of each infection is drawn from the
 * shape of the reported curve. Vaccinations are seeded the same way but are never scaled, since the fit must
 * not change the vaccination history.
 *
 * Seeding infections with their real dates is what puts the population into the correct immunity state at
 * @p t0, and it is also what initialises PAIS: Person::add_new_infection triggers PAIS::init_or_refresh, so
 * every seeded infection that recovered before @p t0 can carry a PAIS into the simulation.
 */
void seed_history(abm::Model& model, const DailySeries* cases, const DailySeries* vaccinations, abm::TimePoint t0,
                  Date start_date, double dark_figure, int lookback_days, RandomNumberGenerator& rng)
{
    // Group the persons by age group once, so both histories can sample from the same lists.
    std::vector<std::vector<abm::PersonId>> persons_by_age(num_age_groups);
    for (auto& person : model.get_persons()) {
        persons_by_age[person.get_age().get()].push_back(person.get_id());
    }

    // Draw `count` day offsets in [-lookback_days, 0), distributed like the given daily curve.
    auto sample_days = [&](const DailySeries& series, size_t age, size_t count) {
        std::vector<double> weights(static_cast<size_t>(lookback_days), 0.0);
        for (int d = 0; d < lookback_days; ++d) {
            const auto date   = offset_date_by_days(start_date, -lookback_days + d);
            const auto offset = get_offset_in_days(date, series.first_date);
            if (offset >= 0 && static_cast<size_t>(offset) < series.values[age].size()) {
                weights[static_cast<size_t>(d)] = std::max(0.0, series.values[age][static_cast<size_t>(offset)]);
            }
        }
        std::vector<int> days(count);
        if (std::accumulate(weights.begin(), weights.end(), 0.0) <= 0.0) {
            return days; // No data in the window, so nothing can be seeded for this age group.
        }
        auto& discrete = DiscreteDistribution<size_t>::get_instance();
        for (size_t i = 0; i < count; ++i) {
            days[i] = -lookback_days + static_cast<int>(discrete(rng, weights));
        }
        return days;
    };

    auto seed = [&](const DailySeries& series, double scale, bool is_infection) {
        for (size_t age = 0; age < num_age_groups; ++age) {
            double total = 0.0;
            for (int d = 0; d < lookback_days; ++d) {
                const auto date   = offset_date_by_days(start_date, -lookback_days + d);
                const auto offset = get_offset_in_days(date, series.first_date);
                if (offset >= 0 && static_cast<size_t>(offset) < series.values[age].size()) {
                    total += std::max(0.0, series.values[age][static_cast<size_t>(offset)]);
                }
            }
            const size_t requested = static_cast<size_t>(total * scale);
            const size_t count     = std::min(requested, persons_by_age[age].size());
            if (requested > persons_by_age[age].size()) {
                // Loudly, because a clamp here silently removes every susceptible of this age group and
                // the epidemic then cannot run at all. Shorten history_lookback_days if this appears.
                log_warning("Seeding {} in age group {} saturated: {} requested for {} persons. The whole age "
                            "group is seeded and no susceptible is left.",
                            is_infection ? "infections" : "vaccinations", age, requested, persons_by_age[age].size());
            }
            if (count == 0) {
                continue;
            }
            const auto chosen = sample_without_replacement(rng, persons_by_age[age].size(), count);
            const auto days   = sample_days(series, age, count);
            for (size_t i = 0; i < count; ++i) {
                const auto pid  = persons_by_age[age][chosen[i]];
                const auto date = t0 + abm::days(days[i]);
                auto& person    = model.get_person(pid);
                auto prng       = abm::PersonalRandomNumberGenerator(model.get_rng(), person);
                if (is_infection) {
                    person.add_new_infection(
                        abm::Infection(prng, abm::VirusVariant::Wildtype, person.get_age(), model.parameters, date),
                        prng, date, model.parameters);
                }
                else {
                    person.add_new_vaccination(abm::ProtectionType::GenericVaccine, date);
                }
            }
            log_info("Seeded {} {} in age group {}.", count, is_infection ? "infections" : "vaccinations", age);
        }
    };

    if (cases != nullptr) {
        seed(*cases, dark_figure, true);
    }
    if (vaccinations != nullptr) {
        seed(*vaccinations, 1.0, false);
    }
}

/// @brief The age bins of the cohort file, in the order used by Cohort::members.
const std::array<std::string, 7> cohort_age_bins{"18-20", "20-29", "30-39", "40-49", "50-59", "60-69", "70+"};

/**
 * @brief Which cohort age bins stand in for each model age group, weighted by the share of the group they cover.
 *
 * A bin that straddles two model age groups is used by both, with the number of years it contributes to each
 * as weight: 30-39 is half 30-34 and half 35-39, and 70+ covers 70-79 of the 60-79 group and all of 80+.
 * Weighting by years is a proxy for weighting by population; it keeps every row of a straddling bin available
 * to both groups, instead of splitting the rows of an already small bin at random.
 * Age groups 0-4 and 5-14 have no stand-in, since the cohort has no one under 18.
 */
const std::array<std::vector<std::pair<size_t, double>>, num_age_groups> cohort_donor_bins{{
    {}, // 0-4
    {}, // 5-14
    {{0, 2.0}, {1, 10.0}, {2, 5.0}}, // 15-34: 18-19, 20-29, 30-34. 15-17 is not covered.
    {{2, 5.0}, {3, 10.0}, {4, 10.0}}, // 35-59: 35-39, 40-49, 50-59
    {{5, 10.0}, {6, 10.0}}, // 60-79: 60-69, 70-79
    {{6, 1.0}}, // 80+
}};

/// @brief The individual history of one cohort member, restricted to before the start of the simulation.
struct CohortMember {
    std::vector<Date> infections{}; ///< Sorted, without duplicates.
    std::vector<Date> vaccinations{}; ///< Sorted, without duplicates.
};

/// @brief The cohort, grouped by the age bins of cohort_age_bins.
struct Cohort {
    std::array<std::vector<CohortMember>, cohort_age_bins.size()> members{};
};

/**
 * @brief Read the cohort history file.
 *
 * Only events strictly before @p start_date are kept, since everything from @p start_date on is simulated.
 * Dates that appear twice in the same row are kept once: the same infection or vaccination recorded twice
 * would otherwise be replayed twice. An unknown age bin or an unparseable date is an error rather than a
 * silent default, so that a changed file format cannot quietly misplace persons.
 */
IOResult<Cohort> read_cohort_file(const std::string& filename, Date start_date)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, "Could not open cohort history file " + filename);
    }
    std::string line;
    if (!std::getline(file, line)) {
        return failure(StatusCode::InvalidFileFormat, "Cohort history file " + filename + " is empty.");
    }
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

    std::map<std::string, size_t> column;
    const auto header = split_csv_line(line);
    for (size_t i = 0; i < header.size(); ++i) {
        column[header[i]] = i;
    }
    const std::vector<std::string> infection_columns{"ih_infection_1_date", "ih_infection_2_date",
                                                     "ih_infection_3_date", "ih_infection_4_date"};
    const std::vector<std::string> vaccination_columns{"ih_vaccine_1_date", "ih_vaccine_2_date",
                                                       "ih_vaccine_3_date"};
    std::vector<std::string> required{"age_group"};
    required.insert(required.end(), infection_columns.begin(), infection_columns.end());
    required.insert(required.end(), vaccination_columns.begin(), vaccination_columns.end());
    for (const auto& name : required) {
        if (column.count(name) == 0) {
            return failure(StatusCode::InvalidFileFormat,
                           "Cohort history file " + filename + " has no column \"" + name + "\".");
        }
    }

    // Parse the given columns of a row into the sorted, distinct dates before start_date.
    auto read_dates = [&](const std::vector<std::string>& values, const std::vector<std::string>& columns,
                          size_t line_number) -> IOResult<std::vector<Date>> {
        std::vector<Date> dates;
        for (const auto& name : columns) {
            const auto& value = values[column.at(name)];
            if (value == "NA" || value.empty()) {
                continue;
            }
            auto date = parse_date(value);
            if (!date) {
                return failure(StatusCode::InvalidValue, "Cohort history file " + filename + " has invalid date \"" +
                                                             value + "\" in line " + std::to_string(line_number) +
                                                             ".");
            }
            if (get_offset_in_days(date.value(), start_date) < 0) {
                dates.push_back(date.value());
            }
        }
        std::sort(dates.begin(), dates.end(), [](const Date& a, const Date& b) {
            return get_offset_in_days(a, b) < 0;
        });
        dates.erase(std::unique(dates.begin(), dates.end()), dates.end());
        return success(dates);
    };

    Cohort cohort;
    size_t line_number = 1;
    size_t num_members = 0;
    while (std::getline(file, line)) {
        ++line_number;
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty()) {
            continue;
        }
        const auto values = split_csv_line(line);
        if (values.size() < header.size()) {
            return failure(StatusCode::InvalidFileFormat, "Cohort history file " + filename +
                                                              " has too few values in line " +
                                                              std::to_string(line_number) + ".");
        }
        const auto& bin_name = values[column["age_group"]];
        const auto bin       = std::find(cohort_age_bins.begin(), cohort_age_bins.end(), bin_name);
        if (bin == cohort_age_bins.end()) {
            return failure(StatusCode::InvalidValue, "Cohort history file " + filename + " has unknown age group \"" +
                                                         bin_name + "\" in line " + std::to_string(line_number) +
                                                         ".");
        }
        CohortMember member;
        BOOST_OUTCOME_TRY(auto&& infections, read_dates(values, infection_columns, line_number));
        BOOST_OUTCOME_TRY(auto&& vaccinations, read_dates(values, vaccination_columns, line_number));
        member.infections   = infections;
        member.vaccinations = vaccinations;
        cohort.members[static_cast<size_t>(bin - cohort_age_bins.begin())].push_back(member);
        ++num_members;
    }
    if (num_members == 0) {
        return failure(StatusCode::InvalidFileFormat, "Cohort history file " + filename + " contains no rows.");
    }
    return success(cohort);
}

/**
 * @brief Seed the infection and vaccination history of the population from the cohort.
 *
 * Every Person is given the full history of one cohort member of their age, drawn with replacement, as in
 * abm_aims_halle.cpp. Differences to that setup, each of which changes the state at @p t0:
 * - Infections and vaccinations are replayed in the order they happened. PAIS::init_or_refresh draws the
 *   PAIS probability from the number of vaccinations the Person has at the time of the infection, so adding
 *   all infections before all vaccinations would treat every historical infection as unvaccinated.
 * - Each infection is drawn with the protection of the Person's latest infection or vaccination before it,
 *   as during the simulation, so prior immunity reduces its severity.
 * - Infection courses ending in death are redrawn. The cohort consists of persons who survived all of their
 *   recorded infections, and the population file is the population after them, so a death here would remove
 *   a Person twice.
 * - An infection that starts while the previous one has not ended is skipped, since the model keeps only one
 *   active infection per Person.
 * - PAIS::update_severity only runs during the simulation, so a PAIS acquired in the history would otherwise
 *   stay in its initial state until @p t0. It is advanced in daily steps from its onset up to every later
 *   infection, which may refresh it, and then up to @p t0.
 */
void seed_cohort_history(abm::Model& model, const Cohort& cohort, abm::TimePoint t0, Date start_date)
{
    constexpr int max_redraws = 100;
    const auto& params        = model.parameters;
    auto to_time_point        = [&](const Date& date) {
        return t0 + abm::days(get_offset_in_days(date, start_date));
    };

    // Advance the PAIS of a Person from `from` to `to` in daily steps. Persons without a PAIS by `to` are
    // skipped, which is most of them.
    auto advance_pais = [&](abm::Person& person, abm::PersonalRandomNumberGenerator& prng, abm::TimePoint from,
                            abm::TimePoint to) {
        if (person.get_pais_state(to) == abm::PAISState::Count) {
            return;
        }
        for (auto t = from + abm::days(1); t <= to; t += abm::days(1)) {
            person.update_pais(params, prng, t, abm::days(1));
        }
    };

    std::array<size_t, num_age_groups> num_unseeded{};
    size_t num_infections = 0, num_vaccinations = 0, num_overlapping = 0, num_redrawn = 0, num_pais = 0;
    auto& discrete    = DiscreteDistribution<size_t>::get_instance();
    auto& uniform_int = UniformIntDistribution<size_t>::get_instance();

    for (auto& person : model.get_persons()) {
        const auto age = person.get_age().get();
        std::vector<size_t> bins;
        std::vector<double> weights;
        for (const auto& [bin, weight] : cohort_donor_bins[age]) {
            if (!cohort.members[bin].empty()) {
                bins.push_back(bin);
                weights.push_back(weight);
            }
        }
        if (bins.empty()) {
            ++num_unseeded[age];
            continue;
        }
        // The Person's own generator, so the drawn history does not depend on the order Persons are visited.
        auto prng          = abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        const auto& donors = cohort.members[bins[discrete(prng, weights)]];
        const auto& member = donors[uniform_int(prng, size_t(0), donors.size() - 1)];

        // Merge both histories into one chronological list. On the same day, the infection comes first,
        // since a vaccination cannot have acted yet.
        std::vector<std::pair<abm::TimePoint, bool>> events; // (time, is_infection)
        for (const auto& date : member.infections) {
            events.emplace_back(to_time_point(date), true);
        }
        for (const auto& date : member.vaccinations) {
            events.emplace_back(to_time_point(date), false);
        }
        std::stable_sort(events.begin(), events.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });
        if (events.empty()) {
            continue;
        }

        auto pais_clock = events.front().first;
        for (const auto& [t, is_infection] : events) {
            advance_pais(person, prng, pais_clock, t);
            pais_clock = t;
            if (!is_infection) {
                person.add_new_vaccination(abm::ProtectionType::GenericVaccine, t);
                ++num_vaccinations;
                continue;
            }
            const auto state = person.get_infection_state(t);
            if (state != abm::InfectionState::Susceptible && state != abm::InfectionState::Recovered) {
                ++num_overlapping;
                continue;
            }
            auto draw = [&]() {
                return abm::Infection(prng, abm::VirusVariant::Wildtype, person.get_age(), params, t,
                                      abm::InfectionState::Exposed, person.get_latest_protection(t), false);
            };
            auto infection = draw();
            for (int i = 0; i < max_redraws && infection.get_highest_infection_state().second ==
                                                    abm::InfectionState::Dead;
                 ++i) {
                infection = draw();
                ++num_redrawn;
            }
            person.add_new_infection(std::move(infection), prng, t, params);
            ++num_infections;
        }
        advance_pais(person, prng, pais_clock, t0);
        if (person.has_active_pais(t0)) {
            ++num_pais;
        }
    }

    for (size_t age = 0; age < num_age_groups; ++age) {
        if (num_unseeded[age] > 0) {
            log_warning("The cohort has no member of age group {}, so its {} persons start without infection or "
                        "vaccination history, i.e. fully susceptible and without PAIS.",
                        age, num_unseeded[age]);
        }
    }
    log_info("Seeded {} infections and {} vaccinations from the cohort; skipped {} infections overlapping an "
             "earlier one and redrew {} fatal courses. {} persons have an active PAIS at t0.",
             num_infections, num_vaccinations, num_overlapping, num_redrawn, num_pais);
}

/// @brief Set the infection parameters that are not fitted. Values follow the Halle setup of abm_aims_halle.cpp.
void set_fixed_infection_parameters(abm::Parameters& params)
{
    set_transition_time<abm::TimeExposedToNoSymptoms>(params, 4.5, 1.5);
    set_transition_time<abm::TimeInfectedNoSymptomsToSymptoms>(params, 1.1, 0.9);
    set_transition_time<abm::TimeInfectedNoSymptomsToRecovered>(params, 8.0, 2.0);
    set_transition_time<abm::TimeInfectedSymptomsToSevere>(params, 6.6, 4.9);
    set_transition_time<abm::TimeInfectedSymptomsToRecovered>(params, 18.1, 6.3);
    set_transition_time<abm::TimeInfectedSevereToCritical>(params, 1.5, 2.0);
    set_transition_time<abm::TimeInfectedSevereToRecovered>(params, 18.1, 6.3);
    set_transition_time<abm::TimeInfectedSevereToDead>(params, 10.7, 4.8);
    set_transition_time<abm::TimeInfectedCriticalToDead>(params, 10.7, 4.8);
    set_transition_time<abm::TimeInfectedCriticalToRecovered>(params, 18.1, 6.3);

    // Protection conferred by a vaccination or an earlier infection. Without these, the seeded vaccination
    // history has no effect at all on the epidemic. Values follow the ABM paper, which cites
    // https://doi.org/10.1016/j.ebiom.2023.104734 for the severity protection and 10.1016/j.vaccine.2023.03.069
    // for the (negligible) protection against infection itself.
    params.get<abm::SeverityProtectionFactor>() =
        TimeSeriesFunctor<ScalarType>{TimeSeriesFunctorType::LinearInterpolation, {{0, 0.8}, {150, 0.8}}};
    params.get<abm::InfectionProtectionFactor>() =
        TimeSeriesFunctor<ScalarType>{TimeSeriesFunctorType::LinearInterpolation, {{0, 0.0}, {150, 0.0}}};
    params.get<abm::HighViralLoadProtectionFactor>() =
        TimeSeriesFunctor<ScalarType>{TimeSeriesFunctorType::LinearInterpolation, {{0, 0.0}, {150, 0.0}}};

    // Transmission is modelled through contacts only, as in the ABM paper.
    params.get<abm::AerosolTransmissionRates>() = 0.0;

    params.get<abm::TestData>()[abm::TestType::PCR]     = {0.9, 0.995, abm::hours(24), abm::TestType::PCR};
    params.get<abm::TestData>()[abm::TestType::Antigen] = {0.9, 0.9999, abm::minutes(15), abm::TestType::Antigen};
    params.get<abm::TestData>()[abm::TestType::Generic] = {0.7, 0.95, abm::minutes(30), abm::TestType::Generic};

    // Age dependent transition probabilities, indexed by the six age groups.
    const std::array<double, num_age_groups> symptoms_per_no_symptoms{0.50, 0.55, 0.60, 0.70, 0.83, 0.90};
    const std::array<double, num_age_groups> severe_per_symptoms{0.02, 0.03, 0.04, 0.07, 0.17, 0.24};
    const std::array<double, num_age_groups> critical_per_severe{0.10, 0.11, 0.12, 0.14, 0.33, 0.62};
    const std::array<double, num_age_groups> deaths_per_critical{0.12, 0.13, 0.15, 0.26, 0.40, 0.48};
    for (size_t i = 0; i < num_age_groups; ++i) {
        const auto age = AgeGroup(i);
        params.get<abm::SymptomsPerInfectedNoSymptoms>()[{abm::VirusVariant::Wildtype, age}] =
            symptoms_per_no_symptoms[i];
        params.get<abm::SeverePerInfectedSymptoms>()[{abm::VirusVariant::Wildtype, age}] = severe_per_symptoms[i];
        params.get<abm::CriticalPerInfectedSevere>()[{abm::VirusVariant::Wildtype, age}] = critical_per_severe[i];
        params.get<abm::DeathsPerInfectedCritical>()[{abm::VirusVariant::Wildtype, age}] = deaths_per_critical[i];
    }
}

/**
 * @brief Set the PAIS parameters, which are not fitted.
 *
 * Values are those of abm_aims_halle.cpp. They are held fixed because PAIS is a forward projection here:
 * with no PAIS data in the objective, the likelihood would be flat in every one of them.
 */
void set_pais_parameters(abm::Parameters& params)
{
    // Probability of developing a PAIS after an infection, by age group and sex, for the three vaccination
    // classes. Age groups 0-4 and 5-14 stay at the default of zero. Rows are the four oldest age groups,
    // columns are {female, male}.
    const std::array<std::array<double, 2>, 4> probability_unvaccinated{
        {{0.108, 0.051}, {0.136, 0.105}, {0.188, 0.108}, {0.192, 0.175}}};
    const std::array<std::array<double, 2>, 4> probability_one_or_two{
        {{0.167, 0.098}, {0.189, 0.133}, {0.166, 0.109}, {0.137, 0.078}}};
    const std::array<std::array<double, 2>, 4> probability_three_or_more{
        {{0.156, 0.107}, {0.193, 0.122}, {0.163, 0.106}, {0.154, 0.101}}};

    const std::array<std::pair<abm::VaccinationClass, const std::array<std::array<double, 2>, 4>*>, 3> by_class{
        {{abm::VaccinationClass::Zero, &probability_unvaccinated},
         {abm::VaccinationClass::OneOrTwo, &probability_one_or_two},
         {abm::VaccinationClass::ThreeOrMore, &probability_three_or_more}}};

    for (const auto& [vaccination_class, table] : by_class) {
        for (size_t row = 0; row < table->size(); ++row) {
            const auto age = AgeGroup(row + 2); // The table starts at age group 2 (15-34).
            params
                .get<abm::PAISProbability>()[{abm::VirusVariant::Wildtype, age, abm::Sex::Female, vaccination_class}] =
                (*table)[row][0];
            params.get<abm::PAISProbability>()[{abm::VirusVariant::Wildtype, age, abm::Sex::Male, vaccination_class}] =
                (*table)[row][1];
        }
    }

    const std::array<double, 3> severity_factor{0.267, 0.247, 0.217};
    const std::array<double, 3> reinfection_protection{0.038, 0.116, 0.118};
    for (size_t i = 0; i < 3; ++i) {
        const auto vaccination_class = static_cast<abm::VaccinationClass>(i);
        params.get<abm::PAISProbabilitySeverityFactor>()[{abm::VirusVariant::Wildtype, vaccination_class}] =
            severity_factor[i];
        params.get<abm::PAISProtectionAtSecondInfection>()[{abm::VirusVariant::Wildtype, vaccination_class}] =
            reinfection_protection[i];
    }

    // Transition rates between PAIS states, in 1/day, corresponding to a medium PAIS lasting about four
    // months and a severe one improving to medium over about six.
    //
    // NOTE: these are not the values of abm_aims_halle.cpp. That matrix is documented as per-day
    // probabilities while PAISTransitionMatrix is consumed as rates, and its Medium row (Medium->Severe
    // 0.99998 per day, Medium->Medium 1e-9) leaves no one in the medium state for more than a day under
    // either reading. The semantics of those numbers need to be settled before they are used here.
    auto rates                              = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(abm::PAISState::Count),
                                                                    static_cast<Eigen::Index>(abm::PAISState::Count))
                                                  .eval();
    const auto medium                       = static_cast<Eigen::Index>(abm::PAISState::Medium);
    const auto severe                       = static_cast<Eigen::Index>(abm::PAISState::Severe);
    const auto healthy                      = static_cast<Eigen::Index>(abm::PAISState::Healthy);
    rates(medium, healthy)                  = 1.0 / 120.0;
    rates(medium, severe)                   = 1.0 / 2000.0;
    rates(severe, medium)                   = 1.0 / 180.0;
    rates(severe, healthy)                  = 1.0 / 365.0;
    params.get<abm::PAISTransitionMatrix>() = rates;
}

} // namespace

const std::vector<FitParameter>& fit_parameters(HistorySource source)
{
    // Bounds taken from the grid search of the ABM paper. No contact reduction is fitted, since testing
    // is the only measure in this setup. The PAIS parameters are not fitted either: they act only on the
    // PAIS output channels, so without PAIS data the likelihood would be flat in all of them.
    // viral_shedding_rate must stay the first entry, since make_model reads it as theta[0].
    static const std::vector<FitParameter> reported_cases{
        {"viral_shedding_rate", 0.1, 20.0},
        {"dark_figure", 2.5, 10.0},
    };
    static const std::vector<FitParameter> cohort{
        {"viral_shedding_rate", 0.1, 20.0},
    };
    return source == HistorySource::Cohort ? cohort : reported_cases;
}

const std::vector<std::string>& observable_channels()
{
    // Only deaths are compared against real data. The PAIS channels are recorded as a forward projection,
    // to be able to inspect them for a parameter vector drawn from the posterior.
    static const std::vector<std::string> channels{"deaths", "pais_medium", "pais_severe"};
    return channels;
}

std::vector<double> sample_prior(RandomNumberGenerator& rng, HistorySource source)
{
    const auto& priors = fit_parameters(source);
    std::vector<double> theta(priors.size());
    auto& uniform = UniformDistribution<double>::get_instance();
    for (size_t i = 0; i < priors.size(); ++i) {
        theta[i] = uniform(rng, priors[i].lower, priors[i].upper);
    }
    return theta;
}

bool is_in_prior_support(const std::vector<double>& theta, HistorySource source)
{
    const auto& priors = fit_parameters(source);
    if (theta.size() != priors.size()) {
        return false;
    }
    for (size_t i = 0; i < priors.size(); ++i) {
        if (theta[i] < priors[i].lower || theta[i] > priors[i].upper) {
            return false;
        }
    }
    return true;
}

IOResult<abm::Model> make_model(const ModelSetup& setup, const std::vector<double>& theta, Date start_date,
                                abm::TimePoint t0, abm::TimePoint tmax, const RandomNumberGenerator& rng)
{
    const auto source = setup.history_source();
    if (!is_in_prior_support(theta, source)) {
        return failure(StatusCode::InvalidValue, "Parameter vector is outside the prior support.");
    }
    if (setup.person_file.empty()) {
        return failure(StatusCode::InvalidValue, "No population file given.");
    }
    if (source == HistorySource::Cohort && (!setup.cases_file.empty() || !setup.vaccinations_file.empty())) {
        return failure(StatusCode::InvalidValue, "A cohort history file cannot be combined with a cases or a "
                                                 "vaccinations file. Give either the cohort or the reported data.");
    }
    const bool has_history =
        source == HistorySource::Cohort || (!setup.cases_file.empty() && !setup.vaccinations_file.empty());
    if (!has_history && !setup.allow_missing_history) {
        return failure(StatusCode::InvalidValue,
                       "No infection and vaccination history given. Pass a cohort history file, or both a cases and "
                       "a vaccinations file, or set allow_missing_history for a smoke test.");
    }

    auto model      = abm::Model(num_age_groups, halle_county_id);
    model.get_rng() = rng;

    // Age group 1 (5-14) goes to school, age groups 2 and 3 (15-34 and 35-59) go to work.
    model.parameters.get<abm::AgeGroupGotoSchool>()              = false;
    model.parameters.get<abm::AgeGroupGotoSchool>()[AgeGroup(1)] = true;
    model.parameters.get<abm::AgeGroupGotoWork>()                = false;
    model.parameters.get<abm::AgeGroupGotoWork>().set_multiple({AgeGroup(2), AgeGroup(3)}, true);

    set_fixed_infection_parameters(model.parameters);
    set_pais_parameters(model.parameters);
    // viral_shedding_rate is the only fitted parameter acting on the model parameters directly;
    // dark_figure, if fitted, is applied when seeding the reported cases history.
    model.parameters.get<abm::InfectionRateFromViralShed>()[abm::VirusVariant::Wildtype] = theta[0];

    BOOST_OUTCOME_TRY(auto&& rows, read_population_file(setup.person_file));
    add_persons_and_locations(model, rows);
    add_infrastructure(model);
    BOOST_OUTCOME_TRY(set_contact_rates(model, setup.contact_dir));
    // No testing and no other intervention is applied for now. tmax is kept in the signature because it
    // is what a testing scheme would need as its validity period once testing is reinstated.
    unused(tmax);

    if (source == HistorySource::Cohort) {
        BOOST_OUTCOME_TRY(auto&& cohort, read_cohort_file(setup.history_file, start_date));
        seed_cohort_history(model, cohort, t0, start_date);
    }
    else if (has_history) {
        BOOST_OUTCOME_TRY(auto&& cases, read_daily_series(setup.cases_file, "new_cases"));
        BOOST_OUTCOME_TRY(auto&& vaccinations, read_daily_series(setup.vaccinations_file, "new_doses"));
        seed_history(model, &cases, &vaccinations, t0, start_date, theta[1], setup.history_lookback_days,
                     model.get_rng());
    }
    else {
        log_warning("Building the model without infection and vaccination history. The population has no immunity "
                    "and no pre-existing PAIS at the start of the simulation.");
    }

    log_info("Halle model built with {} persons and {} locations.", model.get_persons().size(),
             model.get_locations().size());
    return success(std::move(model));
}

} // namespace halle
} // namespace mio
