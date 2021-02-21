#include <epidemiology_io/secir_parameters_io.h>
#include <epidemiology_io/secir_result_io.h>
#include <epidemiology_io/io.h>
#include <epidemiology/utils/memory.h>
#include <epidemiology/utils/uncertain_value.h>
#include <epidemiology/utils/stl_util.h>
#include <epidemiology/utils/graph.h>
#include <epidemiology/migration/migration.h>
#include <epidemiology/secir/damping.h>
#include <epidemiology/model/populations.h>
#include <epidemiology/secir/uncertain_matrix.h>
#include <epidemiology/secir/secir.h>
#include <epidemiology/utils/compiler_diagnostics.h>
#include <epidemiology/utils/date.h>

#include <tixi.h>

#include <json/json.h>
#include <json/value.h>

#include <boost/filesystem.hpp>

#include <numeric>
#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <fstream>

namespace epi
{

namespace details
{
    Date get_max_date(Json::Value root)
    {
        Date max_date = Date(0, 1, 1);
        for (unsigned int i = 0; i < root.size(); i++) {
            auto date_temp = parse_date(root[i]["Date"].asString());
            if (date_temp > max_date) {
                max_date = date_temp;
            }
        }

        return max_date;
    }
    void interpolate_ages(const std::vector<double>& age_ranges, std::vector<std::vector<double>>& interpolation,
                          std::vector<bool>& carry_over)
    {
        std::vector<double> param_ranges = {5., 10., 20., 25., 20., 20.};

        //counter for parameter age groups
        size_t counter = 0;

        //residual of param age groups
        double res = 0.0;
        for (size_t i = 0; i < age_ranges.size(); i++) {

            // if current param age group didn't fit into previous rki age group, transfer residual to current age group
            if (res < 0) {
                interpolation[i].push_back(std::min(-res / age_ranges[i], 1.0));
            }

            if (counter < param_ranges.size() - 1) {
                res += age_ranges[i];
                if (std::abs(res) < age_ranges[i]) {
                    counter++;
                }
                // iterate over param age groups while there is still room in the current rki age group
                while (res > 0) {
                    res -= param_ranges[counter];
                    interpolation[i].push_back((param_ranges[counter] + std::min(res, 0.0)) / age_ranges[i]);
                    if (res >= 0) {
                        counter++;
                    }
                }
                if (res < 0) {
                    carry_over.push_back(true);
                }
                else if (res == 0) {
                    carry_over.push_back(false);
                }
            }
            // if last param age group is reached
            else {
                interpolation[i].push_back((age_ranges[i] + res) / age_ranges[i]);
                if (res < 0 || counter == 0) {
                    carry_over.push_back(true);
                }
                else if (res == 0) {
                    carry_over.push_back(false);
                }
                res = 0;
            }
        }
        // last entries for "unknown" age group
        interpolation.push_back({1.0});
        carry_over.push_back(true);
    }

    void read_rki_data(std::string const& path, const std::string& id_name, std::vector<int> const& vregion, Date date,
                       std::vector<std::vector<double>>& vnum_exp, std::vector<std::vector<double>>& vnum_car,
                       std::vector<std::vector<double>>& vnum_inf, std::vector<std::vector<double>>& vnum_hosp,
                       std::vector<std::vector<double>>& vnum_icu, std::vector<std::vector<double>>& vnum_death,
                       std::vector<std::vector<double>>& vnum_rec, const std::vector<std::vector<int>>& vt_car_to_rec,
                       const std::vector<std::vector<int>>& vt_car_to_inf,
                       const std::vector<std::vector<int>>& vt_exp_to_car,
                       const std::vector<std::vector<int>>& vt_inf_to_rec,
                       const std::vector<std::vector<int>>& vt_inf_to_hosp,
                       const std::vector<std::vector<int>>& vt_hosp_to_rec,
                       const std::vector<std::vector<int>>& vt_hosp_to_icu,
                       const std::vector<std::vector<int>>& vt_icu_to_dead,
                       const std::vector<std::vector<double>>& vmu_C_R, const std::vector<std::vector<double>>& vmu_I_H,
                       const std::vector<std::vector<double>>& vmu_H_U, const std::vector<double>& scaling_factor_inf)
    {
        if (!boost::filesystem::exists(path)) {
            log_error("RKI data file missing: {}.", path);
            return;
        }

        Json::Reader reader;
        Json::Value root;

        std::ifstream rki(path);
        if (!reader.parse(rki, root)) {
            log_error(reader.getFormattedErrorMessages());
            return;
        }

        Date max_date = get_max_date(root);
        assert(max_date != Date(0, 1, 1) && "File is epmty");
        if (max_date == Date(0, 1, 1)) {
            log_error("File is empty");
            return;
        }
        assert(max_date >= date && "Specified date does not exist in divi data");
        if (max_date < date) {
            log_error("Specified date does not exist in divi data");
            return;
        }
        auto days_surplus = get_offset_in_days(max_date, date) - 6;

        if (days_surplus > 0) {
            days_surplus = 0;
        }

        std::vector<std::string> age_names = {"A00-A04", "A05-A14", "A15-A34", "A35-A59", "A60-A79", "A80+", "unknown"};
        std::vector<double> age_ranges     = {5., 10., 20., 25., 20., 20.};

        for (unsigned int i = 0; i < root.size(); i++) {
            auto it = std::find_if(vregion.begin(), vregion.end(), [&root, i, &id_name](auto r) {
                return r == 0 || root[i][id_name].asInt() == r;
            });
            if (it != vregion.end()) {
                auto region_idx = size_t(it - vregion.begin());

                auto& t_exp_to_car  = vt_exp_to_car[region_idx];
                auto& t_car_to_rec  = vt_car_to_rec[region_idx];
                auto& t_car_to_inf  = vt_car_to_inf[region_idx];
                auto& t_inf_to_rec  = vt_inf_to_rec[region_idx];
                auto& t_inf_to_hosp = vt_inf_to_hosp[region_idx];
                auto& t_hosp_to_rec = vt_hosp_to_rec[region_idx];
                auto& t_hosp_to_icu = vt_hosp_to_icu[region_idx];
                auto& t_icu_to_dead = vt_icu_to_dead[region_idx];

                auto& num_car   = vnum_car[region_idx];
                auto& num_inf   = vnum_inf[region_idx];
                auto& num_rec   = vnum_rec[region_idx];
                auto& num_exp   = vnum_exp[region_idx];
                auto& num_hosp  = vnum_hosp[region_idx];
                auto& num_death = vnum_death[region_idx];
                auto& num_icu   = vnum_icu[region_idx];

                auto& mu_C_R = vmu_C_R[region_idx];
                auto& mu_I_H = vmu_I_H[region_idx];
                auto& mu_H_U = vmu_H_U[region_idx];

                auto date_df = parse_date(root[i]["Date"].asString());

                auto it_age = std::find(age_names.begin(), age_names.end() - 1, root[i]["Age_RKI"].asString());
                if (it_age != age_names.end() - 1) {
                    auto age = size_t(it_age - age_names.begin());

                    bool read_icu = false; //params.populations.get({age, SecirCompartments::U}) == 0;

                    if (date_df == offset_date_by_days(date, 0)) {
                        num_inf[age] += (1 - mu_C_R[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        num_rec[age] += root[i]["Confirmed"].asDouble();
                    }
                    if (date_df == offset_date_by_days(date, days_surplus)) {
                        num_car[age] +=
                            (2 * mu_C_R[age] - 1) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // -R9
                    if (date_df == offset_date_by_days(date, -t_car_to_rec[age] + days_surplus)) {
                        num_car[age] -= mu_C_R[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // +R2
                    if (date_df == offset_date_by_days(date, +t_exp_to_car[age] + days_surplus)) {
                        num_exp[age] += mu_C_R[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // +R3
                    if (date_df == offset_date_by_days(date, +t_car_to_inf[age] + days_surplus)) {
                        num_car[age] += (1 - mu_C_R[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        num_exp[age] -= (1 - mu_C_R[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // R2 - R9
                    if (date_df == offset_date_by_days(date, t_exp_to_car[age] - t_car_to_rec[age] + days_surplus)) {
                        num_exp[age] -= mu_C_R[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // R2 + R3
                    if (date_df == offset_date_by_days(date, t_exp_to_car[age] + t_car_to_inf[age] + days_surplus)) {
                        num_exp[age] += (1 - mu_C_R[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // -R4
                    if (date_df == offset_date_by_days(date, -t_inf_to_rec[age])) {
                        num_inf[age] -= (1 - mu_C_R[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // -R6
                    if (date_df == offset_date_by_days(date, -t_inf_to_hosp[age])) {
                        num_inf[age] -= mu_I_H[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        num_hosp[age] += mu_I_H[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // -R6 - R7
                    if (date_df == offset_date_by_days(date, -t_inf_to_hosp[age] - t_hosp_to_icu[age])) {
                        num_inf[age] +=
                            mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        num_hosp[age] -=
                            mu_I_H[age] * mu_H_U[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        if (read_icu) {
                            num_icu[age] +=
                                mu_H_U[age] * mu_I_H[age] * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        }
                    }
                    // -R6 - R5
                    if (date_df == offset_date_by_days(date, -t_inf_to_hosp[age] - t_hosp_to_rec[age])) {
                        num_inf[age] +=
                            mu_I_H[age] * (1 - mu_H_U[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                        num_hosp[age] -=
                            mu_I_H[age] * (1 - mu_H_U[age]) * scaling_factor_inf[age] * root[i]["Confirmed"].asDouble();
                    }
                    // -R10 - R6 - R7
                    if (date_df ==
                        offset_date_by_days(date, -t_icu_to_dead[age] - t_inf_to_hosp[age] - t_hosp_to_icu[age])) {
                        num_death[age] += root[i]["Deaths"].asDouble();
                    }
                    if (read_icu) {
                        // -R6 - R7 - R7
                        if (date_df == offset_date_by_days(date, -t_inf_to_hosp[age] - 2 * t_hosp_to_icu[age])) {
                            num_icu[age] -= mu_I_H[age] * mu_H_U[age] * mu_H_U[age] * scaling_factor_inf[age] *
                                            root[i]["Confirmed"].asDouble();
                        }
                        // -R6 - R5 - R7
                        if (date_df == offset_date_by_days(date, -t_inf_to_hosp[age] - t_hosp_to_icu[age])) {
                            num_icu[age] -= mu_I_H[age] * mu_H_U[age] * (1 - mu_H_U[age]) * scaling_factor_inf[age] *
                                            root[i]["Confirmed"].asDouble();
                        }
                    }
                }
            }
        }

        for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
            auto region = vregion[region_idx];

            auto& num_car   = vnum_car[region_idx];
            auto& num_inf   = vnum_inf[region_idx];
            auto& num_rec   = vnum_rec[region_idx];
            auto& num_exp   = vnum_exp[region_idx];
            auto& num_hosp  = vnum_hosp[region_idx];
            auto& num_death = vnum_death[region_idx];
            auto& num_icu   = vnum_icu[region_idx];

            size_t num_groups = age_ranges.size();
            for (size_t i = 0; i < num_groups; i++) {
                auto try_fix_constraints = [region, &age_names, i](double& value, double error, auto str) {
                    if (value < error) {
                        log_error("{:s} for age group {:s} is {:.4f} for region {:d}, exceeds expected negative value.",
                                  str, age_names[i], value, region);
                        value = 0.0;
                    }
                    else if (value < 0) {
                        log_info("{:s} for age group {:s} is {:.4f} for region {:d}, automatically corrected", str,
                                 age_names[i], value, region);
                        value = 0.0;
                    }
                };

                try_fix_constraints(num_inf[i], -5, "Infected");
                try_fix_constraints(num_car[i], -5, "Carrier");
                try_fix_constraints(num_exp[i], -5, "Exposed");
                try_fix_constraints(num_hosp[i], -5, "Hospitalized");
                try_fix_constraints(num_death[i], -5, "Dead");
                try_fix_constraints(num_icu[i], -5, "ICU");
                try_fix_constraints(num_rec[i], -20, "Recovered");
            }
        }
    }

    void read_divi_data(const std::string& path, const std::string& id_name, const std::vector<int>& vregion, Date date,
                        std::vector<double>& vnum_icu)
    {
        if (!boost::filesystem::exists(path)) {
            log_error("DIVI data file missing: {}.", path);
            return;
        }

        Json::Reader reader;
        Json::Value root;

        std::ifstream divi(path);
        reader.parse(divi, root);

        Date max_date = get_max_date(root);
        assert(max_date != Date(0, 1, 1) && "File is epmty");
        if (max_date == Date(0, 1, 1)) {
            log_error("File is empty");
            return;
        }
        assert(max_date >= date && "Specified date does not exist in divi data");
        if (max_date < date) {
            log_error("Specified date does not exist in divi data");
            return;
        }

        for (unsigned int i = 0; i < root.size(); i++) {
            auto it      = std::find_if(vregion.begin(), vregion.end(), [&root, i, &id_name](auto r) {
                return r == 0 || r == root[i][id_name].asInt();
            });
            auto date_df = parse_date(root[i]["Date"].asString());
            if (it != vregion.end() && date_df == date) {
                auto region_idx = size_t(it - vregion.begin());
                auto& num_icu   = vnum_icu[region_idx];
                num_icu         = root[i]["ICU"].asDouble();
            }
        }
    }

    std::vector<std::vector<double>> read_population_data(const std::string& path, const std::string& id_name,
                                                          const std::vector<int>& vregion)
    {
        if (!boost::filesystem::exists(path)) {
            log_error("Population data file missing: {}.", path);
            std::vector<std::vector<double>> num_population;
            return num_population;
        }

        Json::Reader reader;
        Json::Value root;

        std::ifstream census(path);
        reader.parse(census, root);

        std::vector<std::string> age_names = {"<3 years",    "3-5 years",   "6-14 years",  "15-17 years",
                                              "18-24 years", "25-29 years", "30-39 years", "40-49 years",
                                              "50-64 years", "65-74 years", ">74 years"};
        std::vector<double> age_ranges     = {3., 3., 9., 3., 7., 5., 10., 10., 15., 10., 25.};

        std::vector<std::vector<double>> interpolation(age_names.size());
        std::vector<bool> carry_over;

        interpolate_ages(age_ranges, interpolation, carry_over);

        std::vector<std::vector<double>> vnum_population(vregion.size(), std::vector<double>(age_names.size(), 0.0));

        for (unsigned int i = 0; i < root.size(); i++) {
            auto it = std::find_if(vregion.begin(), vregion.end(), [&root, i, &id_name](auto r) {
                return r == 0 || (int)root[i][id_name].asDouble() / 1000 == r || root[i][id_name] == r;
            });
            if (it != vregion.end()) {
                auto region_idx      = size_t(it - vregion.begin());
                auto& num_population = vnum_population[region_idx];
                for (size_t age = 0; age < age_names.size(); age++) {
                    num_population[age] += root[i][age_names[age]].asDouble();
                }
            }
        }

        std::vector<std::vector<double>> interpol_population(vregion.size(),
                                                             std::vector<double>(age_names.size(), 0.0));
        for (size_t region_idx = 0; region_idx < vregion.size(); ++region_idx) {
            auto& num_population = vnum_population[region_idx];

            int counter = 0;
            for (size_t i = 0; i < interpolation.size() - 1; i++) {
                for (size_t j = 0; j < interpolation[i].size(); j++) {
                    interpol_population[region_idx][counter] += interpolation[i][j] * num_population[i];
                    if (j < interpolation[i].size() - 1 || !carry_over[i]) {
                        counter++;
                    }
                }
            }
        }

        return interpol_population;
    }
} // namespace details

void write_element(const TixiDocumentHandle& handle, const std::string& path, const std::string& element_name,
                   const UncertainValue& element, int io_mode, int num_runs)
{
    auto element_path = path_join(path, element_name);

    if (io_mode == 0) {
        tixiAddDoubleElement(handle, path.c_str(), element_name.c_str(), (double)element, "%g");
    }
    else if (io_mode == 1 || io_mode == 2 || io_mode == 3) {
        assert(element.get_distribution().get() && ("No Distribution detected for " + element_name +
                                                    ". Either define a distribution or choose a different io_mode.")
                                                       .c_str());
        auto distribution = element.get_distribution().get();
        write_distribution(handle, path, element_name, *distribution);
        if (io_mode == 2) {
            tixiAddDoubleElement(handle, element_path.c_str(), "Value", (double)element, "%g");
        }
    }
    else {
        // TODO error handling
        epi::log_error("Wrong input for io_mode.");
    }

    if (io_mode == 3) {
        std::vector<double> predef_sample;
        for (int run = 0; run < num_runs; run++) {
            predef_sample.push_back((double)element);
        }
        write_predef_sample(handle, element_path, predef_sample);
    }
}

void write_distribution(const TixiDocumentHandle& handle, const std::string& path, const std::string& element_name,
                        const ParameterDistribution& distribution)
{

    struct WriteDistVisitor : public ConstParameterDistributionVisitor {
        WriteDistVisitor(const std::string& xml_path, TixiDocumentHandle tixi_handle)
            : handle(tixi_handle)
            , element_path(xml_path)
        {
        }

        void visit(const ParameterDistributionNormal& normal_distribution) override
        {
            tixiAddTextElement(handle, element_path.c_str(), "Distribution", "Normal");
            tixiAddDoubleElement(handle, element_path.c_str(), "Mean", normal_distribution.get_mean(), "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Deviation", normal_distribution.get_standard_dev(),
                                 "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Min", normal_distribution.get_lower_bound(), "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Max", normal_distribution.get_upper_bound(), "%g");
        }

        void visit(const ParameterDistributionUniform& uniform_distribution) override
        {
            tixiAddTextElement(handle, element_path.c_str(), "Distribution", "Uniform");
            tixiAddDoubleElement(handle, element_path.c_str(), "Min", uniform_distribution.get_lower_bound(), "%g");
            tixiAddDoubleElement(handle, element_path.c_str(), "Max", uniform_distribution.get_upper_bound(), "%g");
        }

        TixiDocumentHandle handle;
        std::string element_path;
    };

    tixiCreateElement(handle, path.c_str(), element_name.c_str());
    auto element_path = path_join(path, element_name);

    WriteDistVisitor visitor(element_path, handle);
    distribution.accept(visitor);

    tixiAddFloatVector(handle, element_path.c_str(), "PredefinedSamples", distribution.get_predefined_samples().data(),
                       static_cast<int>(distribution.get_predefined_samples().size()), "%g");
}

std::unique_ptr<UncertainValue> read_element(TixiDocumentHandle handle, const std::string& path, int io_mode)
{
    std::unique_ptr<UncertainValue> value;
    ReturnCode status;
    unused(status);

    if (io_mode == 0) {
        double read_buffer;
        status = tixiGetDoubleElement(handle, path.c_str(), &read_buffer);
        assert(status == SUCCESS && ("Failed to read value at " + path).c_str());
        value = std::make_unique<UncertainValue>(read_buffer);
    }
    else if (io_mode == 1 || io_mode == 2 || io_mode == 3) {
        std::unique_ptr<ParameterDistribution> distribution = read_distribution(handle, path);

        if (io_mode == 2) {
            double read_buffer;
            status = tixiGetDoubleElement(handle, path_join(path, "Value").c_str(), &read_buffer);
            assert(status == SUCCESS && ("Failed to read value at " + path).c_str());
            value = std::make_unique<UncertainValue>(read_buffer);
        }
        value->set_distribution(*distribution.get());
    }
    else {
        // TODO error handling
        epi::log_error("Wrong input for io_mode.");
    }
    return value;
}

std::unique_ptr<ParameterDistribution> read_distribution(TixiDocumentHandle handle, const std::string& path)
{
    ReturnCode status;
    unused(status);
    std::unique_ptr<ParameterDistribution> distribution;

    char* distri_str;
    status = tixiGetTextElement(handle, path_join(path, "Distribution").c_str(), &distri_str);
    assert(status == SUCCESS && ("Failed to read distribution type at " + path).c_str());
    if (strcmp("Normal", distri_str) == 0) {
        double mean;
        double dev;
        double min;
        double max;
        status = tixiGetDoubleElement(handle, path_join(path, "Mean").c_str(), &mean);
        assert(status == SUCCESS && ("Failed to read mean at " + path).c_str());

        status = tixiGetDoubleElement(handle, path_join(path, "Deviation").c_str(), &dev);
        assert(status == SUCCESS && ("Failed to read deviation at " + path).c_str());

        status = tixiGetDoubleElement(handle, path_join(path, "Min").c_str(), &min);
        assert(status == SUCCESS && ("Failed to read min value at " + path).c_str());

        status = tixiGetDoubleElement(handle, path_join(path, "Max").c_str(), &max);
        assert(status == SUCCESS && ("Failed to read max value at " + path).c_str());

        distribution = std::make_unique<ParameterDistributionNormal>(min, max, mean, dev);
    }
    else if (strcmp("Uniform", distri_str) == 0) {
        double min;
        double max;
        status = tixiGetDoubleElement(handle, path_join(path, "Min").c_str(), &min);
        assert(status == SUCCESS && ("Failed to read min value at " + path).c_str());

        status = tixiGetDoubleElement(handle, path_join(path, "Max").c_str(), &max);
        assert(status == SUCCESS && ("Failed to read max value at " + path).c_str());

        distribution = std::make_unique<ParameterDistributionUniform>(min, max);
    }
    else {
        // TODO: true error handling
        epi::log_error("Unknown distribution.");
        assert(false && "Unknown distribution.");
    }

    auto predef_path = path_join(path, "PredefinedSamples");
    int n_predef;
    tixiGetVectorSize(handle, predef_path.c_str(), &n_predef);

    double* predef = nullptr;
    tixiGetFloatVector(handle, predef_path.c_str(), &predef, n_predef);

    for (int i = 0; i < n_predef; i++) {
        distribution->add_predefined_sample(predef[i]);
    }

    return distribution;
}

void write_predef_sample(TixiDocumentHandle handle, const std::string& path, const std::vector<double>& samples)
{
    tixiRemoveElement(handle, path_join(path, "PredefinedSamples").c_str());
    tixiAddFloatVector(handle, path.c_str(), "PredefinedSamples", samples.data(), static_cast<int>(samples.size()),
                       "%g");
}

void write_damping_sampling(TixiDocumentHandle handle, const std::string& path, const std::string& name,
                            const DampingSampling& ds, int io_mode)
{
    auto ds_path = path_join(path, name);
    tixiCreateElement(handle, path.c_str(), name.c_str());
    write_element(handle, ds_path, "Value", ds.get_value(), io_mode, 0);
    tixiAddIntegerElement(handle, ds_path.c_str(), "Level", int(ds.get_level()), "%d");
    tixiAddIntegerElement(handle, ds_path.c_str(), "Type", int(ds.get_type()), "%d");
    tixiAddDoubleElement(handle, ds_path.c_str(), "Time", double(ds.get_time()), "%.18f");
    std::vector<double> matrices_double;
    std::transform(ds.get_matrix_indices().begin(), ds.get_matrix_indices().end(), std::back_inserter(matrices_double), [](auto& i) {
        return double(i);
    });
    tixiAddFloatVector(handle, ds_path.c_str(), "MatrixIndices", matrices_double.data(), int(matrices_double.size()), "%.0f");
    write_matrix(handle, ds_path, "GroupWeights", ds.get_group_weights());
}

DampingSampling read_damping_sampling(TixiDocumentHandle handle, const std::string& path, int io_mode)
{
    ReturnCode status; 
    unused(status);

    auto value = read_element(handle, path_join(path, "Value").c_str(), io_mode);
    int level, type;
    double time;
    status = tixiGetIntegerElement(handle, path_join(path, "Level").c_str(), &level);
    assert(status == SUCCESS && "Failed to read DampingSampling level.");
    status = tixiGetIntegerElement(handle, path_join(path, "Type").c_str(), &type);
    assert(status == SUCCESS && "Failed to read DampingSampling type.");
    tixiGetDoubleElement(handle, path_join(path, "Time").c_str(), &time);
    assert(status == SUCCESS && "Failed to read DampingSampling time.");
    int num_matrices;
    status = tixiGetVectorSize(handle, path_join(path, "MatrixIndices").c_str(), &num_matrices);
    assert(status == SUCCESS && "Failed to read DampingSampling matrix indices.");
    double* matrices_double;
    status = tixiGetFloatVector(handle, path_join(path, "MatrixIndices").c_str(), &matrices_double, num_matrices);
    assert(status == SUCCESS && "Failed to read DampingSampling matrix indices.");
    std::vector<size_t> matrices;
    std::transform(matrices_double, matrices_double + num_matrices, std::back_inserter(matrices), [](auto& d) {
        return int(d);
    });
    auto groups = read_matrix<Eigen::VectorXd>(handle, path_join(path, "GroupWeights"));

    return {*value, DampingLevel(level), DampingType(type), SimulationTime(time), matrices, groups};
}

void write_contact(TixiDocumentHandle handle, const std::string& path, const UncertainContactMatrix& contact_pattern,
                   int io_mode)
{
    tixiCreateElement(handle, path.c_str(), "ContactFreq");
    auto contact_path = path_join(path, "ContactFreq");

    write_damping_matrix_expression_collection(handle, contact_path, contact_pattern.get_cont_freq_mat());

    if (io_mode == 1 || io_mode == 2 || io_mode == 3) {
        tixiAddIntegerElement(handle, contact_path.c_str(), "NumDampings", int(contact_pattern.get_dampings().size()), "%d");
        for (size_t i = 0; i < contact_pattern.get_dampings().size(); ++i)
        {
            auto damping_name = "Damping" + std::to_string(i + 1);
            write_damping_sampling(handle, contact_path, damping_name, contact_pattern.get_dampings()[i], io_mode);
        }
    }
}

UncertainContactMatrix read_contact(TixiDocumentHandle handle, const std::string& path, int io_mode)
{
    ReturnCode status;
    unused(status);

    UncertainContactMatrix contact_patterns{read_damping_matrix_expression_collection(handle, path)};

    if (io_mode == 1 || io_mode == 2 || io_mode == 3) {
        int num_dampings;
        status = tixiGetIntegerElement(handle, path_join(path, "NumDampings").c_str(), &num_dampings);
        assert(status == SUCCESS && "Failed to read dampings.");
        for (size_t i = 0; i < size_t(num_dampings); ++i)
        {
            auto damping_name = "Damping" + std::to_string(i + 1);
            contact_patterns.get_dampings().push_back(read_damping_sampling(handle, path_join(path, damping_name).c_str(), io_mode));
        }
    }
    return contact_patterns;
}

std::vector<int> get_county_ids(const std::string& path)
{
    Json::Reader reader;
    Json::Value root;

    std::vector<int> id;

    std::ifstream census(path_join(path, "county_current_population.json"));
    if (!reader.parse(census, root)) {
        log_error(reader.getFormattedErrorMessages());
        return id;
    }

    for (unsigned int i = 0; i < root.size(); i++) {
        id.push_back(root[i]["ID_County"].asInt());
    }

    return id;
}

} // namespace epi
