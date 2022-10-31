/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
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

#include "abm/output_results.h"

namespace mio
{

namespace abm
{

void OutputResults::store_result_at(const TimePoint t, const World& world)
{
    if (m_print_results) {
        m_result.add_time_point(t.days());
        m_result.get_last_value().setZero();
        for (auto&& locations : world.get_locations()) {
            for (auto& location : locations) {
                m_result.get_last_value() += location.get_subpopulations().cast<double>();
            }
        }
    }

    if (m_print_location_results) {
        for (auto it = m_location_result.begin(); it != m_location_result.end(); ++it) {
            it->second.add_time_point(t.days());
            it->second.get_last_value().setZero();
        }
        for (auto&& locations : world.get_locations()) {
            for (auto& location : locations) {
                m_location_result.at((std::uint32_t)location.get_type()).get_last_value() +=
                    location.get_subpopulations().cast<double>();
            }
        }
    }
}

void OutputResults::set_print_results(const bool print_results, const bool print_location_results)
{
    m_print_results          = print_results;
    m_print_location_results = print_location_results;
}

void OutputResults::print_result_to_file(const std::string& name_of_file)
{
    if (m_print_results) {
        auto f_abm = fopen(name_of_file.c_str(), "w");
        fprintf(f_abm, "# t S E C I I_s I_c R_C R_I D\n");
        for (auto i = 0; i < m_result.get_num_time_points(); ++i) {
            fprintf(f_abm, "%f ", m_result.get_time(i));
            auto v = m_result.get_value(i);
            for (auto j = 0; j < v.size(); ++j) {
                fprintf(f_abm, "%f", v[j]);
                if (j < v.size() - 1) {
                    fprintf(f_abm, " ");
                }
            }
            if (i < m_result.get_num_time_points() - 1) {
                fprintf(f_abm, "\n");
            }
        }
        fclose(f_abm);
    }
    else {
        std::cout << "No results were stored. Set m_print_results to true to store results." << std::endl;
    }

    if (m_print_location_results) {
        for (auto it = m_location_result.begin(); it != m_location_result.end(); ++it) {
            char filepath[256];
            snprintf(filepath, sizeof(filepath), "location_%s.txt",
                     location_type_to_string((mio::abm::LocationType)it->first));
            auto f_abm_loc = fopen(filepath, "w");
            fprintf(f_abm_loc, "# t S E C I I_s I_c R_C R_I D\n");
            for (auto i = 0; i < it->second.get_num_time_points(); ++i) {
                fprintf(f_abm_loc, "%f ", it->second.get_time(i));
                auto v = it->second.get_value(i);
                for (auto j = 0; j < v.size(); ++j) {
                    fprintf(f_abm_loc, "%f", v[j]);
                    if (j < v.size() - 1) {
                        fprintf(f_abm_loc, " ");
                    }
                }
                if (i < it->second.get_num_time_points() - 1) {
                    fprintf(f_abm_loc, "\n");
                }
            }
            fclose(f_abm_loc);
        }
    }
    else {
        std::cout << "No location results were stored. Set print_location_results to true to store results."
                  << std::endl;
    }
}

const char* OutputResults::location_type_to_string(const mio::abm::LocationType l)
{
    switch (l) {
    case mio::abm::LocationType::Home:
        return "Home";
    case mio::abm::LocationType::School:
        return "School";
    case mio::abm::LocationType::Work:
        return "Work";
    case mio::abm::LocationType::SocialEvent:
        return "SocialEvent";
    case mio::abm::LocationType::BasicsShop:
        return "BasicsShop";
    case mio::abm::LocationType::Hospital:
        return "Hospital";
    case mio::abm::LocationType::Car:
        return "Car";
    case mio::abm::LocationType::ICU:
        return "ICU";
    case mio::abm::LocationType::PublicTransport:
        return "PublicTransport";
    case mio::abm::LocationType::TransportWithoutContact:
        return "TransportWithoutContact";
    default:
        return "Not_known_node";
    }
}

} // namespace abm
} // namespace mio
