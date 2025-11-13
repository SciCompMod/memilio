#pragma once

#include "boost/filesystem.hpp"

#include "memilio/epidemiology/damping.h"

#include "control_parameters/control_parameters.h"

#include <vector>
#include <string>
#include <set>

enum class Intervention
{
    Home,
    SchoolClosure,
    HomeOffice,
    GatheringBanFacilitiesClosure,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Count
};

enum class InterventionLevel
{
    Main,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Holidays,
    Count
};

enum class ContactLocation
{
    Home,
    School,
    Work,
    Other,
    Count
};

std::vector<DampingControlParameter> create_control_parameters(const boost::filesystem::path& interventions_file,
                                                               size_t num_age_groups)
{
    std::vector<DampingControlParameter> control_parameters;

    try {
        boost::property_tree::ptree pt;
        boost::property_tree::read_json(interventions_file.string(), pt);

        std::cout << "Available interventions from '" << interventions_file.string() << "':\n";
        std::set<std::string> intervention_names;
        for (const auto& item : pt) {
            try {
                std::string name = item.second.get<std::string>("name");
                intervention_names.insert(name);
                std::cout << " - " << name << "\n";
            }
            catch (...) {
                // ignore malformed entries
            }
        }

        // Helper lambda to reduce duplication
        auto add_control_parameter = [&](const std::string& name, const std::pair<double, double>& allowed_range,
                                         double effectiveness_value, double cost_value, InterventionLevel level_enum,
                                         Intervention type_enum, const std::vector<ContactLocation>& locations_vec) {
            if (intervention_names.find(name) != intervention_names.end()) {
                mio::DampingLevel level(static_cast<size_t>(level_enum));
                mio::DampingType type(static_cast<size_t>(type_enum));

                std::vector<size_t> locations;
                for (auto loc : locations_vec) {
                    locations.push_back(static_cast<size_t>(loc));
                }

                Eigen::VectorX<double> group_weights = Eigen::VectorX<double>::Constant(num_age_groups, 1.0);

                control_parameters.emplace_back(name, allowed_range, effectiveness_value, cost_value, level, type,
                                                locations, group_weights);

                std::cout << "Added control parameter for " << name << " intervention.\n";
            }
        };

        // clang-format off
        // Add control parameters here
        // name, range, effectiveness, cost, intervention level, intervention type, locations
        add_control_parameter("School closure", {0.0, 1.0}, 1.0, 1.0,
            InterventionLevel::Main, Intervention::SchoolClosure, {ContactLocation::School});

        add_control_parameter("Face masks & social distancing School", {0.0, 1.0}, 0.25, 1.0,
            InterventionLevel::PhysicalDistanceAndMasks, Intervention::PhysicalDistanceAndMasks, {ContactLocation::School});

        add_control_parameter("Face masks & social distancing Work", {0.0, 1.0}, 0.25, 1.0,
            InterventionLevel::PhysicalDistanceAndMasks, Intervention::PhysicalDistanceAndMasks, {ContactLocation::Work});

        add_control_parameter("Face masks & social distancing Other", {0.0, 1.0}, 0.35, 1.0,
            InterventionLevel::PhysicalDistanceAndMasks, Intervention::PhysicalDistanceAndMasks, {ContactLocation::Other});

        add_control_parameter("Remote work", {0.0, 1.0}, 0.25, 1.0, 
            InterventionLevel::Main, Intervention::HomeOffice, {ContactLocation::Work});
        // clang-format on
    }
    catch (const std::exception& error) {
        std::cerr << "Failed to read interventions file '" << interventions_file.string() << "': " << error.what()
                  << "\n";
    }

    return control_parameters;
}
