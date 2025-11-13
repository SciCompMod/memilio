#pragma once

#include "boost/filesystem.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "constraints/constraints.h"

#include <vector>
#include <string>
#include <optional>
#include <iostream>

inline void load_constraints_from_file(const boost::filesystem::path& constraints_file,
                                       std::vector<Constraint>& path_constraints,
                                       std::vector<Constraint>& terminal_constraints)
{
    try {
        boost::property_tree::ptree root;
        boost::property_tree::read_json(constraints_file.string(), root);

        // Helper lambda for reading constraint arrays
        auto load_constraints = [](const boost::property_tree::ptree& array) {
            std::vector<Constraint> constraints;
            for (const auto& item : array) {
                const auto& obj = item.second;

                std::string name = obj.get<std::string>("name");
                double min_val   = obj.get<double>("min");
                double max_val   = obj.get<double>("max");

                // Handle optional or null node_id
                std::optional<size_t> node_id_opt = std::nullopt;
                auto node_id_val                  = obj.get_optional<size_t>("node_id");
                if (node_id_val.has_value()) {
                    node_id_opt = std::optional<size_t>(*node_id_val);
                }

                constraints.emplace_back(name, std::make_pair(min_val, max_val), node_id_opt);
            }
            return constraints;
        };

        // Load path_constraints and terminal_constraints (if present)
        if (auto pc = root.get_child_optional("path_constraints")) {
            path_constraints = load_constraints(*pc);
        }

        if (auto tc = root.get_child_optional("terminal_constraints")) {
            terminal_constraints = load_constraints(*tc);
        }

        // Optional: Print debug info
        std::cout << "Loaded " << path_constraints.size() << " path constraints.\n";
        for (const auto& c : path_constraints) {
            std::cout << " - " << c.name() << " [" << c.min() << ", " << c.max() << "]";
            if (c.node_id())
                std::cout << " node_id=" << *c.node_id();
            std::cout << "\n";
        }

        std::cout << "Loaded " << terminal_constraints.size() << " terminal constraints.\n";
        for (const auto& c : terminal_constraints) {
            std::cout << " - " << c.name() << " [" << c.min() << ", " << c.max() << "]";
            if (c.node_id())
                std::cout << " node_id=" << *c.node_id();
            std::cout << "\n";
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error reading constraints from '" << constraints_file.string() << "': " << e.what() << std::endl;
    }
}
