/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: René Schmieding
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
#include "memilio/io/cli.h"

#ifdef MEMILIO_HAS_JSONCPP

#include <ostream>
#include <sstream>
#include <string>
#include <vector>

/// @brief Part of write_help implementation. Writes the header with usage information and default options.
void write_help_preamble(const std::string& executable_name, const std::vector<std::string>& default_options,
                         std::ostream& os)
{
    os << "Usage: " << executable_name;
    if (default_options.size() > 0) {
        os << " [";
        size_t i = 0;
        for (; i < default_options.size() - 1; i++) {
            os << default_options[i] << " ";
        }
        os << default_options[i] << "]";
    }
    os << " <option> <value> ...\n";
    if (default_options.size() > 0) {
        os << "Values of parameter options listed in [brackets] can be optionally entered in that order.\n";
    }
    os << "All values must be entered as json values, i.e. the expression to the right of \"Name : \".\n"
       << "Note that, when entering values, quotation marks may have to be escaped (\\\").\n";
}

/// @brief Part of write_help implementation. Writes out a singular parameter option.
void write_help_parameter(const mio::cli::details::DatalessParameter& parameter, std::ostream& os)
{
    // Max. space and offsets used to print everything nicely.
    constexpr size_t name_space         = 16; // reserved space for name
    constexpr size_t alias_space        = 4; // reserved space for alias
    constexpr size_t name_indent        = 4; // size of "  --"
    constexpr size_t alias_indent       = 3; // size of "  -"
    constexpr size_t description_indent = 2; // size of "  "
    // Write name with "--" prefix and "  " indent
    os << "  --" << parameter.name;
    if (parameter.description.size() > 0 || parameter.alias.size() > 0) {
        if (parameter.name.size() <= name_space) {
            os << std::string(name_space - parameter.name.size(), ' ');
        }
        else {
            os << "\n" << std::string(name_space + name_indent, ' ');
        }
    }
    // Write alias (if available) with "-" prefix and "  " indent
    size_t space = parameter.alias.size() < alias_space ? alias_space - parameter.alias.size() : 0;
    if (parameter.alias.size() > 0) {
        os << "  -" << parameter.alias;
    }
    else {
        space += alias_indent;
    }
    // Write description (if available) and end line indentation
    if (parameter.description.size() > 0) {
        os << std::string(space + description_indent, ' ') << parameter.description;
    }
    os << "\n";
}

void mio::cli::details::write_help(const std::string& executable_name, const AbstractSet& set,
                                   const std::vector<std::string>& default_options, std::ostream& os)
{
    write_help_preamble(executable_name, default_options, os);
    os << "Options:\n";
    for (const auto& parameter : PresetOptions::all_presets) {
        write_help_parameter(parameter, os);
    }
    os << "Parameter options:\n";
    for (const auto& parameter : set.parameters()) {
        write_help_parameter(parameter, os);
    }
}

mio::IOResult<void> mio::cli::details::write_abstract_set_to_file(mio::cli::details::AbstractSet& set,
                                                                  const std::string& filepath)
{
    Json::Value output;
    for (auto& parameter : set.parameters()) {
        BOOST_OUTCOME_TRY(output[parameter.name()], parameter.get());
    }
    return mio::write_json(filepath, output);
}

mio::IOResult<void> mio::cli::details::read_abstract_set_from_file(mio::cli::details::AbstractSet& set,
                                                                   const std::string& filepath)
{
    // read file into json value
    auto json_result = mio::read_json(filepath);
    if (!json_result) {
        return mio::failure(json_result.error());
    }
    // set each parameter manually
    for (auto itr = json_result.value().begin(); itr != json_result.value().end(); itr++) {
        BOOST_OUTCOME_TRY(set.set_param(mio::cli::details::Identifier::make_raw(itr.name()), *itr));
    }
    return mio::success();
}

mio::IOResult<void> mio::cli::details::command_line_interface(const std::string& executable_name,
                                                              const std::span<char*>& argv,
                                                              cli::details::AbstractSet& set,
                                                              const std::vector<std::string>& default_options)
{
    assert(set.parameters().size() > 0 && "At least one parameter is required!");
    // this function glues all functionalities of the cli together. it may repeatedly iterate through all values of
    // argv (starting at 1).
    using namespace mio::cli;
    using namespace mio::cli::details;
    // verify that all default_options are parameter names
    for (const auto& option : default_options) {
        if (!set.contains(Identifier::make_raw(option))) {
            return failure(mio::StatusCode::KeyNotFound, "Default option \"" + option + "\" is not a parameter name.");
        }
    }
    // pre-scan all argumemts before doing anything with them to deal with help and print_option
    // this avoids returning e.g. parsing errors instead of the help dialogue
    for (auto arg_itr = argv.begin() + 1; arg_itr != argv.end(); ++arg_itr) {
        auto id_result = Identifier::parse(*arg_itr);
        // skip non-option arguments
        if (!id_result) {
            continue;
        }
        const auto& id = id_result.value();
        // handle help option
        if (id.matches_parameter(PresetOptions::help)) {
            // print the help dialogue and exit
            std::stringstream ss;
            write_help(executable_name, set, default_options, ss);
            return mio::failure(StatusCode::OK, std::move(ss.str()));
        }
        // handle print_option option
        else if (id.matches_parameter(PresetOptions::print_option)) {
            ++arg_itr; // skip the PrintOption argument
            std::stringstream ss;
            for (; arg_itr != argv.end() && !Identifier::is_option(*arg_itr); ++arg_itr) {
                // try to get the parameter's json value
                BOOST_OUTCOME_TRY(auto&& value, set.get_param(Identifier::make_raw(*arg_itr)));
                // print the name (or alias) and value
                ss << "Option " << *arg_itr << ":\n" << value << "\n";
            }
            // return after all values are printed
            return mio::failure(StatusCode::OK, ss.str());
        }
    }
    // main pass over all args to set options
    auto arg_itr = argv.begin() + 1;
    auto def_itr = default_options.begin();
    // handle parameter options that require values iteratively. assign given values or return an error
    while (arg_itr != argv.end()) {
        const auto id_result = Identifier::parse(*arg_itr);
        // try to parse the first default_options.size() as arguments; afterwards, require an identifier
        if (!id_result) {
            // checking #defaults suffices, as non-option arguments are greedily collected into "arguments" below
            if (def_itr != default_options.end()) {
                const auto& param_name = Identifier::make_raw(*def_itr);
                BOOST_OUTCOME_TRY(set.set_param(param_name, std::string(*arg_itr)));
                ++arg_itr;
                ++def_itr;
                continue;
            }
            else {
                return id_result.error();
            }
        }
        const Identifier current_option(id_result.value());
        ++arg_itr; // go to first argument
        // assert that the first argument is not an identifier (i.e. name or alias)
        if (arg_itr == argv.end() || Identifier::is_option(*arg_itr)) {
            return mio::failure(mio::StatusCode::OutOfRange,
                                "Missing value for option \"" + current_option.string + "\".");
        }
        // collect all argv's that are not identifiers and set i to the position of the next identifier
        std::string arguments(*arg_itr);
        ++arg_itr;
        for (; (arg_itr != argv.end()) && !Identifier::is_option(*arg_itr); ++arg_itr) {
            // here space separated args are joined together. maybe a better way is to make users use 'ticks' to group
            // their input.
            arguments.append(" ").append(*arg_itr);
        }
        // handle built-in options
        if (current_option.matches_parameter(PresetOptions::read_from_json)) {
            BOOST_OUTCOME_TRY(read_abstract_set_from_file(set, arguments));
        }
        else if (current_option.matches_parameter(PresetOptions::write_to_json)) {
            BOOST_OUTCOME_TRY(write_abstract_set_to_file(set, arguments));
        }
        // (try to) set the parameter, to the value given by arguments
        else {
            BOOST_OUTCOME_TRY(set.set_param(current_option, arguments));
        }
    }
    // check if required parameters were set, return an error if not
    if (std::ranges::any_of(set.parameters(), [](auto&& p) {
            return p.is_required();
        })) {
        std::stringstream ss;
        ss << "Missing values for required parameter(s):\n";
        for (const auto& p : set.parameters()) {
            if (p.is_required() == true) {
                ss << "  " << p.name();
            }
        }
        ss << "\n"
           << "Use \"" << executable_name << " --help\" for more info.\n";
        return mio::failure(mio::StatusCode::InvalidValue, ss.str());
    }
    return mio::success();
}

#endif // MEMILIO_HAS_JSONCPP
