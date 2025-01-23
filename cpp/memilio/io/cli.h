/* 
* Copyright (C) 2020-2025 MEmilio
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
#ifndef IO_CLI_H_
#define IO_CLI_H_

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "memilio/io/io.h"
#include "memilio/io/json_serializer.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/utils/parameter_set.h"

#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <functional>

namespace mio
{

namespace details
{

namespace cli
{

// an empty struct used to pass and iterate through types
template <class...>
struct typelist {
};

/// @brief a "parameter" (no Type) for the help option
struct Help {
    static const std::string name()
    {
        return "help";
    }
    static const std::string alias()
    {
        return "h";
    }
    static const std::string description()
    {
        return "Show this dialogue and exit.";
    }
};

/// @brief a "parameter" (no Type) for the print option option
struct PrintOption {
    static const std::string name()
    {
        return "print_option";
    }
    static const std::string description()
    {
        return "Use with parameter option name(s) without \"--\" as value(s). "
               "Prints the current values of specified options in their correct json format, "
               "then exits. Note that quotation marks may have to be escaped (\\\")";
    }
};

/// @brief a "parameter" (no Type) for the json import option
struct ReadFromJson {
    static const std::string name()
    {
        return "read_from_json";
    }
    static const std::string description()
    {
        return "Takes a filepath as value. Reads and assigns parameter option values from the specified json file.";
    }
};

/// @brief a "parameter" (no Type) for the json export option
struct WriteToJson {
    static const std::string name()
    {
        return "write_to_json";
    }
    static const std::string description()
    {
        return "Takes a filepath as value. Writes current values of all parameter options to the specified json file.";
    }
};

/// @brief Type list of options included by default (i.e. non-parameter options).
template <class... Others>
struct PresetOptions : public typelist<Help, PrintOption, ReadFromJson, WriteToJson, Others...> {
};

/// @brief A Field that gets the name of a Parameter
template <class Parameter>
struct Name {
    template <class T = Parameter>
    static const std::string get()
    {
        return T::name();
    }
};

/// @brief A Field that gets the alias of a Parameter, or an empty string if there is none
template <class Parameter>
struct Alias {
    template <class T>
    using alias_expr = decltype(T::alias());

    template <class T>
    constexpr static bool has_alias_v = mio::is_expression_valid<alias_expr, T>::value;

    template <class T = Parameter>
    static const std::string get()
    {
        if constexpr (has_alias_v<T>) {
            return T::alias();
        }
        else {
            return "";
        }
    }
};

/// @brief A Field that gets the description of a Parameter, or an empty string if there is none
template <class Parameter>
struct Description {
    template <class T>
    using description_expr = decltype(T::description());

    template <class T>
    constexpr static bool has_description_v = mio::is_expression_valid<description_expr, T>::value;

    template <class T = Parameter>
    static const std::string get()
    {
        if constexpr (has_description_v<T>) {
            return T::description();
        }
        else {
            return "";
        }
    }
};

/// @brief A Field that matches either name or alias when used as `NameOrAliasMatcher<T>::get() == a_string`.
template <class Parameter>
struct NameOrAliasMatcher {
    static NameOrAliasMatcher get()
    {
        return NameOrAliasMatcher<Parameter>{};
    }
    bool operator==(const std::string& other) const
    {
        return *this == other.data();
    }
    bool operator==(const char* other) const
    {
        return Name<Parameter>::get() == other || (Alias<Parameter>::get() != "" && Alias<Parameter>::get() == other);
    }
};

/// @brief An Identifier is a string used to identify an option. Defines what is considered a name or alias.
class Identifier : public std::string
{
public:
    using std::string::string;

    /// @brief check if the identifier is a name option
    inline bool is_name() const
    {
        return (size() > 2) && (substr(0, 2) == "--");
    }

    /// @brief view only the name of the identifier
    inline const char* get_name() const
    {
        assert(is_name());
        return data() + 2;
    }

    /// @brief check if the identifier is an alias option
    inline bool is_alias() const
    {
        return (size() > 1) && (data()[0] == '-') && (data()[1] != '-');
    }

    /// @brief view only the alias of the identifier
    inline const char* get_alias() const
    {
        assert(is_alias());
        return data() + 1;
    }

    /// @brief check if the identifier is either a name or an alias
    inline bool is_option() const
    {
        return is_name() || is_alias();
    }

    /// @brief check if the identifier matches the Parameter
    template <class Parameter>
    inline bool matches_parameter() const
    {
        if (is_name()) { // match name without "--"
            return Name<Parameter>::get() == substr(2);
        }
        else if (is_alias()) { // match alias without "-"
            return Alias<Parameter>::get() == substr(1);
        }
        else { // not an option
            return false;
        }
    }
};

/// @brief check if the given non-identifier string matches the Parameter by name or alias
template <class Parameter>
inline bool name_matches_parameter(const std::string& name)
{
    return (Name<Parameter>::get() == name) || (name.size() > 0 && Alias<Parameter>::get() == name);
}

/// @brief Helper struct to make write_help_impl write a seperator between presets and other options.
struct Seperator {
    static const std::string name()
    {
        return "\rParameter Options:";
    }
};

/// @brief End recursion
inline void write_help_impl(const typelist<>, std::ostream&)
{
}

/// @brief Recursively write a line of the --help output for Parameter T
template <class T, class... Ts>
void write_help_impl(const typelist<T, Ts...>, std::ostream& os)
{
    // Max. space and offsets used to print everything nicely.
    const size_t name_space         = 16; // reserved space for name
    const size_t alias_space        = 4; // reserved space for alias
    const size_t name_indent        = 4; // size of "  --"
    const size_t alias_indent       = 3; // size of "  -"
    const size_t description_indent = 2; // size of "  "
    // Get strings. Can assume name is not empty
    const std::string name        = Name<T>::get();
    const std::string description = Description<T>::get();
    const std::string alias       = Alias<T>::get();
    // Write name with "--" prefix and "  " indent
    os << "  --" << name;
    if (description.size() > 0 || alias.size() > 0) {
        if (name.size() <= name_space) {
            os << std::string(name_space - name.size(), ' ');
        }
        else {
            os << "\n" << std::string(name_space + name_indent, ' ');
        }
    }
    // Write alias (if available) with "-" prefix and "  " indent
    size_t space = std::max<size_t>(alias_space - alias.size(), 0);
    if (alias.size() > 0) {
        os << "  -" << alias;
    }
    else {
        space += alias_indent;
    }
    // Write description (if available) and end line indentation
    if (description.size() > 0) {
        os << std::string(space + description_indent, ' ') << description;
    }
    os << "\n";
    // Write next entry
    write_help_impl(typelist<Ts...>(), os);
}

/**
 * @brief Writes the output for the --help option.
 * @param[in] executable_name Name of the executable. Usually argv[0] is a good choice.
 * @param[in] parameters A set of parameters.
 * @param[in] os Any instance of std::ostream.
 * @tparam T List of all parameters
 * @tparam Set A parameter set.
 */
template <class... T, template <class...> class Set>
void write_help(const std::string& executable_name, const Set<T...>& /* parameters */, std::ostream& os)
{
    // help text preamble
    os << "Usage: " << executable_name << " <option> <value> ...\n"
       << "Values must be entered as json values, i.e. the expression to the right of \"Name : \".\n"
       << "Options:\n";
    // print options recursively
    write_help_impl(PresetOptions<Seperator, T...>(), os);
}

/// @brief End recursion
template <template <class> /*Field */ class, class Set>
mio::IOResult<void> set_param_impl(const typelist<>, Set&, const char* name, const Json::Value& /* args */)
{
    return mio::failure(mio::StatusCode::KeyNotFound, "No such option \"" + std::string(name) + "\".");
    // end recusrion
}

/// @brief Set the parameter given by name (or alias) to the json value given by args.
template <template <class> class Field, class T, class... Ts, class Set>
mio::IOResult<void> set_param_impl(const typelist<T, Ts...>, Set& parameters, const char* name, const Json::Value& args)
{
    if (Field<T>::get() == name) {
        // deserialize the json value to the parameter's type, then assign it
        auto result = mio::deserialize_json(args, mio::Tag<typename T::Type>());
        if (result) {
            // assign the result to the parameter
            parameters.template get<T>() = result.value();
            return mio::success();
        }
        else { // deserialize failed
            // insert more information to the error message
            std::string msg = "While setting \"" + std::string(name) + "\": " + result.error().message();
            return mio::failure(result.error().code(), msg);
        }
    }
    else {
        // try next parameter
        return set_param_impl<Field>(typelist<Ts...>(), parameters, name, args);
    }
}

/**
 * @brief Set the parameter given by identifer to the json value given by args.
 *
 * The identifier must start with either "--" if it is a parameter's name, or with "-" if it is an alias.
 * The function then recursively tries to match each parameter by this Field, and sets it accordningly.
 *
 * @param[in, out] parameters The set containing the parameter to set.
 * @param[in] identifier The identifier of the parameter.
 * @param[in] args Json value to be convertet to the matching parameter's type.
 * @tparam T List of all parameters.
 * @tparam Set A parameter set.
 * @return Nothing if successfull, an error code otherwise.
 */
template <class... T, template <class...> class Set>
mio::IOResult<void> set_param(Set<T...>& parameters, const Identifier& identifier, const Json::Value& args)
{
    // recursively match parameter either by name or alias and (try to) set it to args
    if (identifier.is_name()) {
        return set_param_impl<Name>(typelist<T...>(), parameters, identifier.get_name(), args);
    }
    else if (identifier.is_alias()) {
        return set_param_impl<Alias>(typelist<T...>(), parameters, identifier.get_alias(), args);
    }
    else {
        return mio::failure(mio::StatusCode::KeyNotFound, "Expected an option, got \"" + identifier + "\".");
    }
}

/**
 * @brief Set the parameter given by identifier to the json value given by args.
 *
 * The identifier must start with either "--" if it is a parameter's name, or with "-" if it is an alias.
 * The function then recursively tries to match each parameter by this Field, and sets it accordningly.
 *
 * @param[in, out] parameters The set containing the parameter to set.
 * @param[in] identifier The identifier of the parameter.
 * @param[in] args String containing a json value to be converted to the matching parameter's type.
 * @tparam T List of all parameters.
 * @tparam Set A parameter set.
 * @return Nothing if successfull, an error code otherwise.
 */
template <class... T, template <class...> class Set>
mio::IOResult<void> set_param(Set<T...>& parameters, const Identifier& identifier, const std::string& args)
{
    // read json value from args, then pass it to set_param(..., Json::Value)
    Json::Value js;
    std::string errors;
    Json::CharReaderBuilder builder;
    const std::unique_ptr<Json::CharReader> parser(builder.newCharReader());
    parser->parse(args.c_str(), args.c_str() + args.size(), &js, &errors);
    return set_param(parameters, identifier, js);
}

/// @brief End recursion.
template <class Set>
mio::IOResult<Json::Value> get_param_impl(typelist<>, Set&, const std::string& name)
{
    return mio::failure(mio::StatusCode::KeyNotFound, "No such option \"" + name + "\".");
}

/// @brief Get the value of the parameter specified by name.
template <class T, class... Ts, class Set>
mio::IOResult<Json::Value> get_param_impl(typelist<T, Ts...>, Set& parameters, const std::string& name)
{
    if (name_matches_parameter<T>(name)) {
        return mio::serialize_json(parameters.template get<T>());
    }
    else {
        return get_param_impl(typelist<Ts...>(), parameters, name);
    }
}

/**
 * @brief Get the parameter given by name.
 * @param[in] parameters The set containing the parameter to get.
 * @param[in] name The name (or alias) of the parameter.
 * @tparam T List of all parameters.
 * @tparam Set A parameter set.
 * @return The json value of the parameter if successfull, an error code otherwise.
 */
template <class... T, template <class...> class Set>
mio::IOResult<Json::Value> get_param(Set<T...>& parameters, const std::string& name)
{
    return get_param_impl(typelist<T...>(), parameters, name);
}

/// @brief Helper struct to provide the correct version of verify() at compile time.
struct OptionVerifier {
    template <template <class T> class Field, class T>
    using is_required = std::is_same<Field<T>, Name<T>>;

    template <template <class> class Field, class OptionA, class OptionB>
    inline static void verify()
    {
        const auto field_a = Field<OptionA>::get();
        const auto field_b = Field<OptionB>::get();
        if constexpr (std::is_same<OptionA, OptionB>::value) {
            assert((!is_required<Field, OptionA>::value || field_a != "") && "Option is missing required field.");
        }
        else {
            if constexpr (is_required<Field, OptionA>::value) {
                assert((field_a != field_b) && "Options may not have duplicate fields. (field required)");
            }
            else {
                assert((field_a == "" || field_b == "" || field_a != field_b) &&
                       "Options may not have duplicate fields. (field optional)");
            }
        }
    }
};

/// @brief End recursion
template <template <class> class Field, class... Ts>
inline void verify_options_impl(typelist<Ts...>, typelist<>, typelist<Ts...>)
{
}

/// @brief Restart recursion for the next A. Skips symmetric matches by reducing Ts.
template <template <class> class Field, class T, class... Ts, class A, class... As>
inline void verify_options_impl(typelist<T, Ts...>, typelist<A, As...>, typelist<>)
{
    // reduce both T- and A-list by the first element, reset B-list
    verify_options_impl<Field>(typelist<Ts...>(), typelist<As...>(), typelist<Ts...>());
}

/**
 * @brief Recursively assert that the parameters do not contain duplicate or empty required Fields.
 * @tparam Field The field of the option to compare.
 * @tparam Ts List of options. Used to reset the list B, Bs if empty.
 * @tparam A Option to verify.
 * @tparam As Queue of options to be verified recursively.
 * @tparam B Second option to verify A against.
 * @tparam Bs Queue of other options to verify against.
 */
template <template <class> class Field, class... Ts, class A, class... As, class B, class... Bs>
inline void verify_options_impl(typelist<Ts...>, typelist<A, As...>, typelist<B, Bs...>)
{
    // TODO: make this a compile time check, e.g. with c++20's constexpr c_str
    OptionVerifier::verify<Field, A, B>();
    // reduce B-list by the first element
    verify_options_impl<Field>(typelist<Ts...>(), typelist<A, As...>(), typelist<Bs...>());
}

/**
 * @brief Assert that no two Names or Aliases in the Set are the same, and that no name is empty.
 * @param[in] parameters Parameter set to verify.
 * @tparam T Lsit of all parameters.
 * @tparam Set A parameter set.
 */
template <class... T, template <class...> class Set>
void verify_options(Set<T...> /* parameters */)
{
    verify_options_impl<Name>(PresetOptions<T...>(), PresetOptions<T...>(), PresetOptions<T...>());
    verify_options_impl<Alias>(PresetOptions<T...>(), PresetOptions<T...>(), PresetOptions<T...>());
}

} // namespace cli

} // namespace details

/**
 * @brief Write a parameter set to the specified filepath.
 * @param[in] parameters An instance of Set.
 * @param[in] filepath The file to write the parameters into.
 * @tparam Set A parameter set.
 */
template <class Set>
IOResult<void> write_parameters_to_file(const Set& parameters, const std::string& filepath)
{
    return mio::write_json(filepath, parameters);
}

/**
 * @brief Read parameters from the specified file into the given parameter set.
 * @param[in, out] parameters An instance of Set<Parameters...>.
 * @param[in] filepath The file to read the parameters from.
 * @tparam Parameters A list of parameter types. 
 * @tparam Set A parameter set template. Will be used as Set<Parameters...>.
 */
template <class... Parameters, template <class...> class Set>
IOResult<void> read_parameters_from_file(Set<Parameters...>& parameters, const std::string& filepath)
{
    using namespace details::cli;
    // read file into json value
    auto json_result = mio::read_json(filepath);
    if (!json_result) {
        return mio::failure(json_result.error());
    }
    // set each parameter manually
    for (auto itr = json_result.value().begin(); itr != json_result.value().end(); itr++) {
        BOOST_OUTCOME_TRY(set_param_impl<details::cli::NameOrAliasMatcher>(typelist<Parameters...>{}, parameters,
                                                                           itr.name().c_str(), *itr));
    }
    return mio::success();
}

/**
 * @brief A cli that takes json values and stores them in a parameter set.
 * 
 * Note that the first element of argv will always be skipped, assuming it is the name of the executable
 * as called from the command line. executable_name simply allows to "clean up" this name, as it may contain a path.
 *
 * The Set template argument is expected to work like mio::ParameterSet.
 * In particular, if not using mio::ParameterSet, the Set is required to
 * a) have a list of (parameter) types as template arguments itself, e.g. ParameterSet<TypeA, TypeB,...>,
 * b) where each type has a member `Type` and `static const std::string name()`,
 * c) and Set has a member function `template <class T> T::Type& get()`.
 * d) Set has to implement a serialize function (see [README](cpp/memilio/io/README.md)).
 * Optionally, each type may have an "alias" or "description" member of the same type as "name".
 * Each name, alias or description should be ASCII only and not contain any control characters.
 *
 * @param[in] executable_name Name of the executable. Usually argv[0] from the main is a good choice.
 * @param[in] argc Argument count, must be the length of argv. Can be directly passed from main.
 * @param[in] argv Argument list for the Programm. Can be directly passed from main.
 * @param[in,out] parameters An instance of the parameter set.
 * @param[in] os Output stream to write to. Default: std::cout.
 * @tparam Set A parameter set.
 * @return Nothing if no errors occured, the error code otherwise.
 */
template <class Set>
mio::IOResult<void> command_line_interface(const std::string& executable_name, const int argc, char** argv,
                                           Set& parameters, std::ostream& os = std::cout)
{
    using namespace mio::details::cli;
    // make sure there are no duplicate names or aliases
    verify_options(parameters);
    // handle help option, scanning all argumemts before doing anything with them
    for (int i = 1; i < argc; i++) {
        if (Identifier(argv[i]).matches_parameter<Help>()) {
            // print the help dialogue and exit
            write_help(executable_name, parameters, os);
            std::exit(0);
        }
    }
    // handle print_option option, scanning all argumemts before doing anything with them
    for (int i = 1; i < argc; i++) {
        if (Identifier(argv[i]).matches_parameter<PrintOption>()) {
            i++; // skip the PrintOption argument
            for (; i < argc && !Identifier(argv[i]).is_option(); i++) {
                // try to get the parameter's json value
                BOOST_OUTCOME_TRY(auto&& value, get_param(parameters, argv[i]));
                // print the name (or alias) and value
                os << "Option " << argv[i] << ":\n" << value << "\n";
            }
            // exit after all values are printed
            std::exit(0);
        }
    }
    // handle parameter options that require values iteratively. assign given values or return an error
    int i = 1;
    while (i < argc) {
        const Identifier identifier(argv[i]);
        i++; // go to first argument
        // assert that the first argument is not an identifier (i.e. name or alias)
        if (i == argc || Identifier(argv[i]).is_option()) {
            return mio::failure(mio::StatusCode::OutOfRange, "Missing value for option \"" + identifier + "\".");
        }
        // collect all argv's that are not identifiers and set i to the position of the next identifier
        std::string arguments(argv[i]);
        i++;
        for (; (i < argc) && !Identifier(argv[i]).is_option(); i++) {
            arguments.append(" ").append(argv[i]);
        }
        // handle built-in options
        if (identifier.matches_parameter<ReadFromJson>()) {
            BOOST_OUTCOME_TRY(read_parameters_from_file(parameters, arguments));
        }
        else if (identifier.matches_parameter<WriteToJson>()) {
            BOOST_OUTCOME_TRY(write_parameters_to_file(parameters, arguments));
        }
        // (try to) set the parameter, to the value given by arguments
        else {
            BOOST_OUTCOME_TRY(set_param(parameters, identifier, arguments));
        }
    }
    return mio::success();
}

/**
 * @brief A cli that takes json values and stores them in a parameter set.
 * 
 * Note that the first element of argv will always be skipped, assuming it is the name of the executable
 * as called from the command line. executable_name simply allows to "clean up" this name, as it may contain a path.
 *
 * The Set template argument is expected to work like mio::ParameterSet.
 * In particular, if not using mio::ParameterSet, the Set is required to
 * a) have a list of (parameter) types as template arguments itself, e.g. ParameterSet<TypeA, TypeB,...>,
 * b) where each type has a member `Type` and `static const std::string name()`,
 * c) and Set has a member function `template <class T> T::Type& get()`.
 * d) Set has to implement a serialize function (see [README](cpp/memilio/io/README.md)).
 * Optionally, each type may have an "alias" or "description" member of the same type as "name".
 * Each name, alias or description should be ASCII only and not contain any control characters.
 *
 * @param[in] executable_name Name of the executable. Usually argv[0] from the main is a good choice.
 * @param[in] argc Argument count, must be the length of argv. Can be directly passed from main.
 * @param[in] argv Argument list for the Programm. Can be directly passed from main.
 * @param[in] os Output stream to write to. Default: std::cout.
 * @tparam Parameters A list of parameter types. 
 * @tparam Set A parameter set template. Will be used as Set<Parameters...>. Default: mio::ParameterSet.
 * @return An instance of Set<Parameters...> if no errors occured, the error code otherwise.
 */
template <class... Parameters, template <class...> class Set = ParameterSet>
mio::IOResult<Set<Parameters...>> command_line_interface(const std::string& executable_name, const int argc,
                                                         char** argv, std::ostream& os = std::cout)
{
    static_assert(sizeof...(Parameters) != 0, "At least one Parameter is required.");
    // create a new parameter set, and pass it and all other arguments to the main cli function
    Set<Parameters...> parameters{};
    auto result = command_line_interface(executable_name, argc, argv, parameters, os);
    // check the result, return parameters if appropriate
    if (result) {
        return mio::IOResult<Set<Parameters...>>(mio::success(std::move(parameters)));
    }
    else {
        return mio::IOResult<Set<Parameters...>>(mio::failure(result.error()));
    }
}

} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // IO_CLI_H_