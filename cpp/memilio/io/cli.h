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
#ifndef MIO_IO_CLI_H
#define MIO_IO_CLI_H

#include "memilio/config.h" // needed for defining MEMILIO_HAS_JSONCPP

#ifdef MEMILIO_HAS_JSONCPP

#include "memilio/io/io.h"
#include "memilio/io/json_serializer.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/string_literal.h"
#include "memilio/utils/type_list.h"

#include <cassert>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace mio
{

namespace cli
{

namespace details
{

/// @brief A Field that gets the name of a Parameter.
struct Name {
    template <class Parameter>
    static const std::string get(const Parameter& p)
    {
        return p.name();
    }
};

/// @brief A Field that gets the alias of a Parameter, or an empty string if there is none.
struct Alias {
    template <class T>
    using alias_expr = decltype(std::declval<T>().alias());

    template <class T>
    constexpr static bool has_alias_v = mio::is_expression_valid<alias_expr, T>::value;

    template <class Parameter>
    static const std::string get(const Parameter& p)
    {
        if constexpr (has_alias_v<Parameter>) {
            return p.alias();
        }
        else {
            return "";
        }
    }
};

/// @brief A Field that gets the description of a Parameter, or an empty string if there is none.
struct Description {
    template <class T>
    using description_expr = decltype(std::declval<T>().description());

    template <class T>
    constexpr static bool has_description_v = mio::is_expression_valid<description_expr, T>::value;

    template <class Parameter>
    static const std::string get(const Parameter& p)
    {
        if constexpr (has_description_v<Parameter>) {
            return p.description();
        }
        else {
            return "";
        }
    }
};

/// @brief A Field that gets the required flag of a Parameter, or false if there is none.
struct IsRequired {
    template <class T>
    using is_required_expr = decltype(std::declval<T>().is_required());

    template <class T>
    constexpr static bool has_is_required_v = is_expression_valid<is_required_expr, T>::value;

    template <class Parameter>
    constexpr static bool get(const Parameter& p)
    {
        if constexpr (has_is_required_v<Parameter>) {
            return p.is_required();
        }
        else {
            return false;
        }
    }
};

/// @brief Struct containing all non-data members of a parameter. Serves as base for AbstractParameter.
struct DatalessParameter {
    std::string name, alias = "", description = "";
    bool is_required = false;
};

/// @brief Subset of DatalessParameter, used to get user input in ParameterSetBuilder::add.
struct OptionalFields {
    std::string alias = "", description = "";
    bool is_required = false;
};

/// @brief Different types of identifiers, determines the content and handling of the Identifier's string.
enum class IdentifierType
{
    Name, // e.g. "--help"
    Alias, // e.g. "-h"
    Raw // either a name or alias without dashes, e.g. "help"
};

/**
 * @brief An Identifier is a string used to identify a parameter option.
 * This class defines what strings are considered a name or alias.
 * To create an Identifier, either parse a command line argument or make a raw identifier from a known-good string.
 */
class Identifier
{
private:
    /// @brief Basic constructor. Use via parse or make_raw.
    Identifier(std::string str, IdentifierType t)
        : string(std::move(str))
        , type(t)
    {
    }

public:
    /**
     * @brief Create an Identifier from a command line argument.
     * @param[in] str A string that is expected to be an option.
     * @return An Identifier if str is a valid option, an error if not.
     */
    static IOResult<Identifier> parse(const std::string& str)
    {
        if (is_name(str)) {
            return mio::success(Identifier(str, IdentifierType::Name));
        }
        else if (is_alias(str)) {
            return mio::success(Identifier(str, IdentifierType::Alias));
        }
        else {
            return mio::failure(mio::StatusCode::InvalidValue, "Expected an option, got \"" + str + "\".");
        }
    }

    /**
     * @brief Create an Identifier from a known-good option.
     * @param[in] str A string that is either a name or alias.
     * @return A raw Identifier of the given string.
     */
    static Identifier make_raw(const std::string& raw_name)
    {
        return Identifier(raw_name, IdentifierType::Raw);
    }

    /// @brief Check if a string is a name option.
    static bool is_name(std::string_view unverified_id)
    {
        // check for two dashes and at least one more character
        return (unverified_id.size() > 2) && (unverified_id[0] == '-' && unverified_id[1] == '-');
    }

    /// @brief Check if a string is an alias option.
    static bool is_alias(std::string_view unverified_id)
    {
        // check for exactly one dash and at least one more non-dash character
        return (unverified_id.size() > 1) && (unverified_id[0] == '-') && (unverified_id[1] != '-');
    }

    /// @brief Check if a string is either a name or alias option.
    static bool is_option(std::string_view unverified_id)
    {
        return is_name(unverified_id) || is_alias(unverified_id);
    }

    /// @brief Get a view of the underlying string, skipping any dashes.
    inline const std::string_view strip() const
    {
        switch (type) {
        case IdentifierType::Name:
            return std::string_view{string.data() + 2, string.size() - 2};
        case IdentifierType::Alias:
            return std::string_view{string.data() + 1, string.size() - 1};
        case IdentifierType::Raw:
            return std::string_view{string};
        }
        assert(false && "Invalid IdentifierType");
        // this return is practically inaccessible, but needed to avoid a "reached end of non-void function" warning
        return std::string_view{};
    }

    /// @brief Test whether the identifier matches a given parameter by name or alias.
    inline bool matches_parameter(const DatalessParameter& p) const
    {
        switch (type) {
        case IdentifierType::Name:
            return p.name == strip();
        case IdentifierType::Alias:
            return p.alias == strip();
        case IdentifierType::Raw:
            return (p.name == string) || (string.size() > 0 && p.alias == string);
        }
        assert(false && "Invalid IdentifierType");
        // this return is practically inaccessible, but needed to avoid a "reached end of non-void function" warning
        return false;
    }

    std::string string; ///< The underlying string, usually containing a command line argument.
    IdentifierType type; ///< The type set by parse or make_raw.
};

/// @brief Static container holding all preset options. Used for matching input arguments and writing help text.
struct PresetOptions {
    const inline static DatalessParameter help{"help", "h", "Show this dialogue and exit.", false};

    const inline static DatalessParameter print_option{
        "print_option", "",
        "Use with parameter option name(s) without \"--\" as value(s). Prints the current values of specified "
        "options in their correct json format, then exits.",
        false};

    const inline static DatalessParameter read_from_json{
        "read_from_json", "",
        "Takes a filepath as value. Reads and assigns parameter option values from the specified json file.", false};

    const inline static DatalessParameter write_to_json{
        "write_to_json", "",
        "Takes a filepath as value. Writes current values of all parameter options to the specified json file.", false};

    const inline static std::vector<DatalessParameter> all_presets{help, print_option, read_from_json, write_to_json};
};

/**
 * @brief A class to make both compile-time and run-time parameters usable in the CLI.
 * Data ownership is defined at runtime, so that it can be used as an interface for external parameter sets or directly
 * as a data class. 
 * Uses type erasure, so type information must be stored elsewhere for direct data access. 
 */
class AbstractParameter : public DatalessParameter
{
public:
    /**
     * @brief Create a parameter with type-erased data.
     * @tparam Type The original type of the pointer given by "value".
     * @param p Defines name, alias, description and is_required members of the parameter.
     * @param value A shared pointer to the parameter value. The deleter defines whether this is owning or non-owning.
     * @
     */
    template <class Type>
    AbstractParameter(mio::Tag<Type>, const DatalessParameter& p, std::shared_ptr<void>&& value)
        : DatalessParameter(p)
        , m_data(value)
        , m_serialize([](const std::shared_ptr<void>& param) -> IOResult<Json::Value> {
            return mio::serialize_json(*static_cast<Type*>(param.get()));
        })
        , m_set([](std::shared_ptr<void>& param, const std::string& name, const Json::Value& json) -> IOResult<void> {
            mio::unused(param); // prevent a "set, but unused" compiler warning
            // deserialize the json value to the parameter's type, then assign it
            auto param_result = mio::deserialize_json(json, mio::Tag<Type>());
            if (param_result) {
                // assign the json to the parameter
                *(static_cast<Type*>(param.get())) = std::move(param_result.value());
                return mio::success();
            }
            else { // deserialize failed
                // insert more information to the error message
                std::string msg = "While setting \"" + name + "\": " + param_result.error().message();
                return mio::failure(param_result.error().code(), msg);
            }
        })
    {
    }

    /**
     * @brief Create a parameter with non-owning data pointer.
     * Used when building an AbstractSet from mio::ParameterSet.
     */
    template <class Param>
    AbstractParameter(mio::Tag<Param>, typename Param::Type& value)
        : AbstractParameter(mio::Tag<typename Param::Type>{},
                            DatalessParameter{Name::get(Param{}), Alias::get(Param{}), Description::get(Param{}),
                                              IsRequired::get(Param{})},
                            std::shared_ptr<void>(static_cast<void*>(&value), [](void*) {}))
    {
    }

    /**
     * @brief Copy another AbstractParameter.
     * Used when building an AbstractSet from mio::cli::ParameterSet.
     */
    template <class Param>
    AbstractParameter(mio::Tag<Param>, AbstractParameter& other)
        : AbstractParameter(other)
    {
    }

    /**
     * @brief Retrieve the value of the parameter.
     * Note that this function only fails if serialization fails, which is very unlikely.
     * @return A Json representation of the parameter's value if successful, an error otherwise.
     */
    IOResult<Json::Value> get() const
    {
        auto json_result = m_serialize(m_data);
        if (json_result) {
            return json_result;
        }
        else {
            std::string msg = "While getting \"" + name() + "\": " + json_result.error().message();
            return mio::failure(json_result.error().code(), msg);
        }
    }

    /**
     * @brief Set the value of the parameter.
     * @param[in] value A Json representation of the new value.
     * @return Nothing if successful, an error otherwise.
     */
    IOResult<void> set(const Json::Value& value)
    {
        return m_set(m_data, name(), value);
    }

    /**
     * @brief Access to the DatalessParameter members mirroring a mio::ParameterSet parameter.
     * @{
     */
    const std::string& name() const
    {
        return DatalessParameter::name;
    }
    const std::string& alias() const
    {
        return DatalessParameter::alias;
    }
    const std::string& description() const
    {
        return DatalessParameter::description;
    }
    bool is_required() const
    {
        return DatalessParameter::is_required;
    }
    /** @} */

    /// @brief Access the raw data. Must be casted to the original type for safe usage.
    void* data()
    {
        return m_data.get();
    }

private:
    std::shared_ptr<void> m_data; ///< An owning or non-owning pointer to the parameter value.
    IOResult<Json::Value> (*m_serialize)(
        const std::shared_ptr<void>&); ///< Function pointer with type info needed for get implementation.
    IOResult<void> (*m_set)(std::shared_ptr<void>&, const std::string&,
                            const Json::Value&); ///< Function pointer with type info needed for set implementation.
};

class AbstractSet
{
public:
    using MapType = std::unordered_map<std::string_view, AbstractParameter&>;

    /**
     * @brief Create a new AbstractSet from a given set of parameters.
     * Keep in mind that the AbstractSet usually uses the original set's storage, only acting as a type erased
     * interface. Therefore the AbstractSet must not outlive the original set.
     * This method verifies, that the set does not contain parameters with empty names, or duplicate names or aliases. 
     * @tparam Set Either a class similar to mio::ParameterSet, or a vector of AbstractParameters.
     * @param[in] parameter_set The original set of parameters. Mind its lifetime!
     * @return A new AbstractSet if the parameters are valid, an error otherwise.
     */
    template <class Set>
    static IOResult<AbstractSet> build(Set& parameter_set)
    {
        AbstractSet set(parameter_set);
        BOOST_OUTCOME_TRY(fill_maps(set));
        return mio::success(std::move(set));
    }

    /**
     * @brief Get a parameter's value.
     * @param[in] id An Identifier matching the parameter to get.
     * @return The parameter's value in Json representation if successful, an error otherwise.
     */
    IOResult<Json::Value> get_param(const Identifier& id)
    {
        auto param = find(id);
        if (!param) {
            return IOResult<Json::Value>(param.error().code(), "Could not get parameter: " + param.error().message());
        }
        else {
            return param.value()->second.get();
        }
    }

    /**
     * @brief Set a parameter's value.
     * @param[in] id An Identifier matching the parameter to be set.
     * @param[in] args A string containing a Json representation of the new value.
     * @return Nothing if successful, an error otherwise.
     */
    IOResult<void> set_param(const Identifier& id, const std::string& args)
    {
        Json::Value js;
        std::string errors;
        Json::CharReaderBuilder builder;
        const std::unique_ptr<Json::CharReader> parser(builder.newCharReader());
        parser->parse(args.c_str(), args.c_str() + args.size(), &js, &errors);
        // do not directly raise errors, to avoid hiding e.g. a "parameter not found"
        return set_param(id, js, errors);
    }

    /**
     * @brief Set a parameter's value.
     * @param[in] id An Identifier matching the parameter to be set.
     * @param[in] value The new value for the parameter in Json representation.
     * @param[in] parse_errors Optional argument with errors from parsing the Json value.
     * @return Nothing if successful, an error otherwise.
     */
    IOResult<void> set_param(const Identifier& id, const Json::Value& value, const std::string& parse_errors = "")
    {
        auto param = find(id);
        if (!param) {
            return IOResult<void>(param.error().code(), "Could not set parameter: " + param.error().message());
        }
        else {
            param.value()->second.DatalessParameter::is_required = false; // mark as set
            // try to set the value. append parsing errors if not successful
            IOResult<void> result = param.value()->second.set(value);
            if (result || parse_errors.empty()) {
                return result;
            }
            else {
                return mio::failure(result.error().code(), result.error().message() + "\n" + parse_errors);
            }
        }
    }

    /// @brief Test if the set contains an item matching the given id.
    bool contains(const Identifier& id)
    {
        return static_cast<bool>(find(id));
    }

    /**
     * @brief Access all AbstractParameters.
     * @{
     */
    std::vector<AbstractParameter>& parameters()
    {
        return m_parameters;
    }
    const std::vector<AbstractParameter>& parameters() const
    {
        return m_parameters;
    }
    /** @} */

private:
    /// @brief Do not allow creating sets without build method.
    AbstractSet() = default;

    /// @brief Construct AbstractSet with non-owning Parameters from mio::ParameterSet or similar.
    template <class... T, template <class...> class Set>
    AbstractSet(Set<T...>& parameter_set)
        : m_parameters(std::vector{AbstractParameter(mio::Tag<T>{}, parameter_set.template get<T>())...})
    {
    }

    /// @brief Construct a set directly from a vector of AbstractParameter%s.
    AbstractSet(std::vector<AbstractParameter>& parameters)
        : m_parameters(parameters)
    {
    }

    /**
     * @brief Obtain a reference to a parameter matching the given Identifier.
     * Uses maps instead of linear search. Checks the potentially smaller alias map first, then the name map.
     * @param[in] id An Identifier.
     * @return A map iterator if a matching parameter was found, an error otherwise.
     */
    IOResult<MapType::iterator> find(const Identifier& id)
    {
        MapType::iterator param_itr = m_map_by_alias.find(id.strip());
        if (param_itr != m_map_by_alias.end()) {
            return mio::success(param_itr);
        }
        param_itr = m_map_by_name.find(id.strip());
        if (param_itr != m_map_by_name.end()) {
            return mio::success(param_itr);
        }
        return mio::failure(mio::StatusCode::KeyNotFound, "No such option \"" + id.string + "\".");
    }

    /**
     * @brief Initialize the maps of an abstract set used for finding parameters.
     * Also verifies that there are no empty names, or duplicate names or aliases.
     * @param[in] set A newly constructed set with all parameters and default-initialized (empty) maps.
     * @return Nothing if parameters are good, error otherwise.
     */
    static IOResult<void> fill_maps(AbstractSet& set)
    {
        for (AbstractParameter& p : set.m_parameters) {
            if (p.name().empty()) {
                return mio::failure(mio::StatusCode::InvalidValue, "An option is missing a name.");
            }
            if (!set.m_map_by_name.insert({p.name(), p}).second) {
                return mio::failure(mio::StatusCode::InvalidValue, "Options may not have duplicate names. "
                                                                   "Found two instances of \"" +
                                                                       p.name() + "\".");
            }
            if (!p.alias().empty() && !set.m_map_by_alias.insert({p.alias(), p}).second) {
                return mio::failure(mio::StatusCode::InvalidValue, "Options may not have duplicate aliases. "
                                                                   "Found two instances of \"" +
                                                                       p.alias() + "\".");
            }
        }
        return mio::success();
    }

    std::vector<AbstractParameter> m_parameters;
    MapType m_map_by_name;
    MapType m_map_by_alias;
};

/**
 * @brief Writes the output for the --help option.
 * @param[in] executable_name Name of the executable. Usually argv[0] is a good choice.
 * @param[in] set A set of parameters.
 * @param[in] os Any instance of std::ostream.
 * @tparam T List of all parameters
 * @tparam Set A parameter set.
 */
void write_help(const std::string& executable_name, const AbstractSet& set,
                const std::vector<std::string>& default_options, std::ostream& os);

/// @brief Implementation of the CLI. See the main function below for details.
mio::IOResult<void> command_line_interface(const std::string& executable_name, const int argc, char** argv,
                                           cli::details::AbstractSet& parameters,
                                           const std::vector<std::string>& default_options);

/// @brief Helper to allow storing a StringLiteral in a TypeList.
template <StringLiteral>
struct NamedType;

/**
 * @brief A minimalistic class similar to mio::ParameterSet intended for use with the CLI.
 * Uses StringLiteral%s to restore type information from AbstractParameter%s for outside access.
 */
template <class... Tags>
class ParameterSet
{
    using Names = TypeList<type_at_index_t<0, Tags>...>;
    using Types = TypeList<type_at_index_t<1, Tags>...>;

public:
    /**
     * @brief Create a ParameterSet. Use the ParameterSetBuilder instead of calling this directly.
     * @param[in] parameters A vector of AbstractParameter%s holding their own data.
     */
    ParameterSet(std::vector<details::AbstractParameter>&& parameters)
        : m_parameters(parameters)
    {
    }

    /**
     * @brief Access a parameter value by its name.
     * @tparam Name The name of the parameter.
     * @return A reference to the parameter's value.
     */
    template <mio::StringLiteral Name>
    inline auto& get()
    {
        constexpr size_t parameter_index = index_of_type_v<details::NamedType<Name>, Names>;
        using ReturnType                 = type_at_index_t<parameter_index, Types>;
        return *static_cast<ReturnType*>(m_parameters[parameter_index].data());
    }

    /**
     * @brief Direct access to a parameter onject. Only intended for building an AbstractSet in the CLI.
     * @tparam Tag The Tag corresponding to a parameter.
     * @return A reference to a parameter.
     */
    template <class Tag>
    inline details::AbstractParameter& get()
    {
        constexpr size_t parameter_index = index_of_type_v<Tag, Tags...>;
        return m_parameters[parameter_index];
    }

private:
    std::vector<details::AbstractParameter> m_parameters; ///< Vector of parameters with ownership over their values.
};

/**
 * @brief Write a parameter set to the specified filepath.
 * @param[in] parameters An instance of Set.
 * @param[in] filepath The file to write the parameters into. 
 */
mio::IOResult<void> write_abstract_set_to_file(mio::cli::details::AbstractSet& set, const std::string& filepath);

/**
 * @brief Read parameters from the specified file into the given parameter set.
 * @param[in, out] parameters An instance of Set<Parameters...>.
 * @param[in] filepath The file to read the parameters from.
 */
mio::IOResult<void> read_abstract_set_from_file(mio::cli::details::AbstractSet& set, const std::string& filepath);

} // namespace details

/// @brief A builder for a parameter set to be used with the CLI, as direct construction of such a set is inconvenient.
template <class... Params>
class ParameterSetBuilder
{
    // make private constructor accessible independent of template arguments
    template <class...>
    friend class ParameterSetBuilder;

public:
    /// @brief Set up a builder for a parameter set that can be used with the CLI.
    ParameterSetBuilder() = default;

    /**
     * @brief Add a new parameter to the builder. 
     * Call this directly on a new builder or another add() call. Finalize with the build() method.
     * Do not store intermediate results.
     * Usage example:
     * ```
     * auto parameters = ParameterSetBuilder().add<"A">(1).add<"B">(2).build();
     * ```
     * @tparam Name The name of the new parameter as a StringLiteral. Must be unique in this set.
     * @tparam Type The type of the new parameter value. Use simple types for better CLI usage.
     * @param[in] initial_value An initial value for the parameter.
     * @param[in] optionals Can be used to set alias, description and is_required. 
     * @return Returns a builder object. Use add() on it to add more parameters, or use build() to get a parameter set.
     */
    template <mio::StringLiteral Name, class Type>
    [[nodiscard]] inline auto add(Type&& initial_value, details::OptionalFields&& optionals = {}) &&
    {
        using ValueType = std::decay_t<Type>; // get base type in case Type was deduced and is e.g. const or &
        // since we get *this as rvalue, we can move the parameters
        auto new_params = std::move(m_parameters);
        // create a new owning data pointer, stored as void*
        std::shared_ptr<void> data(new ValueType(initial_value), std::default_delete<ValueType>{});
        // create a new abstract parameter, then move all parameters to a new builder
        new_params.emplace_back(Tag<ValueType>{},
                                details::DatalessParameter{std::string(Name), optionals.alias, optionals.description,
                                                           optionals.is_required},
                                std::move(data));
        return ParameterSetBuilder<Params..., TypeList<details::NamedType<Name>, ValueType>>{std::move(new_params)};
    }

    /// @brief Finalize the builder and create a parameter set.
    [[nodiscard]] inline mio::cli::details::ParameterSet<Params...> build() &&
    {
        return mio::cli::details::ParameterSet<Params...>(std::move(m_parameters));
    }

private:
    /// @brief Private constructor used by "add".
    ParameterSetBuilder(std::vector<details::AbstractParameter>&& parameters)
        : m_parameters(std::move(parameters))
    {
    }

    std::vector<details::AbstractParameter> m_parameters; ///< Vector of parameters, holding initial values.
};

} // namespace cli

/**
 * @brief Write a parameter set to the specified filepath.
 * @param[in] parameters An instance of Set.
 * @param[in] filepath The file to write the parameters into.
 * @tparam Set A parameter set.
 */
template <class Set>
IOResult<void> write_parameters_to_file(Set& parameters, const std::string& filepath)
{
    BOOST_OUTCOME_TRY(auto&& set, cli::details::AbstractSet::build(parameters));
    return write_abstract_set_to_file(set, filepath);
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
    BOOST_OUTCOME_TRY(auto&& set, cli::details::AbstractSet::build(parameters));
    return read_abstract_set_from_file(set, filepath);
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
 * b) where each type has a member `static const std::string name()`,
 * c) and Set has a member function `template <class T> T::Type& get()`.
 * d) Each member stored in Set must be serializable and deserializable (see documentation).
 * Optionally, each parameter may have an "alias" or "description" member of the same type as "name".
 * Each name, alias or description should be ASCII only and not contain any control characters. Names and aliases must
 * be unique.
 * Optionally, each parameter may have "is_required" boolean member function. If any required parameter is not set by
 * the given arguments, an error is returned.
 *
 * Important: Always check the result as follows, since it is used for all CLI output:
 * ```
 * auto result = mio::command_line_interface( ... );
 * if (!result) {
 *     std::cout << result.error().message();
 *     return result.error().code().value();
 * }
 * ```
 *
 * @param[in] executable_name Name of the executable. Usually argv[0] from the main is a good choice.
 * @param[in] argc Argument count, must be the length of argv. Can be directly passed from main.
 * @param[in] argv Argument list for the programm. Can be directly passed from main.
 * @param[in,out] parameters An instance of the parameter set.
 * @param[in] default_options Parameter names that allow setting parameter values with leading non-option arguments.
 * @tparam Set A parameter set.
 * @return Nothing if no errors occured, the error code and message (including the help dialogue!) otherwise.
 */
template <class Set>
mio::IOResult<void> command_line_interface(const std::string& executable_name, const int argc, char** argv,
                                           Set& parameters, const std::vector<std::string>& default_options = {})
{
    // parse the parameters into an iterable format
    BOOST_OUTCOME_TRY(auto&& set, cli::details::AbstractSet::build(parameters));
    return cli::details::command_line_interface(executable_name, argc, argv, set, default_options);
}

/**
 * @brief A cli that takes json values and stores them in a parameter set.
 * 
 * This overload creates a new parameter set instead of using an existing one.
 * See the other overload for more details.
 *
 * @param[in] executable_name Name of the executable. Usually argv[0] from the main is a good choice.
 * @param[in] argc Argument count, must be the length of argv. Can be directly passed from main.
 * @param[in] argv Argument list for the programm. Can be directly passed from main.
 * @tparam Parameters A list of parameter tags. 
 * @tparam Set A parameter set template. Will be used as Set<Parameters...>. Default: mio::ParameterSet.
 * @return An instance of Set<Parameters...> if no errors occured, the error code and message otherwise.
 */
template <class... Parameters, template <class...> class Set = mio::ParameterSet>
mio::IOResult<Set<Parameters...>> command_line_interface(const std::string& executable_name, const int argc,
                                                         char** argv,
                                                         const std::vector<std::string>& default_options = {})
{
    static_assert(sizeof...(Parameters) != 0, "At least one Parameter is required.");
    // create a new parameter set, and pass it and all other arguments to the main cli function
    Set<Parameters...> parameters{};
    auto result = command_line_interface(executable_name, argc, argv, parameters, default_options);
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

#endif // MIO_IO_CLI_H
