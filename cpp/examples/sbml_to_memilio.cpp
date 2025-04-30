#include "memilio/config_internal.h"
#include "memilio/utils/logging.h"
#include "memilio/io/io.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <sbml/SBMLTypes.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

/**
*    @brief Get the filename
*
*    @param filepath A path including the filename
*
*    Uses `boost::filesystem::path` to get the path to the input file, then uses `stem()` to get the name of the input 
*    file without file ending. 
*    This function is regularly used to get the namespace used for the generated code.
**/
std::string get_filename(const char* filepath)
{
    boost::filesystem::path p(filepath);
    return p.stem().string();
}

/**
*    @brief Creates a folder with the given name to store the generated model files
*
*    @param filename The name of the folder to be created
*
*    Extracts the filename using :cpp:func:`get_filename(const char* filename)`, converts it to lower case and creates 
*    a folder with the resulting name using `boost::filesystem`.
**/
void create_folder(const std::string& filename)
{
    std::string lowercase_name = boost::to_lower_copy<std::string>(filename);
    boost::filesystem::path p(lowercase_name);
    boost::filesystem::create_directory(p);
    mio::log_info("Creating folder at ./{}", p.string());
}

/**
*    @brief Formats a formula for the example.cpp file
*
*    @param model The model where the formula stems from
*    @param math_string The formula to be processed, preprocessed as char*
*    @param name_namespace The namespace of the model
*
*    This function goes through the parts of the formula as identified by spaces. It first strips leading/trailing 
*    braces and leading minus sings. Then it checks if the part is a parameter, species or compartment. If it is a 
*    parameter, it is replaced by the corresponding parameter in the model. If it is a species, it is replaced by the 
*    corresponding species in the model. If it is a compartment, it is replaced by the size of the compartment. If it is
*    pi, it is replaced by M_PI. If it is time, it is replaced by 0.0. The leading and trailing braces are then 
*    re-added to the part and the parts are joined together.
**/
std::string format_main_formula(Model* model, char* math_string, std::string* name_namespace)
{
    std::string return_string = "";
    std::vector<std::string> formula_parts;
    boost::split(formula_parts, math_string, boost::is_any_of(" "), boost::algorithm::token_compress_on);
    for (size_t i = 0; i < formula_parts.size(); i++) {
        std::string leading_bracket = "", trailing_bracket = "";
        while (formula_parts[i].find_first_of('(') != std::string::npos) {
            size_t first_index = formula_parts[i].find_first_of('(');
            leading_bracket    = leading_bracket + formula_parts[i].substr(0, first_index + 1);
            formula_parts[i].erase(0, first_index + 1);
        }
        if (formula_parts[i][0] == '-') {
            leading_bracket = leading_bracket + "-";
            formula_parts[i].erase(0, 1);
        }
        while (formula_parts[i].back() == ')' || formula_parts[i].back() == ',') {
            trailing_bracket = trailing_bracket + formula_parts[i].back();
            formula_parts[i].pop_back();
        }
        if (model->getListOfParameters()->getElementBySId(formula_parts[i]) != NULL) {
            formula_parts[i] = "params.template get<mio::" + *name_namespace + "::" + formula_parts[i] + "<double>>()";
        }
        else if (model->getListOfSpecies()->getElementBySId(formula_parts[i]) != NULL) {
            formula_parts[i] =
                "model.populations[{mio::" + *name_namespace + "::InfectionState::" + formula_parts[i] + "}]";
        }
        else if (model->getListOfCompartments()->getElementBySId(formula_parts[i]) != NULL) {
            double size =
                Compartment_getSize((Compartment*)model->getListOfCompartments()->getElementBySId(formula_parts[i]));
            formula_parts[i] = std::to_string(size);
        }
        if (formula_parts[i] == "pi") {
            formula_parts[i] = "M_PI";
        }
        if (formula_parts[i] == "time") {
            formula_parts[i] = "0.0";
        }
        if (leading_bracket.size() > 0) {
            formula_parts[i] = leading_bracket + formula_parts[i];
        }
        if (trailing_bracket.size() > 0) {
            formula_parts[i] = formula_parts[i] + trailing_bracket;
        }
    }
    return boost::algorithm::join(formula_parts, " ");
}

/**
*    @brief Returns the initial assignment of a species for the example.cpp file
*
*    @param species The species to be processed
*    @param model The model where the species stems from
*    @param name_namespace The namespace of the model
*
*    First it checks whether an initial amount is given for the species. Then it checks the other possibilites for 
*    initial assignements:
*    - If an initial concentration is given, it is used
*    - If an initial assignment is given, it is used. This could be
*        - a values
*        - a other parameter
*        - a formula, which would then be parsed using :cpp:func:`format_main_formula(Model * model, char* math_string, std::string * name_namespace)`
*
*    If none of these are given, it returns "NAN", otherwise the identified value.
**/
std::string get_initial_assignment(Species* species, Model* model, std::string* name_namespace)
{
    double base_value = species->getInitialAmount();
    if (!isnan(species->getInitialConcentration())) {
        base_value = species->getInitialConcentration();
    }
    auto initial_assignment = model->getInitialAssignmentBySymbol(species->getId());
    if (initial_assignment != NULL) {
        auto formula = initial_assignment->getMath();
        if (formula->isNumber()) {
            base_value = initial_assignment->getMath()->getValue();
        }
        else if (formula->isName()) {
            return "model.parameters.get<mio::" + *name_namespace + "::" + formula->getName() + "<double>>()";
        }
        else {
            return format_main_formula(model, SBML_formulaToL3String(formula), name_namespace);
        }
    }
    if (isnan(base_value)) {
        return "NAN";
    }
    return std::to_string(base_value);
}

/**
*    @brief Returns the string to change variables inside an event for the example.cpp file
*
*    @param formula The formula to be processed
*    @param model The model where the formula stems from
*    @param name_namespace The namespace of the model
*
*    Checks whether the formula is a parameter, species or compartment. In the first two cases it returns the code to
*    get the corresponding values, in the last case it prints an error message. 
**/
mio::IOResult<std::string> format_event_variable(std::string formula, Model* model, std::string* name_namespace)
{
    std::string return_string = "";

    if (model->getListOfParameters()->getElementBySId(formula) != NULL) {
        return_string = "sim.get_model().parameters.get<mio::" + *name_namespace + "::" + formula + "<ScalarType>>()";
    }
    else if (model->getListOfSpecies()->getElementBySId(formula) != NULL) {
        return_string =
            "sim.get_result().get_last_value()[sim.get_model().populations.get_flat_index(mio::" + *name_namespace +
            "::InfectionState::" + formula + ")]";
    }
    else if (model->getListOfCompartments()->getElementBySId(formula) != NULL) {
        mio::log_error("Unfortunately compartments can not be changed at the moment.");
        return mio::failure(mio::StatusCode::InvalidValue);
    }
    return mio::success(return_string);
}

/**
*    @brief Formats the formula for an event for use in the example.cpp file
*
*    @param formula The formula to be processed
*    @param model The model where the formula stems from
*    @param name_namespace The namespace of the model
*
*    Goes through the parts of the formula as identified by spaces. It first strips leading/trailing braces and leading
*    minus signs. Then it checks if the part is a parameter, species or compartment. If it is a parameter, it is 
*    replaced by the corresponding parameter in the model. If it is a species, it is replaced by the corresponding 
*    species in the model. If it is a compartment, it is replaced by the size of the compartment. If it is pi, it is 
*    replaced by M_PI. If it is time, it is replaced by t. The leading and trailing braces are then re-added to the 
*    part and the parts are joined together.
**/
mio::IOResult<std::string> format_event_formulas(std::string formula, Model* model, std::string* name_namespace)
{
    std::string return_string = "";
    std::vector<std::string> formula_parts;
    boost::split(formula_parts, formula, boost::is_any_of(" "), boost::algorithm::token_compress_on);
    for (size_t i = 0; i < formula_parts.size(); i++) {
        std::string leading_bracket = "", trailing_bracket = "";
        while (formula_parts[i].find_first_of('(') != std::string::npos) {
            size_t first_index = formula_parts[i].find_first_of('(');
            leading_bracket    = leading_bracket + formula_parts[i].substr(0, first_index + 1);
            formula_parts[i].erase(0, first_index + 1);
        }
        if (formula_parts[i][0] == '-') {
            leading_bracket = leading_bracket + "-";
            formula_parts[i].erase(0, 1);
        }
        while (formula_parts[i].back() == ')' || formula_parts[i].back() == ',') {
            trailing_bracket = trailing_bracket + formula_parts[i].back();
            formula_parts[i].pop_back();
        }

        if (model->getListOfParameters()->getElementBySId(formula_parts[i]) != NULL) {
            formula_parts[i] =
                "sim.get_model().parameters.get<mio::" + *name_namespace + "::" + formula_parts[i] + "<ScalarType>>()";
        }
        else if (model->getListOfSpecies()->getElementBySId(formula_parts[i]) != NULL) {
            formula_parts[i] =
                "sim.get_result().get_last_value()[sim.get_model().populations.get_flat_index(mio::" + *name_namespace +
                "::InfectionState::" + formula_parts[i] + ")]";
        }
        else if (model->getListOfCompartments()->getElementBySId(formula_parts[i]) != NULL) {
            double size =
                Compartment_getSize((Compartment*)model->getListOfCompartments()->getElementBySId(formula_parts[i]));
            formula_parts[i] = std::to_string(size);
        }
        if (formula_parts[i] == "pi") {
            formula_parts[i] = "M_PI";
        }
        if (formula_parts[i] == "time") {
            formula_parts[i] = "t";
        }
        if (leading_bracket.size() > 0) {
            formula_parts[i] = leading_bracket + formula_parts[i];
        }
        if (trailing_bracket.size() > 0) {
            formula_parts[i] = formula_parts[i] + trailing_bracket;
        }
    }
    return mio::success(boost::algorithm::join(formula_parts, " "));
}

/**
*    @brief Formats the trigger formula for an event for use in the example.cpp file
*
*    @param formula The formula to be processed
*    @param model The model where the formula stems from
*    @param name_namespace The namespace of the model
*
*    Goes through the parts of the formula as identified by spaces.
*    - It first strips leading/trailing braces and leading minus sings. 
*    - Then it checks if the part is a parameter, species or compartment. 
*        - If it is a parameter, it is replaced by the corresponding parameter in the model. 
*        - If it is a species, it prints an error. 
*        - If it is a compartment, it is replaced by the size of the compartment. 
*    - If it is pi, it is replaced by M_PI. 
*    - If it is time, it prints an error. 
*    The leading and trailing braces are then re-added to the part and the parts are joined together.
**/
mio::IOResult<std::string> format_event_trigger(std::string formula, Model* model, std::string* name_namespace)
{
    std::string return_string = "";
    std::vector<std::string> formula_parts;
    boost::split(formula_parts, formula, boost::is_any_of(" "), boost::algorithm::token_compress_on);
    for (size_t i = 0; i < formula_parts.size(); i++) {
        std::string leading_bracket = "", trailing_bracket = "";
        while (formula_parts[i].find_first_of('(') != std::string::npos) {
            size_t first_index = formula_parts[i].find_first_of('(');
            leading_bracket    = leading_bracket + formula_parts[i].substr(0, first_index + 1);
            formula_parts[i].erase(0, first_index + 1);
        }
        if (formula_parts[i][0] == '-') {
            leading_bracket = leading_bracket + "-";
            formula_parts[i].erase(0, 1);
        }
        while (formula_parts[i].back() == ')' || formula_parts[i].back() == ',') {
            trailing_bracket = trailing_bracket + formula_parts[i].back();
            formula_parts[i].pop_back();
        }

        if (model->getListOfParameters()->getElementBySId(formula_parts[i]) != NULL) {
            formula_parts[i] =
                "params.template get<mio::" + *name_namespace + "::" + formula_parts[i] + "<ScalarType>>()";
        }
        else if (model->getListOfSpecies()->getElementBySId(formula_parts[i]) != NULL) {
            return mio::failure(mio::StatusCode::InvalidValue);
        }
        else if (model->getListOfCompartments()->getElementBySId(formula_parts[i]) != NULL) {
            double size =
                Compartment_getSize((Compartment*)model->getListOfCompartments()->getElementBySId(formula_parts[i]));
            formula_parts[i] = std::to_string(size);
        }
        if (formula_parts[i] == "pi") {
            formula_parts[i] = "M_PI";
        }
        if (formula_parts[i] == "time") {
            return mio::failure(mio::StatusCode::InvalidValue);
        }
        if (leading_bracket.size() > 0) {
            formula_parts[i] = leading_bracket + formula_parts[i];
        }
        if (trailing_bracket.size() > 0) {
            formula_parts[i] = formula_parts[i] + trailing_bracket;
        }
    }
    return mio::success(boost::algorithm::join(formula_parts, " "));
}

/**
*    @brief Returns the maximal time mentioned in a formula or an error code
*
*    @param trigger The formula to be processed as ASTNode*
*    @param model The model where the formula stems from
*    @param name_namespace The namespace of the model
*
*    This function assumes a very specific layout of the formula. It needs to be a comparison between exaclty two 
*    values. If that is not the case, it prints an error. Then it expects one of the two nodes to be a time node. The 
*    other node has the comparison value and is returned. It is always assumed, but never checked, that we have a 
*    positive comparison, i.e. the time is indeed the maximum time.
*    This function is used to find the maximum time until which an event runs. It is used for the generation of the     
*    example.cpp file.
**/
mio::IOResult<std::string> find_tmax(const ASTNode* trigger, Model* model, std::string* name_namespace)
{
    auto trigger_type                      = trigger->getType();
    std::vector<ASTNodeType_t> comparisons = {AST_RELATIONAL_GEQ, AST_RELATIONAL_GT, AST_RELATIONAL_LEQ,
                                              AST_RELATIONAL_LT};
    bool is_comparison                     = false;
    for (auto i : comparisons) {
        if (i == trigger_type) {
            is_comparison = true;
        }
    }
    if (!is_comparison) {
        return mio::failure(mio::StatusCode::InvalidValue);
    }
    if (trigger->getNumChildren() > 2) {
        return mio::failure(mio::StatusCode::InvalidValue);
    }
    auto left_child  = trigger->getLeftChild();
    auto right_child = trigger->getRightChild();
    if (left_child->getType() != AST_NAME_TIME && right_child->getType() != AST_NAME_TIME) {
        return mio::failure(mio::StatusCode::InvalidValue);
    }
    if (left_child->getType() == AST_NAME_TIME) {
        return format_event_trigger(SBML_formulaToL3String(right_child), model, name_namespace);
    }
    else {
        return format_event_trigger(SBML_formulaToL3String(left_child), model, name_namespace);
    }
    mio::log_debug("No tmax as number was found in the formula.");
    return mio::failure(mio::StatusCode::InvalidValue);
}

/**
*    @brief Checks that a model is suitable, i.e. only uses one compartment
*
*    @param model The model to be checked
*
*    This function checks whether the model has exactly one compartment. If it does not, it returns false and prints an 
*    error. 
**/
bool verify_model_suitability(Model* model)
{
    if (model->getNumCompartments() != 1) {
        mio::log_error("Currently only models with exactly 1 compartment are supported!");
        return false;
    }
    return true;
}

/**
*    @brief Creates the infection_states.h file
*
*    @param model The model to be processed
*    @param filename The name of the file to be processed
*
*    Extracts the filename using :cpp:func:`get_filename(const char* filename)`, converts it to lower case and creates 
*    a file with the resulting name using `boost::filesystem`. Then it writes the species in the model as infection 
*    states to the file.
**/
bool create_infection_state(Model* model, const std::string& filename)
{
    std::string lowercase_name = boost::to_lower_copy<std::string>(filename);
    size_t number_species      = model->getListOfSpecies()->size();
    std::string uppercase_name = boost::to_upper_copy<std::string>(filename);

    std::ofstream infection_state;
    infection_state.open(lowercase_name + "/infection_state.h", std::ios::out);

    infection_state << "#ifndef " << uppercase_name << "_INFECTIONSTATE_H" << std::endl;
    infection_state << "#define " << uppercase_name << "_INFECTIONSTATE_H" << std::endl;
    infection_state << std::endl << std::endl;
    infection_state << "namespace mio" << std::endl;
    infection_state << "{" << std::endl;
    infection_state << "namespace " << lowercase_name << std::endl;
    infection_state << "{" << std::endl;
    infection_state << "enum class InfectionState" << std::endl;
    infection_state << "{" << std::endl;
    for (size_t i = 0; i < number_species; i++) {
        Species curr_species = *(Species*)model->getListOfSpecies()->get(i);
        infection_state << "  " << curr_species.getId() << "," << std::endl;
    }
    infection_state << "  Count\n};\n}\n}\n\n#endif" << std::endl;
    infection_state.close();

    return true;
}

/**
*    @brief Creates the parameters.h file
*
*    @param model The model to be processed
*    @param filename The name of the file to be processed
*
*    Extracts the filename using :cpp:func:`get_filename(const char* filename)`, converts it to lower case and creates 
*    a file with the resulting name using `boost::filesystem`. Then it creates one struct for every parameter in the 
*    model. It uses the value as returned by libsbml as default value. (This may be overwritten in the example.cpp file.)
**/
bool create_parameters(Model* model, const std::string& filename)
{
    std::string lowercase_name = boost::to_lower_copy<std::string>(filename);
    std::string uppercase_name = boost::to_upper_copy<std::string>(filename);

    size_t number_parameters = model->getListOfParameters()->size();

    std::ofstream parameters;
    parameters.open(lowercase_name + "/parameters.h", std::ios::out);

    parameters << "#ifndef " << uppercase_name << "_PARAMETERS_H" << std::endl;
    parameters << "#define " << uppercase_name << "_PARAMETERS_H" << std::endl;
    parameters
        << "\n#include \"memilio/epidemiology/age_group.h\"\n#include "
           "\"memilio/epidemiology/uncertain_matrix.h\"\n#include \"memilio/utils/custom_index_array.h\"\n#include "
           "\"memilio/utils/parameter_set.h\"\n#include \"memilio/utils/uncertain_value.h\"\n"
        << std::endl;
    parameters << "namespace mio\n{" << std::endl;
    parameters << "namespace " << lowercase_name << "\n{" << std::endl;

    std::string parameterset_initializer = "using ParametersBase = ParameterSet<";

    if (number_parameters == 0) {
        parameters << "template <typename FP = ScalarType>\nstruct DummyParam {\n    using Type = double;\n    static "
                      "Type get_default()\n    {\n        return 0.0;\n    }\n    static std::string name()\n    {\n   "
                      "     return \"DummyParameter\";\n    }\n};\n"
                   << std::endl;
        parameterset_initializer = parameterset_initializer + "DummyParam<FP>, ";
    }
    else {
        for (size_t i = 0; i < number_parameters; i++) {
            if (i != 0) {
                parameterset_initializer += ", ";
            }
            Parameter param  = *(Parameter*)model->getListOfParameters()->get(i);
            double value     = param.getValue();
            std::string name = param.getId();
            parameters << "template <typename FP = ScalarType>" << std::endl;
            parameters << "struct " << name << " {\n    using Type = double;\n  static Type get_default()\n{"
                       << std::endl;
            parameters << "return " << value << ";\n}" << std::endl;
            parameters << "static std::string name()\n{\n   return \"" << name << "\";\n}\n};\n" << std::endl;
            parameterset_initializer = parameterset_initializer + name + "<FP>";
        }
    }
    parameters << "template <typename FP = ScalarType>" << std::endl;
    parameters << parameterset_initializer << ">;\n " << std::endl;
    parameters << "template <typename FP = ScalarType>\nclass Parameters : public ParametersBase<FP>\n{\npublic:\n    "
                  "Parameters()\n        : ParametersBase<FP>()\n    {\n    }\n\npublic:\n    template <class "
                  "IOContext>\n    static IOResult<Parameters> deserialize(IOContext& io)\n    {\n        "
                  "BOOST_OUTCOME_TRY(auto&& base, ParametersBase<FP>::deserialize(io));\n        return "
                  "success(Parameters(std::move(base)));\n    }\n\n};\n}\n}\n#endif\n"
               << std::endl;

    parameters.close();

    return true;
}

/**
*    @brief Create the model.cpp file
*
*    @param filename The name of the file to be processed
*
*    Creates the generic model.cpp file. The file location is generated using the filename.
**/
bool create_model_cpp(const std::string& filename)
{

    std::string lowercase_name = boost::to_lower_copy<std::string>(filename);

    std::ofstream model_cpp;
    model_cpp.open(lowercase_name + "/model.cpp", std::ios::out);

    model_cpp << "#include \"" << lowercase_name << "/model.h\"\n\nnamespace mio\n{" << std::endl;
    model_cpp << "namespace " << lowercase_name << "\n{\n\n}\n}" << std::endl;
    model_cpp.close();

    return true;
}

/**
*    @brief Create the model.h file
*
*    @param model The model to be processed
*    @param filename The name of the file to be processed
*
*    Creates the model.h file based on the filename. First some generic input is written to the file. Then the 
*   `get_derivatives`-function is generated. It
*    - creates an index variable for every species
*    - creates a lambda function for every function definition in the model to have it at hand
*    - goes through every reaction in the model and generates the corresponding code. This is done in local scopes to make sure local parameters don't cause errors.
*    - checks every rule whether it is a RateRule and generates code for it if it is.
*    
*    In the end it appends a generic serialization function.
**/
bool create_model_h(Model* model, const std::string& filename)
{
    std::string lowercase_name = boost::to_lower_copy<std::string>(filename);
    std::string uppercase_name = boost::to_upper_copy<std::string>(filename);
    size_t number_species      = model->getListOfSpecies()->size();
    size_t number_reactions    = model->getListOfReactions()->size();
    std::ofstream model_h;
    model_h.open(lowercase_name + "/model.h", std::ios::out);

    model_h << "#ifndef " << uppercase_name << "_MODEL_H\n#define " << uppercase_name << "_MODEL_H" << std::endl;
    model_h << "\n#include <cmath>\n#include \"memilio/compartments/compartmentalmodel.h\"\n#include "
               "\"memilio/epidemiology/age_group.h\"\n#include \"memilio/epidemiology/populations.h\"\n#include "
               "\"memilio/epidemiology/contact_matrix.h\""
            << std::endl;
    model_h << "#include \"" << lowercase_name << "/infection_state.h\"\n#include \"" << lowercase_name
            << "/parameters.h\"" << std::endl;
    model_h << "namespace mio\n{" << std::endl;
    model_h << "namespace " << lowercase_name << std::endl;
    model_h << "{\ntemplate <typename FP = ScalarType>\nclass Model\n    : public mio::CompartmentalModel<FP, "
               "InfectionState, mio::Populations<FP, InfectionState>, Parameters<FP>>\n{\n    using Base =\n        "
               "mio::CompartmentalModel<FP, InfectionState, mio::Populations<FP, InfectionState>, "
               "Parameters<FP>>;\n\npublic:\n    using typename Base::ParameterSet;\n    using typename "
               "Base::Populations;\n\n    Model(const Populations& pop, const ParameterSet& params)\n        : "
               "Base(pop, params)\n    {\n    }\n\n    Model()\n        : Base(Populations(InfectionState::Count), "
               "ParameterSet())\n    {\n    }"
            << std::endl;
    model_h << "void get_derivatives(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, "
               "FP t, Eigen::Ref<Eigen::VectorX<FP>> dydt) const override\n{"
            << std::endl;
    model_h << "auto& params = this->parameters;" << std::endl;
    for (size_t i = 0; i < number_species; i++) {
        Species curr_species = *(Species*)model->getListOfSpecies()->get(i);
        model_h << "    size_t " << curr_species.getId()
                << "i = this->populations.get_flat_index({InfectionState::" << curr_species.getId() << "});"
                << std::endl;
    }
    for (size_t i = 0; i < model->getListOfFunctionDefinitions()->size(); i++) {
        FunctionDefinition function_def = *(FunctionDefinition*)model->getListOfFunctionDefinitions()->get(i);
        std::string lambda_string       = SBML_formulaToL3String(function_def.getMath());
        std::string variable_string     = "";
        std::string formula_string      = "";
        if (lambda_string.find_last_of(',') != std::string::npos) {
            size_t delim = lambda_string.find_last_of(',');
            std::vector<std::string> variables;
            boost::split(variables, lambda_string.substr(0, delim + 1), boost::is_any_of(", ()"));
            for (
                size_t variable_index = 1; variable_index < variables.size();
                variable_index++) { //We need to skip the leading "lambda" but can't filter for this as it might also be a variable.
                if (variables[variable_index] != "") {
                    variable_string = variable_string + " double " + variables[variable_index] + ",";
                }
            }
            variable_string.pop_back();
            formula_string = lambda_string.substr(delim + 1, lambda_string.size() - delim - 2);
        }
        model_h << "auto " << function_def.getId() << " = [](" << variable_string << "){" << std::endl;
        model_h << "return " << formula_string << ";" << std::endl;
        model_h << "};" << std::endl;
    }
    model_h << std::endl;
    for (size_t i = 0; i < number_reactions; i++) {
        model_h << "{" << std::endl;
        auto curr_reaction = *(Reaction*)model->getListOfReactions()->get(i);
        auto products      = curr_reaction.getListOfProducts();
        auto educts        = curr_reaction.getListOfReactants();
        auto modifiers     = curr_reaction.getListOfModifiers();
        auto formula       = curr_reaction.getKineticLaw();
        if (formula->getListOfLocalParameters()->size() != 0) { // for SBML Level 3
            for (size_t j = 0; j < formula->getListOfLocalParameters()->size(); j++) {
                auto param = *(Parameter*)formula->getListOfLocalParameters()->get(j);
                model_h << "double " << param.getId() << " = " << param.getValue() << ";" << std::endl;
            }
        }
        if (formula->getListOfParameters()->size() != 0) { // for outdated SBML
            for (size_t j = 0; j < formula->getListOfParameters()->size(); j++) {
                auto param = *(Parameter*)formula->getListOfParameters()->get(j);
                model_h << "double " << param.getId() << " = " << param.getValue() << ";" << std::endl;
            }
        }
        auto math_string = SBML_formulaToL3String(formula->getMath());
        std::vector<std::string> formula_parts;
        boost::split(formula_parts, math_string, boost::is_any_of(" "), boost::algorithm::token_compress_on);
        for (size_t j = 0; j < formula_parts.size(); j++) {
            std::string leading_bracket = "", trailing_bracket = "";
            while (formula_parts[j].find_first_of('(') != std::string::npos) {
                size_t first_index = formula_parts[j].find_first_of('(');
                leading_bracket    = leading_bracket + formula_parts[j].substr(0, first_index + 1);
                formula_parts[j].erase(0, first_index + 1);
            }
            if (formula_parts[j][0] == '-') {
                leading_bracket = leading_bracket + "-";
                formula_parts[j].erase(0, 1);
            }
            while (formula_parts[j].back() == ')' || formula_parts[j].back() == ',') {
                trailing_bracket = trailing_bracket + formula_parts[j].back();
                formula_parts[j].pop_back();
            }
            if (formula->getListOfLocalParameters()->getElementBySId(formula_parts[j]) != NULL) {
            }
            else if (model->getListOfParameters()->getElementBySId(formula_parts[j]) != NULL) {
                formula_parts[j] = "params.template get<" + formula_parts[j] + "<FP>>()";
            }
            else if (model->getListOfSpecies()->getElementBySId(formula_parts[j]) != NULL) {
                bool prod = false, ed = false, mod = false;
                for (size_t k = 0; k < products->size(); k++) {
                    if (products->get(k)->getSpecies() == formula_parts[k]) {
                        prod = true;
                    }
                }
                for (size_t k = 0; k < educts->size(); k++) {
                    if (educts->get(k)->getSpecies() == formula_parts[k]) {
                        ed = true;
                    }
                }
                for (size_t k = 0; k < modifiers->size(); k++) {
                    if (modifiers->get(k)->getSpecies() == formula_parts[k]) {
                        mod = true;
                    }
                }
                if (prod) {
                    formula_parts[j] = "pop[" + formula_parts[j] + "j]";
                }
                else if (ed) {
                    formula_parts[j] = "y[" + formula_parts[j] + "j]";
                }
                else if (mod) {
                    formula_parts[j] = "pop[" + formula_parts[j] + "j]";
                }
            }
            else if (model->getListOfCompartments()->getElementBySId(formula_parts[j]) != NULL) {
                double size = Compartment_getSize(
                    (Compartment*)model->getListOfCompartments()->getElementBySId(formula_parts[j]));
                formula_parts[j] = std::to_string(size);
            }
            if (formula_parts[j] == "pi") {
                formula_parts[j] = "M_PI";
            }
            if (formula_parts[j] == "time") {
                formula_parts[j] = "t";
            }
            if (leading_bracket.size() > 0) {
                formula_parts[j] = leading_bracket + formula_parts[j];
            }
            if (trailing_bracket.size() > 0) {
                formula_parts[j] = formula_parts[j] + trailing_bracket;
            }
        }

        std::string output_formula = boost::algorithm::join(formula_parts, " ");
        for (size_t j = 0; j < educts->size(); j++) {
            auto educt = *(SpeciesReference*)educts->get(j);
            model_h << "dydt[" << educt.getSpecies() << "j] -= " << output_formula << ";" << std::endl;
        }
        for (size_t j = 0; j < products->size(); j++) {
            auto product = *(SpeciesReference*)products->get(j);
            model_h << "dydt[" << product.getSpecies() << "i] += " << output_formula << ";" << std::endl;
        }
        model_h << "}" << std::endl;
    }
    for (size_t i = 0; i < model->getListOfRules()->size(); i++) {
        auto rule = model->getListOfRules()->get(i);
        if (rule->getType() == RULE_TYPE_RATE) {
            auto math_string = SBML_formulaToL3String(rule->getMath());
            std::vector<std::string> formula_parts;
            boost::split(formula_parts, math_string, boost::is_any_of(" "), boost::algorithm::token_compress_on);
            for (size_t j = 0; j < formula_parts.size(); j++) {
                std::string leading_bracket = "", trailing_bracket = "";
                while (formula_parts[j].find_first_of('(') != std::string::npos) {
                    size_t first_index = formula_parts[j].find_first_of('(');
                    leading_bracket    = leading_bracket + formula_parts[j].substr(0, first_index + 1);
                    formula_parts[j].erase(0, first_index + 1);
                }
                if (formula_parts[j][0] == '-') {
                    leading_bracket = leading_bracket + "-";
                    formula_parts[j].erase(0, 1);
                }
                while (formula_parts[j].back() == ')' || formula_parts[j].back() == ',') {
                    trailing_bracket = trailing_bracket + formula_parts[j].back();
                    formula_parts[j].pop_back();
                }
                if (model->getListOfParameters()->getElementBySId(formula_parts[j]) != NULL) {
                    formula_parts[j] = "params.template get<" + formula_parts[j] + "<FP>>()";
                }
                else if (model->getListOfSpecies()->getElementBySId(formula_parts[j]) != NULL) {
                    formula_parts[j] = "y[" + formula_parts[j] + "i]";
                }
                else if (model->getListOfCompartments()->getElementBySId(formula_parts[j]) != NULL) {
                    double size = Compartment_getSize(
                        (Compartment*)model->getListOfCompartments()->getElementBySId(formula_parts[j]));
                    formula_parts[j] = std::to_string(size);
                }
                if (formula_parts[j] == "pi") {
                    formula_parts[j] = "M_PI";
                }
                if (formula_parts[j] == "time") {
                    formula_parts[j] = "t";
                }
                if (leading_bracket.size() > 0) {
                    formula_parts[j] = leading_bracket + formula_parts[j];
                }
                if (trailing_bracket.size() > 0) {
                    formula_parts[j] = formula_parts[j] + trailing_bracket;
                }
            }
            std::string output_formula = boost::algorithm::join(formula_parts, " ");
            model_h << "dydt[" << rule->getVariable() << "i] += " << output_formula << ";" << std::endl;
        }
    }
    model_h << "}" << std::endl;
    model_h << "template <class IOContext>\n    void serialize(IOContext& io) const\n    {\n        auto obj = "
               "io.create_object(\"Model\");\n        obj.add_element(\"Parameters\", this->parameters);\n        "
               "obj.add_element(\"Populations\", this->populations);\n    }\n\n    template <class IOContext>\n    "
               "static IOResult<Model> deserialize(IOContext& io)\n    {\n        auto obj = "
               "io.expect_object(\"Model\");\n        auto par = obj.expect_element(\"Parameters\", "
               "Tag<ParameterSet>{});\n        auto pop = obj.expect_element(\"Populations\", Tag<Populations>{});\n   "
               "     return apply(\n            io,\n            [](auto&& par_, auto&& pop_) {\n                "
               "return Model{pop_, par_};\n            },\n            par, pop);\n    }\n};\n\n}\n}\n\n#endif"
            << std::endl;
    model_h.close();

    return true;
}

/**
*    @brief Create the CMakeLists.txt file for the model folder
*
*    @param filename The name of the file to be processed
*
*    Creates the generic CMakeLists.txt file for the model folder, assuming that only the files generated by this 
*    program are needed. The file location is generated using the filename and the file is stored in the generated 
*    folder.
**/
bool create_cmake(const std::string& filename)
{
    std::string lowercase_name = boost::to_lower_copy<std::string>(filename);

    std::ofstream cmakelists;
    cmakelists.open(lowercase_name + "/CMakeLists.txt", std::ios::out);
    cmakelists << "add_library(" << lowercase_name << std::endl;
    cmakelists << "infection_state.h\nmodel.h\nmodel.cpp\nparameters.h\n)" << std::endl;
    cmakelists << "target_link_libraries(" << lowercase_name << " PUBLIC memilio)" << std::endl;
    cmakelists << "target_include_directories(" << lowercase_name
               << " PUBLIC\n$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/"
                  "..>\n$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>\n)"
               << std::endl;
    cmakelists << "target_compile_options(" << lowercase_name << " PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})"
               << std::endl;
    cmakelists.close();

    return true;
}

/**
*    @brief Creates the corresponding file to example.cpp
*
*    @param model The model to be processed
*    @param filename The name of the file to be processed
*
*    The filename is based on the given filename. 
*    The file first contains generic input, then the model, it's parameters and specis are initialized. 
*    The model is simulated in steps to add events that are triggered by a specific time. In the end the model is 
*    simulated for another 50 days. The results are printed to the console and saved in a file as a table.
**/
bool create_example_cpp(Model* model, const std::string& filename)
{
    std::string lowercase_name = boost::to_lower_copy<std::string>(filename);
    size_t number_species      = model->getListOfSpecies()->size();

    std::ofstream example;
    example.open(lowercase_name + ".cpp", std::ios::out);
    example << "#include \"memilio/compartments/simulation.h\"\n#include \"memilio/config.h\"\n#include "
               "\"memilio/math/euler.h\"\n#include \"memilio/math/integrator.h\"\n#include "
               "\"memilio/utils/logging.h\"\n#include <memory>\n"
            << std::endl;
    example << "#include \"" << lowercase_name << "/infection_state.h\"" << std::endl;
    example << "#include \"" << lowercase_name << "/model.h\"" << std::endl;
    example << "#include \"" << lowercase_name << "/parameters.h\"" << std::endl;
    example << "\n int main()\n{\n  mio::set_log_level(mio::LogLevel::warn);" << std::endl;

    example << "ScalarType t0 = " << 0.0 << ";" << std::endl; //read out real start time
    example << "ScalarType tmax = " << 150.0 << ";" << std::endl; // read out real end time
    example << "ScalarType dt = " << 0.1 << ";" << std::endl; // is this a parameter?
    example << "mio::log_info(\"Simulating model " << lowercase_name << "; t={} ... {} with dt = {}.\", t0, tmax, dt);"
            << std::endl;
    example << "mio::" << lowercase_name << "::Model<ScalarType> model{};" << std::endl;
    example << "auto& params = model.parameters;" << std::endl;
    std::string table = "";
    for (size_t i = 0; i < number_species; i++) {
        Species curr_species = *(Species*)model->getListOfSpecies()->get(i);
        std::string number   = get_initial_assignment(&curr_species, model, &lowercase_name);
        example << "model.populations[{mio::" << lowercase_name << "::InfectionState::" << curr_species.getId()
                << "}] = " << number << ";" << std::endl;
        table = table + "\"" + curr_species.getId() + "\",";
    }
    table.pop_back();

    //Add initial assignments for Parameters...
    for (size_t i = 0; i < model->getNumInitialAssignments(); i++) {
        auto assignment = model->getListOfInitialAssignments()->get(i);
        for (size_t j = 0; j < model->getNumParameters(); j++) {
            if (model->getParameter(j)->getIdAttribute() == assignment->getSymbol())
                example << "params.get<mio::" << lowercase_name << "::" << assignment->getSymbol()
                        << "<ScalarType>>()  = "
                        << format_main_formula(model, SBML_formulaToL3String(assignment->getMath()), &lowercase_name)
                        << ";" << std::endl;
        }
    }

    example << "std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator = "
               "std::make_shared<mio::EulerIntegratorCore<ScalarType>>();"
            << std::endl;
    example << "auto sim = mio::Simulation<ScalarType, mio::" << lowercase_name
            << "::Model<ScalarType>>(model, t0, dt);" << std::endl;
    example << "sim.set_integrator(integrator);" << std::endl;
    for (size_t i = 0; i < model->getListOfEvents()->size(); i++) {
        auto event                      = model->getListOfEvents()->get(i);
        auto trigger                    = event->getTrigger()->getMath();
        mio::IOResult<std::string> tmax = find_tmax(trigger, model, &lowercase_name);
        if (!tmax) {
            mio::log_info("There is an event that I can not implement, I will ignore it.");
        }
        else {
            example << "tmax = " << tmax.value() << ";" << std::endl;
            example << "sim.advance(tmax);" << std::endl;
            for (size_t j = 0; j < event->getListOfEventAssignments()->size(); j++) {
                auto assignment                               = event->getListOfEventAssignments()->get(j);
                auto variable                                 = assignment->getVariable();
                mio::IOResult<std::string> formatted_variable = format_event_variable(variable, model, &lowercase_name);
                if (!formatted_variable) {
                    return false;
                }
                std::string formula                          = SBML_formulaToL3String(assignment->getMath());
                mio::IOResult<std::string> formatted_formula = format_event_formulas(formula, model, &lowercase_name);
                example << formatted_variable.value() << " = " << formatted_formula.value() << ";" << std::endl;
            }
        }
    }
    example << "sim.advance(tmax + 50); \nauto sir = sim.get_result();" << std::endl;

    example << "sir.print_table({" << table << "});" << std::endl;
    example << "std::cout << \"\\nnumber total: \" << sir.get_last_value().sum() << \"\\n\";" << std::endl;

    example << "const std::string file_name = \"" << lowercase_name << ".dat\";" << std::endl;
    example << "std::ofstream file(file_name);\nstd::cout << \"Writing output to \" << file_name << std::endl;"
            << std::endl;
    example << "sir.print_table({" << table << "}, 10, 3, file);" << std::endl;
    example << "file.close();" << std::endl;

    example << "}";
    example.close();

    return true;
}

/**
*    @brief Create the additional lines for the CmakeLists.txt file
*
*    @param filename The name of the file to be processed
*
*    Extracts the filename using :cpp:func:`get_filename(const char* filename)`. Then it writes the generic compilation 
*    commands to a file called "CMakeListsAddition.txt" that is stored in the current working directory.
**/
bool modify_cmakelists(const std::string& filename)
{
    std::string lowercase_name = boost::to_lower_copy<std::string>(filename);
    std::ofstream modifications;
    modifications.open("CMakeListsAddition.txt", std::ios::app);
    modifications << "add_executable(ex_" << lowercase_name << " " << lowercase_name << ".cpp)" << std::endl;
    modifications << "target_link_libraries(ex_" << lowercase_name << " PRIVATE memilio " << lowercase_name << ")"
                  << std::endl;
    modifications << "target_compile_options(ex_" << lowercase_name
                  << " PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})" << std::endl;
    modifications.close();
    modifications.open("CMakeListsFolderNames.txt", std::ios::app);
    modifications << "add_subdirectory(models/" << lowercase_name << ")" << std::endl;
    modifications.close();
    return true;
}

/** 
*    @brief Calls clang-format to format the generated files
*
*    @param filename The name of the file to be processed
*
*    Calls clang-format using the `system` command. The command is generated using the filename. Dpending on the return 
*    code, an error or info message is printed.
**/
void format_files(const std::string& filename)
{
    std::string lowercase_name = boost::to_lower_copy<std::string>(filename);
    std::string command =
        std::string("clang-format -i --files= ") + lowercase_name + ".cpp " + lowercase_name + "/**.[hc]*";
    int status = system(command.c_str());
    mio::log_debug("Return status: {}", status);
    if (status != 0) {
        mio::log_error("Error while trying to format the files, files are functional but hard to read.");
    }
    else {
        mio::log_info("Formatted all files.");
    }
}

/** 
*    @brief Runs the program
*    @param argc The number of command line arguments
*    @param argv The command line arguments, here the name of the file that should be processed
*
*    This function expects exaclty one parameter which should be the path to the file that should be processed. If more 
*    or less parameters are given, it raises an error. If the correct number of parameters is given, it treats the file 
*    as a SBML file and generates the MEmilio code for it.
**/
int main(int argc, char* argv[])
{
    mio::set_log_level(mio::LogLevel::debug);

    if (argc != 2) {
        mio::log_error("Please provide a SBML file at startup!");
        return 1;
    }
    const char* filename = argv[1];
    SBMLDocument* document;
    SBMLReader reader;

    document = reader.readSBML(filename);

    if (SBMLDocument_getNumErrors(document) > 0) {
        if (XMLError_isFatal(SBMLDocument_getError(document, 0)) ||
            XMLError_isError(SBMLDocument_getError(document, 0))) {
            mio::log_error("Fatal error while trying to read the input file!");
            return 3;
        }
    }

    auto model = document->getModel();

    if (!verify_model_suitability(model)) {
        return 4;
    }

    std::string core_filename = get_filename(filename);

    create_folder(core_filename);

    if (!create_infection_state(model, core_filename)) {
        return 5;
    }

    if (!create_parameters(model, core_filename)) {
        return 6;
    }

    if (!create_model_cpp(core_filename)) {
        return 7;
    }

    if (!create_model_h(model, core_filename)) {
        return 8;
    }

    if (!create_cmake(core_filename)) {
        return 9;
    }

    if (!create_example_cpp(model, core_filename)) {
        return 10;
    }

    if (!modify_cmakelists(core_filename)) {
        return 11;
    }
    mio::log_info("Created all files.");
    mio::log_info("Formatting files.");
    format_files(core_filename);
}