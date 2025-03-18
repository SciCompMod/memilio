#ifndef OSIR_DOC_H
#define OSIR_DOC_H
#include <unordered_map>
#include <string>

namespace pymio
{

namespace details
{

    std::unordered_map<std::string, std::string> docstrings = {
        {"interpolate_simulation_result", ""},
        {"interpolate_ensemble_results", ""}
    };

    // void initialize_docs() {
    //     docs["my_function"] = "This function does something useful.\n\nArgs:\n    x (int): An input value\n\nReturns:\n    int: The result.";
    //     docs["MyClass"] = "This is a class that does things.\n\nMethods:\n    my_method(): Performs an action.";
    // }

} // namespace details
} // namespace pymio

#endif //OSIR_DOC_H