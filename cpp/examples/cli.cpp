/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding
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

int main(int argc, char** argv)
{
    // Print a message if no arguments were given.
    if (argc == 1) {
        std::cout << "This is a small example on how to use the command line interface. "
                     "Use \"-h\" to show the help dialogue.\n";
    }
    // Create a parameter set for the CLI using the builder. This defines the (parameter) options that the user can set
    // through the command line.
    //
    // To add() a parameter, you need to specify the name and type as template arguments, then pass an initial value.
    // The type can sometimes be deduced from the initial value, so it can potentially be omitted.
    // After the initial value you can set some optional fields of the parameter:
    // - alias, which allows you to add a shorthand for setting values
    // - description, which contains details on, e.g., what the parameter does and what values are accepted.
    // - is_required, which makes the CLI check whether the parameter was set. If not, it exits with an error.
    //
    // As a general rule, use simple types! The more complicated the type, the more complex is the Json representation
    // that the user has to input.
    //
    // Instead of using the builder, you can also define and pass a mio::ParameterSet as parameters.
    // The main difference (for the CLI) is that the mio::ParameterSet uses struct names to "get" parameters, while
    // the mio::cli::ParameterSet uses StringLiteral%s.
    auto parameters = mio::cli::ParameterSetBuilder()
                          .add<"Name", std::vector<std::string>>({"FirstName", "LastName"},
                                                                 {"n", "Enter your name as list of strings.", false})
                          .add<"Age">(0, {"a", "Enter your age."})
                          .add<"Greeting">(std::string("Hello World!"),
                                           {.description = "Enter a custom greeting.", .is_required = false})
                          .build();
    // Define some default options. This is an optional feature, that allows users to set some options in the given
    // order as the first arguments, without specifying their name or alias.
    auto default_options = std::vector<std::string>{"Name", "Age"};
    // Parse command line arguments and/or set parameters. This next line as well as the following check on its result
    // are required to use the CLI.
    auto result = mio::command_line_interface(argv[0], argc, argv, parameters, default_options);
    // Catch and print help output, printed options, and errors.
    if (!result) {
        std::cout << result.error().message(); // Do not use formatted_message().
        return result.error().code().value(); // Use exit here when not used in main().
    }
    // Now, do something with the parameters!
    // Note that the CLI only verifies that the user input is parsable, not plausible. If a parameter has certain value
    // requirements, like "Age > 0", you must check this yourself.
    std::cout << parameters.get<"Greeting">() << "\n"
              << "Name: ";
    for (auto& name : parameters.get<"Name">()) {
        std::cout << name << " ";
    }
    std::cout << "\n";
    if (parameters.get<"Age">() > 0) {
        std::cout << "Age: " << parameters.get<"Age">() << "\n";
    }

    return 0;
}
