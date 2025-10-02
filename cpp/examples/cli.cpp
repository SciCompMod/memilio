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

#include <vector>

struct Name {
    using Type = std::vector<std::string>;

    static Type get_default()
    {
        return Type{"FirstName", "LastName"};
    }
    const static std::string name()
    {
        return "Name";
    }
    const static std::string alias()
    {
        return "n";
    }
    const static std::string description()
    {
        return "Enter your name as list of strings.";
    }
};

struct Age {
    using Type = int;
    const static std::string name()
    {
        return "Age";
    }
    const static std::string alias()
    {
        return "a";
    }
    const static std::string description()
    {
        return "Enter your age.";
    }
};

struct Greeting {
    using Type = std::string;

    static Type get_default()
    {
        return Type{"Hello World!"};
    }
    const static std::string name()
    {
        return "Greeting";
    }
    const static std::string description()
    {
        return "Enter a custom greeting.";
    }
};

int main(int argc, char** argv)
{
    if (argc == 1) { // Print this if no arguments were given
        std::cout << "This is a small example on how to use the command line interface. "
                     "Use \"-h\" to show the help dialogue.\n";
    }
    // create parameter set
    auto parameters = mio::ParameterSet<Name, Age, Greeting>{};
    // get command line options
    auto result = mio::command_line_interface("cli_example", argc, argv, parameters);
    // catch errors
    if (!result) {
        std::cout << result.error().formatted_message();
        return result.error().code().value();
    }
    // do something with the parameters
    std::cout << parameters.get<Greeting>() << "\n"
              << "Name: ";
    for (auto& name : parameters.get<Name>()) {
        std::cout << name << " ";
    }
    std::cout << "\n";
    if (parameters.get<Age>() > 0) {
        std::cout << "Age: " << parameters.get<Age>() << "\n";
    }
}
