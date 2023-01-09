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
    // create parameter set
    auto parameter = mio::ParameterSet<Name, Age, Greeting>{};
    // get command line options
    auto result = mio::command_line_interface("cli", argc, argv, parameter);
    // catch errors
    if (!result) {
        std::cout << result.error().formatted_message();
        return result.error().code().value();
    }
    // do something with the parameters
    std::cout << parameter.get<Greeting>() << "\n"
              << "Name: ";
    for (auto& name : parameter.get<Name>()) {
        std::cout << name << " ";
    }
    std::cout << "\n";
    if (parameter.get<Age>() > 0) {
        std::cout << "Age: " << parameter.get<Age>() << "\n";
    }
}
