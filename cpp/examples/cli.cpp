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
    auto p = mio::ParameterSet<Name, Age, Greeting>{};
    // get command line options
    auto r = mio::command_line_interface("cli", argc, argv, p);
    // catch errors
    if (!r) {
        std::cout << r.error().formatted_message();
        return r.error().code().value();
    }
    // do something with the parameters
    std::cout << p.get<Greeting>() << "\n"
              << "Name: ";
    for (auto& name : p.get<Name>()) {
        std::cout << name << " ";
    }
    std::cout << "\n";
    if (p.get<Age>() > 0) {
        std::cout << "Age: " << p.get<Age>() << "\n";
    }
}