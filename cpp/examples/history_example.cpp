#include "history.h"
#include <iostream>
#include <vector>

// Define a Logger
struct MyLogger {
    using Type = int;

    Type log(const int& input)
    {
        return input * 2;
    }

    bool should_log(const int& input)
    {
        return input % 2 == 0;
    }
};

// Define a Writer
struct MyWriter {
    using Data = std::vector<int>;

    template <class Logger>
    static void add_record(const typename Logger::Type& t, Data& data)
    {
        data.push_back(t);
    }
};

int main()
{
    // Create a History object
    History<MyWriter, MyLogger> history;

    // Log some data
    for (int i = 0; i < 10; ++i) {
        history.log(i);
    }

    // Get the log data
    auto logData = history.get_log<MyLogger>();

    // Print the log data
    for (const auto& data : logData) {
        std::cout << data << std::endl;
    }

    return 0;
}