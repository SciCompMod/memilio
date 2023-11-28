#include "memilio/io/history.h"
#include <iostream>
#include <vector>

// Define a Logger
struct MyLogger : mio::LogAlways {
    using Type = int;
    static Type log(const int& input)
    {
        return input * 2;
    }
};

int main()
{
    // Create a History object
    mio::History<mio::DataWriterToMemory, MyLogger> history;

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