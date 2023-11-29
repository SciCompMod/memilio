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
    // Create a History object with an Writer and a Logger.
    // The Writer is used to store the data in the History object and the Logger is used to specify which data is stored, as described in the Logger section.
    mio::History<mio::DataWriterToMemory, MyLogger> history;

    // Log some data
    for (int i = 0; i < 10; ++i) {
        history.log(i); // The data specified by the logger is stored in the History object using the Writer.
    }

    // Get the log data. One has to specify which from which logger the data should be retrieved.
    std::vector<int> logData = std::get<0>(history.get_log()); // The first logger is MyLogger

    // Print the log data
    for (const auto& data : logData) {
        std::cout << data << std::endl;
    }

    return 0;
}