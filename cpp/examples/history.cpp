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
    // Create a History object with the predefined Writer DataWriterToMemory and our above defined Logger MyLogger.
    // The Writer DataWriterToMemory uses an std::vector to store the data.
    mio::History<mio::DataWriterToMemory, MyLogger> history;

    // Log some data.
    for (int i = 0; i < 10; ++i) {
        history.log(i); // The data specified by the Logger is stored in the History object using the Writer.
    }

    // Get the log data.
    // For this one has to specify which from which Logger the data should be retrieved.
    // Also the data type of this is an std::vector<int> as the Writer has a tuple of std::vector<x> where x is specified by the Type of the logger, in this case int.
    std::vector<int> logData = std::get<0>(history.get_log()); // The first logger is MyLogger

    // Print the logged data.
    std::cout << "Logged data from MyLogger:" << std::endl;
    for (const auto& data : logData) {
        std::cout << data << std::endl;
    }

    return 0;
}