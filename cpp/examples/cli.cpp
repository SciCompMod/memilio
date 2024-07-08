#include "memilio/io/cli.h"
#include "memilio/io/io.h"
#include "memilio/io/json_serializer.h"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <system_error>
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

using HourlyContactMatrix = std::array<Eigen::MatrixXd, 24>;

mio::IOResult<HourlyContactMatrix> read_hourly_contact_matrix(const std::string& csv_file)
{
    std::ifstream file(csv_file);
    if (!file.good()) {
        return mio::failure(mio::StatusCode::FileNotFound, "Could not open " + csv_file + ".");
    }

    std::string reader;

    std::vector<double> matrix_entries;
    int current_hour = 0;
    int t, row, col;
    double val;

    HourlyContactMatrix hcm;

    // skip first line
    std::getline(file, reader);
    // possible EOF here
    // read csv
    while (std::getline(file, reader)) {
        int status = sscanf(reader.c_str(), "%i,%i,%i,%lf\n", &t, &row, &col, &val);

        if (status != 4) {
            return mio::failure(mio::StatusCode::InvalidFileFormat,
                                "Unexpected format while reading " + csv_file + ". Line reads \"" + reader + "\"");
        }

        printf("%i,%i,%i,%lf\n", t, row, col, val);

        if (t > current_hour) {
            size_t n          = std::round(std::sqrt(matrix_entries.size()));
            hcm[current_hour] = Eigen::MatrixXd(n, n);
            for (size_t i = 0; i < matrix_entries.size(); i++) {
                hcm[current_hour].data()[i] = matrix_entries[i];
            }
            matrix_entries.clear();
            ++current_hour;
        }

        matrix_entries.push_back(val);
    }

    size_t n          = std::round(std::sqrt(matrix_entries.size()));
    hcm[current_hour] = Eigen::MatrixXd(n, n);
    for (size_t i = 0; i < matrix_entries.size(); i++) {
        hcm[current_hour].data()[i] = matrix_entries[i];
    }

    return mio::success(hcm);
}

int main(int argc, char** argv)
{
    std::string filename = "/home/schm_r6/Documents/24h_networks_csv/office_5_5.csv";

    auto res = read_hourly_contact_matrix(filename);

    if (!res) {
        std::cout << res.error().formatted_message();
        return res.error().code().value();
    }

    std::cout << mio::serialize_json(res.value()).value() << "\n";

    return 0;

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
