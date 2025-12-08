#include "memilio/io/cli.h"

#include <iostream>
#include <fstream>
#include <string>


void run(const std::string& output_file_path) {
    std::ofstream out(output_file_path);
    if (!out) {
        std::cerr << "Failed to open file.\n";
        return;
    }

    out <<
        "Time, SchoolClosure, HomeOffice, PhysicalDistancingSchool, PhysicalDistancingWork, PhysicalDistancingOther\n"
        "0.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000\n"
        "7.0000000000, 1.0000000000, 0.0000663063, 0.0000000000, 0.0000663063, 1.0000000000\n"
        "14.0000000000, 1.0000000000, 0.0017841277, 0.0000000000, 0.0017841277, 1.0000000000\n"
        "21.0000000000, 1.0000000000, 0.0018360008, 0.0000000000, 0.0018360008, 1.0000000000\n"
        "28.0000000000, 1.0000000000, 0.0007428344, 0.0000000000, 0.0007428344, 1.0000000000\n"
        "35.0000000000, 1.0000000000, 0.0000460844, 0.0000000000, 0.0000460844, 1.0000000000\n"
        "42.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000\n"
        "49.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000\n";

    out.close();
}
int main(int argc, char** argv)
{
    auto cli_parameters = mio::cli::ParameterSetBuilder()
                            .add<"DataDirectory">(std::string{})
                            .add<"OutputFileName">(std::string{"OptimalControl_output.txt"})
                            .add<"ConstraintInfectedCases">(1e5)
                            .build();

    auto cli_result = mio::command_line_interface(argv[0], argc, argv, cli_parameters);
    if (!cli_result) {
        std::cout << cli_result.error().message();  
        return cli_result.error().code().value();  
    }

    auto output_file_path = cli_parameters.get<"DataDirectory">() + "/" + cli_parameters.get<"OutputFileName">();
    run(output_file_path);
    return 0;
}