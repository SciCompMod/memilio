#include <epidemiology/seir.h>
#include <parameter_studies/parameter_studies.h>

int main (int argc, char *argv[])
{
    std::string input_filename;

    if (argc > 1) {
        // If provided, the first argument is the input file
        input_filename = argv[1];
    }
    else {
        // If not provided, we use a sample input file
        input_filename = "parameter_studies_example_input.txt";
    }

    // Create parameter study
    parameter_study_t<simulate_seir> parameter_study(input_filename);

    // Run parameter study
    parameter_study.run();

    // Gather output
    parameter_study.gather_output();

    return 0;
}