#include "sbml2memilio.h"

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cout << "Please provide a SBML file at startup!" << std::endl << std::endl;
        return 1;
    }
    const char* filename = argv[1];
    SBMLDocument* document;
    SBMLReader reader;

    document = reader.readSBML(filename);

    auto model = document->getModel();

    if (!verify_model_suitability(model)) {
        return 1;
    }

    create_folder(filename);

    if (!create_infection_state(model, filename)) {
        return 1;
    }

    if (!create_parameters(model, filename)) {
        return 1;
    }

    if (!create_model_cpp(model, filename)) {
        return 1;
    }

    if (!create_model_h(model, filename)) {
        return 1;
    }

    if (!create_cmake(model, filename)) {
        return 1;
    }

    if (!create_example_cpp(model, filename)) {
        return 1;
    }

    if (!modify_cmakelists(model, filename)) {
        return 1;
    }

    format_files(filename);

    std::cout << "Das war ein voller Erfolg!" << std::endl;
    return 0;
}