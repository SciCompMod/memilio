#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <cstddef>

#include "boost/outcome/try.hpp"
#include "boost/outcome/result.hpp"
#include "boost/optional.hpp"
#include "boost/filesystem.hpp"

#include "memilio/io/epi_data.h"
#include "memilio/utils/logging.h"


#include "optimization_model/optimization_model.h"


// Example: ./../build/bin/secirvvs_ESID /home/jli/Memilio-Branches/memilio/data/Germany

mio::IOResult<void> read_json_files(const boost::filesystem::path& data_directory)
{
    BOOST_OUTCOME_TRY(auto&& intervention_list, mio::read_json((data_directory / "data" / "intervention_list.json").string()));
    BOOST_OUTCOME_TRY(auto&& parameter_list, mio::read_json((data_directory / "data" / "parameter_list.json").string()));
    BOOST_OUTCOME_TRY(auto&& scenarios, mio::read_json((data_directory / "data" / "scenarios.json").string()));
    BOOST_OUTCOME_TRY(auto&& scenario_data_run, mio::read_json((data_directory / "data" / "scenario_data_run.json").string()));
    mio::unused(intervention_list, parameter_list, scenarios, scenario_data_run);
    return mio::success();
}

int main(int argc, char* argv[])
{

    // ------------------------------ //
    // --- Parse 'data_directory' --- //
    // ------------------------------ //
    boost::filesystem::path data_directory;
    if (argc > 2) {
        std::cerr << "Usage: " << argv[0] << " [data_directory]\n";
        return 1;
    }
    if (argc == 2) {
        data_directory = argv[1];
    }
    else {
        data_directory = boost::filesystem::current_path();
    }

    mio::set_log_level(mio::LogLevel::warn);

    // double INF = 1e19;

    // ---------------------------- //
    // --- Create model problem --- //
    // ---------------------------- //
    double t0   = 0.0;
    double tmax = 56.0; // 8 weeks in days
    
    auto result = read_json_files(data_directory);
    OptimizationModel model(data_directory, t0, tmax, 6);
    mio::osecirvvs::Model<double> compmodel = model.create_model<double>();

    std::cout << compmodel.parameters.template get<mio::osecirvvs::ReducTimeInfectedMild<double>>()[mio::AgeGroup(0)] << "\n";
    return 0;
}
