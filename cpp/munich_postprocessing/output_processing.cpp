#include <algorithm>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cstddef>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <fstream>
#include <chrono>
#include <map>
#include "memilio/io/hdf5_cpp.h"
#include "models/abm/location_type.h"
#include "memilio/utils/logging.h"
#include "boost/filesystem.hpp"
#include "output_processing.h"

int main()
{
    // int num_sims                    = 1;
    // int start_sim                   = 1;
    // std::string sim_result_folder   = "/hpc_data/bick_ju/INSIDeMunichPaper/IScienceRebuttal/OnlyMovement/";
    // std::string output_file_loctype = sim_result_folder + std::to_string(start_sim) + "_output_v5.h5";
    // std::string output_file_ww_area = sim_result_folder + std::to_string(start_sim) + "_output_v4.h5";
    // std::string save_file_loctype   = sim_result_folder + std::to_string(start_sim) + "_num_agents_loc_type.txt";
    // std::string save_file_area      = sim_result_folder + std::to_string(start_sim) + "_num_agents_area.txt";
    // std::string save_file_area2 =
    //     sim_result_folder + std::to_string(start_sim) + "_num_agents_area_inhabitants_commuters.txt";
    // calculate_agents_per_area_inhabitants_commuters(output_file_ww_area, save_file_area2);
    // calculate_agents_per_quantity(sim_result_folder, save_file_loctype, num_sims, start_sim, "locType", "v5");
    // calculate_agents_per_quantity(sim_result_folder, save_file_area, num_sims, start_sim, "area", "v4");
    // std::string output_dir_infections           = "/hpc_data/bick_ju/1000000/";
    // std::string save_file_loctype_infections    = output_dir_infections + "num_agents_infections_loctype.txt";
    // std::string save_file_ww_area_infections    = output_dir_infections + "num_agents_infections_area.txt";
    // std::string save_file_hh_size_ag_infections = output_dir_infections + "num_agents_infections_hh_size_ag.txt";
    // std::string person_file = "/hpc_data/bick_ju/INSIDeMunichPaper/IScienceRebuttal/Input/persons_scaled.csv";
    // //calculate_infections_per_hh_size(output_dir_infections, person_file, save_file_hh_size_ag_infections, num_sims);
    // // calculate_infections_per_quantity(output_dir_infections, save_file_loctype_infections, "locType", "v5", num_sims);
    // // calculate_infections_per_quantity(output_dir_infections, save_file_ww_area_infections, "area", "v4", num_sims);
    // std::string save_file_loctype_ag = sim_result_folder + std::to_string(start_sim) + "_num_agents_loc_type_ag.txt";
    // std::string save_file_area_ag    = sim_result_folder + std::to_string(start_sim) + "_num_agents_area_ag.txt";
    // calculate_agents_per_quantity_age_groups(sim_result_folder, save_file_loctype_ag, "locType", person_file, num_sims,
    //                                          start_sim, "v5");
    // calculate_agents_per_quantity_age_groups(sim_result_folder, save_file_area_ag, "area", person_file, num_sims,
    //                                          start_sim, "v4");

    std::string area_to_station_file =
        "/hpc_data/bick_ju/INSIDeMunichPaper/IScienceRebuttal/Input/tandler_upstream_gebiete.json";
    std::string infection_paths_file =
        "/hpc_data/bick_ju/INSIDeMunichPaper/IScienceRebuttal/example_output/1_infection_paths.txt";
    std::string hdf5_output_file = "/hpc_data/bick_ju/INSIDeMunichPaper/IScienceRebuttal/example_output/1_output_v4.h5";
    std::string save_file        = "/hpc_data/bick_ju/INSIDeMunichPaper/IScienceRebuttal/example_output/";
    int timepoint                = 15 * 24;
    std::string station_name     = "S_M4";
    save_file += "infection_age_shedding_" + std::to_string(timepoint) + "_" + station_name + ".csv";

    calculate_infected_per_measurement_station(area_to_station_file, infection_paths_file, hdf5_output_file, save_file,
                                               timepoint, station_name);

    return 0;
}
