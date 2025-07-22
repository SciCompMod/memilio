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

void calculate_agents_per_area_inhabitants_commuters(std::string output_file, std::string save_file)
{
    auto start = std::chrono::high_resolution_clock::now();
    //Load H5 file
    mio::H5File sim_output{H5Fopen(output_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)};
    std::string group_name = "data";
    //Open group
    mio::H5Group h5group{H5Gopen(sim_output.id, group_name.c_str(), H5P_DEFAULT)};
    //Open dataset
    mio::H5DataSet dset_loc_ids{H5Dopen(h5group.id, "loc_ids", H5P_DEFAULT)};
    //Open dataspace
    mio::H5DataSpace dataspace{H5Dget_space(dset_loc_ids.id)};
    //Get dimensions
    int ndims     = H5Sget_simple_extent_ndims(dataspace.id);
    hsize_t* dims = (hsize_t*)malloc(ndims * sizeof(hsize_t));
    H5Sget_simple_extent_dims(dataspace.id, dims, NULL);

    hsize_t rows = dims[0];
    hsize_t cols = dims[1];
    std::cout << "Num agents: " << rows << std::endl;
    int* data = (int*)malloc(rows * cols * sizeof(int));
    std::cout << "Starting reading dataset...\n";
    herr_t status = H5Dread(dset_loc_ids.id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    std::cout << "Finished reading dataset!\n";
    std::vector<int> data1(data, data + dims[0] * dims[1]);
    std::cout << "Starting getting max area...\n";
    int max_area = *std::max_element(data1.begin(), data1.end());
    std::cout << "Got max area!\n";
    //Save data in matrix
    // Matrix has tuple with (NumAgents, NumInhabitants, NumCommuters)
    std::vector<std::vector<std::tuple<int, int, int>>> matrix(
        dims[1], std::vector<std::tuple<int, int, int>>(max_area + 1, {0, 0, 0}));
    for (size_t i = 0; i < dims[0]; ++i) {
        for (size_t j = 0; j < dims[1]; ++j) {
            int area = data[i * dims[1] + j];
            // First tuple element is the number of agents
            std::get<0>(matrix[j][area]) += 1;
            if (area == data[i * dims[1]]) {
                std::get<1>(matrix[j][area]) += 1;
            }
            else {
                std::get<2>(matrix[j][area]) += 1;
            }
        }
    }
    std::cout << "Created matrix!\n";

    auto file = fopen(save_file.c_str(), "w");
    if (file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", save_file);
    }
    else {
        fprintf(file, "t,Area,NumAgents, NumInhabitants, NumCommuters");
        fprintf(file, "\n");
        for (size_t i = 0; i < matrix.size(); ++i) {
            for (size_t j = 0; j < matrix[0].size(); ++j) {
                fprintf(file, "%d,", static_cast<int>(i));
                fprintf(file, "%d,", static_cast<int>(j));
                fprintf(file, "%d,", std::get<0>(matrix[i][j]));
                fprintf(file, "%d,", std::get<1>(matrix[i][j]));
                fprintf(file, "%d,", std::get<2>(matrix[i][j]));
                fprintf(file, "\n");
            }
        }
        fclose(file);
    }
    auto stop     = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time to calculate agents per quantity: " << duration.count() << "[s]" << std::endl;
}

int main()
{
    int sim_num                     = 0;
    std::string output_file_loctype = "/hpc_data/bick_ju/OnlyMovement/" + std::to_string(sim_num) + "_output_v5.h5";
    std::string output_file_ww_area = "/hpc_data/bick_ju/OnlyMovement/" + std::to_string(sim_num) + "_output_v4.h5";
    std::string save_file_loctype =
        "/hpc_data/bick_ju/OnlyMovement/" + std::to_string(sim_num) + "_num_agents_loc_type.txt";
    std::string save_file_area = "/hpc_data/bick_ju/OnlyMovement/" + std::to_string(sim_num) + "_num_agents_area.txt";
    std::string save_file_area2 =
        "/hpc_data/bick_ju/OnlyMovement/" + std::to_string(sim_num) + "_num_agents_area_inhabitants_commuters.txt";
    //calculate_agents_per_area_inhabitants_commuters(output_file_ww_area, save_file_area2);
    // calculate_agents_per_quantity(output_file_loctype, save_file_loctype, "locType");
    // calculate_agents_per_quantity(output_file_ww_area, save_file_area, "area");
    std::string output_dir_infections           = "/hpc_data/bick_ju/1000000/";
    std::string save_file_loctype_infections    = output_dir_infections + "num_agents_infections_loctype.txt";
    std::string save_file_ww_area_infections    = output_dir_infections + "num_agents_infections_area.txt";
    std::string save_file_hh_size_ag_infections = output_dir_infections + "num_agents_infections_hh_size_ag.txt";
    std::string person_file                     = "/home/bick_ju/Documents/INSIDe/persons/persons_final_corr.csv";
    int num_sims                                = 10;
    //calculate_infections_per_hh_size(output_dir_infections, person_file, save_file_hh_size_ag_infections, num_sims);
    // calculate_infections_per_quantity(output_dir_infections, save_file_loctype_infections, "locType", "v5", num_sims);
    // calculate_infections_per_quantity(output_dir_infections, save_file_ww_area_infections, "area", "v4", num_sims);
    std::string save_file_loctype_ag =
        "/hpc_data/bick_ju/OnlyMovement/new_060725/" + std::to_string(sim_num) + "_num_agents_loc_type_ag.txt";
    std::string save_file_area_ag =
        "/hpc_data/bick_ju/OnlyMovement/new_060725/" + std::to_string(sim_num) + "_num_agents_area_ag.txt";
    //std::string person_file = "/home/bick_ju/Documents/INSIDe/persons/persons_scaled.csv";
    // calculate_agents_per_quantity_age_groups(output_file_loctype, save_file_loctype_ag, "locType", person_file);
    // calculate_agents_per_quantity_age_groups(output_file_ww_area, save_file_area_ag, "area", person_file);
    return 0;
}
