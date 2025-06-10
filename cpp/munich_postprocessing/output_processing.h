#include <string>
#include <chrono>
#include <vector>
#include <fstream>
#include "memilio/io/hdf5_cpp.h"
#include "memilio/utils/logging.h"
#include "boost/filesystem.hpp"

void calculate_infections_per_quantity(std::string output_path, std::string save_file, std::string area_or_type,
                                       std::string v4_or_v5, int num_sims)
{
    auto start    = std::chrono::high_resolution_clock::now();
    auto txt_file = fopen(save_file.c_str(), "w");
    if (txt_file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", save_file);
    }
    fprintf(txt_file, "Sim,t,%s,NumInfected,NewInfections,NumAgents", area_or_type.c_str());
    fprintf(txt_file, "\n");

    for (int sim = 0; sim < num_sims; ++sim) {
        auto start_sim = std::chrono::high_resolution_clock::now();
        //Location dataset
        std::string file = output_path + std::to_string(sim) + "_output_" + v4_or_v5 + ".h5";
        //Load H5 file
        mio::H5File sim_output{H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)};
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
        int* data    = (int*)malloc(rows * cols * sizeof(int));
        //Read data to vector
        herr_t status = H5Dread(dset_loc_ids.id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        std::vector<int> data1(data, data + dims[0] * dims[1]);
        //Get the maximum area/location type; is needed as dimension for result matrix
        int max_area = *std::max_element(data1.begin(), data1.end());

        //Transmission recovery tp dataset
        //Open dataset
        mio::H5DataSet dset_transm_rec_tp{H5Dopen(h5group.id, "transm_recovery_tp", H5P_DEFAULT)};
        //Open dataspace
        mio::H5DataSpace dataspace2{H5Dget_space(dset_transm_rec_tp.id)};
        //Get dimensions
        int ndims2     = H5Sget_simple_extent_ndims(dataspace2.id);
        hsize_t* dims2 = (hsize_t*)malloc(ndims2 * sizeof(hsize_t));
        H5Sget_simple_extent_dims(dataspace2.id, dims2, NULL);
        hsize_t rows2 = dims2[0];
        hsize_t cols2 = dims2[1];
        int* data2    = (int*)malloc(rows2 * cols2 * sizeof(int));
        //Read data to vector
        herr_t status2 = H5Dread(dset_transm_rec_tp.id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2);

        //Save data in matrix
        //Number of new infections for every time point per area/location type and
        std::vector<std::vector<int>> matrix_new_infections(dims[1], std::vector<int>(max_area + 1));
        //Number of infected agents for every time point per area/location type
        std::vector<std::vector<int>> matrix_infected(dims[1], std::vector<int>(max_area + 1));
        //Number of total agents for every time point per area/location type
        std::vector<std::vector<int>> matrix(dims[1], std::vector<int>(max_area + 1));
        //Iterate over all agents
        for (size_t agent = 0; agent < dims[0]; ++agent) {
            for (size_t t = 0; t < dims[1]; ++t) {
                int loc_type = data[agent * dims[1] + t];
                matrix[t][loc_type] += 1;
            }
            int t_transm = data2[agent * dims2[1]];
            //If time point of transmission is bigger than last time point, the agent has not been infected during the simulation
            if (t_transm > int(dims[1])) {
                continue;
            }
            //Get recovery time point or simulation end time point if agent does not recover during the simulation time frame
            int t_rec = std::min(data2[agent * dims2[1] + 1], static_cast<int>(dims[1]));
            if (t_transm >= 0) { //Agent got infected during the simulation
                int loc_type_at_infection = data[agent * dims[1] + t_transm];
                matrix_new_infections[t_transm][loc_type_at_infection] += 1;
            }
            for (int j = std::max(0, t_transm); j < t_rec; ++j) {
                int loc_type = data[agent * dims[1] + j];
                matrix_infected[j][loc_type] += 1;
            }
        }

        for (size_t t = 0; t < matrix_infected.size(); ++t) {
            for (size_t loc = 0; loc < matrix_infected[0].size(); ++loc) {
                fprintf(txt_file, "%d,", sim);
                fprintf(txt_file, "%d,", static_cast<int>(t));
                fprintf(txt_file, "%d,", static_cast<int>(loc));
                fprintf(txt_file, "%d,", matrix_infected[t][loc]);
                fprintf(txt_file, "%d,", matrix_new_infections[t][loc]);
                fprintf(txt_file, "%d,", matrix[t][loc]);
                fprintf(txt_file, "\n");
            }
        }
        auto stop_sim = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop_sim - start_sim);
        std::cout << "Sim " << std::to_string(sim)
                  << ": Time to calculate infections per quantity: " << duration.count() << "[s]" << std::endl;
    }
    fclose(txt_file);
    auto stop     = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Total time to calculate infections per quantity: " << duration.count() << "[s]" << std::endl;
}
