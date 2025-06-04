#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <chrono>
#include "memilio/io/hdf5_cpp.h"
#include "models/abm/location_type.h"
#include "memilio/utils/logging.h"

void calculate_agents_per_quantity(std::string output_file, std::string save_file, std::string area_or_type)
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
    std::vector<std::vector<int>> matrix(dims[1], std::vector<int>(max_area + 1));
    for (size_t i = 0; i < dims[0]; ++i) {
        for (size_t j = 0; j < dims[1]; ++j) {
            int loc_type = data[i * dims[1] + j];
            matrix[j][loc_type] += 1;
        }
    }
    std::cout << "Created matrix!\n";

    auto file = fopen(save_file.c_str(), "w");
    if (file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", save_file);
    }
    else {
        fprintf(file, "t,%s,NumAgents", area_or_type.c_str());
        fprintf(file, "\n");
        for (size_t i = 0; i < matrix.size(); ++i) {
            for (size_t j = 0; j < matrix[0].size(); ++j) {
                fprintf(file, "%d,", static_cast<int>(i));
                fprintf(file, "%d,", static_cast<int>(j));
                fprintf(file, "%d,", matrix[i][j]);
                fprintf(file, "\n");
            }
        }
        fclose(file);
    }
    auto stop     = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time to calculate agents per quantity: " << duration.count() << std::endl;
}

void calculate_infections_per_quantity(std::string output_path, std::string save_file, std::string area_or_type,
                                       std::string v4_or_v5, int num_sims)
{
    auto start    = std::chrono::high_resolution_clock::now();
    auto txt_file = fopen(save_file.c_str(), "w");
    if (txt_file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", save_file);
    }
    fprintf(txt_file, "Sim,t,%s,NumInfected,NewInfections", area_or_type.c_str());
    fprintf(txt_file, "\n");

    for (int sim = 0; sim < num_sims; ++sim) {
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

        hsize_t rows  = dims[0];
        hsize_t cols  = dims[1];
        int* data     = (int*)malloc(rows * cols * sizeof(int));
        herr_t status = H5Dread(dset_loc_ids.id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        std::vector<int> data1(data, data + dims[0] * dims[1]);
        int max_area = *std::max_element(data1.begin(), data1.end());

        //transmission recovery tp dataset
        //Open dataset
        mio::H5DataSet dset_transm_rec_tp{H5Dopen(h5group.id, "transm_recovery_tp", H5P_DEFAULT)};
        //Open dataspace
        mio::H5DataSpace dataspace2{H5Dget_space(dset_transm_rec_tp.id)};
        //Get dimensions
        int ndims2     = H5Sget_simple_extent_ndims(dataspace2.id);
        hsize_t* dims2 = (hsize_t*)malloc(ndims2 * sizeof(hsize_t));
        H5Sget_simple_extent_dims(dataspace2.id, dims2, NULL);

        hsize_t rows2  = dims2[0];
        hsize_t cols2  = dims2[1];
        int* data2     = (int*)malloc(rows2 * cols2 * sizeof(int));
        herr_t status2 = H5Dread(dset_transm_rec_tp.id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2);

        //Save data in matrix
        std::vector<std::vector<int>> matrix_new_infections(dims[1], std::vector<int>(max_area + 1));
        std::vector<std::vector<int>> matrix_infected(dims[1], std::vector<int>(max_area + 1));
        for (size_t i = 0; i < dims[0]; ++i) {
            int t_transm = data2[i * dims2[1]];
            if (t_transm > int(dims[1])) {
                continue;
            }
            int t_rec = std::min(data2[i * dims2[1] + 1], static_cast<int>(dims[1]));
            if (t_transm >= 0) {
                int loc_type_at_infection = data[i * dims[1] + t_transm];
                matrix_new_infections[t_transm][loc_type_at_infection] += 1;
            }
            for (int j = std::max(0, t_transm); j < t_rec; ++j) {
                int loc_type = data[i * dims[1] + j];
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
                fprintf(txt_file, "\n");
            }
        }
    }
    fclose(txt_file);
    auto stop     = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time to calculate infections per quantity: " << duration.count() << std::endl;
}

int main()
{
    int sim_num                     = 0;
    std::string output_file_loctype = "V:/bick_ju/test/" + std::to_string(sim_num) + "_output_v5.h5";
    std::string output_file_ww_area = "V:/bick_ju/test/" + std::to_string(sim_num) + "_output_v4.h5";
    std::string save_file_loctype   = "V:/bick_ju/test/" + std::to_string(sim_num) + "_num_agents_loc_type.txt";
    std::string save_file_loctype_infections =
        "V:/bick_ju/test/" + std::to_string(sim_num) + "_num_agents_infections.txt";
    std::string save_file_area = "V:/bick_ju/test/" + std::to_string(sim_num) + "_num_agents_area.txt";
    calculate_agents_per_quantity(output_file_loctype, save_file_loctype, "locType");
    calculate_agents_per_quantity(output_file_ww_area, save_file_area, "area");
    calculate_infections_per_quantity("V:/bick_ju/test/", save_file_loctype_infections, "locType", "v5", 1);
    return 0;
}
