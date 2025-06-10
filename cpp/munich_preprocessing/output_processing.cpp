#include <algorithm>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <chrono>
#include <map>
#include "memilio/io/hdf5_cpp.h"
#include "models/abm/location_type.h"
#include "memilio/utils/logging.h"
#include "boost/filesystem.hpp"

std::string determine_age_group(uint32_t age)
{
    if (age <= 4) {
        return "0-4";
    }
    else if (age <= 15) {
        return "5-15";
    }
    else if (age <= 34) {
        return "16-34";
    }
    else if (age <= 59) {
        return "35-59";
    }
    else if (age <= 79) {
        return "60-79";
    }
    else if (age > 79) {
        return "80+";
    }
    else {
        return "NaN";
    }
}

int stringToMinutes(const std::string& input)
{
    size_t colonPos = input.find(":");
    if (colonPos == std::string::npos) {
        // Handle invalid input (no colon found)
        return -1; // You can choose a suitable error code here.
    }

    std::string xStr = input.substr(0, colonPos);
    std::string yStr = input.substr(colonPos + 1);

    int x = std::stoi(xStr);
    int y = std::stoi(yStr);
    return x * 60 + y;
}

int longLatToInt(const std::string& input)
{
    double y = std::stod(input) * 1e+5; //we want the 5 numbers after digit
    return (int)y;
}

void split_line(std::string string, std::vector<int32_t>* row)
{
    std::vector<std::string> strings;
    boost::split(strings, string, boost::is_any_of(","));
    std::transform(strings.begin(), strings.end(), std::back_inserter(*row), [&](std::string s) {
        if (s.find(":") != std::string::npos) {
            return stringToMinutes(s);
        }
        else if (s.find(".") != std::string::npos) {
            return longLatToInt(s);
        }
        else {
            return std::stoi(s);
        }
    });
}

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
    std::cout << "Time to calculate agents per quantity: " << duration.count() << "[s]" << std::endl;
}

void calculate_agents_per_quantity_age_groups(std::string output_file, std::string save_file, std::string area_or_type,
                                              std::string person_file)
{
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<int> ages;
    // Read in persons
    const boost::filesystem::path p = person_file;
    if (!boost::filesystem::exists(p)) {
        mio::log_error("Cannot read in data. File does not exist.");
    }
    // File pointer
    std::fstream fin;

    // Open an existing file
    fin.open(person_file, std::ios::in);
    std::vector<int32_t> row;
    std::vector<std::string> row_string;
    std::string line;

    // Read the Titles from the Data file
    std::getline(fin, line);
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    std::vector<std::string> titles;
    boost::split(titles, line, boost::is_any_of(","));
    uint32_t count_of_titles              = 0;
    std::map<std::string, uint32_t> index = {};
    for (auto const& title : titles) {
        index.insert({title, count_of_titles});
        row_string.push_back(title);
        count_of_titles++;
    }

    while (std::getline(fin, line)) {
        row.clear();

        // read columns in this row
        split_line(line, &row);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        uint32_t age = row[index["age"]];
        ages.push_back(age);
    }

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

    hsize_t rows  = dims[0];
    hsize_t cols  = dims[1];
    int* data     = (int*)malloc(rows * cols * sizeof(int));
    herr_t status = H5Dread(dset_loc_ids.id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    std::vector<int> data1(data, data + dims[0] * dims[1]);
    int max_area = *std::max_element(data1.begin(), data1.end());

    //Age groups
    std::map<std::string, int> age_groups{{"0-4", 0},   {"5-15", 1},  {"16-34", 2},
                                          {"35-59", 3}, {"60-79", 4}, {"80+", 5}};
    //Age groups
    std::map<int, std::string> index_to_age_groups{{0, "0-4"},   {1, "5-15"},  {2, "16-34"},
                                                   {3, "35-59"}, {4, "60-79"}, {5, "80+"}};
    //Save data in matrix
    std::vector<std::vector<std::vector<int>>> matrix(dims[1],
                                                      std::vector<std::vector<int>>(max_area + 1, std::vector<int>(6)));
    for (size_t agent = 0; agent < dims[0]; ++agent) {
        for (size_t j = 0; j < dims[1]; ++j) {
            int loc_type   = data[agent * dims[1] + j];
            std::string ag = determine_age_group(ages[agent]);
            if (ag == "NaN") {
                mio::log(mio::LogLevel::err, "Age group does not exists.}");
            }
            matrix[j][loc_type][age_groups[ag]] += 1;
        }
    }
    std::cout << "Created matrix!\n";

    auto file = fopen(save_file.c_str(), "w");
    if (file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", save_file);
    }
    else {
        fprintf(file, "t,%s,AgeGroup,NumAgents", area_or_type.c_str());
        fprintf(file, "\n");
        for (size_t agent = 0; agent < matrix.size(); ++agent) {
            for (size_t quantity = 0; quantity < matrix[0].size(); ++quantity) {
                for (size_t age_group = 0; age_group < matrix[0][0].size(); ++age_group) {
                    fprintf(file, "%d,", static_cast<int>(agent));
                    fprintf(file, "%d,", static_cast<int>(quantity));
                    fprintf(file, "%s,", index_to_age_groups[age_group].c_str());
                    fprintf(file, "%d,", matrix[agent][quantity][age_group]);
                    fprintf(file, "\n");
                }
            }
        }
        fclose(file);
    }
    auto stop     = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time to calculate agents per quantity age groups: " << duration.count() << "[s]" << std::endl;
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
        //Iterate over all agents
        for (size_t agent = 0; agent < dims[0]; ++agent) {
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

int main()
{
    // int sim_num                     = 0;
    // std::string output_file_loctype = "/hpc_data/bick_ju/OnlyMovement/" + std::to_string(sim_num) + "_output_v5.h5";
    // std::string output_file_ww_area = "/hpc_data/bick_ju/OnlyMovement/" + std::to_string(sim_num) + "_output_v4.h5";
    // std::string save_file_loctype =
    //     "/hpc_data/bick_ju/OnlyMovement/" + std::to_string(sim_num) + "_num_agents_loc_type.txt";
    // std::string save_file_area = "/hpc_data/bick_ju/OnlyMovement/" + std::to_string(sim_num) + "_num_agents_area.txt";
    // calculate_agents_per_quantity(output_file_loctype, save_file_loctype, "locType");
    // calculate_agents_per_quantity(output_file_ww_area, save_file_area, "area");
    std::string output_dir_infections        = "";
    std::string save_file_loctype_infections = output_dir_infections + "num_agents_infections_loctype.txt";
    std::string save_file_ww_area_infections = output_dir_infections + "num_agents_infections_area.txt";
    int num_sims                             = 10;
    calculate_infections_per_quantity(output_dir_infections, save_file_loctype_infections, "locType", "v5", num_sims);
    calculate_infections_per_quantity(output_dir_infections, save_file_ww_area_infections, "area", "v4", num_sims);
    // std::string save_file_loctype_ag =
    //     "/hpc_data/bick_ju/OnlyMovement/" + std::to_string(sim_num) + "_num_agents_loc_type_ag.txt";
    // std::string save_file_area_ag =
    //     "/hpc_data/bick_ju/OnlyMovement/" + std::to_string(sim_num) + "_num_agents_area_ag.txt";
    // std::string person_file = "/home/bick_ju/Documents/INSIDe/persons/persons_final_corr.csv";
    // calculate_agents_per_quantity_age_groups(output_file_loctype, save_file_loctype_ag, "locType", person_file);
    // calculate_agents_per_quantity_age_groups(output_file_ww_area, save_file_area_ag, "area", person_file);
    return 0;
}
