#include <cstddef>
#include <cstdio>
#include <string>
#include <chrono>
#include <vector>
#include <fstream>
#include "memilio/io/hdf5_cpp.h"
#include "memilio/utils/logging.h"
#include "boost/filesystem.hpp"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

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

int determine_age_group_index(uint32_t age)
{
    if (age <= 4) {
        return 0;
    }
    else if (age <= 15) {
        return 1;
    }
    else if (age <= 34) {
        return 2;
    }
    else if (age <= 59) {
        return 3;
    }
    else if (age <= 79) {
        return 4;
    }
    else if (age > 79) {
        return 5;
    }
    else {
        return -1;
    }
}

void calculate_infections_per_hh_size(std::string output_path, std::string person_file, std::string save_file,
                                      int num_sims)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<int> homes;
    std::map<int, int> home_to_hh_size;
    std::vector<int> age_groups;

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
    int p_index = -1;
    while (std::getline(fin, line)) {
        row.clear();
        p_index++;

        // read columns in this row
        split_line(line, &row);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        int home_id = row[index["home_id"]];
        homes.push_back(home_id);
        auto iter_home = home_to_hh_size.find(home_id);
        if (iter_home == home_to_hh_size.end()) {
            home_to_hh_size.insert({home_id, 1});
        }
        else {
            home_to_hh_size[home_id] += 1;
        }
        int age_group = determine_age_group_index(row[index["age"]]);
        age_groups.push_back(age_group);
    }

    // Get maximum household size and maximum age group
    int max_hh_size =
        std::max_element(home_to_hh_size.begin(), home_to_hh_size.end(), [](const auto& a, const auto& b) {
            return a.second < b.second;
        })->second;

    int max_ag = *std::max_element(age_groups.begin(), age_groups.end());

    auto txt_file = fopen(save_file.c_str(), "w");
    if (txt_file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", save_file);
    }
    fprintf(txt_file, "Sim,t,");
    for (int i = 0; i < max_hh_size; ++i) {
        fprintf(txt_file, "NumAgents%dHH,NumInfected%dHH,", i + 1, i + 1);
    }
    for (int i = 0; i < max_ag + 1; ++i) {
        fprintf(txt_file, "NumAgentsAg%d,NumInfectedAg%d,", i, i);
    }
    fprintf(txt_file, "\n");

    for (int sim = 0; sim < num_sims; ++sim) {
        auto start_sim = std::chrono::high_resolution_clock::now();
        //Location dataset
        std::string file = output_path + std::to_string(sim) + "_output_v4.h5";
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
        //Number of agents per hh size
        std::vector<std::vector<int>> matrix_num_agents_hh_size(dims[1], std::vector<int>(max_hh_size));
        //Number of infected per hh size
        std::vector<std::vector<int>> matrix_num_infected_hh_size(dims[1], std::vector<int>(max_hh_size));
        //Number of agents per age group
        std::vector<std::vector<int>> matrix_num_agents_ag(dims[1], std::vector<int>(max_ag + 1));
        //Number of infected per age group
        std::vector<std::vector<int>> matrix_num_infected_ag(dims[1], std::vector<int>(max_ag + 1));

        //Iterate over all agents
        for (size_t agent = 0; agent < dims[0]; ++agent) {
            for (size_t t = 0; t < dims[1]; ++t) {
                matrix_num_agents_hh_size[t][home_to_hh_size[homes[agent]] - 1] += 1;
                matrix_num_agents_ag[t][age_groups[agent]] += 1;
            }
            int t_transm = data2[agent * dims2[1]];
            //If time point of transmission is bigger than last time point, the agent has not been infected during the simulation
            if (t_transm > int(dims[1])) {
                continue;
            }
            //Get recovery time point or simulation end time point if agent does not recover during the simulation time frame
            int t_rec = std::min(data2[agent * dims2[1] + 1], static_cast<int>(dims[1]));
            for (int j = std::max(0, t_transm); j < t_rec; ++j) {
                matrix_num_infected_hh_size[j][home_to_hh_size[homes[agent]] - 1] += 1;
                matrix_num_infected_ag[j][age_groups[agent]] += 1;
            }
        }

        for (size_t t = 0; t < matrix_num_agents_hh_size.size(); ++t) {
            fprintf(txt_file, "%d,", sim);
            fprintf(txt_file, "%d,", static_cast<int>(t));
            for (size_t hh_size = 0; hh_size < matrix_num_agents_hh_size[0].size(); ++hh_size) {
                fprintf(txt_file, "%d,", matrix_num_agents_hh_size[t][hh_size]);
                fprintf(txt_file, "%d,", matrix_num_infected_hh_size[t][hh_size]);
            }
            for (size_t ag = 0; ag < matrix_num_agents_ag[0].size(); ++ag) {
                fprintf(txt_file, "%d,", matrix_num_agents_ag[t][ag]);
                fprintf(txt_file, "%d,", matrix_num_infected_ag[t][ag]);
            }
            fprintf(txt_file, "\n");
        }
        auto stop_sim = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop_sim - start_sim);
        std::cout << "Sim " << std::to_string(sim)
                  << ": Time to calculate infections hh size and age group: " << duration.count() << "[s]" << std::endl;
    }
    fclose(txt_file);
    auto stop     = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Total time to calculate infections hh size and age group: " << duration.count() << "[s]" << std::endl;
}
