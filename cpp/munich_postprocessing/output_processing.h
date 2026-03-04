#include "json/reader.h"
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <string>
#include <chrono>
#include <vector>
#include <fstream>
#include "memilio/io/hdf5_cpp.h"
#include "memilio/utils/logging.h"
#include "boost/filesystem.hpp"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <json/value.h>

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
                                       std::string v4_or_v5, int num_sims, int start_sim)
{
    auto start    = std::chrono::high_resolution_clock::now();
    auto txt_file = fopen(save_file.c_str(), "w");
    if (txt_file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", save_file);
    }
    fprintf(txt_file, "Sim,t,%s,NumInfected,NewInfections,NumAgents", area_or_type.c_str());
    fprintf(txt_file, "\n");

    for (int sim = start_sim; sim < num_sims + start_sim; ++sim) {
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
                                      int num_sims, int start_sim)
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

    for (int sim = start_sim; sim < num_sims + start_sim; ++sim) {
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

void calculate_agents_per_quantity(std::string output_path, std::string save_file, int num_sims, int start_sim,
                                   std::string area_or_type, std::string v4_or_v5)
{
    auto file = fopen(save_file.c_str(), "w");
    if (file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", save_file);
    }
    fprintf(file, "Sim,t,%s,NumAgents", area_or_type.c_str());
    fprintf(file, "\n");
    for (int sim = start_sim; sim < num_sims + start_sim; ++sim) {
        auto start              = std::chrono::high_resolution_clock::now();
        std::string output_file = output_path + std::to_string(sim) + "_output_" + v4_or_v5 + ".h5";
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
        //Save data in matrix
        std::vector<std::vector<int>> matrix(dims[1], std::vector<int>(max_area + 1));
        for (size_t i = 0; i < dims[0]; ++i) {
            for (size_t j = 0; j < dims[1]; ++j) {
                int loc_type = data[i * dims[1] + j];
                matrix[j][loc_type] += 1;
            }
        }

        for (size_t i = 0; i < matrix.size(); ++i) {
            for (size_t j = 0; j < matrix[0].size(); ++j) {
                fprintf(file, "%d,", sim);
                fprintf(file, "%d,", static_cast<int>(i));
                fprintf(file, "%d,", static_cast<int>(j));
                fprintf(file, "%d,", matrix[i][j]);
                fprintf(file, "\n");
            }
        }
        auto stop     = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        std::cout << "Sim " << std::to_string(sim) << ": Time to calculate agents per quantity: " << duration.count()
                  << "[s]" << std::endl;
    }
    fclose(file);
}

std::string get_age_group_string(uint32_t age)
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

void calculate_agents_per_quantity_age_groups(std::string output_path, std::string save_file, std::string area_or_type,
                                              std::string person_file, int num_sims, int start_sim,
                                              std::string v4_or_v5)
{
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

    auto file = fopen(save_file.c_str(), "w");
    if (file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", save_file);
    }
    fprintf(file, "Sim,t,%s,AgeGroup,NumAgents", area_or_type.c_str());
    fprintf(file, "\n");
    for (int sim = start_sim; sim < num_sims + start_sim; ++sim) {
        auto start              = std::chrono::high_resolution_clock::now();
        std::string output_file = output_path + std::to_string(sim) + "_output_" + v4_or_v5 + ".h5";
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
        std::vector<std::vector<std::vector<int>>> matrix(
            dims[1], std::vector<std::vector<int>>(max_area + 1, std::vector<int>(6)));
        for (size_t agent = 0; agent < dims[0]; ++agent) {
            for (size_t j = 0; j < dims[1]; ++j) {
                int loc_type   = data[agent * dims[1] + j];
                std::string ag = get_age_group_string(ages[agent]);
                if (ag == "NaN") {
                    mio::log(mio::LogLevel::err, "Age group does not exists.}");
                }
                matrix[j][loc_type][age_groups[ag]] += 1;
            }
        }

        for (size_t agent = 0; agent < matrix.size(); ++agent) {
            for (size_t quantity = 0; quantity < matrix[0].size(); ++quantity) {
                for (size_t age_group = 0; age_group < matrix[0][0].size(); ++age_group) {
                    fprintf(file, "%d,", sim);
                    fprintf(file, "%d,", static_cast<int>(agent));
                    fprintf(file, "%d,", static_cast<int>(quantity));
                    fprintf(file, "%s,", index_to_age_groups[age_group].c_str());
                    fprintf(file, "%d,", matrix[agent][quantity][age_group]);
                    fprintf(file, "\n");
                }
            }
        }
        auto stop     = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        std::cout << "Sim " << std::to_string(sim)
                  << ": Time to calculate agents per quantity age groups: " << duration.count() << "[s]" << std::endl;
    }
    fclose(file);
}

double get_viral_load(int infection_age, double time_infected_in_days, double T_E, double T_INS, double T_symptoms)
{
    double peak           = 8.1;
    double t_peak_in_days = (T_E + T_INS) / 24.;
    double incline        = std::abs(peak / t_peak_in_days);
    double end_date       = time_infected_in_days;
    if (T_symptoms < 1) {
        t_peak_in_days = 0.5 * time_infected_in_days;
        peak           = incline * t_peak_in_days;
    }
    double decline = -peak / (end_date - t_peak_in_days);
    if (infection_age > 0 && infection_age / 24. <= end_date) {
        if (infection_age / 24. <= peak / incline) {
            return incline * (infection_age / 24.);
        }
        else {
            return peak + decline * (infection_age / 24. - peak / incline);
        }
    }
    else {
        return 0.;
    }
}

double get_infectivity(double T_E, double T_INS, double T_ISy, double T_ISev, double T_ICri, int infection_age)
{
    double time_shift            = 0.6 * T_E;
    double time_infected_in_days = (T_E + T_INS + T_ISy + T_ISev + T_ICri) / 24.;
    if (time_shift >= infection_age)
        return 0;
    ScalarType infectivity = 1 / (1 + exp(-(-7 + 1 * get_viral_load(infection_age - time_shift, time_infected_in_days,
                                                                    T_E, T_INS, T_ISy + T_ISev + T_ICri))));
    if ((infectivity < 0) || (infectivity > (1.0 / (1 + exp(-(-7 + 1 * 8.1)))))) {
        mio::log_error("Error: Infectivity smaller zero or above peak value.");
    }
    return 0.75 * infectivity;
}

void calculate_infected_per_measurement_station(std::string area_to_station_file, std::string infection_paths_file,
                                                std::string hdf5_output_file, std::string save_file, int timepoint,
                                                std::string station_name)
{
    // Get areas corresponding to measurement stations
    std::ifstream area_to_station_stream(area_to_station_file);
    if (!area_to_station_stream) {
        std::cerr << "Failed to open file\n";
        return;
    }

    Json::Value area_to_station;
    Json::CharReaderBuilder builder;
    JSONCPP_STRING errs;

    bool ok = Json::parseFromStream(builder, area_to_station_stream, &area_to_station, &errs);

    std::vector<int> areas;

    if (area_to_station.isMember(station_name)) {
        const Json::Value& arr = area_to_station[station_name];

        for (const auto& v : arr) {
            areas.push_back(std::stoi(v.asString()));
        }
    }

    // Get infected agent ids influencing the measurement station
    //Load H5 file
    mio::H5File sim_output{H5Fopen(hdf5_output_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)};
    std::string group_name = "data";
    //Open group
    mio::H5Group h5group{H5Gopen(sim_output.id, group_name.c_str(), H5P_DEFAULT)};

    //WW areas dataset
    //Open dataset
    mio::H5DataSet dset_loc_ids{H5Dopen(h5group.id, "loc_ids", H5P_DEFAULT)};
    //Open dataspace
    mio::H5DataSpace dataspace{H5Dget_space(dset_loc_ids.id)};
    //Get dimensions
    int ndims1     = H5Sget_simple_extent_ndims(dataspace.id);
    hsize_t* dims1 = (hsize_t*)malloc(ndims1 * sizeof(hsize_t));
    H5Sget_simple_extent_dims(dataspace.id, dims1, NULL);
    hsize_t rows1 = dims1[0];
    hsize_t cols1 = dims1[1];
    int* data1    = (int*)malloc(rows1 * cols1 * sizeof(int));
    //Read data to vector
    herr_t status1 = H5Dread(dset_loc_ids.id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1);

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

    std::vector<size_t> infected_agent_ids_at_station;
    std::vector<int> infection_ages;

    for (size_t agent = 0; agent < dims1[0]; ++agent) {
        int ww_area_at_t = data1[agent * dims1[1] + timepoint];
        if (std::find(areas.begin(), areas.end(), ww_area_at_t) != areas.end()) {
            int t_transm   = data2[agent * dims2[1]];
            int t_recovery = data2[agent * dims2[1] + 1];
            if (t_transm <= timepoint && t_recovery > timepoint) {
                infected_agent_ids_at_station.push_back(agent);
                infection_ages.push_back(timepoint - t_transm);
            }
        }
    }

    std::vector<double> shedding_values;

    // Read in infection paths
    int agent_id;
    double S, E, I_ns, I_sy, I_sev, I_cri, R, D;

    std::ifstream file(infection_paths_file);
    std::string header;
    std::getline(file, header);
    std::string line;

    size_t current_id = 0;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (iss >> agent_id >> S >> E >> I_ns >> I_sy >> I_sev >> I_cri >> R >> D) {
            if (static_cast<size_t>(agent_id) == infected_agent_ids_at_station[current_id]) {
                shedding_values.push_back(get_infectivity(E, I_ns, I_sy, I_sev, I_cri, infection_ages[current_id]));
                current_id += 1;
            }
        }
    }

    if (infected_agent_ids_at_station.size() != infection_ages.size() ||
        infected_agent_ids_at_station.size() != shedding_values.size()) {
        std::cerr << "Vector sizes do not match\n";
        return;
    }

    // Write results to file
    std::ofstream out(save_file);
    if (!out) {
        std::cerr << "Failed to open output file\n";
        return;
    }

    out << "Id,infection_age,shedding\n";

    for (size_t i = 0; i < infected_agent_ids_at_station.size(); ++i) {
        out << infected_agent_ids_at_station[i] << "," << infection_ages[i] << "," << shedding_values[i] << "\n";
    }

    std::cout << "\n";
}

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
