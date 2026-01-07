#include "../include/multi_run_simulator.h"
#include "../include/file_utils.h"
#include "abm/time.h"
#include "memilio/io/io.h"
#include "abm/analyze_result.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <hdf5.h>
#include "city_builder.h"

std::string location_type_to_string(mio::abm::LocationType type)
{
    switch (type) {
    case mio::abm::LocationType::Home:
        return "Home";
    case mio::abm::LocationType::Work:
        return "Work";
    case mio::abm::LocationType::School:
        return "School";
    case mio::abm::LocationType::SocialEvent:
        return "SocialEvent";
    case mio::abm::LocationType::BasicsShop:
        return "BasicsShop";
    case mio::abm::LocationType::Hospital:
        return "Hospital";
    case mio::abm::LocationType::ICU:
        return "ICU";
    case mio::abm::LocationType::EventPanvadere:
        return "EventPanvadere";
    default:
        return "Unknown";
    }
}

mio::IOResult<MultiRunResults> MultiRunSimulator::run_multi_simulation(const MultiRunConfig& config,
                                                                       mio::RandomNumberGenerator rng)
{
    MultiRunResults results;
    results.event_type      = config.event_config.type;
    results.simulation_type = config.simulation_type;

    for (int run = 0; run < config.num_runs; ++run) {
        if (run % 10 == 0 || run == config.num_runs - 1) {
            std::cout << "Run " << run + 1 << "/" << config.num_runs << std::endl;
        }

        // Step 1: Build city
        // std::cout << "Building city..." << std::endl;
        auto run_rng_counter =
            mio::rng_totalsequence_counter<uint64_t>(static_cast<uint32_t>(0), mio::Counter<uint32_t>(0));
        rng.set_counter(run_rng_counter);
        auto base_world = CityBuilder::build_world(config.city_config, rng);
        // CityBuilder::print_city_summary(config.city_config);

        for (auto& person : base_world.get_persons()) {
            // Reset the infection state for each person
            person.get_rng_counter() =
                mio::Counter<uint32_t>(run * 1000000 + run); // Ensure unique counter for each run
        }
        run_rng_counter =
            mio::rng_totalsequence_counter<uint64_t>(static_cast<uint32_t>(run), mio::Counter<uint32_t>(0));
        rng.set_counter(run_rng_counter);
        // We could reset it to zero if we dont want random assignment of infections.

        // Step 2: Get map from specific event to ids of persons in simulation
        // std::cout << "Mapping events to person IDs..." << std::endl;
        BOOST_OUTCOME_TRY(auto event_map, EventSimulator::map_events_to_persons(base_world, config.event_config.type));

        // Step 3: Get initial infections
        std::vector<uint32_t> initial_infections;
        if (config.simulation_type == SimType::Panvadere) {
            // std::cout << "Loading infections from Panvadere..." << std::endl;
            BOOST_OUTCOME_TRY(auto infections,
                              EventSimulator::initialize_from_panvadere(config.event_config.type, event_map));
            initial_infections = infections;
        }
        else if (config.simulation_type == SimType::Memilio) {
            // std::cout << "Simulating event infections..." << std::endl;
            BOOST_OUTCOME_TRY(auto infections,
                              EventSimulator::initialize_from_event_simulation(config.event_config.type, event_map));
            initial_infections = infections;
        }

        // Step 4: Run multiple simulations
        // std::cout << "Running " << config.num_runs << " simulations..." << std::endl;
        results.all_runs.reserve(config.num_runs);

                auto single_result = run_single_simulation_with_infections(
            std::move(base_world), initial_infections, config.infection_parameter_k, config.simulation_days);

        if (single_result.has_value()) {
            results.all_runs.push_back(single_result.value());
            results.infection_parameter_k = config.infection_parameter_k;
            results.successful_runs++;
        }
        else {
            std::cerr << "Run " << run << " failed: " << single_result.error().message() << std::endl;
        }
    }
    return mio::success(results);
}

std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>>
calculate_contact_hours(const std::vector<std::vector<std::vector<std::tuple<uint32_t, uint32_t>>>>& simulation_data)
{
    // Map to store contact hours: (person1_id, person2_id) -> total_hours
    // We use a map with pair as key to accumulate hours across all runs and timesteps
    std::map<std::tuple<uint32_t, uint32_t, uint32_t>, uint32_t> contact_hours;

    // Iterate through all simulation runs
    for (const auto& run : simulation_data) {
        // Iterate through all time steps in this run
        for (const auto& timestep : run) {
            // Group persons by location for this timestep
            std::unordered_map<uint32_t, std::vector<uint32_t>> location_to_persons;

            // Build location -> persons mapping for current timestep
            for (const auto& person_location : timestep) {
                uint32_t person_id   = std::get<0>(person_location);
                uint32_t location_id = std::get<1>(person_location);
                location_to_persons[location_id].push_back(person_id);
            }

            // For each location, count contacts between all pairs of persons
            for (const auto& location_entry : location_to_persons) {
                uint32_t location_id            = location_entry.first;
                const auto& persons_at_location = location_entry.second;

                // If there are at least 2 persons at this location
                if (persons_at_location.size() >= 2) {
                    // Generate all pairs of persons at this location
                    for (size_t i = 0; i < persons_at_location.size(); ++i) {
                        for (size_t j = i + 1; j < persons_at_location.size(); ++j) {
                            uint32_t person1 = persons_at_location[i];
                            uint32_t person2 = persons_at_location[j];

                            // Ensure consistent ordering (smaller ID first)
                            if (person1 > person2) {
                                std::swap(person1, person2);
                            }

                            // Increment contact hours for this pair
                            contact_hours[{person1, person2, location_id}]++;
                        }
                    }
                }
            }
        }
    }

    // Convert map to vector of tuples
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> result;
    result.reserve(contact_hours.size());

    for (const auto& entry : contact_hours) {
        result.emplace_back(std::get<0>(entry.first), std::get<1>(entry.first), std::get<2>(entry.first), entry.second);
    }

    return result;
}

mio::IOResult<void> save_contact_intenseness(
    const std::vector<std::vector<std::vector<std::tuple<uint32_t, uint32_t>>>>& history_person_and_location_id,
    const std::string& base_dir)
{
    // we  have std::vector<std::vector<std::tuple<uint32_t, uint32_t>>> which describes a vector where we have:
    // 1. Person ID
    // 2. Location ID
    // inner vector the persons
    // middle vector represents the time steps
    // outer vector represents the different simulation runs

    auto contact_map = calculate_contact_hours(history_person_and_location_id);

    // Save contact intensiveness to a CSV file
    std::string contact_intensiveness_file = base_dir + "/contact_intensiveness.csv";
    std::ofstream contact_file(contact_intensiveness_file);
    if (contact_file.is_open()) {
        contact_file << "Person 1, Person 2, Location ID, Contact Hours\n";
        for (const auto& entry : contact_map) {
            contact_file << std::get<0>(entry) << ", " << std::get<1>(entry) << ", " << std::get<2>(entry) << ", "
                         << std::get<3>(entry) << "\n";
        }
        contact_file.close();
    }

    return mio::success();
}

mio::IOResult<void> save_amount_infected_per_id(
    const std::vector<std::vector<std::vector<std::tuple<uint32_t, uint32_t, mio::abm::LocationType>>>>&
        detailed_infection,
    const std::string& base_dir)
{
    // We count the number of infections per person ID
    std::unordered_map<uint32_t, uint32_t> infection_count;

    for (const auto& time_series : detailed_infection) {
        for (const auto& entry : time_series) {
            for (size_t i = 0; i < entry.size(); ++i) {
                if (entry.size() < 2) {
                    continue; // Skip entries that do not have enough data
                }
                uint32_t person_id = std::get<0>(entry[i]);
                infection_count[person_id] += 1;
            }
        }
    }

    // Save the infection counts to a CSV file
    std::string infection_count_file = base_dir + "/infection_count.csv";
    std::ofstream count_file(infection_count_file);
    if (count_file.is_open()) {
        count_file << "Person ID, Infected Amount\n";
        for (const auto& entry : infection_count) {
            count_file << entry.first << ", " << entry.second << "\n";
        }
        count_file.close();
    }

    return mio::success();
}

mio::IOResult<void> save_detailed_infection_for_all_runs(
    const std::vector<std::vector<std::vector<std::tuple<uint32_t, uint32_t, mio::abm::LocationType>>>>&
        detailed_infection,
    const std::string& base_dir)
{
    // Create directory for all runs
    std::string all_runs_dir = base_dir + "/all_runs_detailed_infections";
    std::filesystem::create_directories(all_runs_dir);

    std::cout << "Saving detailed infection data for all " << detailed_infection.size() << " runs..." << std::endl;

    // Save detailed infection data for all runs
    for (size_t run_idx = 0; run_idx < detailed_infection.size(); ++run_idx) {
        std::vector<uint32_t> persons_at_event;
        bool patient_zero = true;

        std::string run_file = all_runs_dir + "/run_" + std::to_string(run_idx) + "_detailed_infection.csv";
        std::ofstream run_output(run_file);

        if (run_output.is_open()) {
            run_output << "Timestep,Person_ID,Location_ID,Location_Type\n";

            for (int timestep = 0; timestep < (int)detailed_infection[run_idx].size(); ++timestep) {
                for (const auto& infection_entry : detailed_infection[run_idx][timestep]) {
                    auto use_timestep = timestep;
                    if (std::get<2>(infection_entry) == mio::abm::LocationType::EventPanvadere) {
                        persons_at_event.push_back(std::get<0>(infection_entry));
                        use_timestep = timestep - 48;
                        if (patient_zero) {
                            use_timestep = timestep - 96;
                            patient_zero = false;
                        }
                    }
                    run_output << use_timestep << "," << std::get<0>(infection_entry) << ","
                               << std::get<1>(infection_entry) << ","
                               << location_type_to_string(std::get<2>(infection_entry)) << "\n";
                }
            }
            run_output.close();
        }

        if ((run_idx + 1) % 10 == 0 || run_idx == detailed_infection.size() - 1) {
            std::cout << "Saved detailed infection data for run " << run_idx + 1 << "/" << detailed_infection.size()
                      << std::endl;
        }
    }

    std::cout << "All detailed infection files saved to: " << all_runs_dir << std::endl;
    return mio::success();
}

mio::IOResult<void> save_detailed_infection_and_contact_for_best_run(
    const std::vector<std::vector<std::vector<std::tuple<uint32_t, uint32_t, mio::abm::LocationType>>>>&
        detailed_infection,
    const std::string& dir_h5_50_percentile_run,
    std::vector<std::vector<mio::TimeSeries<ScalarType>>>& history_infected_amount,
    std::vector<std::vector<std::vector<std::tuple<uint32_t, uint32_t>>>>& history_contact_hours, EventType event_type,
    std::vector<std::vector<std::vector<std::tuple<uint32_t, bool>>>>& history_infected_status)
{
    // First we need to find the run, which has the lowest discrete L2 norm in relation to the median value of all runs, saved in dir_h5_50_percentile_run(saved as h5)
    // First we calculate a vector of all runs and the amount of infected at each timestep from history_infected_amount
    std::vector<std::vector<ScalarType>> infected_amount_per_run;
    for (const auto& time_series : history_infected_amount) {
        std::vector<ScalarType> run_data;
        for (int i = 0; i < time_series[0].get_num_time_points(); ++i) {
            run_data.push_back(static_cast<ScalarType>(time_series[0].get_value(i)(0)));
        }
        infected_amount_per_run.push_back(run_data);
    }

    // Now we read the h5 median file and save it as a vector
    auto median_file = dir_h5_50_percentile_run + "/Results.h5";
    std::vector<ScalarType> median_data;

    // Open the HDF5 file using official HDF5 C API
    hid_t file_id = H5Fopen(median_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        std::cout << "Failed to open HDF5 file: " << median_file << std::endl;
        return mio::success();
    }
    std::cout << "Successfully opened: " << median_file << std::endl;

    // Based on the Python code, the data is stored under group "0"
    hid_t group_id = H5Gopen2(file_id, "0", H5P_DEFAULT);
    if (group_id < 0) {
        std::cout << "Failed to open group '0'" << std::endl;
        H5Fclose(file_id);
        return mio::success();
    }
    std::cout << "Successfully opened group '0'" << std::endl;

    // First, let's see what objects are in this group
    hsize_t num_objs;
    H5Gget_num_objs(group_id, &num_objs);
    std::cout << "Number of objects in group '0': " << num_objs << std::endl;

    for (hsize_t i = 0; i < num_objs; i++) {
        char name[256];
        H5Gget_objname_by_idx(group_id, i, name, sizeof(name));
        std::cout << "Object " << i << ": " << name << std::endl;
    }

    // Try to open the "Total" dataset (which contains the infection data we need)
    hid_t dataset_id = H5Dopen2(group_id, "Total", H5P_DEFAULT);
    if (dataset_id < 0) {
        std::cout << "Failed to open dataset 'Total'" << std::endl;
        H5Gclose(group_id);
        H5Fclose(file_id);
        return mio::success();
    }
    std::cout << "Successfully opened dataset 'Total'" << std::endl;

    // Get the dataspace of the dataset
    hid_t dataspace_id = H5Dget_space(dataset_id);
    if (dataspace_id < 0) {
        std::cout << "Failed to get dataspace" << std::endl;
        H5Dclose(dataset_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return mio::success();
    }

    // Get the number of dimensions and their sizes
    int ndims = H5Sget_simple_extent_ndims(dataspace_id);
    std::vector<hsize_t> dims(ndims);
    H5Sget_simple_extent_dims(dataspace_id, dims.data(), nullptr);

    std::cout << "Dataset 'Total' dimensions: ";
    for (int i = 0; i < ndims; ++i) {
        std::cout << dims[i];
        if (i < ndims - 1)
            std::cout << " x ";
    }
    std::cout << std::endl;

    // Calculate total number of elements
    hsize_t total_elements = 1;
    for (int i = 0; i < ndims; ++i) {
        total_elements *= dims[i];
    }
    std::cout << "Total elements: " << total_elements << std::endl;

    // Read the data - try different data types
    std::vector<double> temp_data(total_elements);
    herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data.data());

    if (status < 0) {
        std::cout << "Failed to read as double, trying float..." << std::endl;
        std::vector<float> temp_float_data(total_elements);
        status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_float_data.data());
        if (status >= 0) {
            for (size_t i = 0; i < total_elements; ++i) {
                temp_data[i] = static_cast<double>(temp_float_data[i]);
            }
            std::cout << "Successfully read as float" << std::endl;
        }
        else {
            std::cout << "Failed to read as float, trying int..." << std::endl;
            std::vector<int> temp_int_data(total_elements);
            status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_int_data.data());
            if (status >= 0) {
                for (size_t i = 0; i < total_elements; ++i) {
                    temp_data[i] = static_cast<double>(temp_int_data[i]);
                }
                std::cout << "Successfully read as int" << std::endl;
            }
        }
    }
    else {
        std::cout << "Successfully read as double" << std::endl;
    }

    if (status < 0) {
        std::cout << "Failed to read data from HDF5 file" << std::endl;
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return mio::success();
    }

    // std::cout << "Data read successfully, checking values..." << std::endl;
    // std::cout << "First 1000 values: ";
    // for (size_t i = 0; i < std::min(size_t(1000), temp_data.size()); ++i) {
    //     std::cout << temp_data[i] << " ";
    // }
    // std::cout << std::endl;

    // Convert to ScalarType and handle different dimensions
    if (ndims == 1) {
        // 1D dataset - use directly
        median_data.reserve(dims[0]);
        for (hsize_t i = 0; i < dims[0]; ++i) {
            median_data.push_back(static_cast<ScalarType>(temp_data[i]));
        }
    }
    else if (ndims == 2) {
        // 2D dataset - sum across all infection states for each time point
        median_data.reserve(dims[0]);
        for (hsize_t i = 0; i < dims[0]; ++i) {
            // Sum across all infection states for each time point
            ScalarType sum = 0.0;
            for (hsize_t j = 0; j < dims[1]; ++j) {
                sum += static_cast<ScalarType>(temp_data[i * dims[1] + j]);
            }
            median_data.push_back(sum);
        }
    }
    else {
        // Higher dimensional - flatten and take first elements
        median_data.reserve(dims[0]);
        hsize_t step = total_elements / dims[0];
        for (hsize_t i = 0; i < dims[0]; ++i) {
            median_data.push_back(static_cast<ScalarType>(temp_data[i * step]));
        }
    }

    // Clean up HDF5 resources
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Gclose(group_id);
    H5Fclose(file_id);

    // Now find the run with the lowest L2 norm relative to median
    size_t best_run_index  = 0;
    ScalarType min_l2_norm = std::numeric_limits<ScalarType>::max();

    for (size_t run_idx = 0; run_idx < infected_amount_per_run.size(); ++run_idx) {
        const auto& run_data = infected_amount_per_run[run_idx];

        // Calculate L2 norm (squared differences)
        ScalarType l2_norm = 0.0;
        size_t min_size    = std::min(run_data.size(), median_data.size());

        for (size_t i = 0; i < min_size; ++i) {
            ScalarType diff = run_data[i] - median_data[i];
            l2_norm += diff * diff;
        }

        if (l2_norm < min_l2_norm) {
            min_l2_norm    = l2_norm;
            best_run_index = run_idx;
        }
        std::cout << "Run " << run_idx << ": L2 norm = " << l2_norm << std::endl;
    }

    // Save the detailed infection data for the best run
    // we also save the infected
    std::vector<uint32_t> persons_at_event;
    bool patient_zero = true;
    if (best_run_index < detailed_infection.size()) {
        std::string best_run_file = dir_h5_50_percentile_run + "/../../best_run_detailed_infection.csv";
        std::ofstream best_run_output(best_run_file);

        if (best_run_output.is_open()) {
            best_run_output << "Timestep,Person_ID,Location_ID,Location_Type\n";

            for (int timestep = 0; timestep < (int)detailed_infection[best_run_index].size(); ++timestep) {
                for (const auto& infection_entry : detailed_infection[best_run_index][timestep]) {
                    auto use_timestep = timestep;
                    if (std::get<2>(infection_entry) == mio::abm::LocationType::EventPanvadere) {
                        persons_at_event.push_back(std::get<0>(infection_entry));
                        use_timestep = timestep - 48;
                        if (patient_zero) {
                            use_timestep = timestep - 96;
                            patient_zero = false;
                        }
                    }
                    best_run_output << use_timestep << "," << std::get<0>(infection_entry) << ","
                                    << std::get<1>(infection_entry) << ","
                                    << location_type_to_string(std::get<2>(infection_entry)) << "\n";
                }
            }
            best_run_output.close();
        }
    }

    // We also want to save who was where at which timepoint
    std::string contact_hours_file = dir_h5_50_percentile_run + "/../../best_run_contact_data.csv";
    std::ofstream contact_hours_output(contact_hours_file);
    auto event_time = 0;
    if (event_type == EventType::WorkMeeting_Baseline_Meetings || event_type == EventType::WorkMeeting_Many_Meetings) {
        event_time = 8;
    }
    else {
        event_time = 2; // For other events, we assume 2 hours
    }
    if (contact_hours_output.is_open()) {
        contact_hours_output << "Timestep,Person_ID,Location_ID\n";

        // We need to add the time at the Event where persons infected each other, we asusme it was 2 days before so -48hours and then for the work event 8 hours and for the restaurant event 2 hours.
        for (auto&& entry : persons_at_event) {
            for (int hours = 0; hours < event_time; ++hours) {
                contact_hours_output << -48 + hours << "," << entry << "," << mio::abm::INVALID_LOCATION_INDEX
                                     << "\n"; // Assuming location ID 0 for the event
            }
        }

        for (size_t timestep = 0; timestep < history_contact_hours[best_run_index].size(); ++timestep) {
            for (const auto& contact_entry : history_contact_hours[best_run_index][timestep]) {
                contact_hours_output << timestep << "," << std::get<0>(contact_entry) << ","
                                     << std::get<1>(contact_entry) << "\n";
            }
        }
        contact_hours_output.close();
    }

    // Save infected status
    auto infected_status_file = dir_h5_50_percentile_run + "/../../best_run_infected_status.csv";
    std::ofstream infected_status_output(infected_status_file);
    if (infected_status_output.is_open()) {
        infected_status_output << "Timestep,Person_ID,Infected\n";

        for (int hours = 0; hours < event_time; ++hours) {
            contact_hours_output << -48 + hours << "," << persons_at_event[0] << "," << "true"
                                 << "\n"; // Assuming infected
        }

        for (size_t timestep = 0; timestep < history_infected_status[best_run_index].size(); ++timestep) {
            for (const auto& infected_entry : history_infected_status[best_run_index][timestep]) {
                infected_status_output << timestep << "," << std::get<0>(infected_entry) << ","
                                       << std::get<1>(infected_entry) << "\n";
            }
        }
        infected_status_output.close();
    }

    return mio::success();
}

mio::IOResult<void> MultiRunSimulator::save_multi_run_results(const MultiRunResults& results,
                                                              const std::string& base_dir, const MultiRunConfig& config)
{
    BOOST_OUTCOME_TRY(create_result_folders(base_dir));

    // Save percentile results
    std::string summary_file = base_dir + "/summary.txt";
    std::ofstream summary(summary_file);
    if (summary.is_open()) {
        summary << "Multi-Run Simulation Summary\n";
        summary << "Event Type: " << EventSimulator::event_type_to_string(results.event_type) << "\n";
        summary << "Initialization: " << EventSimulator::simulation_type_to_string(results.simulation_type) << "\n";
        summary << "Infection Parameter K: " << results.infection_parameter_k << "\n";
        summary << "Total Runs: " << results.all_runs.size() << "\n";
        summary.close();
    }
    
    // Save average household size of initially infected persons
    std::string household_size_file = base_dir + "/average_household_size_initial_infections.csv";
    std::ofstream household_size_stream(household_size_file);
    if (household_size_stream.is_open()) {
        household_size_stream << "Run,Average_Household_Size,Average_Persons_Above_Age_Group_0\n";
        for (size_t run_idx = 0; run_idx < results.all_runs.size(); ++run_idx) {
            household_size_stream << run_idx << "," 
                                << results.all_runs[run_idx].average_household_size_of_initial_infections << ","
                                << results.all_runs[run_idx].average_persons_above_age_group_0_in_initial_households << "\n";
        }
        household_size_stream.close();
    }

    // Save household id
    auto household_id_file = base_dir + "/household_id.csv";
    std::ofstream household_file(household_id_file);
    if (household_file.is_open()) {
        household_file << "Person_ID,Household_ID\n";
        for (const auto& person_info : results.all_runs[0].history_household_id) {
            household_file << std::get<0>(person_info) << "," << std::get<1>(person_info) << "\n";
        }
        household_file.close();
    }
    // Save work id
    auto work_id_file = base_dir + "/work_id.csv";
    std::ofstream work_file(work_id_file);
    if (work_file.is_open()) {
        work_file << "Person_ID,Work_ID\n";
        for (const auto& person_info : results.all_runs[0].history_work_id) {
            work_file << std::get<0>(person_info) << "," << std::get<1>(person_info) << "\n";
        }
        work_file.close();
    }

    //Location id and person id
    // std::string location_type_and_id_file = base_dir + "/location_type_and_id.txt";
    // std::ofstream location_file(location_type_and_id_file);
    // if (location_file.is_open()) {
    //     location_file << "Person ID, Location Type\n";
    //     for (const auto& timestep : results.all_runs[0].infection_per_location_type_and_id) {
    //         for (const auto& person_info : timestep) {
    //             location_file << std::get<0>(person_info) << ", " << location_type_to_string(std::get<1>(person_info))
    //                           << "\n";
    //         }
    //         location_file << "\n"; // Separate time steps with a newline
    //     }
    //     location_file.close();
    // }
    // Save id to location_type
    auto loc_id_and_type_file = base_dir + "/location_id_and_type.csv";
    std::ofstream loc_id_and_type_stream(loc_id_and_type_file);
    if (loc_id_and_type_stream.is_open()) {
        loc_id_and_type_stream << "Location_ID,Location_Type\n";
        for (const auto& loc_info : results.all_runs[0].loc_id_and_type) {
            loc_id_and_type_stream << std::get<0>(loc_info) << "," << location_type_to_string(std::get<1>(loc_info))
                                   << "\n";
        }
        loc_id_and_type_stream.close();
    }

    CityBuilder::save_city_to_file(config.city_config, base_dir + "/city_config.txt");

    //Save percentile results
    auto ensembl_inf_per_loc_type      = std::vector<std::vector<mio::TimeSeries<ScalarType>>>{};
    auto ensembl_inf_per_age_group     = std::vector<std::vector<mio::TimeSeries<ScalarType>>>{};
    auto ensemble_amount_of_infections = std::vector<std::vector<mio::TimeSeries<ScalarType>>>{};

    auto ensemble_contact_hours = std::vector<std::vector<std::vector<std::tuple<uint32_t, uint32_t>>>>{};

    auto ensemble_detailed_infection =
        std::vector<std::vector<std::vector<std::tuple<uint32_t, uint32_t, mio::abm::LocationType>>>>{};
    auto ensemble_infection_history = std::vector<std::vector<std::vector<std::tuple<uint32_t, bool>>>>{};

    auto ensemble_params              = std::vector<std::vector<mio::abm::World>>{};
    auto ensemble_params_no_agegroups = std::vector<std::vector<mio::abm::World>>{};

    for (const auto& run : results.all_runs) {
        ensembl_inf_per_loc_type.push_back(run.infection_per_loc_type);
        ensembl_inf_per_age_group.push_back(run.infection_state_per_age_group);
        ensemble_amount_of_infections.push_back(run.history_infected_amount);
        ensemble_params.push_back(run.ensemble_params);
        ensemble_params_no_agegroups.push_back(run.ensemble_params_no_agegroups);
        ensemble_contact_hours.push_back(run.history_person_and_location_id);
        ensemble_detailed_infection.push_back(run.history_detailed_infection);
        ensemble_infection_history.push_back(run.history_infected_status);
    }

    // BOOST_OUTCOME_TRY(save_results(ensembl_inf_per_loc_type, ensemble_params, {0},
    //                                base_dir + "/infection_per_location_type_per_age_group", false, true));
    // BOOST_OUTCOME_TRY(save_results(ensembl_inf_per_age_group, ensemble_params, {0},
    //                                base_dir + "/infection_state_per_age_group", false, true));

    BOOST_OUTCOME_TRY(save_results(ensemble_amount_of_infections, ensemble_params_no_agegroups, {0},
                                   base_dir + "/amount_of_infections", false, true));

    // BOOST_OUTCOME_TRY(save_contact_intenseness(ensemble_contact_hours, base_dir));

    // BOOST_OUTCOME_TRY(save_amount_infected_per_id(ensemble_detailed_infection, base_dir));

    // BOOST_OUTCOME_TRY(save_detailed_infection_for_all_runs(ensemble_detailed_infection, base_dir));

    // BOOST_OUTCOME_TRY(save_detailed_infection_and_contact_for_best_run(
    //     ensemble_detailed_infection, base_dir + "/amount_of_infections/p50", ensemble_amount_of_infections,
    //     ensemble_contact_hours, results.event_type, ensemble_infection_history));

    // Short sanity check if detailed infections final infection count is the same as the amount of infected agent counted otherwise
    // Just print the run number and an x vs y print
    // for (size_t run_idx = 0; run_idx < ensemble_detailed_infection.size(); ++run_idx) {
    //     size_t amount_infected = 0;
    //     for (size_t t = 0; t < ensemble_detailed_infection[run_idx].size(); ++t) {
    //         amount_infected += ensemble_detailed_infection[run_idx][t].size();
    //     }
    //     std::cout << "Run " << run_idx + 1 << ": " << amount_infected
    //               << " infected agents from detailed infection data\n";
    //     const size_t amount_infected_non_detailed =
    //         static_cast<size_t>(ensemble_amount_of_infections[run_idx][0].get_value(
    //             ensemble_amount_of_infections[run_idx][0].get_num_time_points() - 1)(0));
    //     std::cout << "Run " << run_idx + 1 << ": " << amount_infected_non_detailed
    //               << " infected agents from non-detailed infection data\n";
    // }

    return mio::success();
}

void assign_infection_state(mio::abm::World& world, const std::vector<uint32_t>& initial_infections,
                            mio::abm::TimePoint t)
{
    for (const auto& person_id : initial_infections) {
        auto rng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), world.get_person(person_id));
        world.get_person(person_id).add_new_infection(
            mio::abm::Infection(rng, mio::abm::VirusVariant::Alpha, world.get_person(person_id).get_age(),
                                world.parameters, t, mio::abm::InfectionState::InfectedNoSymptoms));
    }

    // std::cout << "Household map created with " << household_map.size() << " entries." << std::endl;
    // for (const auto& [person_id_pre, person_id_sim] : household_map) {
    //     std::cout << "Person ID Pre: " << person_id_pre << " Person ID Sim: " << person_id_sim << " HouseholdID: "
    //               << city.get_person(person_id_sim).get_assigned_location_index(mio::abm::LocationType::Home)
    //               << "House has "
    //               << city
    //                      .get_individualized_location(mio::abm::LocationId{
    //                          city.get_person(person_id_sim).get_assigned_location_index(mio::abm::LocationType::Home),
    //                          mio::abm::LocationType::Home})
    //                      .get_persons()
    //                      .size()
    //               << " persons." << std::endl;
    // }

    // print initial infections / also print household and workplace id
    // std::cout << "Initial infections assigned to persons: \n";
    // for (const auto& person_id : initial_infections) {
    //     std::cout << person_id << " ";
    //     auto& person = world.get_person(person_id);
    //     std::cout << "(Household ID: " << person.get_assigned_location_index(mio::abm::LocationType::Home)
    //               << ", Work ID: " << person.get_assigned_location_index(mio::abm::LocationType::Work) << ") ";
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;
}

mio::IOResult<SimulationResults>
MultiRunSimulator::run_single_simulation_with_infections(mio::abm::World&& base_world,
                                                         const std::vector<uint32_t>& initial_infections,
                                                         double k_parameter, int simulation_days)
{
    // TODO: Implement single simulation run with initial infections
    // 1. Apply initial infections to specified people/households
    // 2. Set infection parameter K
    // 3. Run simulation for specified days
    // 4. Return results
    SimulationResults results;
    auto t0   = mio::abm::TimePoint(0); // Start time per simulation
    auto tmax = mio::abm::TimePoint(0) + mio::abm::days(simulation_days); // End time per simulation

    // Apply initial infections
    assign_infection_state(base_world, initial_infections, t0);
    
    // Calculate average household size of initially infected persons
    double total_household_size = 0.0;
    double total_persons_above_age_group_0 = 0.0;
    for (const auto& person_id : initial_infections) {
        auto& person = base_world.get_person(person_id);
        auto household_index = person.get_assigned_location_index(mio::abm::LocationType::Home);
        auto& household = base_world.get_individualized_location(
            mio::abm::LocationId{household_index, mio::abm::LocationType::Home});
        total_household_size += household.get_persons().size();
        
        // Count persons with age group > 0
        for (const auto& household_person_id : household.get_persons()) {
            auto household_person = household_person_id->get_age();
            if (household_person != mio::AgeGroup(0)) {
                total_persons_above_age_group_0 += 1.0;
            }
        }
    }
    double average_household_size = initial_infections.empty() ? 0.0 : total_household_size / initial_infections.size();
    double average_persons_above_age_group_0 = initial_infections.empty() ? 0.0 : total_persons_above_age_group_0 / initial_infections.size();
    
    // Set infection parameter K
    base_world.parameters.get<mio::abm::InfectionRateFromViralShed>() = k_parameter;

    auto sim = mio::abm::Simulation(t0, std::move(base_world));

    // DEBUG
    mio::History<mio::abm::TimeSeriesWriter, LogInfectionPerLocationTypePerAgeGroup> historyInfectionPerLocationType{
        Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups())};
    mio::History<mio::abm::TimeSeriesWriter, LogInfectionStatePerAgeGroup> historyInfectionStatePerAgeGroup{
        Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups())};
    mio::History<mio::abm::DataWriterToMemory, LogPersonIdAbndLocationTyoe> historyLocationTypeAndId;

    // FIGURES
    mio::History<mio::abm::TimeSeriesWriter, LogAmountOfInfections> historyAmountInfected{
        Eigen::Index(1)}; // 1st figure
    mio::History<mio::abm::DataWriterToMemory, LogHouseholdId, LogLocationIdAndPersonId, LogWhoInfected, LogWorkId>
        historyContactNetwork; // 2nd figure
    mio::History<mio::abm::DataWriterToMemory, LogInfectionDetailed> historyInfectionDetailed; // 3rd figure
    mio::History<mio::abm::DataWriterToMemory, LogLocationIdAndType> historyLocationIdAndType; // 4th figure
    mio::History<mio::abm::DataWriterToMemory, LogIsInfected> historyIsInfected; // 5th figure

    sim.advance(tmax, historyInfectionPerLocationType, historyInfectionStatePerAgeGroup, historyAmountInfected,
                historyContactNetwork, historyInfectionDetailed, historyLocationIdAndType, historyLocationTypeAndId,
                historyIsInfected);

    // DEBUG
    results.infection_per_loc_type =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionPerLocationType.get_log())};
    results.infection_state_per_age_group =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionStatePerAgeGroup.get_log())};
    results.ensemble_params                    = std::vector<mio::abm::World>{sim.get_world()};
    results.ensemble_params_no_agegroups       = std::vector<mio::abm::World>{mio::abm::World(1)};
    results.infection_per_location_type_and_id = std::get<0>(historyLocationTypeAndId.get_log());

    // FIGURES
    results.history_infected_amount =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyAmountInfected.get_log())};
    results.history_household_id =
        std::vector<std::tuple<uint32_t, uint32_t>>{std::get<0>(historyContactNetwork.get_log())[0]};
    results.history_work_id =
        std::vector<std::tuple<uint32_t, uint32_t>>{std::get<3>(historyContactNetwork.get_log())[0]};
    results.history_person_and_location_id =
        std::vector<std::vector<std::tuple<uint32_t, uint32_t>>>{std::get<1>(historyContactNetwork.get_log())};
    results.history_detailed_infection =
        std::vector<std::vector<std::tuple<uint32_t, uint32_t, mio::abm::LocationType>>>{
            std::get<0>(historyInfectionDetailed.get_log())};
    results.loc_id_and_type =
        std::vector<std::tuple<uint32_t, mio::abm::LocationType>>{std::get<0>(historyLocationIdAndType.get_log())[0]};
    results.history_infected_status =
        std::vector<std::vector<std::tuple<uint32_t, bool>>>{std::get<0>(historyIsInfected.get_log())};
    results.average_household_size_of_initial_infections = average_household_size;
    results.average_persons_above_age_group_0_in_initial_households = average_persons_above_age_group_0;

    // Placeholder implementation
    return mio::success(results);
}
