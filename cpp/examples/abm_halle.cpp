#include "abm/lockdown_rules.h"
#include "abm/model.h"
#include "abm/common_abm_loggers.h"
#include "memilio/io/directories.h"

#include "H5Tpublic.h"
#include "H5public.h"
#include "abm/infection_state.h"
#include "abm/location_id.h"
#include "abm/location_type.h"
#include "abm/time.h"
#include "memilio/io/io.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"
// #include "pybind_util.h"
// #include "utils/custom_index_array.h"
// #include "utils/parameter_set.h"
// #include "utils/index.h"
#include "abm/simulation.h"
#include "abm/household.h"
#include "abm/personal_rng.h"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "memilio/io/hdf5_cpp.h"
// #include "munich_postprocessing/output_processing.h"


// #include "pybind11/attr.h"
// #include "pybind11/cast.h"
// #include "pybind11/pybind11.h"
// #include "pybind11/operators.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <filesystem>
#include <iostream>

#include <cstring> // for strerror
#include <cerrno>

#include <chrono>


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

void write_mapping_to_file(std::string filename, std::map<int, std::vector<std::string>>& mapping)
{
    auto file = fopen(filename.c_str(), "w");
    if (file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", filename);
    }
    else {
        for (auto it = mapping.begin(); it != mapping.end(); it++) {
            fprintf(file, "%d", it->first);
            for (auto s : it->second) {
                fprintf(file, " %s", s.c_str());
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }
}

mio::AgeGroup determine_age_group(uint32_t age)
{
    if (age <= 4) {
        return mio::AgeGroup(0);
    }
    else if (age <= 6) {
        return mio::AgeGroup(1);
    }
    else if (age <= 15) {
        return mio::AgeGroup(2);
    }
    else if (age <= 18) {
        return mio::AgeGroup(3);
    }
    else if (age <= 21) {
        return mio::AgeGroup(4);
    }
    else if (age <= 35) {
        return mio::AgeGroup(5);
    }
    else if (age <= 60) {
        return mio::AgeGroup(6);
    }
    else if (age <= 65) {
        return mio::AgeGroup(7);
    }
    else if (age <= 75) {
        return mio::AgeGroup(8);
    }
    else if (age <= 80) {
        return mio::AgeGroup(9);
    }
    else {
        return mio::AgeGroup(10);
    }
}

void initialize_model(mio::abm::Model& model, std::string person_file,
    size_t max_work_size, size_t max_school_size){
// void initialize_model(mio::abm::Model& model, std::string person_file, std::string outfile){
    // Mapping of ABM locations to traffic areas/cells
    // - each traffic area is mapped to a vector containing strings with LocationType and LocationId
    std::map<int, std::vector<std::string>> loc_area_mapping;
    // Mapping of traffic data location ids to ABM location ids
    std::map<int, mio::abm::LocationId> home_locations;
    std::map<int, mio::abm::LocationId> shop_locations;
    std::map<int, mio::abm::LocationId> event_locations;
    std::map<int, mio::abm::LocationId> school_locations;
    std::map<int, mio::abm::LocationId> work_locations;
    std::vector<mio::abm::LocationId> hospitals;
    std::vector<mio::abm::LocationId> icus;
    std::map<std::pair<mio::abm::LocationType, mio::abm::LocationId>, mio::abm::LocationId> hosp_to_icu;
    std::vector<double> hospital_weights;
    std::vector<double> icu_weights;
    // Mapping of assigned agents to school and work locations
    std::map<int, std::map<mio::abm::LocationId, size_t>> school_sizes;
    std::map<int, std::map<mio::abm::LocationId, size_t>> work_sizes;



    // File pointer
    std::fstream fin;
    // Open an existing file
    fin.open(person_file, std::ios::in);

    if (!fin.is_open()) {
        std::cerr << "Error: could not open file " << person_file << "\n";
        return ;  // or handle error appropriately
    }
    else{std::cout << "file open\n";}

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
        // splitline(line, &row);
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        int beds          = row[index["beds"]];
        int icu_beds      = row[index["icu_beds"]];
        int hospital_zone = row[index["hospital_zone"]];
        auto hosp         = model.add_location(mio::abm::LocationType::Hospital);
        hospitals.push_back(hosp);
        hospital_weights.push_back(beds);
        std::string loc =
            "0" + std::to_string(static_cast<int>(mio::abm::LocationType::Hospital)) + std::to_string(hosp.get());
        auto zone_iter = loc_area_mapping.find(hospital_zone);
        if (zone_iter == loc_area_mapping.end()) {
            loc_area_mapping.insert({hospital_zone, {loc}});
        }
        else {
            loc_area_mapping[hospital_zone].push_back(loc);
        }
        // Add icu if there is one
        if (icu_beds > 0) {
            auto icu = model.add_location(mio::abm::LocationType::ICU);
            icus.push_back(icu);
            icu_weights.push_back(icu_beds);
            hosp_to_icu.insert({std::make_pair(mio::abm::LocationType::Hospital, hosp), icu});
            std::string loc_icu =
                "0" + std::to_string(static_cast<int>(mio::abm::LocationType::ICU)) + std::to_string(icu.get());
            zone_iter = loc_area_mapping.find(hospital_zone);
            if (zone_iter == loc_area_mapping.end()) {
                loc_area_mapping.insert({hospital_zone, {loc_icu}});
            }
            else {
                loc_area_mapping[hospital_zone].push_back(loc_icu);
            }
        }

        uint32_t age = row[index["age"]];

        int home_id   = row[index["home_id"]];
        int home_zone = row[index["home_zone"]];

        mio::abm::LocationId home;

        auto iter_home = home_locations.find(home_id);
        // check whether home location already exists in model
        if (iter_home == home_locations.end()) {
            // if home location does not exists yet, create new location and insert it to mapping
            home = model.add_location(mio::abm::LocationType::Home);
            home_locations.insert({home_id, home});
            std::string loc_home =
                "0" + std::to_string(static_cast<int>(mio::abm::LocationType::Home)) + std::to_string(home.get());
            auto zone_iter_home = loc_area_mapping.find(home_zone);
            if (zone_iter_home == loc_area_mapping.end()) {
                loc_area_mapping.insert({home_zone, {loc_home}});
            }
            else {
                loc_area_mapping[home_zone].push_back(loc_home);
            }
        }
        else {
            home = home_locations[home_id];
        }
        // Add person to model and assign home location to it
        auto pid     = model.add_person(home, determine_age_group(age));
        auto& person = model.get_person(pid);
        person.set_assigned_location(mio::abm::LocationType::Home, home,model.get_id());

        int shop_id   = row[index["shop_id"]];
        int shop_zone = row[index["shop_zone"]];

        mio::abm::LocationId shop;

        auto iter_shop = shop_locations.find(shop_id);
        // Check whether shop location already exists in model
        if (iter_shop == shop_locations.end()) {
            // Create shop location and add it to mapping
            shop = model.add_location(mio::abm::LocationType::BasicsShop);
            // Shops with ids -1 are individual locations each
            if (shop_id != -1) {
                shop_locations.insert({shop_id, shop});
            }
            std::string loc_shop =
                "0" + std::to_string(static_cast<int>(mio::abm::LocationType::BasicsShop)) + std::to_string(shop.get());
            auto zone_iter_shop = loc_area_mapping.find(shop_zone);
            if (zone_iter_shop == loc_area_mapping.end()) {
                loc_area_mapping.insert({shop_zone, {loc_shop}});
            }
            else {
                loc_area_mapping[shop_zone].push_back(loc_shop);
            }
        }
        else {
            shop = shop_locations[shop_id];
        }
        // Assign shop to person
        person.set_assigned_location(mio::abm::LocationType::BasicsShop, shop,model.get_id());
        // model.get_location(shop).increase_size();

        int event_id   = row[index["event_id"]];
        int event_zone = row[index["event_zone"]];

        mio::abm::LocationId event;

        auto iter_event = event_locations.find(event_id);
        // Check whether event location already exists in model
        if (iter_event == event_locations.end()) {
            //Create event location and add it to mapping
            event = model.add_location(mio::abm::LocationType::SocialEvent);
            // Events with id -1 are individual locations each
            if (event_id != -1) {
                event_locations.insert({event_id, event});
            }
            std::string loc_event = "0" + std::to_string(static_cast<int>(mio::abm::LocationType::SocialEvent)) +
                              std::to_string(event.get());
            auto zone_iter_event = loc_area_mapping.find(event_zone);
            if (zone_iter_event == loc_area_mapping.end()) {
                loc_area_mapping.insert({event_zone, {loc_event}});
            }
            else {
                loc_area_mapping[event_zone].push_back(loc_event);
            }
        }
        else {
            event = event_locations[event_id];
        }
        // Assign event location to person
        person.set_assigned_location(mio::abm::LocationType::SocialEvent, event,model.get_id());
        // model.get_location(event).increase_size();

        // Check if person is school-aged
        if (person.get_age() == mio::AgeGroup(1) || person.get_age() == mio::AgeGroup(2) || person.get_age() == mio::AgeGroup(3)) {
            int school_id   = row[index["school_id"]];
            int school_zone = row[index["school_zone"]];

            mio::abm::LocationId school;

            auto iter_school = school_locations.find(school_id);
            // Check whether school location is already in model
            if (iter_school == school_locations.end()) {
                // Add schools locations to model and insert it in mapping
                school = model.add_location(mio::abm::LocationType::School);
                // schools with id -1 are individual locations each
                if (school_id != -1) {
                    school_locations.insert({school_id, school});
                    // Add school to size map
                    school_sizes[school_id].insert({school, 1});
                }
                std::string loc_school = "0" + std::to_string(static_cast<int>(mio::abm::LocationType::School)) +
                                  std::to_string(school.get());
                auto zone_iter_school = loc_area_mapping.find(school_zone);
                if (zone_iter_school == loc_area_mapping.end()) {
                    loc_area_mapping.insert({school_zone, {loc_school}});
                }
                else {
                    loc_area_mapping[school_zone].push_back(loc_school);
                }
            }
            else {
                school = school_locations[school_id];
                if (school_sizes[school_id][school] == max_school_size) {
                    // Check if a new school has to be open or if there still is a school that has capacity
                    bool found = false;
                    for (auto const& [key, val] : school_sizes[school_id]) {
                        if (val < max_school_size) {
                            found = true;
                            school_sizes[school_id][key] += 1;
                            school = key;
                            break;
                        }
                    }
                    if (!found) {
                        // Create new school
                        school = model.add_location(mio::abm::LocationType::School);
                        school_locations.insert({school_id, school});
                        // Add school to size map
                        school_sizes[school_id].insert({school, 1});
                        std::string loc_school = "0" + std::to_string(static_cast<int>(mio::abm::LocationType::School)) +
                                          std::to_string(school.get());
                        auto zone_iter_school = loc_area_mapping.find(school_zone);
                        if (zone_iter_school == loc_area_mapping.end()) {
                            loc_area_mapping.insert({school_zone, {loc_school}});
                        }
                        else {
                            loc_area_mapping[school_zone].push_back(loc_school);
                        }
                    }
                }
                else {
                    school_sizes[school_id][school] += 1;
                }
            }
            // Assign school location to person
            person.set_assigned_location(mio::abm::LocationType::School, school,model.get_id());
            // model.get_location(school).increase_size();
        }
        // Check if person is work-aged
        if (person.get_age() == mio::AgeGroup(4) || person.get_age() == mio::AgeGroup(5) || person.get_age() == mio::AgeGroup(6) || person.get_age() == mio::AgeGroup(7)) {
            int work_id   = row[index["work_id"]];
            int work_zone = row[index["work_zone"]];

            if (work_zone == -2) {
                mio::log_error("Person with id {} has work age but no work zone", row[index["puid"]]);
            }

            mio::abm::LocationId work;

            auto iter_work = work_locations.find(work_id);
            // Check whether work location already exists in model
            if (iter_work == work_locations.end()) {
                // Add work location to model and insert it in mapping
                work = model.add_location(mio::abm::LocationType::Work);
                // Locations with id -1 are individual locations each
                if (work_id != -1) {
                    work_locations.insert({work_id, work});
                    // Add work to size map
                    work_sizes[work_id].insert({work, 1});
                }
                std::string loc_id =
                    "0" + std::to_string(static_cast<int>(mio::abm::LocationType::Work)) + std::to_string(work.get());
                auto zone_iter_id = loc_area_mapping.find(work_zone);
                if (zone_iter_id == loc_area_mapping.end()) {
                    loc_area_mapping.insert({work_zone, {loc_id}});
                }
                else {
                    loc_area_mapping[work_zone].push_back(loc_id);
                }
            }
            else {
                work = work_locations[work_id];
                if (work_sizes[work_id][work] == max_work_size) {
                    // Check if a new work has to be opened or if there still is a school that has capacity
                    bool found = false;
                    for (auto const& [key, val] : work_sizes[work_id]) {
                        if (val < max_work_size) {
                            found = true;
                            work_sizes[work_id][key] += 1;
                            work = key;
                            break;
                        }
                    }
                    if (!found) {
                        // Create new work
                        work = model.add_location(mio::abm::LocationType::Work);
                        work_locations.insert({work_id, work});
                        // Add work to size map
                        work_sizes[work_id].insert({work, 1});
                        std::string loc_work = "0" + std::to_string(static_cast<int>(mio::abm::LocationType::Work)) +
                                          std::to_string(work.get());
                        auto zone_iter_work = loc_area_mapping.find(work_zone);
                        if (zone_iter_work == loc_area_mapping.end()) {
                            loc_area_mapping.insert({work_zone, {loc_work}});
                        }
                        else {
                            loc_area_mapping[work_zone].push_back(loc_work);
                        }
                    }
                }
                else {
                    work_sizes[work_id][work] += 1;
                }
            }
            // Assign work location to person
            person.set_assigned_location(mio::abm::LocationType::Work, work,model.get_id());
            // model.get_location(work).increase_size();
        }
        // Assign Hospital and ICU
        size_t hosp_int = mio::DiscreteDistribution<size_t>::get_instance()(model.get_rng(), hospital_weights);
        person.set_assigned_location(mio::abm::LocationType::Hospital, hospitals[hosp_int],model.get_id());
        // model.get_location(hospitals[hosp_int]).increase_size();
        if (hosp_to_icu.count(std::make_pair(mio::abm::LocationType::Hospital, hospitals[hosp_int])) > 0) {
            person.set_assigned_location(
                mio::abm::LocationType::ICU,
                hosp_to_icu[std::make_pair(mio::abm::LocationType::Hospital, hospitals[hosp_int])],model.get_id());
            // model.get_location(hosp_to_icu[std::make_pair(mio::abm::LocationType::Hospital, hospitals[hosp_int])])
            //     .increase_size();
        }
        else {
            size_t icu = mio::DiscreteDistribution<size_t>::get_instance()(model.get_rng(), icu_weights);
            person.set_assigned_location(mio::abm::LocationType::ICU, icus[icu],model.get_id());
            // model.get_location(icus[icu]).increase_size();
        }

        // write_mapping_to_file(outfile, loc_area_mapping);
    }
}
int main()
{
    std::cout << "Start" << "\n";

    auto start = std::chrono::high_resolution_clock::now();
    // ------------------------------------------------
    // ------------ Ages, Fam, Child/Adult ------------
    // ------------------------------------------------

    mio::set_log_level(mio::LogLevel::warn);
    size_t num_age_groups            = 11;
    // const auto age_group_0_to_3      = mio::AgeGroup(0);
    // const auto age_group_4_to_6      = mio::AgeGroup(1);
    // const auto age_group_7_to_15     = mio::AgeGroup(2);
    // const auto age_group_16_to_18    = mio::AgeGroup(3);
    // const auto age_group_19_to_21    = mio::AgeGroup(4);
    // const auto age_group_22_to_35    = mio::AgeGroup(5);
    // const auto age_group_36_to_60    = mio::AgeGroup(6);
    // const auto age_group_61_to_65    = mio::AgeGroup(7);
    // const auto age_group_66_to_75    = mio::AgeGroup(8);
    // const auto age_group_76_to_80    = mio::AgeGroup(9);
    // const auto age_group_81_and_over = mio::AgeGroup(10);

    // Create the model with 11 age groups.
    std::string path = "/home/wulf_ka/home/abm/memilio/cpp/examples/df_abm_short.csv";
    std::string out  = "/home/wulf_ka/home/abm/memilio/cpp/examples/out";
    auto model  = mio::abm::Model(num_age_groups);
    // Set same infection parameter for all age groups. For example, the incubation period is log normally distributed with parameters 4 and 1.
    model.parameters.get<mio::abm::TimeExposedToNoSymptoms>() = mio::ParameterDistributionLogNormal(4., 1.);

    // Test init
    initialize_model(model, path, 50,50);

    auto end_init = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_init = end_init - start;

    std::cout << "Time for init: " << elapsed_init.count() << " seconds\n";

    //     // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    // model.parameters.get<mio::abm::AgeGroupGotoSchool>()                 = false;
    // model.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_0_to_3]   = true;
    // model.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_4_to_6]  = true;
    // model.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_7_to_15] = true;
    // model.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_16_to_18] = true;
    // // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 and 35-59)
    // model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple(
    //     {age_group_19_to_21, age_group_22_to_35, age_group_36_to_60, age_group_61_to_65}, true);

    // // zu Testzwecken weil sich das ABM sonst beschwert dass die nicht genutzt werden
    // model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple(
    //     {age_group_66_to_75, age_group_76_to_80, age_group_81_and_over}, true);

    // Check if the parameters satisfy their contraints.
    model.parameters.check_constraints();

    // ------------------------------------
    // ------------ Households ------------
    // ------------------------------------


    // --------------------------------
    // ------------ Events ------------
    // --------------------------------



    // -----------------------------------------
    // ------------ Infection Param ------------
    // -----------------------------------------

    // Increase aerosol transmission for all locations
    model.parameters.get<mio::abm::AerosolTransmissionRates>() = 5.0;
    // Increase contact rate for all people between 15 and 34 (i.e. people meet more often in the same location)
    // model.get_location(work).get_infection_parameters().get<mio::abm::ContactRates>().get_baseline()(
    //     age_group_19_to_21.get(), age_group_19_to_21.get()) = 10.0;

    // People can get tested at work (and do this with 0.5 probability) from time point 0 to day 10.
    auto validity_period       = mio::abm::days(1);
    auto probability           = 0.5;
    auto start_date            = mio::abm::TimePoint(0);
    auto end_date              = mio::abm::TimePoint(0) + mio::abm::days(10);
    auto test_type             = mio::abm::TestType::Antigen;
    auto test_parameters       = model.parameters.get<mio::abm::TestData>()[test_type];
    auto testing_criteria_work = mio::abm::TestingCriteria();
    auto testing_scheme_work = mio::abm::TestingScheme(testing_criteria_work, validity_period, start_date, end_date,
                                                        test_parameters, probability);
    model.get_testing_strategy().add_scheme(mio::abm::LocationType::Work, testing_scheme_work);

    // Assign infection state to each person.
    // The infection states are chosen randomly with the following distribution
    std::vector<ScalarType> infection_distribution{0.5, 0.3, 0.05, 0.05, 0.05, 0.05, 0.0, 0.0};
    for (auto& person : model.get_persons()) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), infection_distribution));
        auto rng = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                            model.parameters, start_date, infection_state));
        }
    }

    // ------------------------------------------
    // ------------ Events to People ------------
    // ------------------------------------------


    // Set start and end time for the simulation.
    auto t0   = mio::abm::TimePoint(0);
    auto tmax = t0 + mio::abm::days(10);
    auto sim  = mio::abm::Simulation(t0, std::move(model));

    // Create a history object to store the time series of the infection states.
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)};

    // Run the simulation until tmax with the history object.
    sim.advance(tmax, historyTimeSeries);

    // Write results to a file. Also print the filepath to make it easier to find
    auto outpath = mio::create_directories_or_exit(mio::example_results_dir("abm_minimal")) / "history.txt";
    std::ofstream outfile(outpath);
    std::get<0>(historyTimeSeries.get_log())
        .print_table(outfile, {"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4);
    std::cout << "Results written to " << outpath << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Time total: " << elapsed.count() << " seconds\n";

    return 0;
}
