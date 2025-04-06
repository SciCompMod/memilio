#include "abm/abm.h"
#include "abm/household.h"
#include <algorithm>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <ostream>
#include <string>
#include <iostream>
#include <utility>
#include <vector>
#include "abm/common_abm_loggers.h"
#include "abm/location.h"
#include "abm/location_type.h"
#include "abm/parameters.h"
#include "abm/person.h"
#include "abm/simulation.h"
#include "abm/time.h"
#include "abm/world.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/io/history.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"

// #define ADD_OVERFLOW_LOCATIONS

constexpr std::string_view loc_type_name(mio::abm::LocationType t)
{
    using mio::abm::LocationType;
    switch (t) {
    case LocationType::Home:
        return "Home";
    case LocationType::BasicsShop:
        return "BasicsShop";
    case LocationType::Car:
        return "Car";
    case LocationType::Cemetery:
        return "Cemetery";
    case LocationType::Hospital:
        return "Hospital";
    case LocationType::ICU:
        return "ICU";
    case LocationType::PublicTransport:
        return "PublicTransport";
    case LocationType::TransportWithoutContact:
        return "TransportWithoutContact";
    case LocationType::School:
        return "School";
    case LocationType::Work:
        return "Work";
    case LocationType::SocialEvent:
        return "SocialEvent";
    case LocationType::Count:
        return "InvalidValueCount";
    default:
        return "InvalidValue";
    }
}

const std::string_view age_group_name(const mio::AgeGroup& a)
{
    switch (a.get()) {
    case 0:
        return "age_group_0_to_4";
    case 1:
        return "age_group_5_to_24";
    case 2:
        return "age_group_25_to_64";
    case 3:
        return "age_group_65_plus";
    default:
        return "InvalidValue";
    }
}

struct LocationPerPersonLogger : public mio::LogAlways {
    using Type = std::tuple<mio::abm::TimePoint, std::vector<uint32_t>, std::vector<uint32_t>>;

    Type log(const mio::abm::Simulation& sim)
    {
        location_per_person.resize(sim.get_world().get_persons().size());
        std::transform(sim.get_world().get_persons().begin(), sim.get_world().get_persons().end(),
                       location_per_person.begin(), [](const mio::abm::Person& p) {
                           return p.get_location().get_index();
                       });
        local_person_index.resize(sim.get_world().get_persons().size());
        std::transform(sim.get_world().get_persons().begin(), sim.get_world().get_persons().end(),
                       local_person_index.begin(), [](const mio::abm::Person& p) {
                           const auto& loc = p.get_location();
                           if (loc.is_contact_matrix_assignment_dynamic()) {
                               const auto p_loc =
                                   std::find_if(loc.get_persons().begin(), loc.get_persons().end(), [&p](auto&& val) {
                                       return val->get_person_id() == p.get_person_id();
                                   });
                               assert(p_loc != loc.get_persons().end() && "person is not assigned to location!");
                               return std::distance(loc.get_persons().begin(), p_loc);
                           }
                           else {
                               const auto p_loc = std::find(loc.get_assigned_persons().begin(),
                                                            loc.get_assigned_persons().end(), p.get_person_id());
                               assert(p_loc != loc.get_assigned_persons().end() &&
                                      "person is not assigned to location!");
                               return std::distance(loc.get_assigned_persons().begin(), p_loc);
                           }
                       });
        return std::tuple(sim.get_time(), location_per_person, local_person_index);
    }

    std::vector<uint32_t> location_per_person;
    std::vector<uint32_t> local_person_index;
};

struct LocationSizesLogger : public mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, std::vector<uint32_t>>;

    Type log(const mio::abm::Simulation& sim)
    {
        persons_per_location.resize(sim.get_world().get_locations().size());
        std::transform(sim.get_world().get_locations().begin(), sim.get_world().get_locations().end(),
                       persons_per_location.begin(), [](const mio::abm::Location& l) {
                           return l.get_number_persons();
                       });
        return std::pair(sim.get_time() - mio::abm::hours(1), persons_per_location);
    }

    std::vector<uint32_t> persons_per_location;
};

void write_contacts_flat(std::ostream& out, const mio::abm::World& world,
                         const std::vector<LocationPerPersonLogger::Type>& loclog)
{
    std::vector<std::vector<uint32_t>> assigned_persons_per_location(world.get_locations().size());
    for (auto& person : world.get_persons()) {
        for (auto& assigned_loc : person.get_assigned_locations()) {
            if (assigned_loc != mio::abm::INVALID_LOCATION_INDEX) {
                assigned_persons_per_location[assigned_loc].push_back(person.get_person_id());
            }
        }
    }

    std::bitset<(uint32_t)mio::abm::LocationType::Count> missing_loc_types;
    out << "# Hour, PersonId1, PersonId2, Intensity, LocationId, LocationType\n";
    for (size_t t = 1; t < loclog.size(); t++) { // iterate over time points
        const auto& log_time            = std::get<0>(loclog[t]);
        const auto& location_per_person = std::get<1>(loclog[t]);
        const auto& local_person_index  = std::get<2>(loclog[t]);
        for (size_t p = 0; p < location_per_person.size(); p++) { // iterate over population
            assert(location_per_person[p] != mio::abm::INVALID_LOCATION_INDEX);
            const auto& loc = world.get_locations()[location_per_person[p]]; // p's location at time t
            // skip locations without assigments
            if (loc.get_assigned_persons().size() == 0) {
                missing_loc_types[(uint32_t)loc.get_type()] = true;
                continue;
            }
            // write out all contacts with a positive intensity
            for (auto contact_id : assigned_persons_per_location[loc.get_index()]) {
                if (contact_id == mio::abm::INVALID_PERSON_ID || contact_id == p) {
                    // skip loops (contact with oneself) and filler entries (invalid ids)
                    continue;
                }
                if (location_per_person[contact_id] != loc.get_index()) {
                    // skip false contacts, i.e. ones that are not at p's location
                    continue;
                }
                const auto& t_matrix = loc.get_contact_matrices()[log_time.hour_of_day()];
                const auto intensity = t_matrix(local_person_index[p], local_person_index[contact_id]);
                if (intensity > 0) {
                    out << log_time.hours() << ", " << p << ", " << contact_id << ", " << intensity << ", "
                        << loc.get_index() << ", " << loc_type_name(loc.get_type()) << "\n";
                }
            }
        }
    }
    if (missing_loc_types.any()) {
        for (uint32_t i = 0; i < (uint32_t)mio::abm::LocationType::Count; i++) {
            if (missing_loc_types[i]) {
                mio::log_warning("Missing contact matrix for a Location of type: {}",
                                 loc_type_name((mio::abm::LocationType)i));
            }
        }
    }
}

mio::IOResult<mio::abm::HourlyContactMatrix> read_hourly_contact_matrix(const std::string& csv_file)
{
    std::ifstream file(csv_file);
    if (!file.good()) {
        return mio::failure(mio::StatusCode::FileNotFound, "Could not open " + csv_file + ".");
    }

    std::string reader;

    std::vector<double> matrix_entries;
    int current_hour = 0;
    int t, row, col;
    double val;

    mio::abm::HourlyContactMatrix hcm;

    // skip first line
    std::getline(file, reader);
    // possible EOF here
    // read csv
    while (std::getline(file, reader)) {
        const auto fmt_string = "%i,%i,%i,%lf\n";
        const int status      = sscanf(reader.c_str(), fmt_string, &t, &row, &col, &val);

        if (status != 4) {
            return mio::failure(mio::StatusCode::InvalidFileFormat, "Unexpected format while reading " + csv_file +
                                                                        ". Expexted \"" + fmt_string +
                                                                        "\". Line reads \"" + reader + "\"\n");
        }

        // collect matrix_entries until t hits the next hour
        if (t > current_hour) {
            const double exact_sqrt = std::sqrt((double)matrix_entries.size());
            if (std::ceil(exact_sqrt) != std::floor(exact_sqrt)) {
                return mio::failure(mio::StatusCode::InvalidFileFormat, "Unexpected format while reading " + csv_file +
                                                                            ". Matrix at t=" + std::to_string(t) +
                                                                            " is not square.\n");
            }
            const size_t n = (size_t)exact_sqrt;
            if (current_hour > 1 && (size_t)hcm[current_hour - 1].cols() != n) {
                return mio::failure(mio::StatusCode::InvalidFileFormat,
                                    "Unexpected format while reading " + csv_file + ". Matrix at t=" +
                                        std::to_string(current_hour) + " differs in size from prior matrix.\n");
            }
            hcm[current_hour] = Eigen::MatrixXd(n, n);
            for (size_t i = 0; i < matrix_entries.size(); i++) {
                hcm[current_hour].data()[i] = matrix_entries[i];
            }
            matrix_entries.clear();
            ++current_hour;
        }

        matrix_entries.push_back(val);
    }

    const size_t n    = std::round(std::sqrt(matrix_entries.size()));
    hcm[current_hour] = Eigen::MatrixXd(n, n);
    for (size_t i = 0; i < matrix_entries.size(); i++) {
        hcm[current_hour].data()[i] = matrix_entries[i];
    }

    return mio::success(hcm);
}

void assign_contact_matrix(mio::abm::World& world, mio::abm::LocationId location, std::string filename,
                           unsigned location_capacity, bool dynamically_assigned)
{
    auto res = read_hourly_contact_matrix(filename);
    if (!res) {
        std::cout << res.error().formatted_message();
        exit(1);
    }

    std::vector<uint32_t> assigned_persons_for_location;
    for (auto& person : world.get_persons()) {
        if (person.get_assigned_locations()[(uint32_t)location.type] == location.index) {
            assigned_persons_for_location.push_back(person.get_person_id());
        }
    }
    assert((dynamically_assigned || (location_capacity >= assigned_persons_for_location.size())) &&
           "Insufficiently large contact network for this location.");
    if (location_capacity > assigned_persons_for_location.size()) {
        mio::log_info("Fewer persons assigned to location {} than expected: Expected number assigned ({}) {} network "
                      "capacity ({}). Unassigned network nodes are ignored.",
                      location.index, dynamically_assigned ? ">=" : "==", location_capacity,
                      assigned_persons_for_location.size());
    }
    // fill the vector with invalid ids to avoid dealing with differing vector sizes during computation.
    assigned_persons_for_location.resize(location_capacity, mio::abm::INVALID_PERSON_ID);

    world.get_individualized_location(location).assign_contact_matrices(res.value(), assigned_persons_for_location,
                                                                        dynamically_assigned);
}

void assign_random_contact_matrix(mio::abm::World& world, mio::abm::LocationId location, std::string base_path,
                                  std::string type_of_location, unsigned network_capacity,
                                  unsigned time_step_in_minutes, bool dynamically_assigned)
{
    // find contact network directory
    const std::filesystem::path file_dir =
        mio::path_join(base_path, type_of_location, "num_nodes=" + std::to_string(network_capacity),
                       "time=" + std::to_string(time_step_in_minutes));

    if (!std::filesystem::exists(file_dir)) {
        mio::log(mio::LogLevel::critical, "Cannot find directory \"{}\"", file_dir.string());
        exit(1);
    }

    // count files so we can draw a random network
    const int num_files =
        (int)std::distance(std::filesystem::directory_iterator(file_dir), std::filesystem::directory_iterator());

    if (num_files == 0) {
        mio::log(mio::LogLevel::critical, "No files in directory \"{}\"", file_dir.string());
        exit(1);
    }

    // make random choice
    int file_choice = mio::UniformIntDistribution<int>::get_instance()(world.get_rng(), 0, num_files - 1);

    auto file_itr = std::filesystem::directory_iterator(file_dir);
    for (int i = 0; i < file_choice; i++) {
        ++file_itr;
    }

    const std::string file_name = file_itr->path().string();

    mio::log_info("Picked file {} for location {} of type {}.", file_name, location.index,
                  loc_type_name(location.type));

    // assign contact network
    assign_contact_matrix(world, location, file_name, network_capacity, dynamically_assigned);

    // make checks and apply fixes:
    // if dynamically_assigned is false, we assume the contact matrix columns and assigned persons should be paired 1:1.
    // otherwise, we allow any number of assigned persons, give up the 1:1 relation in favor of a ad-hoc assignemt of
    // columns to the available persons at any time step, and set the capacity of the location to restricted_capacity.
    if (dynamically_assigned) {
        auto& loc = world.get_individualized_location(location);
        loc.set_capacity(network_capacity, 0);
    }
}

const static Eigen::MatrixXd& full_home_contact_matrix()
{
    // data used for abm paper
    // age groups (paper): 0_to_4, 5_to_14, 15_to_34, 35_to_59, 60_to_79, 80_plus
    const static Eigen::MatrixXd contact_weights = []() {
        Eigen::MatrixXd m(6, 6);
        m << 0.4413, 0.4504, 1.2383, 0.8033, 0.0494, 0.0017, //
            0.0485, 0.7616, 0.6532, 1.1614, 0.0256, 0.0013, //
            0.1800, 0.1795, 0.8806, 0.6413, 0.0429, 0.0032, //
            0.0495, 0.2639, 0.5189, 0.8277, 0.0679, 0.0014, //
            0.0087, 0.0394, 0.1417, 0.3834, 0.7064, 0.0447, //
            0.0292, 0.0648, 0.1248, 0.4179, 0.3497, 0.1544;
        return m;
    }();

    // age weights translate from paper age groups to the ones we use here
    // we want to sum over contact weights proportional to age range overlap
    // age groups: 0_to_4, 5_to_24, 25_to_64, 65_plus
    const static Eigen::MatrixXd age_weights = []() {
        Eigen::MatrixXd m(4, 6);
        m << 1, 0, 0, 0, 0, 0, //
            0, 1, 0.5, 0, 0, 0, //
            0, 0, 0.5, 1, 0.25, 0, //
            0, 0, 0, 0, 0.75, 1;
        return m;
    }();

    // translate from paper age groups to the ones we use here
    const static Eigen::MatrixXd full_contact_matrix = []() {
        Eigen::MatrixXd m = age_weights * contact_weights * age_weights.transpose();
        // normalize so the weights can be used as proportions of an hour
        m.array() /= m.maxCoeff();
        return m;
    }();

    return full_contact_matrix;
}

void assign_home_contact_matrix(mio::abm::World& world, mio::abm::Location& home)
{
    // find out who lives here
    std::vector<uint32_t> assigned_persons;
    for (auto& person : world.get_persons()) {
        if (person.get_assigned_locations()[(uint32_t)home.get_type()] == home.get_index()) {
            assigned_persons.push_back(person.get_person_id());
        }
    }

    const auto noise = []() {
        return 1; // no noise
        // return mio::UniformDistribution<double>::get_instance()(mio::thread_local_rng(), 0.8, 1.0);
    };

    // create contact matrix for the location
    Eigen::MatrixXd home_contact_matrix(assigned_persons.size(), assigned_persons.size());
    for (Eigen::Index i = 0; i < home_contact_matrix.rows(); i++) {
        for (Eigen::Index j = 0; j < home_contact_matrix.cols(); j++) {
            const auto age_i          = (size_t)world.get_persons()[assigned_persons[i]].get_age();
            const auto age_j          = (size_t)world.get_persons()[assigned_persons[j]].get_age();
            home_contact_matrix(i, j) = noise() * full_home_contact_matrix()(age_i, age_j);
        }
    }

    // set a contact matrix for every hour
    mio::abm::HourlyContactMatrix hcm;
    for (auto& matrix : hcm) {
        matrix = home_contact_matrix;
    }
    // add it to the location
    home.assign_contact_matrices(hcm, assigned_persons, false);
}

void create_sampled_world(mio::abm::World& world, const fs::path& input_dir, const mio::abm::TimePoint& t0,
    int max_num_persons, mio::Date start_date_sim)
{
    //Set global infection parameters (similar to infection parameters in SECIR model) and initialize the world

    set_parameters(world.parameters);
    set_local_parameters(world);
    restart_timer(timer, "time taken for setting up parameters and local parameters");

    // Create the world object from statistical data.
    create_world_from_data(world, (input_dir / "mobility/braunschweig_result_ffa8_modified2.csv").generic_string(),
        max_num_persons);
    world.use_migration_rules(false);
    restart_timer(timer, "time taken for braunschweig trip input");

    // Assign an infection state to each person.
    assign_infection_state(world, t0);
    restart_timer(timer, "time taken for assigning infection state");

    // Assign vaccination status to each person.
    assign_vaccination_state(world, start_date_sim);
    restart_timer(timer, "time taken for assigning vaccination state");

}


int main()
{
    mio::set_log_level(mio::LogLevel::warn);

    size_t num_age_groups         = 4;
    const auto age_group_0_to_4   = mio::AgeGroup(0);
    const auto age_group_5_to_24  = mio::AgeGroup(1);
    const auto age_group_25_to_64 = mio::AgeGroup(2);
    const auto age_group_65_plus  = mio::AgeGroup(3);

    const unsigned time_step_in_minutes = 60;
    std::string contacts_path =
        // "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/contacts/microcontacts/24h_networks_csv";
        "/home/schm_r6/Documents/memilio_contact_networks";
    // "/localdata2/dial_mo/Graphs/24h_networks_csv";

    std::cout << "Base home contact matrix:\n" << full_home_contact_matrix() << "\n";

    // Create the world with 4 age groups.
    auto world = mio::abm::World(num_age_groups);

    world.parameters.get<mio::abm::AgeGroupGotoSchool>()                        = false;
    world.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_24]     = true;
    world.parameters.get<mio::abm::AgeGroupGotoWork>()                          = false;
    world.parameters.get<mio::abm::AgeGroupGotoWork>()[age_group_25_to_64]      = true;

    world.parameters.get<mio::abm::BasicShoppingRate>()[age_group_0_to_4]   = 0.1;
    world.parameters.get<mio::abm::BasicShoppingRate>()[age_group_5_to_24]  = 0.5;
    world.parameters.get<mio::abm::BasicShoppingRate>()[age_group_25_to_64] = 0.5;
    world.parameters.get<mio::abm::BasicShoppingRate>()[age_group_65_plus]  = 0.3;

    world.parameters.get<mio::abm::SocialEventRate>().get_baseline().array() = 0.5;

    // Setup the world

    //1. Households
    int n_households = 500;
    auto child       = mio::abm::HouseholdMember(num_age_groups); // A child is 25/75% 0-4 or 5-24.
    child.set_age_weight(age_group_0_to_4, 1);
    child.set_age_weight(age_group_5_to_24, 3);

    auto parent = mio::abm::HouseholdMember(num_age_groups); // A parent is 100% 25-64.
    parent.set_age_weight(age_group_25_to_64, 1);

    auto grandparent = mio::abm::HouseholdMember(num_age_groups); // A grandparent is 100% 65+.
    grandparent.set_age_weight(age_group_65_plus, 1);

    // One-person household
    auto onePersonHousehold_group = mio::abm::HouseholdGroup();
    auto onePersonHousehold_full  = mio::abm::Household();
    onePersonHousehold_full.add_members(parent, 1);
    // onePersonHousehold_full.add_members(grandparent, 1);
    onePersonHousehold_group.add_households(onePersonHousehold_full, n_households * 0.4);

    // Two-person household with two parents.
    auto twoPersonHousehold_group = mio::abm::HouseholdGroup();
    auto twoPersonHousehold_full  = mio::abm::Household();
    twoPersonHousehold_full.add_members(child, 1);
    twoPersonHousehold_full.add_members(parent, 1);
    twoPersonHousehold_group.add_households(twoPersonHousehold_full, n_households * 0.3);

    // Three-person household with two parents and one child.
    auto threePersonHousehold_group = mio::abm::HouseholdGroup();
    auto threePersonHousehold_full  = mio::abm::Household();
    threePersonHousehold_full.add_members(child, 1);
    threePersonHousehold_full.add_members(parent, 2);
    threePersonHousehold_group.add_households(threePersonHousehold_full, n_households * 0.2);

    // Four-person household with two parents and two children.
    auto fourPersonHousehold_group = mio::abm::HouseholdGroup();
    auto fourPersonHousehold_full  = mio::abm::Household();
    fourPersonHousehold_full.add_members(child, 2);
    fourPersonHousehold_full.add_members(parent, 2);
    fourPersonHousehold_group.add_households(fourPersonHousehold_full, n_households * 0.1);

    add_household_group_to_world(world, onePersonHousehold_group);
    add_household_group_to_world(world, twoPersonHousehold_group);
    add_household_group_to_world(world, threePersonHousehold_group);
    add_household_group_to_world(world, fourPersonHousehold_group);

    //2. Schools
    std::vector<mio::abm::LocationId> schools;
    const std::vector<int> school_network_capacities{100, 200};
    for (auto size : school_network_capacities) {
        schools.push_back(world.add_location(mio::abm::LocationType::School));
        world.get_individualized_location(schools.back())
            .get_infection_parameters()
            .set<mio::abm::MaximumContacts>(size);
    }

    //3. Workplaces
    std::vector<mio::abm::LocationId> works;
    const std::vector<int> work_network_capacities{10, 10, 10, 10, 10, 10, 10, 10, 10,  10,  20,
                                                   20, 20, 20, 30, 30, 30, 50, 50, 100, 100, 100};
    for (auto size : work_network_capacities) {
        works.push_back(world.add_location(mio::abm::LocationType::Work));
        world.get_individualized_location(works.back()).get_infection_parameters().set<mio::abm::MaximumContacts>(size);
    }

    //4. Social Events
    std::vector<mio::abm::LocationId> social_events;
    // note that the required event capacity varies, since children may be in agegroup 0, which is excluded from events
    const std::vector<int> social_event_network_capacities    = {90, 90, 90, 90, 90, 15, 30};
    const std::vector<int> social_event_assignment_capacities = {180, 180, 180, 180, 180, 15, 30};

    for (auto size : social_event_network_capacities) {
        social_events.push_back(world.add_location(mio::abm::LocationType::SocialEvent));
        world.get_individualized_location(social_events.back())
            .get_infection_parameters()
            .set<mio::abm::MaximumContacts>(size);
    }

    //5. Shops
    std::vector<mio::abm::LocationId> shops;
    const std::vector<int> shop_network_capacities{20, 20, 20, 20, 20, 40, 40, 40, 40, 40, 40,
                                                   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};
    const std::vector<int> shop_assignment_capacities = {40, 40, 40, 40, 40, 80, 80, 80, 80, 80, 80,
                                                         80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80};
    for (auto size : shop_network_capacities) {
        shops.push_back(world.add_location(mio::abm::LocationType::BasicsShop));
        world.get_individualized_location(shops.back()).get_infection_parameters().set<mio::abm::MaximumContacts>(size);
    }

    //6. Hospitals
    auto hospital = world.add_location(mio::abm::LocationType::Hospital);
    world.get_individualized_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    // Add ICU with 10 maximum contacts.
    auto icu = world.add_location(mio::abm::LocationType::ICU);
    world.get_individualized_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(10);

    // Infection distribution remains unchanged
    std::vector<double> infection_distribution{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (auto& person : world.get_persons()) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), infection_distribution));
        auto rng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         world.parameters, mio::abm::TimePoint(0), infection_state));
        }
    }

    // Assign shared locations to the people
    for (auto& person : world.get_persons()) {
        //assign hospital and ICU
        person.set_assigned_location(hospital);
        person.set_assigned_location(icu);
    }

#ifdef ADD_OVERFLOW_LOCATIONS
    mio::log_warning("Using overflow locations without contact matrix!");

    auto social_event_other = world.add_location(mio::abm::LocationType::SocialEvent);
    auto school_other       = world.add_location(mio::abm::LocationType::School);
    auto work_other         = world.add_location(mio::abm::LocationType::Work);
    auto shop_other         = world.add_location(mio::abm::LocationType::BasicsShop);

    std::cout << "\"other\" location ids:\n"
              << "  " << loc_type_name(social_event_other.type) << ": " << social_event_other.index << "\n"
              << "  " << loc_type_name(school_other.type) << ": " << school_other.index << "\n"
              << "  " << loc_type_name(work_other.type) << ": " << work_other.index << "\n"
              << "  " << loc_type_name(shop_other.type) << ": " << shop_other.index << "\n";

    // Assign people to a randomized locations
    int num_assignees_school_other       = 0;
    int num_assignees_work_other         = 0;
    int num_assignees_shop_other         = 0;
    int num_assignees_social_event_other = 0;
#endif

    std::vector<double> school_dist{school_network_capacities.begin(), school_network_capacities.end()};
    std::vector<double> work_dist{work_network_capacities.begin(), work_network_capacities.end()};
    std::vector<double> shop_dist{shop_assignment_capacities.begin(), shop_assignment_capacities.end()};
    std::vector<double> social_event_dist{social_event_assignment_capacities.begin(),
                                          social_event_assignment_capacities.end()};

    auto rng = world.get_rng();
    for (auto& person : world.get_persons()) {
        //to the schools
        if (person.get_age() == age_group_5_to_24) {
            if (std::accumulate(school_dist.begin(), school_dist.end(), 0) > 0) {
                auto school = mio::DiscreteDistribution<size_t>::get_instance()(rng, school_dist);
                person.set_assigned_location(schools[school]);
                --school_dist[school];
            }
            else {
#ifdef ADD_OVERFLOW_LOCATIONS
                person.set_assigned_location(school_other);
                ++num_assignees_school_other;
#else
                mio::log(mio::LogLevel::critical, "School capacity too low.");
                return 1;
#endif
            }
        }
        //to the workplaces
        if (person.get_age() == age_group_25_to_64) {
            if (std::accumulate(work_dist.begin(), work_dist.end(), 0) > 0) {
                auto work = mio::DiscreteDistribution<size_t>::get_instance()(rng, work_dist);
                person.set_assigned_location(works[work]);
                --work_dist[work];
            }
            else {
#ifdef ADD_OVERFLOW_LOCATIONS
                person.set_assigned_location(work_other);
                ++num_assignees_work_other;
#else
                mio::log(mio::LogLevel::critical, "Work capacity too low.");
                return 1;
#endif
            }
        }

        //to the shops
        if (std::accumulate(shop_dist.begin(), shop_dist.end(), 0) > 0) {
            auto shop = mio::DiscreteDistribution<size_t>::get_instance()(rng, shop_dist);
            person.set_assigned_location(shops[shop]);
            --shop_dist[shop];
        }
        else {
#ifdef ADD_OVERFLOW_LOCATIONS
            person.set_assigned_location(shop_other);
            ++num_assignees_shop_other;
#else
            mio::log(mio::LogLevel::critical, "Shop capacity too low.");
            return 1;
#endif
        }

        //to social events
        if (person.get_age() != age_group_0_to_4) {
            if (std::accumulate(social_event_dist.begin(), social_event_dist.end(), 0) > 0) {
                auto event = mio::DiscreteDistribution<size_t>::get_instance()(rng, social_event_dist);
                person.set_assigned_location(social_events[event]);
                --social_event_dist[event];
            }
            else {
#ifdef ADD_OVERFLOW_LOCATIONS
                person.set_assigned_location(social_event_other);
                ++num_assignees_social_event_other;
#else
                mio::log(mio::LogLevel::critical, "Social event capacity too low.");
                return 1;
#endif
            }
        }
    }

#ifdef ADD_OVERFLOW_LOCATIONS
    std::cout << "Assignees in school_other: " << num_assignees_school_other << "\n"
              << "Assignees in work_other: " << num_assignees_work_other << "\n"
              << "Assignees in shop_other: " << num_assignees_shop_other << "\n"
              << "Assignees in social_event_other: " << num_assignees_social_event_other << "\n";
#endif

    // Load contact matrices
    for (auto& location : world.get_locations()) {
        if (location.get_type() == mio::abm::LocationType::Home) {
            assign_home_contact_matrix(world, location);
        }
    }
    for (size_t i = 0; i < works.size(); i++) {
        assign_random_contact_matrix(world, works[i], contacts_path, "office", work_network_capacities[i],
                                     time_step_in_minutes, false);
    }

    for (size_t i = 0; i < schools.size(); i++) {
        assign_random_contact_matrix(world, schools[i], contacts_path, "highschool", school_network_capacities[i],
                                     time_step_in_minutes, false);
    }

    for (size_t i = 0; i < shops.size(); i++) {
        assign_random_contact_matrix(world, shops[i], contacts_path, "supermarked", shop_network_capacities[i],
                                     time_step_in_minutes, true);
    }

    for (size_t i = 0; i < social_events.size(); i++) {
        assign_random_contact_matrix(world, social_events[i], contacts_path, "gallery",
                                     social_event_network_capacities[i], time_step_in_minutes, true);
    }

    // Run the simulatio
    auto t0   = mio::abm::TimePoint(0);
    auto tmax = t0 + mio::abm::days(10);
    world.parameters.check_constraints();
    auto sim = mio::abm::Simulation(t0, std::move(world));
    // Create a history object to store the time series of the infection states.
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)};
    mio::History<mio::DataWriterToMemory, LocationPerPersonLogger, LocationSizesLogger> loclog;
    // Run the simulation until tmax with the history object.
    sim.advance(tmax, historyTimeSeries, loclog);

    {
        std::ofstream contactfile("micro_abm_contacts.csv");
        write_contacts_flat(contactfile, sim.get_world(), std::get<0>(loclog.get_log()));
        std::cout << "Contacts written to micro_abm_contacts.csv" << std::endl;
    }

    {
        std::ofstream outfile("debug_location-sizes.csv");
        for (auto& loc_sizes : std::get<1>(loclog.get_log())) {
            outfile << "t=" << loc_sizes.first.hours() << ": ";
            auto max = std::max_element(loc_sizes.second.begin(), loc_sizes.second.end());
            // print out maximum as well as ID and type
            outfile << "max=" << *max
                    << " argmax id=" << std::distance(loc_sizes.second.begin(), loc_sizes.second.begin())
                    << " argmax type="
                    << loc_type_name(
                           sim.get_world().get_locations()[std::distance(loc_sizes.second.begin(), max)].get_type())
                    << " all=";
            // print out all locations
            for (auto& loc_size : loc_sizes.second) {
                outfile << loc_size << " ";
            }
            outfile << "\n";
        }
        std::cout << "Debug: Current Location population written to debug_location-sizes.csv\n";
    }

    {
        std::ofstream outfile("metadata_location-props.csv");
        outfile << "# LocationID, LocationType (Value), LocationType (Name)\n";
        for (auto& location : sim.get_world().get_locations()) {
            outfile << location.get_index() << ", " << (uint32_t)location.get_type() << ", "
                    << loc_type_name(location.get_type()) << "\n";
        }
        std::cout << "Metadata: Location info written to metadata_location-props.csv\n";
    }

    {
        std::ofstream outfile("metadata_person-props.csv");
        outfile << "# PersonID, AgeGroup (Value), AgeGroup (Name), Assigned LocationID (Home), Assigned LocationID "
                   "(School), Assigned LocationID (Work), Assigned LocationID (SocialEvent), Assigned LocationID "
                   "(BasicsShop), Assigned LocationID (Hospital), Assigned LocationID (ICU), Assigned LocationID "
                   "(Car), Assigned LocationID (PublicTransport), Assigned LocationID (TransportWithoutContact), "
                   "Assigned LocationID (Cemetery)\n";
        for (auto& person : sim.get_world().get_persons()) {
            outfile << person.get_person_id() << ", " << person.get_age().get() << ", "
                    << age_group_name(person.get_age());
            for (auto id : person.get_assigned_locations()) {
                outfile << ", " << id;
            }
            outfile << "\n";
        }
        std::cout << "Metadata: Person info written to metadata_person-props.csv\n";
    }

    {
        std::ofstream outfile("micro_abm.txt");
        std::get<0>(historyTimeSeries.get_log())
            .print_table({"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4, outfile);
        std::cout << "Results written to micro_abm.txt" << std::endl;
    }
    return 0;
}
