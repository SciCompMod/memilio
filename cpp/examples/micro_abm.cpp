#include "abm/abm.h"
#include "abm/household.h"
#include <algorithm>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <ostream>
#include <string>
#include <iostream>
#include <utility>
#include <vector>
#include "abm/common_abm_loggers.h"
#include "abm/location.h"
#include "abm/location_type.h"
#include "abm/person.h"
#include "abm/simulation.h"
#include "abm/time.h"
#include "abm/world.h"
#include "memilio/io/history.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"

struct LocationPerPersonLogger : public mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, std::vector<uint32_t>>;

    Type log(const mio::abm::Simulation& sim)
    {
        location_per_person.resize(sim.get_world().get_persons().size());
        std::transform(sim.get_world().get_persons().begin(), sim.get_world().get_persons().end(),
                       location_per_person.begin(), [](const mio::abm::Person& p) {
                           return p.get_location().get_index();
                       });
        return std::pair(sim.get_time(), location_per_person);
    }

    std::vector<uint32_t> location_per_person;
};

void write_contacts_flat(std::ostream& out, const mio::abm::World& world,
                         const std::vector<LocationPerPersonLogger::Type>& loclog)
{
    std::bitset<(uint32_t)mio::abm::LocationType::Count> missing_loc_types;
    out << "# Hour, PersonId1, PersonId2, Intensity\n";
    for (size_t t = 1; t < loclog.size(); t++) { // iterate over time points
        for (size_t p = 0; p < loclog[t].second.size(); p++) { // iterate over population
            assert(loclog[t].second[p] != mio::abm::INVALID_LOCATION_INDEX);
            const auto& loc = world.get_locations()[loclog[t].second[p]]; // p's location at time t
            // skip locations without assigments
            if (loc.get_assigned_persons().size() == 0) {
                missing_loc_types[(uint32_t)loc.get_type()] = true;
                continue;
            }
            const auto p_loc = std::find(loc.get_assigned_persons().begin(), loc.get_assigned_persons().end(), p);
            assert(p_loc != loc.get_assigned_persons().end());
            // write out all contacts with a positive intensity
            for (size_t contact = 0; contact < loc.get_assigned_persons().size(); contact++) {
                if (loc.get_assigned_persons()[contact] == mio::abm::INVALID_PERSON_ID ||
                    loc.get_assigned_persons()[contact] == p) {
                    // skip loops (contact with oneself) and filler entries (invalid ids)
                    continue;
                }
                const auto& t_matrix = loc.get_contact_matrices()[loclog[t].first.hour_of_day()];
                const auto intensity = t_matrix(p_loc - loc.get_assigned_persons().begin(), contact);
                if (intensity > 0) {
                    out << loclog[t].first.hours() << ", " << p << ", " << loc.get_assigned_persons()[contact] << ", "
                        << intensity << "\n";
                }
            }
        }
    }
    if (missing_loc_types.any()) {
        for (uint32_t i = 0; i < (uint32_t)mio::abm::LocationType::Count; i++) {
            if (missing_loc_types[i]) {
                mio::log_warning("Missing contact matrix for a Location of type: {}", i);
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
        int status = sscanf(reader.c_str(), "%i,%i,%i,%lf\n", &t, &row, &col, &val);

        if (status != 4) {
            return mio::failure(mio::StatusCode::InvalidFileFormat,
                                "Unexpected format while reading " + csv_file + ". Line reads \"" + reader + "\"");
        }

        if (t > current_hour) {
            size_t n          = std::round(std::sqrt(matrix_entries.size()));
            hcm[current_hour] = Eigen::MatrixXd(n, n);
            for (size_t i = 0; i < matrix_entries.size(); i++) {
                hcm[current_hour].data()[i] = matrix_entries[i];
            }
            matrix_entries.clear();
            ++current_hour;
        }

        matrix_entries.push_back(val);
    }

    size_t n          = std::round(std::sqrt(matrix_entries.size()));
    hcm[current_hour] = Eigen::MatrixXd(n, n);
    for (size_t i = 0; i < matrix_entries.size(); i++) {
        hcm[current_hour].data()[i] = matrix_entries[i];
    }

    return mio::success(hcm);
}

void assign_contact_matrix(mio::abm::World& world, mio::abm::LocationId location, std::string filename,
                           unsigned location_capacity = 0)
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
    if (location_capacity > 0) {
        assert(location_capacity >= assigned_persons_for_location.size());
        assigned_persons_for_location.resize(location_capacity, mio::abm::INVALID_PERSON_ID);
    }

    world.get_individualized_location(location).assign_contact_matrices(res.value(), assigned_persons_for_location);
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
    home.assign_contact_matrices(hcm, assigned_persons);
}

int main()
{
    size_t num_age_groups         = 4;
    const auto age_group_0_to_4   = mio::AgeGroup(0);
    const auto age_group_5_to_24  = mio::AgeGroup(1);
    const auto age_group_25_to_64 = mio::AgeGroup(2);
    const auto age_group_65_plus  = mio::AgeGroup(3);

    std::cout << "Base home contact matrix:\n" << full_home_contact_matrix() << "\n";

    // Create the world with 4 age groups.
    auto world = mio::abm::World(num_age_groups);

    world.parameters.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    world.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_24] = true;
    world.parameters.get<mio::abm::AgeGroupGotoWork>()                      = false;
    world.parameters.get<mio::abm::AgeGroupGotoWork>()[age_group_25_to_64]  = true;

    // Setup the world

    //1. Households
    int n_households = 100; // we do 100 Households -> 50/30/15/5 split for 1/2/3/4 person households
    auto child       = mio::abm::HouseholdMember(num_age_groups); // A child is 50/50% 0-4 or 5-14.
    child.set_age_weight(age_group_0_to_4, 1);
    child.set_age_weight(age_group_5_to_24, 3);

    auto parent = mio::abm::HouseholdMember(num_age_groups); // A parent is 100% 25-64.
    parent.set_age_weight(age_group_25_to_64, 1);

    auto grandparent = mio::abm::HouseholdMember(num_age_groups); // A grandparent is 100% 65+.
    grandparent.set_age_weight(age_group_65_plus, 1);

    // One-person household
    auto onePersonHousehold_group = mio::abm::HouseholdGroup();
    auto onePersonHousehold_full  = mio::abm::Household();
    onePersonHousehold_full.add_members(parent, 3);
    onePersonHousehold_full.add_members(grandparent, 1);
    onePersonHousehold_group.add_households(onePersonHousehold_full, n_households * 0.5);

    // Two-person household with two parents.
    auto twoPersonHousehold_group = mio::abm::HouseholdGroup();
    auto twoPersonHousehold_full  = mio::abm::Household();
    twoPersonHousehold_full.add_members(child, 1);
    twoPersonHousehold_full.add_members(child, 1);
    twoPersonHousehold_group.add_households(twoPersonHousehold_full, n_households * 0.3);

    // Three-person household with two parents and one child.
    auto threePersonHousehold_group = mio::abm::HouseholdGroup();
    auto threePersonHousehold_full  = mio::abm::Household();
    threePersonHousehold_full.add_members(child, 1);
    threePersonHousehold_full.add_members(parent, 2);
    threePersonHousehold_group.add_households(threePersonHousehold_full, n_households * 0.15);

    // Four-person household with two parents and two children.
    auto fourPersonHousehold_group = mio::abm::HouseholdGroup();
    auto fourPersonHousehold_full  = mio::abm::Household();
    fourPersonHousehold_full.add_members(child, 2);
    fourPersonHousehold_full.add_members(parent, 2);
    fourPersonHousehold_group.add_households(fourPersonHousehold_full, n_households * 0.05);

    add_household_group_to_world(world, onePersonHousehold_group);
    add_household_group_to_world(world, twoPersonHousehold_group);
    add_household_group_to_world(world, threePersonHousehold_group);
    add_household_group_to_world(world, fourPersonHousehold_group);

    //2. Schools
    std::vector<mio::abm::LocationId> schools;
    const std::vector<int> school_sizes{60, 60};
    std::vector<std::string> school_files{"highschool_60_60.csv", "primaryschool_60_60.csv"};
    for (auto size : school_sizes) {
        schools.push_back(world.add_location(mio::abm::LocationType::School));
        world.get_individualized_location(schools.back())
            .get_infection_parameters()
            .set<mio::abm::MaximumContacts>(size);
    }
    // auto school_other = world.add_location(mio::abm::LocationType::School);
    // world.get_individualized_location(school_other).get_infection_parameters().set<mio::abm::MaximumContacts>(1000);

    //3. Workplaces
    std::vector<mio::abm::LocationId> works;
    const std::vector<int> work_sizes{20, 20, 15, 15, 15, 5, 100};
    std::vector<std::string> work_files{"office_20_20.csv", "office_20_20.csv", "office_15_15.csv",  "office_15_15.csv",
                                        "office_15_15.csv", "office_5_5.csv",   "office_100_100.csv"};
    for (auto size : work_sizes) {
        works.push_back(world.add_location(mio::abm::LocationType::Work));
        world.get_individualized_location(works.back()).get_infection_parameters().set<mio::abm::MaximumContacts>(size);
    }
    // auto work_other = world.add_location(mio::abm::LocationType::Work);
    // world.get_individualized_location(work_other).get_infection_parameters().set<mio::abm::MaximumContacts>(1000);

    //4. Social Events
    // Maximum contacs limit the number of people that a person can infect while being at this location.
    auto event = world.add_location(mio::abm::LocationType::SocialEvent);
    world.get_individualized_location(event).get_infection_parameters().set<mio::abm::MaximumContacts>(100);

    //5. Shops
    std::vector<mio::abm::LocationId> shops;
    const std::vector<int> shop_sizes{5, 5, 5, 5, 5, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30};
    std::vector<std::string> shop_files(5, "supermarked_5_5.csv");
    shop_files.resize(shop_sizes.size(), "supermarked_30_30.csv");
    for (auto size : shop_sizes) {
        shops.push_back(world.add_location(mio::abm::LocationType::BasicsShop));
        world.get_individualized_location(shops.back()).get_infection_parameters().set<mio::abm::MaximumContacts>(size);
    }
    // auto shop_other = world.add_location(mio::abm::LocationType::BasicsShop);
    // world.get_individualized_location(shop_other).get_infection_parameters().set<mio::abm::MaximumContacts>(1000);

    //6. Hospitals
    auto hospital = world.add_location(mio::abm::LocationType::Hospital);
    world.get_individualized_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    // Add ICU with 5 maximum contacs.
    auto icu = world.add_location(mio::abm::LocationType::ICU);
    world.get_individualized_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);

    std::vector<double> infection_distribution{0.5, 0.3, 0.05, 0.05, 0.05, 0.05, 0.0, 0.0};
    for (auto& person : world.get_persons()) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), infection_distribution));
        auto rng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         world.parameters, mio::abm::TimePoint(0), infection_state));
        }
    }

    // Assign locations to the people
    for (auto& person : world.get_persons()) {
        //assign shop and event
        person.set_assigned_location(event);
        //assign hospital and ICU
        person.set_assigned_location(hospital);
        person.set_assigned_location(icu);
        //assign work/school to people depending on their age
    }

    // Assign people to a random school or workplace, and a shop
    // int num_assignees_school_other = 0;
    // int num_assignees_work_other   = 0;
    // int num_assignees_shop_other   = 0;

    std::vector<double> school_dist{school_sizes.begin(), school_sizes.end()};
    std::vector<double> work_dist{work_sizes.begin(), work_sizes.end()};
    std::vector<double> shop_dist{shop_sizes.begin(), shop_sizes.end()};

    auto rng = world.get_rng();
    for (auto& person : world.get_persons()) {
        //to the schools
        if (person.get_age() == age_group_5_to_24) {
            if (std::accumulate(school_sizes.begin(), school_sizes.end(), 0) > 0) {
                auto school = mio::DiscreteDistribution<size_t>::get_instance()(rng, school_dist);
                person.set_assigned_location(schools[school]);
                --school_dist[school];
            }
            else {
                assert(false && "School capacity too low.");
                // person.set_assigned_location(school_other);
                // ++num_assignees_school_other;
            }
        }
        //to the workplaces
        if (person.get_age() == age_group_25_to_64) {
            if (std::accumulate(work_sizes.begin(), work_sizes.end(), 0) > 0) {
                auto work = mio::DiscreteDistribution<size_t>::get_instance()(rng, work_dist);
                person.set_assigned_location(works[work]);
                --work_dist[work];
            }
            else {
                assert(false && "Work capacity too low.");
                // person.set_assigned_location(work_other);
                // ++num_assignees_work_other;
            }
        }

        //to the shops
        if (std::accumulate(shop_sizes.begin(), shop_sizes.end(), 0) > 0) {
            auto shop = mio::DiscreteDistribution<size_t>::get_instance()(rng, shop_dist);
            person.set_assigned_location(shops[shop]);
            --shop_dist[shop];
        }
        else {
            assert(false && "Shop capacity too low.");
            // person.set_assigned_location(shop_other);
            // ++num_assignees_shop_other;
        }
    }
    // std::cout << "Assignees in school_other: " << num_assignees_school_other << "\n"
    //           << "Assignees in work_other: " << num_assignees_work_other << "\n"
    //           << "Assignees in shop_other: " << num_assignees_shop_other << "\n";

    // Load contact matrices
    std::string contacts_path =
        "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/contacts/microcontacts/24h_networks_csv";

    for (auto& location : world.get_locations()) {
        if (location.get_type() == mio::abm::LocationType::Home) {
            assign_home_contact_matrix(world, location);
        }
    }
    for (size_t i = 0; i < works.size(); i++) {
        assign_contact_matrix(world, works[i], mio::path_join(contacts_path, work_files.at(i)), work_sizes.at(i));
    }
    for (size_t i = 0; i < schools.size(); i++) {
        assign_contact_matrix(world, schools[i], mio::path_join(contacts_path, school_files.at(i)), school_sizes.at(i));
    }
    for (size_t i = 0; i < shops.size(); i++) {
        assign_contact_matrix(world, shops[i], mio::path_join(contacts_path, shop_files.at(i)), shop_sizes.at(i));
    }

    // Run the simulation
    auto t0   = mio::abm::TimePoint(0);
    auto tmax = t0 + mio::abm::days(10);
    world.parameters.check_constraints();
    auto sim = mio::abm::Simulation(t0, std::move(world));
    // Create a history object to store the time series of the infection states.
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)};
    mio::History<mio::DataWriterToMemory, LocationPerPersonLogger> loclog;
    // Run the simulation until tmax with the history object.
    sim.advance(tmax, historyTimeSeries, loclog);

    {
        std::ofstream contactfile("micro_abm_contacts.csv");
        write_contacts_flat(contactfile, sim.get_world(), std::get<0>(loclog.get_log()));
        std::cout << "Contacts written to micro_abm_contacts.csv" << std::endl;
    }

    {
        std::ofstream outfile("micro_abm.txt");
        std::get<0>(historyTimeSeries.get_log())
            .print_table({"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4, outfile);
        std::cout << "Results written to micro_abm.txt" << std::endl;
    }
    return 0;
}
