#include "abm/abm.h"
#include "abm/household.h"
#include <algorithm>
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

void write_contacts(std::ostream& out, const mio::abm::World& world,
                    const std::vector<LocationPerPersonLogger::Type>& loclog)
{
    out << "# Hour, PersonId, ContactId1, Intensity1, ContactId2, Intensity2, ...\n";
    // std::cout << loclog.size() << "\n";
    // std::vector<std::pair<int, uint32_t>> stays(world.get_persons().size(), {0, 0}); // t0, loc
    // for (size_t p = 0; p < stays.size(); p++) {
    //     stays[p].second = loclog[0].second[p];
    //     stays[p].first  = loclog[0].first.seconds();
    // }
    for (size_t t = 1; t < loclog.size(); t++) {
        for (size_t p = 0; p < loclog[t].second.size(); p++) {
            std::string preamble(std::to_string(loclog[t].first.hours()) + ", " + std::to_string(p));
            // write contacts
            assert(loclog[t].second[p] != mio::abm::INVALID_LOCATION_INDEX);
            const auto& loc = world.get_locations()[loclog[t].second[p]];
            for (size_t contact = 0; contact < loc.get_assigned_persons().size(); contact++) {
                auto p_loc = std::find(loc.get_assigned_persons().begin(), loc.get_assigned_persons().end(), p);
                const auto& matrices = loc.get_contact_matrices();
                const auto& t_matrix = matrices[loclog[t].first.hour_of_day()];
                assert(p_loc != loc.get_assigned_persons().end());
                const auto intensity = t_matrix(p_loc - loc.get_assigned_persons().begin(), contact);
                if (intensity > 0) {
                    out << preamble << ", " << loc.get_assigned_persons()[contact] << ", " << intensity;
                    preamble = "";
                }
            }
            // // reset stay
            // stays[p].first  = loclog[t].first.seconds();
            // stays[p].second = loclog[t].second[p];

            if (preamble.empty()) {
                out << "\n";
            }
        }
    }
}

void write_contacts_flat(std::ostream& out, const mio::abm::World& world,
                         const std::vector<LocationPerPersonLogger::Type>& loclog)
{
    out << "# Hour, PersonId1, PersonId2, Intensity\n";
    for (size_t t = 1; t < loclog.size(); t++) {
        for (size_t p = 0; p < loclog[t].second.size(); p++) {
            // write contacts
            assert(loclog[t].second[p] != mio::abm::INVALID_LOCATION_INDEX);
            const auto& loc = world.get_locations()[loclog[t].second[p]];
            for (size_t contact = 0; contact < loc.get_assigned_persons().size(); contact++) {
                auto p_loc = std::find(loc.get_assigned_persons().begin(), loc.get_assigned_persons().end(), p);
                const auto& matrices = loc.get_contact_matrices();
                const auto& t_matrix = matrices[loclog[t].first.hour_of_day()];
                assert(p_loc != loc.get_assigned_persons().end());
                const auto intensity = t_matrix(p_loc - loc.get_assigned_persons().begin(), contact);
                if (intensity > 0) {
                    out << loclog[t].first.hours() << ", " << p << ", " << loc.get_assigned_persons()[contact] << ", "
                        << intensity << "\n";
                }
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

void assign_contact_matrix(mio::abm::World& world, mio::abm::LocationId location, std::string filename)
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
    auto school_class1 = world.add_location(mio::abm::LocationType::School);
    world.get_individualized_location(school_class1).get_infection_parameters().set<mio::abm::MaximumContacts>(60);
    auto school_class2 = world.add_location(mio::abm::LocationType::School);
    world.get_individualized_location(school_class2).get_infection_parameters().set<mio::abm::MaximumContacts>(60);
    auto school_class_other = world.add_location(mio::abm::LocationType::School);
    world.get_individualized_location(school_class_other).get_infection_parameters().set<mio::abm::MaximumContacts>(60);

    //3. Workplaces
    auto work1 = world.add_location(mio::abm::LocationType::Work);
    world.get_individualized_location(work1).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto work2 = world.add_location(mio::abm::LocationType::Work);
    world.get_individualized_location(work2).get_infection_parameters().set<mio::abm::MaximumContacts>(15);
    auto work3 = world.add_location(mio::abm::LocationType::Work);
    world.get_individualized_location(work3).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto work4 = world.add_location(mio::abm::LocationType::Work);
    world.get_individualized_location(work4).get_infection_parameters().set<mio::abm::MaximumContacts>(100);
    auto work_other = world.add_location(mio::abm::LocationType::Work);
    world.get_individualized_location(work_other).get_infection_parameters().set<mio::abm::MaximumContacts>(10000);

    //4. Social Events
    // Maximum contacs limit the number of people that a person can infect while being at this location.
    auto event = world.add_location(mio::abm::LocationType::SocialEvent);
    world.get_individualized_location(event).get_infection_parameters().set<mio::abm::MaximumContacts>(100);

    //5. Shops
    auto shop1 = world.add_location(mio::abm::LocationType::BasicsShop);
    world.get_individualized_location(shop1).get_infection_parameters().set<mio::abm::MaximumContacts>(30);
    auto shop2 = world.add_location(mio::abm::LocationType::BasicsShop);
    world.get_individualized_location(shop2).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto shop3 = world.add_location(mio::abm::LocationType::BasicsShop);
    world.get_individualized_location(shop3).get_infection_parameters().set<mio::abm::MaximumContacts>(10000);

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

    // assign people to the right workplace, randomly to work1, work2 or work3  and to the right school
    int number_of_persons_max_to_work_1     = 20;
    int current_number_of_persons_to_work_1 = 0;

    int number_of_persons_max_to_work_2     = 15;
    int current_number_of_persons_to_work_2 = 0;

    int number_of_persons_max_to_work_3     = 5;
    int current_number_of_persons_to_work_3 = 0;

    int number_of_persons_max_to_work_4     = 100;
    int current_number_of_persons_to_work_4 = 0;

    int number_of_persons_max_to_highschool     = 60;
    int current_number_of_persons_to_highschool = 0;

    int number_of_persons_max_to_primary_school     = 60;
    int current_number_of_persons_to_primary_school = 0;

    int number_of_persons_max_to_supermarket_1     = 5;
    int current_number_of_persons_to_supermarket_1 = 0;

    int number_of_persons_max_to_supermarket_2     = 30;
    int current_number_of_persons_to_supermarket_2 = 0;

    auto rng = world.get_rng();
    for (auto& person : world.get_persons()) {
        //to the 2 different schools
        if (person.get_age() == age_group_5_to_24) {
            std::vector<double> school_distribution{1.0, 1.0};
            auto school = mio::DiscreteDistribution<size_t>::get_instance()(rng, school_distribution);
            if (school == 0 && current_number_of_persons_to_primary_school < number_of_persons_max_to_primary_school) {
                person.set_assigned_location(school_class1);
                current_number_of_persons_to_primary_school++;
            }
            else if (school == 1 && current_number_of_persons_to_highschool < number_of_persons_max_to_highschool) {
                person.set_assigned_location(school_class2);
                current_number_of_persons_to_highschool++;
            }
            else {
                person.set_assigned_location(school_class_other);
            }
        }
        //to the 3 different workplaces
        if (person.get_age() == age_group_25_to_64) {
            std::vector<double> work_distribution{20.0, 15.0, 5.0, 100.0};
            auto work = mio::DiscreteDistribution<size_t>::get_instance()(rng, work_distribution);
            if (work == 0 && current_number_of_persons_to_work_1 < number_of_persons_max_to_work_1) {
                person.set_assigned_location(work1);
                current_number_of_persons_to_work_1++;
            }
            else if (work == 1 && current_number_of_persons_to_work_2 < number_of_persons_max_to_work_2) {
                person.set_assigned_location(work2);
                current_number_of_persons_to_work_2++;
            }
            else if (work == 2 && current_number_of_persons_to_work_3 < number_of_persons_max_to_work_3) {
                person.set_assigned_location(work3);
                current_number_of_persons_to_work_3++;
            }
            else if (work == 3 && current_number_of_persons_to_work_4 < number_of_persons_max_to_work_4) {
                person.set_assigned_location(work4);
                current_number_of_persons_to_work_4++;
            }
            else {
                person.set_assigned_location(work_other);
            }
        }
        //to the 2 different shops
        std::vector<double> shop_distribution{30.0, 5.0};
        auto shop = mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), shop_distribution);
        if (shop == 0 && current_number_of_persons_to_supermarket_1 < number_of_persons_max_to_supermarket_1) {
            person.set_assigned_location(shop1);
            current_number_of_persons_to_supermarket_1++;
        }
        else if (shop == 1 && current_number_of_persons_to_supermarket_2 < number_of_persons_max_to_supermarket_2) {
            person.set_assigned_location(shop2);
            current_number_of_persons_to_supermarket_2++;
        }
        else {
            person.set_assigned_location(shop3);
        }
    }

    std::string contacts_path =
        "/Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/data/contacts/microcontacts/24h_networks_csv";

    for (auto& location : world.get_locations()) {
        if (location.get_type() == mio::abm::LocationType::Home) {
            assign_home_contact_matrix(world, location);
        }
    }

    assign_contact_matrix(world, work1, mio::path_join(contacts_path, "office_20_20.csv"));
    assign_contact_matrix(world, work2, mio::path_join(contacts_path, "office_15_15.csv"));
    assign_contact_matrix(world, work3, mio::path_join(contacts_path, "office_5_5.csv"));
    assign_contact_matrix(world, work4, mio::path_join(contacts_path, "office_100_100.csv"));

    // assign_contact_matrix(world, school_class1, mio::path_join(contacts_path, "highschool_60_60.csv"));
    // assign_contact_matrix(world, school_class2, mio::path_join(contacts_path, "primaryschool_60_60.csv"));

    assign_contact_matrix(world, shop1, mio::path_join(contacts_path, "supermarked_5_5.csv"));
    assign_contact_matrix(world, shop2, mio::path_join(contacts_path, "supermarked_30_30.csv"));

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
    }

    {
        std::ofstream outfile("micro_abm.txt");
        std::get<0>(historyTimeSeries.get_log())
            .print_table({"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4, outfile);
        std::cout << "Results written to micro_abm.txt" << std::endl;
    }
    return 0;
}
