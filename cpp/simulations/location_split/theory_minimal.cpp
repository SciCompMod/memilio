/**
 * @file theory_minimal.cpp
 * @brief Minimal experiments with the MEmilio-ABM in the homogeneous setting of the note
 * "State transitions in the MEmilio-ABM" (Sec. 3.1), to show its one-step results in their pure form and
 * how a fixed group accumulates them over time.
 *
 * Setting: one age group, one contact per day between every pair of agents in a room (1 x 1 contact
 * matrix), identical constant viral shedding of every infectious agent, no mobility. Agents are placed in
 * rooms (Locations of type Home) and stay there unless an experiment moves them. Newly infected agents stay
 * Exposed for 10 days and do not shed, so the number i of infectious agents in a room is constant, as in
 * the propositions of the note. The infection coefficient is chosen such that a = Lambda * dt takes a
 * prescribed value, with dt = 1 hour.
 *
 * Experiments (every file has the closed form of the note next to the simulated mean):
 *  1. one_step_single.csv: one infectious agent with s susceptibles, one step (Eq. (9), Proposition 1).
 *  2. one_step_clustering.csv: x infectious and x susceptible agents in one room or in x rooms of one
 *     infectious and one susceptible agent, one step (Proposition 2).
 *  3. contact_time.csv: schools of 200 agents with 10 infectious agents, as one room, in fixed rooms or in
 *     rooms redrawn every 8 hours, for 160 hours. With i fixed, the infection hazards of a susceptible add up
 *     over the steps, so a fixed room of N agents gives sum_r s_r (1 - exp(-a i_r T / (N_r - 1))) after
 *     T steps: Eq. (9) with a replaced by a T.
 *
 * Usage: location_split_theory [output directory] [number of seeds]
 */
#include "abm/infection.h"
#include "abm/location_type.h"
#include "abm/model.h"
#include "abm/parameters.h"
#include "abm/personal_rng.h"
#include "memilio/io/io.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/parameter_distributions.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace
{
using mio::abm::InfectionState;
using mio::abm::LocationId;
using mio::abm::PersonId;

const auto virus     = mio::abm::VirusVariant::Wildtype;
const auto age       = mio::AgeGroup(0);
const auto t0        = mio::abm::TimePoint(0);
const auto dt        = mio::abm::hours(1);
constexpr double dt_days = 1.0 / 24.0;

/// @brief A Model with the homogeneous setting described in the file comment.
mio::abm::Model make_model(uint32_t seed)
{
    mio::abm::Model model(1);
    auto& params = model.parameters;
    // Newly infected agents stay Exposed (no shedding) for the whole experiment; infectious agents stay
    // infectious. The viral load jumps to its peak and stays there (decline 0.001 per day), and every agent
    // has the same shedding factor, so all infectious agents shed the same amount.
    params.get<mio::abm::TimeExposedToNoSymptoms>()           = mio::ParameterDistributionConstant(10.0);
    params.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()  = mio::ParameterDistributionConstant(1000.0);
    params.get<mio::abm::TimeInfectedNoSymptomsToRecovered>() = mio::ParameterDistributionConstant(1000.0);
    params.get<mio::abm::ViralLoadDistributions>()[{virus, age}] = mio::abm::ViralLoadDistributionsParameters{
        mio::AbstractParameterDistribution(mio::ParameterDistributionConstant(8.1)),
        mio::AbstractParameterDistribution(mio::ParameterDistributionConstant(1e4)),
        mio::AbstractParameterDistribution(mio::ParameterDistributionConstant(-0.001))};
    params.get<mio::abm::ViralShedFactor>()[{virus, age}] =
        mio::AbstractParameterDistribution(mio::ParameterDistributionConstant(1.0));
    params.get<mio::abm::AerosolTransmissionRates>() = 0.0;
    model.use_mobility_rules(false);
    model.reset_rng(std::vector<uint32_t>{seed, 0x7e0e7u, 2026u, 11u, 26u, 1u});
    return model;
}

LocationId add_room(mio::abm::Model& model)
{
    const auto room = model.add_location(mio::abm::LocationType::Home);
    model.get_location(room).get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() =
        Eigen::MatrixXd::Constant(1, 1, 1.0);
    return room;
}

PersonId add_agent(mio::abm::Model& model, LocationId room, bool infectious)
{
    const auto id = model.add_person(room, age);
    model.assign_location(id, room);
    if (infectious) {
        auto& person = model.get_person(id);
        auto rng     = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        person.add_new_infection(
            mio::abm::Infection(rng, virus, age, model.parameters, t0, InfectionState::InfectedNoSymptoms));
    }
    return id;
}

/// @brief Set the infection coefficient so that a = Lambda * dt, using the shedding of one infectious agent.
/// @return The infection coefficient k.
double set_a(mio::abm::Model& model, PersonId infectious, double a)
{
    const double shed = model.get_person(infectious).get_infection().get_viral_shed(t0 + dt / 2);
    // Per day: tau = k * contacts * shed * i / (N - 1); one step: a = k * contacts * shed * dt.
    const double k                                                     = a / (shed * dt_days);
    model.parameters.get<mio::abm::InfectionRateFromViralShed>()[{virus}] = k;
    return k;
}

/// @brief Susceptible, present and summed shedding per room at time t.
struct RoomState {
    int susceptible = 0, present = 0, infectious = 0;
    double shed     = 0.0;
};

std::map<uint32_t, RoomState> room_states(const mio::abm::Model& model, mio::abm::TimePoint t)
{
    std::map<uint32_t, RoomState> rooms;
    for (const auto& person : model.get_persons()) {
        auto& r = rooms[person.get_location().get()];
        r.present += 1;
        const auto state = person.get_infection_state(t);
        if (state == InfectionState::Susceptible) {
            r.susceptible += 1;
        }
        else if (state == InfectionState::InfectedNoSymptoms) {
            r.infectious += 1;
            r.shed += person.get_infection().get_viral_shed(t + dt / 2);
        }
    }
    return rooms;
}

/// @brief Expected new infections in the next step given the current state, as the ABM computes them.
double expected_new(const std::map<uint32_t, RoomState>& rooms, double k)
{
    double expected = 0.0;
    for (const auto& [id, r] : rooms) {
        if (r.susceptible > 0 && r.present > 1) {
            const double rate = k * r.shed / (r.present - 1); // per day, one contact per pair and day
            expected += r.susceptible * (1.0 - std::exp(-rate * dt_days));
        }
    }
    return expected;
}

struct Mean {
    double sum = 0, sum_sq = 0;
    int64_t n  = 0;
    void add(double x)
    {
        sum += x;
        sum_sq += x * x;
        ++n;
    }
    double mean() const
    {
        return n ? sum / n : 0.0;
    }
    double se() const
    {
        if (n < 2) {
            return 0.0;
        }
        const double m = mean();
        return std::sqrt(std::max(0.0, (sum_sq / n - m * m) / (n - 1)));
    }
};

// -------------------------------------------------------------------------------------------------------
// 1. One infectious agent with s susceptibles, one step.
// -------------------------------------------------------------------------------------------------------
void one_step_single(const std::string& dir, int seeds)
{
    std::ofstream out(dir + "/one_step_single.csv");
    out << "a,susceptible,rooms,observed_per_room,se,expected_abm,closed_form\n";
    for (double a : {0.05, 0.5, 1.0, 2.0}) {
        for (int s : {1, 2, 3, 5, 10, 20, 50, 100, 200}) {
            const int copies = std::max(50, 20000 / (s + 1));
            Mean per_room;
            double expected = 0.0;
            for (int seed = 0; seed < seeds; ++seed) {
                auto model = make_model(static_cast<uint32_t>(1000 * s + seed));
                PersonId first_infectious = PersonId::invalid_ID();
                for (int c = 0; c < copies; ++c) {
                    const auto room = add_room(model);
                    const auto id   = add_agent(model, room, true);
                    if (c == 0) {
                        first_infectious = id;
                    }
                    for (int j = 0; j < s; ++j) {
                        add_agent(model, room, false);
                    }
                }
                const double k = set_a(model, first_infectious, a);
                expected += expected_new(room_states(model, t0), k) / copies;
                model.evolve(t0, dt);
                std::map<uint32_t, int> new_per_room;
                for (const auto& person : model.get_persons()) {
                    const auto state = person.get_infection_state(t0 + dt);
                    new_per_room[person.get_location().get()] += state == InfectionState::Exposed;
                }
                for (const auto& [room, n] : new_per_room) {
                    per_room.add(n);
                }
            }
            const double closed = s * (1.0 - std::exp(-a / s));
            out << a << "," << s << "," << per_room.n << "," << per_room.mean() << "," << per_room.se() << ","
                << expected / seeds << "," << closed << "\n";
        }
    }
}

// -------------------------------------------------------------------------------------------------------
// 2. x infectious and x susceptible agents: one room (clustered) or x rooms of 1 + 1 (split), one step.
// -------------------------------------------------------------------------------------------------------
void one_step_clustering(const std::string& dir, int seeds)
{
    std::ofstream out(dir + "/one_step_clustering.csv");
    out << "a,x,arrangement,groups,observed_per_group,se,expected_abm,closed_form\n";
    for (double a : {0.05, 0.5, 1.0}) {
        for (int x : {1, 2, 3, 5, 10, 20}) {
            for (const std::string arrangement : {"clustered", "split"}) {
                const int copies = std::max(50, 20000 / (2 * x));
                Mean per_group;
                double expected = 0.0;
                for (int seed = 0; seed < seeds; ++seed) {
                    auto model = make_model(static_cast<uint32_t>(100000 + 1000 * x + seed));
                    PersonId first_infectious = PersonId::invalid_ID();
                    std::vector<std::vector<LocationId>> groups;
                    for (int c = 0; c < copies; ++c) {
                        std::vector<LocationId> rooms;
                        if (arrangement == "clustered") {
                            rooms.push_back(add_room(model));
                        }
                        for (int j = 0; j < x; ++j) {
                            if (arrangement == "split") {
                                rooms.push_back(add_room(model));
                            }
                            const auto id = add_agent(model, rooms.back(), true);
                            if (c == 0 && j == 0) {
                                first_infectious = id;
                            }
                            add_agent(model, rooms.back(), false);
                        }
                        groups.push_back(rooms);
                    }
                    const double k = set_a(model, first_infectious, a);
                    expected += expected_new(room_states(model, t0), k) / copies;
                    model.evolve(t0, dt);
                    std::map<uint32_t, int> new_per_room;
                    for (const auto& person : model.get_persons()) {
                        new_per_room[person.get_location().get()] +=
                            person.get_infection_state(t0 + dt) == InfectionState::Exposed;
                    }
                    for (const auto& rooms : groups) {
                        int n = 0;
                        for (const auto& room : rooms) {
                            n += new_per_room[room.get()];
                        }
                        per_group.add(n);
                    }
                }
                const double closed = arrangement == "clustered" ? x * (1.0 - std::exp(-a * x / (2.0 * x - 1.0)))
                                                                 : x * (1.0 - std::exp(-a));
                out << a << "," << x << "," << arrangement << "," << per_group.n << "," << per_group.mean() << ","
                    << per_group.se() << "," << expected / seeds << "," << closed << "\n";
            }
        }
    }
}

// -------------------------------------------------------------------------------------------------------
// 3. Schools of 200 agents with 10 infectious agents over 160 hours.
// -------------------------------------------------------------------------------------------------------
void contact_time(const std::string& dir, int seeds)
{
    constexpr int school_size = 200, infectious_per_school = 10, schools = 50, hours = 160;
    std::ofstream out(dir + "/contact_time.csv");
    out << "a,arrangement,room_size,hour,schools,observed_cumulative,se,expected_abm_cumulative,closed_form\n";
    struct Arm {
        std::string name;
        int room_size;
        int redraw_every; ///< 0: fixed rooms
    };
    const std::vector<Arm> arms{{"one room", 200, 0}, {"fixed", 50, 0}, {"fixed", 25, 0}, {"fixed", 10, 0},
                                {"fixed", 5, 0},      {"redrawn", 5, 8}};
    for (double a : {0.0475, 0.5}) {
        for (const auto& arm : arms) {
            std::vector<Mean> observed(hours + 1);
            std::vector<double> expected(hours + 1, 0.0), closed(hours + 1, 0.0);
            for (int seed = 0; seed < seeds; ++seed) {
                auto model = make_model(static_cast<uint32_t>(200000 + 100 * arm.room_size + 10 * arm.redraw_every + seed));
                const int rooms_per_school = school_size / arm.room_size;
                std::vector<std::vector<LocationId>> rooms(schools);
                std::vector<std::vector<PersonId>> members(schools);
                PersonId first_infectious = PersonId::invalid_ID();
                for (int sc = 0; sc < schools; ++sc) {
                    for (int r = 0; r < rooms_per_school; ++r) {
                        rooms[sc].push_back(add_room(model));
                    }
                    // Infectious agents at random places: shuffle who is infectious, then fill the rooms in order.
                    std::vector<bool> infectious(school_size, false);
                    std::fill(infectious.begin(), infectious.begin() + infectious_per_school, true);
                    std::shuffle(infectious.begin(), infectious.end(), model.get_rng());
                    for (int p = 0; p < school_size; ++p) {
                        const auto id = add_agent(model, rooms[sc][p / arm.room_size], infectious[p]);
                        members[sc].push_back(id);
                        if (sc == 0 && infectious[p] && first_infectious == PersonId::invalid_ID()) {
                            first_infectious = id;
                        }
                    }
                }
                const double k = set_a(model, first_infectious, a);

                // Closed form for the rooms as they are at the start: exact for fixed rooms, Eq. (9) with a T.
                const auto start = room_states(model, t0);
                for (int h = 0; h <= hours; ++h) {
                    double c = 0.0;
                    for (const auto& [id, r] : start) {
                        if (r.present > 1) {
                            c += r.susceptible * (1.0 - std::exp(-a * r.infectious * h / (r.present - 1)));
                        }
                    }
                    closed[h] += c / schools;
                }

                double cumulative_expected    = 0.0;
                std::vector<int> base(schools, school_size - infectious_per_school);
                for (int h = 0; h < hours; ++h) {
                    const auto t = t0 + mio::abm::hours(h);
                    if (arm.redraw_every > 0 && h > 0 && h % arm.redraw_every == 0) {
                        for (int sc = 0; sc < schools; ++sc) {
                            std::shuffle(members[sc].begin(), members[sc].end(), model.get_rng());
                            // Assign the new room as well: the infection related mobility rules move agents to
                            // their assigned Home, which would otherwise send them back to their first room.
                            for (int p = 0; p < school_size; ++p) {
                                model.assign_location(members[sc][p], rooms[sc][p / arm.room_size]);
                                model.change_location(members[sc][p], rooms[sc][p / arm.room_size]);
                            }
                        }
                    }
                    cumulative_expected += expected_new(room_states(model, t), k);
                    model.evolve(t, dt);
                    expected[h + 1] += cumulative_expected / schools;
                }
                // Cumulative new infections per school, per hour.
                for (int h = 0; h <= hours; ++h) {
                    const auto t = t0 + mio::abm::hours(h);
                    std::vector<int> susceptible(schools, 0);
                    for (int sc = 0; sc < schools; ++sc) {
                        for (const auto& id : members[sc]) {
                            susceptible[sc] += model.get_person(id).get_infection_state(t) == InfectionState::Susceptible;
                        }
                        observed[h].add(base[sc] - susceptible[sc]);
                    }
                }
            }
            for (int h = 0; h <= hours; ++h) {
                out << a << "," << arm.name << "," << arm.room_size << "," << h << "," << observed[h].n << ","
                    << observed[h].mean() << "," << observed[h].se() << "," << expected[h] / seeds << ","
                    << (arm.redraw_every > 0 ? std::string("") : std::to_string(closed[h] / seeds)) << "\n";
            }
            std::cout << "  contact time: a = " << a << ", " << arm.name << " rooms of " << arm.room_size << " done"
                      << std::endl;
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    const std::string dir = argc > 1 ? argv[1] : "theory_minimal_output";
    const int seeds       = argc > 2 ? std::stoi(argv[2]) : 20;
    mio::set_log_level(mio::LogLevel::warn);
    if (!mio::create_directory(dir)) {
        std::cerr << "Could not create " << dir << "\n";
        return 1;
    }
    std::cout << "one step, single infectious agent" << std::endl;
    one_step_single(dir, seeds);
    std::cout << "one step, clustering" << std::endl;
    one_step_clustering(dir, seeds);
    std::cout << "contact time" << std::endl;
    contact_time(dir, seeds);
    std::cout << "results written to " << dir << std::endl;
    return 0;
}
