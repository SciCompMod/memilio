#include "../include/custom_loggers.h"
#include "../include/constants.h"

LogInfectionStatePerAgeGroup::Type LogInfectionStatePerAgeGroup::log(const mio::abm::Simulation& sim)
{
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(
        Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups()));
    auto curr_time     = sim.get_time();
    const auto persons = sim.get_world().get_persons();

    for (auto i = size_t(0); i < persons.size(); ++i) {
        auto& p    = persons[i];
        auto index = (((size_t)(mio::abm::InfectionState::Count)) * ((uint32_t)p.get_age().get())) +
                     ((uint32_t)p.get_infection_state(curr_time));
        sum[index] += 1;
    }
    return std::make_pair(curr_time, sum);
}

LogInfectionPerLocationTypePerAgeGroup::Type
LogInfectionPerLocationTypePerAgeGroup::log(const mio::abm::Simulation& sim)
{
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(
        Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups()));
    auto curr_time       = sim.get_time();
    const auto locations = sim.get_world().get_locations();

    for (auto i = size_t(0); i < locations.size(); ++i) {
        auto& loc = locations[i];
        for (auto person_id_inf : loc.get_infected_persons()) {
            auto& person = sim.get_world().get_person(person_id_inf);
            auto index   = (((size_t)mio::abm::LocationType::Count) * ((uint32_t)person.get_age().get())) +
                         ((uint32_t)loc.get_type());
            sum[index] += 1;
        }
    }

    return std::make_pair(curr_time, sum);
}

LogLocationTypeAndId::Type LogLocationTypeAndId::log(const mio::abm::Simulation& sim)
{

    Type location_information{};
    location_information.reserve(sim.get_world().get_persons().size());
    for (auto&& person : sim.get_world().get_persons()) {
        location_information.push_back(std::make_tuple(person.get_person_id(), person.get_location().get_type()));
    }
    return location_information;
}

LogWhoInfected::Type LogWhoInfected::log(const mio::abm::Simulation& sim)
{
    Type who_infected{};
    who_infected.reserve(sim.get_world().get_persons().size());
    for (auto&& location : sim.get_world().get_locations()) {
        for (auto&& person_id : location.get_infected_persons()) {
            who_infected.push_back(person_id);
        }
    }
    return who_infected;
}

LogLocationIdAndPersonId::Type LogLocationIdAndPersonId::log(const mio::abm::Simulation& sim)
{

    Type location_information{};
    location_information.reserve(sim.get_world().get_persons().size());
    for (auto&& person : sim.get_world().get_persons()) {
        location_information.push_back(std::make_tuple(person.get_person_id(), person.get_location().get_index()));
    }
    return location_information;
}

LogInfectionDetailed::Type LogInfectionDetailed::log(const mio::abm::Simulation& sim)
{
    Type infected_detailed{};
    infected_detailed.reserve(sim.get_world().get_persons().size());
    for (auto&& location : sim.get_world().get_locations()) {
        for (auto&& person_id : location.get_infected_persons()) {
            infected_detailed.push_back(std::make_tuple(person_id, location.get_index(), location.get_type()));
        }
    }
    return infected_detailed;
}

LogAmountOfInfections::Type LogAmountOfInfections::log(const mio::abm::Simulation& sim)
{
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;

    static Type log(const mio::abm::Simulation& sim)
    {
        Eigen::VectorXd sum = Eigen::VectorXd::Zero(1);
        auto curr_time      = sim.get_time();
        const auto persons  = sim.get_world().get_persons();

        for (auto&& person : sim.get_world().get_persons()) {
            if (person.get_infection_state(curr_time) != mio::abm::InfectionState::Susceptible) {
                sum[0] += 1;
            }
        }
        return std::make_pair(curr_time, sum);
    }
};