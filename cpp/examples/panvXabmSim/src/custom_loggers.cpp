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
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(Eigen::Index(mio::abm::LocationType::Count));
    auto curr_time      = sim.get_time();

    // If there is no interresting person ids to logged, log all persons.
    // Otherwise log accordingly
    for (auto&& loc : sim.get_world().get_locations()) {
        sum[(int)(loc.get_type())] += 1;
    }
    return std::make_pair(curr_time, sum);
}