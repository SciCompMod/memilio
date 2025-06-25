#include "../include/custom_loggers.h"
#include "../include/constants.h"

LogInfectionStatePerAgeGroup::Type LogInfectionStatePerAgeGroup::log(const mio::abm::Simulation& sim)
{
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(
        Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups()));
    auto curr_time     = sim.get_time();
    const auto persons = sim.get_world().get_persons();

    for (auto i = size_t(0); i < persons.size(); ++i) {
        auto& p = persons[i];
        if (p.get_should_be_logged()) {
            auto index = (((size_t)(mio::abm::InfectionState::Count)) * ((uint32_t)p.get_age().get())) +
                         ((uint32_t)p.get_infection_state(curr_time));
            sum[index] += 1;
        }
    }
    return std::make_pair(curr_time, sum);
}

LogInfectionPerLocationTypePerAgeGroup::Type
LogInfectionPerLocationTypePerAgeGroup::log(const mio::abm::Simulation& sim)
{
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(
        Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups()));
    auto curr_time     = sim.get_time();
    auto prev_time     = sim.get_prev_time();
    const auto persons = sim.get_world().get_persons();

    for (auto i = size_t(0); i < persons.size(); ++i) {
        auto& p = persons[i];
        if (p.get_should_be_logged()) {
            if ((p.get_infection_state(prev_time) != mio::abm::InfectionState::Exposed) &&
                (p.get_infection_state(curr_time) == mio::abm::InfectionState::Exposed)) {
                auto index = (((size_t)(mio::abm::LocationType::Count)) * ((uint32_t)p.get_age().get())) +
                             ((uint32_t)p.get_location().get_type());
                sum[index] += 1;
            }
        }
    }
    return std::make_pair(curr_time, sum);
}