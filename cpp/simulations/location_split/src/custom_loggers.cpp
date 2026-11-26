#include "custom_loggers.h"

#include "abm/model.h"

LogInfectionStatePerAgeGroup::Type LogInfectionStatePerAgeGroup::log(const mio::abm::Simulation<>& sim)
{
    const auto num_states = (size_t)mio::abm::InfectionState::Count;
    const auto num_groups = (size_t)sim.get_model().parameters.get_num_groups();

    Eigen::VectorX<ScalarType> sum = Eigen::VectorX<ScalarType>::Zero(Eigen::Index(num_states * num_groups));
    const auto curr_time           = sim.get_time();

    for (const auto& person : sim.get_model().get_persons()) {
        const auto index = num_states * (size_t)person.get_age().get() + (size_t)person.get_infection_state(curr_time);
        sum[Eigen::Index(index)] += 1;
    }
    return std::make_pair(curr_time, sum);
}

LogTotalInfections::Type LogTotalInfections::log(const mio::abm::Simulation<>& sim)
{
    Eigen::VectorX<ScalarType> sum = Eigen::VectorX<ScalarType>::Zero(1);
    const auto curr_time           = sim.get_time();

    for (const auto& person : sim.get_model().get_persons()) {
        if (person.get_infection_state(curr_time) != mio::abm::InfectionState::Susceptible) {
            sum[0] += 1;
        }
    }
    return std::make_pair(curr_time, sum);
}
