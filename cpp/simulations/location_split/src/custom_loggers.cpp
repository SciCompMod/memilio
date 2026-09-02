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

LogPersonStates::Type LogPersonStates::log(const mio::abm::Simulation<>& sim)
{
    Type states{};
    states.reserve(sim.get_model().get_persons().size());
    const auto curr_time = sim.get_time();

    for (const auto& person : sim.get_model().get_persons()) {
        states.push_back(PersonState{person.get_id(), person.get_location(), person.get_location_type(),
                                     person.get_infection_state(curr_time), person.get_age()});
    }
    return states;
}

LogInfectiousPersons::Type LogInfectiousPersons::log(const mio::abm::Simulation<>& sim)
{
    Type infectious{};
    const auto curr_time = sim.get_time();

    for (const auto& person : sim.get_model().get_persons()) {
        if (person.is_infected(curr_time) && person.get_infection().get_viral_shed(curr_time) > 0) {
            infectious.push_back(person.get_id());
        }
    }
    return infectious;
}

namespace
{
/// @brief Collect the Location of the given LocationType assigned to each Person.
std::vector<std::tuple<mio::abm::PersonId, mio::abm::LocationId>>
assigned_locations(const mio::abm::Simulation<>& sim, mio::abm::LocationType type)
{
    std::vector<std::tuple<mio::abm::PersonId, mio::abm::LocationId>> result{};
    result.reserve(sim.get_model().get_persons().size());
    for (const auto& person : sim.get_model().get_persons()) {
        result.emplace_back(person.get_id(), person.get_assigned_location(type));
    }
    return result;
}
} // namespace

LogHouseholdId::Type LogHouseholdId::log(const mio::abm::Simulation<>& sim)
{
    return assigned_locations(sim, mio::abm::LocationType::Home);
}

LogWorkId::Type LogWorkId::log(const mio::abm::Simulation<>& sim)
{
    return assigned_locations(sim, mio::abm::LocationType::Work);
}

LogLocationIdAndType::Type LogLocationIdAndType::log(const mio::abm::Simulation<>& sim)
{
    Type location_information{};
    location_information.reserve(sim.get_model().get_locations().size());
    for (const auto& location : sim.get_model().get_locations()) {
        location_information.emplace_back(location.get_id(), location.get_type());
    }
    return location_information;
}
