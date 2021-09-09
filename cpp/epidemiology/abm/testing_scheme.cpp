#include "epidemiology/abm/testing_scheme.h"
#include "epidemiology/abm/world.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/abm/parameters.h"
#include "epidemiology/utils/random_number_generator.h"

namespace epi
{

TestingScheme::TestingScheme(TimeSpan interval, double probability)
    : m_time_interval(interval)
    , m_probability(probability)
{
}

TestingScheme::TestingScheme()
    : TestingScheme(seconds(std::numeric_limits<int>::max()), 1)
{
}

bool TestingScheme::run_scheme(Person& person, const GlobalTestingParameters& params) const
{
    if (person.get_time_since_test() > m_time_interval) {
        double random = UniformDistribution<double>::get_instance()();
        if (random < m_probability) {
            person.get_tested(params.get<Sensitivity>(), params.get<Specificity>());
        }
    }
    return !person.is_in_quarantine();
}

} // namespace epi