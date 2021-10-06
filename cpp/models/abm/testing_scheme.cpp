#include "abm/testing_scheme.h"
#include "abm/world.h"
#include "abm/location.h"
#include "abm/parameters.h"
#include "memilio/utils/random_number_generator.h"

namespace mio
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
    if (person.get_time_since_negative_test() > m_time_interval) {
        double random = UniformDistribution<double>::get_instance()();
        if (random < m_probability) {
            return !person.get_tested(params.get<AntigenTest>());
        }
    }
    return true;
}

} // namespace mio