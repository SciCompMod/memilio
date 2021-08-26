#ifndef EPI_ABM_TESTING_SCHEME_H
#define EPI_ABM_TESTING_SCHEME_H

#include "epidemiology/abm/time.h"
#include "epidemiology/abm/parameters.h"
#include "epidemiology/abm/time.h"

#include <functional>

namespace epi
{

class Person;

/**
 * Testing Scheme to regular test people
 */
class TestingScheme
{
public:
    /**
     * create a testing scheme.
     * @param interval the interval in which people who go to the location get tested
     * @param probability probability with which a person gets tested
     */
    TestingScheme(TimeSpan interval, double probability);

    /**
     * create a default testing scheme such that no regular testing happens
     */
    TestingScheme();

    /**
     * get the time interval of this testing scheme
     */
    TimeSpan get_interval() const
    {
        return m_time_interval;
    }

    /**
     * get probability of this testing scheme
     */
    double get_probability() const
    {
        return m_probability;
    }

    /**
     * set the time interval of this testing scheme
     */
    void set_interval(TimeSpan t)
    {
        m_time_interval = t;
    }

    /**
     * set probability of this testing scheme
     */
    void set_probability(double p)
    {
        m_probability = p;
    }

    void run_scheme(Person& person, const GlobalTestingParameters& params) const;

private:
    TimeSpan m_time_interval;
    double m_probability;
};

} // namespace epi

#endif
