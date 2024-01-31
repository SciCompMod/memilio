#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

// TODO: find a meaningfull header name

#include "abm/person.h"
#include "abm/location.h"
#include <vector>

namespace mio
{

namespace abm
{

// add the contribution of person to the exposure rates at location
void add_exposure_contribution(Location::AirExposureRates& local_air_exposure,
                               Location::ContactExposureRates& local_contact_exposure, const Person& person,
                               const Location& location, TimePoint t, TimeSpan dt);

// let a person interact with a location for and at some time
void interact(Person& person, const Location& location, const Location::AirExposureRates& local_air_exposure,
              const Location::ContactExposureRates& local_contact_exposure, const TimePoint t, const TimeSpan dt,
              const Parameters& global_parameters, PersonalRandomNumberGenerator& personal_rng);

// interact, but it computes the correct exposures for you
void interact(Person& person, const Location& location, const std::vector<Person>& local_population, const TimePoint t,
              const TimeSpan dt, const Parameters& global_parameters, PersonalRandomNumberGenerator& personal_rng);

// move a person to another location. returns false if the person was at the target location already, true otherwise
bool migrate(Person& person, const Location& destination, const std::vector<uint32_t>& cells = {0},
             const TransportMode mode = TransportMode::Unknown);
} // namespace abm

} // namespace mio

#endif