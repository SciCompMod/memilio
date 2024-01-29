#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

// TODO: find a meaningfull header name

#include "abm/person.h"
#include "abm/location.h"

namespace mio
{

namespace abm
{

// let a person interact with a location for and at some time
void interact(Person& person, Location& location, TimePoint t, TimeSpan dt, const Parameters& global_parameters,
              Person::RandomNumberGenerator& personal_rng);

// move a person to another location. returns false if the person was at the target location already, true otherwise
bool migrate(Person& person, Location& destination, const std::vector<uint32_t>& cells = {0},
             TransportMode mode = TransportMode::Unknown);
} // namespace abm

} // namespace mio

#endif