#ifndef EPI_ABM_WORLD_BUILDER_H
#define EPI_ABM_WORLD_BUILDER_H

#include "epidemiology/abm/age.h"
#include "epidemiology/abm/world.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/utils/eigen.h"

namespace epi
{
Eigen::VectorXd find_optimal_locations(Eigen::VectorXd& num_people_sorted, int num_locs,
                                       Eigen::MatrixXd& contact_matrix);
void create_locations(uint32_t num_locs, LocationType type, World& world, Eigen::MatrixXd& contact_matrix);
} // namespace epi

#endif
