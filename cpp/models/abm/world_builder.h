/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Elisabeth Kluth
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef EPI_ABM_WORLD_BUILDER_H
#define EPI_ABM_WORLD_BUILDER_H

#include "abm/age.h"
#include "abm/world.h"
#include "abm/person.h"
#include "memilio/math/eigen.h"

namespace mio
{
/*
 * Helper function:
 * Compute vector valued function f that has to be optimized
 */
void compute_f(const Eigen::VectorXd& x, Eigen::VectorXd& f, Eigen::VectorXd& y);

/* 
 * Helper function:
 * Compute Jacobian of vector valued function f using Eigen
 */
void compute_df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac, Eigen::VectorXd& y);

/*
 * Compute the quadratic function that is to be minimized
 */
double myfunc(const std::vector<double>& x, std::vector<double>& grad, void* my_func_data);

/*
 * Compute the constraint for the optimization
 */
void constraints(unsigned m, double* result, unsigned x_len, const double* x, double* grad, void* data);

/*
 * Determines how to divide people over different locations such that in average a given contact matrix is met.
 * @param num_locs number of locations
 * @param num_people_sorted number of people of each age group
 * @param contact_matrix in average contacts of people should be given by this matrix 
 */
Eigen::VectorXd find_optimal_locations(Eigen::VectorXd& num_people_sorted, int num_locs,
                                       Eigen::MatrixXd& contact_matrix, Eigen::VectorXd& size_locs);

/*
 * Create and assign new locations such that a given contact matrix is met.
 * @param num_locs number of new locations
 * @param type location type of new locations
 * @param contact_matrix in average contacts of people at the new locations should be given by this matrix.
 */
void create_locations(uint32_t num_locs, LocationType type, World& world, Eigen::MatrixXd& contact_matrix,
                      Eigen::VectorXd& size_locs);

/*
 * Compute the average contact matrix of the given people
 * @param M average contact matrix
 * @param x number of people sorted by age group and location
 * @param num_locs number of locations
 */
void compute_average_contact_matrix(Eigen::MatrixXd& M, Eigen::VectorXd& x, uint32_t num_locs);
} // namespace mio

#endif
