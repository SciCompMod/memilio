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
#include "abm/world_builder.h"
#include "memilio/math/eigen.h"

#include "nlopt.hpp"
#include <cmath>

namespace mio
{

void compute_f(const Eigen::VectorXd& x, Eigen::VectorXd& f, Eigen::VectorXd& y, int l, int n)
{
    // reshape x to get a matrix for easier notation
    Eigen::Map<const Eigen::MatrixXd> M(x.data(), n, l);

    // the average of the contact matrices per location is given
    // go through every entry of the contact matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // sum over enty (i,j) of the contact matrices of every location
            // scale entry of local contact matrix according to effective contacts
            double sum = 0;
            for (int k = 0; k < l; k++) {
                double coef = std::min(x(x.size() - 1) / M.col(k).sum(), 1.0);
                if (i != j) {
                    sum += coef * M(j, k) * M(i, k);
                }
                else {
                    sum += coef * (M(i, k) - 1) * M(i, k);
                }
            }

            //compare average of local contact matrices to original contact matrix
            f(i + n * j) = sum - y(n + l + i + n * j) * M.row(i).sum();
        }
    }
}

// implement Jacobian manually for better runtime
void compute_df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac, Eigen::VectorXd& y, int l, int n)
{
    fjac.setZero();
    // reshape x to get a matrix for easier notation
    Eigen::Map<const Eigen::MatrixXd> M(x.data(), n, l);

    // use the same structure as in the functional
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // compute derivative w.r.t. every entry of x (characterized by location k and age group m)
            for (int k = 0; k < l; k++) {
                double der_eff_contacts;
                for (int m = 0; m < n; m++) {
                    // compute the two terms via product rule
                    // first determine the coefficients that only depend on location
                    double der_term1, der_term2;
                    der_eff_contacts = 0;
                    der_term2        = std::min(x(x.size() - 1) / M.col(k).sum(), 1.0);
                    if (std::min(x(x.size() - 1) / M.col(k).sum(), 1.0) < 1.0) {
                        der_term1        = -x(x.size() - 1) / (M.col(k).sum() * M.col(k).sum());
                        der_eff_contacts = 1 / M.col(k).sum();
                        if (i == j) {
                            der_eff_contacts *= (M(j, k) - 1) * M(i, k);
                        }
                        else {
                            der_eff_contacts *= M(j, k) * M(i, k);
                        }
                    }
                    else {
                        der_term1 = 0;
                    }
                    if (i == j) {
                        der_term1 *= (M(j, k) - 1) * M(i, k);
                    }
                    else {
                        der_term1 *= M(j, k) * M(i, k);
                    }
                    // now compute and save derivative
                    if (i == m && i != j) {
                        der_term2 *= M(j, k);
                    }
                    else if (j == m && i != j) {
                        der_term2 *= M(i, k);
                    }
                    else if (j == m && i == j) {
                        der_term2 *= (2 * M(i, k) - 1);
                    }
                    else {
                        der_term2 = 0;
                    }
                    if (i == m) {
                        der_term2 -= y(n + l + i + n * j);
                    }
                    fjac(i + n * j, m + n * k) = der_term1 + der_term2;
                }
                fjac(i + n * j, n * l) += der_eff_contacts;
            }
        }
    }
}

struct Data {
    int l, n;
    Eigen::VectorXd y;
};

double myfunc(const std::vector<double>& x, std::vector<double>& grad, void* erased_data)
{
    // x in Eigen vector
    const double* ptr = &x[0];
    Eigen::Map<const Eigen::VectorXd> x_eigen(ptr, x.size());
    // remap my_func_data
    Data* data         = reinterpret_cast<Data*>(erased_data);
    Eigen::VectorXd& y = data->y;
    int l              = data->l;
    int n              = data->n;

    // computation of multidimensional functional
    Eigen::VectorXd f(n * n);
    compute_f(x_eigen, f, y, l, n);

    if (!grad.empty()) {
        // compute Jacobi matrix of f
        Eigen::MatrixXd fjac(n * n, n * l + 1);
        compute_df(x_eigen, fjac, y, l, n);
        // compute gradient of norm of f
        for (int i = 0; i < n * n; i++) {
            double sum = 0;
            for (int j = 0; j < n * n; j++) {
                sum += f(j) * fjac(j, i);
            }
            grad[i] = 1 / f.norm() * sum;
        }
    }

    // quadratic minimization problem
    return f.norm();
}

void constraints(unsigned /*m*/, double* result, unsigned /*x_len*/, const double* x, double* grad, void* erased_data)
{
    Data* data         = reinterpret_cast<Data*>(erased_data);
    Eigen::VectorXd& y = data->y;
    int l              = data->l;
    int n              = data->n;
    // reshape x to get a matrix for easier notation
    Eigen::Map<const Eigen::MatrixXd> M(x, n, l);

    Eigen::Map<Eigen::VectorXd> g(result, n + l);
    if (grad != nullptr) {
        // at every location contains a certain number of people
        // at the moment: every location has the same number of people
        g.segment(0, l) = (M.colwise().sum().transpose() - y.segment(0, l));

        // the number of people per age group is given
        g.segment(l, n) = (M.rowwise().sum() - y.segment(l, n));
    }
}

Eigen::VectorXd find_optimal_locations(Eigen::VectorXd& num_people_sorted, int num_locs,
                                       Eigen::MatrixXd& contact_matrix, Eigen::VectorXd& size_locs)
{
    int n = int(contact_matrix.rows());
    int l = num_locs;

    Eigen::VectorXd y((n + 1) * n + l);

    // size of the new locations
    y.segment(0, l) = size_locs;
    // number of people per age group
    y.segment(l, n) = num_people_sorted;
    // entries of contact matrix (entry (i,j) goes to (j*n + i))
    Eigen::Map<Eigen::VectorXd> v(contact_matrix.data(), n * n);
    y.segment(l + n, n * n) = v;
    Data data{l, n, std::move(y)};

    //last entry of x: effective contacts
    // set the following starting values provide a rough fit.
    std::vector<double> x(n * l + 1, num_people_sorted.sum() / (n * l));
    x[n * l] = 100;

    //get instance of optimization problem
    nlopt::opt opt(nlopt::GN_ISRES, n * l + 1);

    //set lower bounds for variables: number of people at a location are nonnegative
    std::vector<double> lb(n * l + 1, 0);
    // effective contacts are at least 1
    lb[n * l] = 1;
    std::vector<double> ub(n * l + 1, num_people_sorted.sum());
    ub[n * l] = 1000;
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    opt.set_min_objective(myfunc, &data);
    const std::vector<double> tol(n * l + 1, 1e-8);

    opt.add_equality_mconstraint(constraints, &data, tol);
    opt.set_xtol_rel(1e-6);
    opt.set_ftol_abs(1e-6);
    double minf;

    try {
        /*nlopt::result result = */ opt.optimize(x, minf);
        //std::cout << "result : " << int(result) << std::endl;
    }
    catch (std::exception& e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    double* ptr = &x[0];
    Eigen::Map<Eigen::VectorXd> x_eigen(ptr, n * l + 1);
    return x_eigen;
}

void create_locations(uint32_t num_locs, LocationType type, World& world, Eigen::MatrixXd& contact_matrix,
                      Eigen::VectorXd& size_locs)
{
    auto num_ages = uint32_t(contact_matrix.rows());

    // sort people according to their age groups
    std::vector<std::vector<uint32_t>> people_sorted(num_ages);
    uint32_t index = 0;
    for (auto& person : world.get_persons()) {
        people_sorted[size_t(person.get_age())].push_back(index);
        ++index;
    }

    // get vector with number of people per age group
    Eigen::VectorXd num_people_sorted(num_ages);
    for (size_t i = 0; i < num_ages; i++) {
        num_people_sorted(i) = double(people_sorted[i].size());
    }

    // find optimal number of people per age group for every new location by minimizing nonlinear optimization problem
    Eigen::VectorXd sol = find_optimal_locations(num_people_sorted, num_locs, contact_matrix, size_locs);

    // create locations and assign them to people
    Eigen::VectorXi current_index = Eigen::VectorXi::Zero(num_ages);
    for (size_t i = 0; i < num_locs; i++) {
        auto loc = world.add_location(type);
        //add enough people of each age group
        for (size_t j = 0; j < num_ages; j++) {
            int num = (int)round(sol(j + num_ages * i));
            for (int counter = 0; counter < num && current_index(j) < num_people_sorted(j); counter++) {
                auto person = world.get_persons().begin() + people_sorted[j][current_index(j)];
                person->set_assigned_location(loc);
                current_index(j)++;
            }
        }
    }
    // add people that have not been assigned before
    uint32_t loc_index = 0;
    for (size_t j = 0; j < num_ages; j++) {
        while (current_index(j) < num_people_sorted(j) - 1) {
            auto person = world.get_persons().begin() + people_sorted[j][current_index(j)];
            person->set_assigned_location({loc_index, type});
            current_index(j)++;
            if (loc_index < world.get_locations()[(size_t)type].size() - 1) {
                loc_index++;
            }
            else {
                loc_index = 0;
            }
        }
    }
}

void compute_average_contact_matrix(Eigen::MatrixXd& average_contacs, Eigen::VectorXd& x, uint32_t num_locs)
{
    auto num_ages = uint32_t(average_contacs.rows());
    // reshape x to get a matrix for easier notation
    Eigen::Map<const Eigen::MatrixXd> M(x.data(), num_ages, num_locs);

    //compute average for every entry of the matrix
    for (uint32_t i = 0; i < num_ages; i++) {
        for (uint32_t j = 0; j < num_ages; j++) {
            //sum over all locations
            double sum = 0;
            for (uint32_t l = 0; l < num_locs; l++) {
                double coef = std::min(x(x.size() - 1) / M.col(l).sum(), 1.0);
                if (i == j) {
                    sum += coef * M(i, l) * (M(i, l) - 1);
                }
                else {
                    sum += coef * M(i, l) * M(j, l);
                }
            }
            if (M.row(i).sum() != 0) {
                average_contacs(i, j) = sum / M.row(i).sum();
            }
        }
    }
}
} // namespace mio
