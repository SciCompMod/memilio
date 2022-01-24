/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Rene Schmieding
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
// functions and operators neccessary for a Contolled Stepper to work with Eigen::VectorXd
// these have to be declared *before* the includes

#include "memilio/math/eigen.h"

namespace std {
Eigen::VectorXd abs(Eigen::VectorXd x) {
    // elementwise operations are defined on arrays within Eigen
    // casts to and from array supposedly cost no runtime when using compiler optimisation 
    return x.array().abs().matrix();
}
}

#include "memilio/math/stepper_wrapper.h"

Eigen::VectorXd operator+ (const double s, const Eigen::VectorXd& v) {
    return (v.array() + s).matrix();
}

Eigen::VectorXd operator/ (const Eigen::VectorXd& v, const Eigen::VectorXd& w) {
    return (v.array() / w.array()).matrix();
}
