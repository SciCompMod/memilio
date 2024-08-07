/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Anna Wendler
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

#ifndef IDE_SECIR_INITIALIZE_FROM_ODE_H
#define IDE_SECIR_INITIALIZE_FROM_ODE_H

#include "ide_secir/model_ide.h"
#include "ide_secir/parameters.h"
#include "ode_secir/model.h"
#include "infection_state.h"
#include "memilio/config.h"
#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"

namespace mio
{
namespace isecir
{
/**
* @brief Takes the resulting compartments of an ODE simulation and computes the respective flows and stores them in a TimeSeries. 
*
* With t_max and t_window it can be determined for which time window the flows will be computed. 
* dt_reference is the time step size of the ODE simulation. By default, we compute the flows with the same time step size.
* It is also possible to compute the corresponding flows for a bigger time step which can be given by dt_comparison. Here we assume, that
* dt_comparison is a multiple of dt_reference. 
*
* @param[in] model_ode ODE model that is used.
* @param[in] compartments TimeSeries containing compartments from ODE simulation
* @param[out] flows TimeSeries where the computed flows will be stored. 
* @param[in] t_max Maximal time for which the flows are computed.
* @param[in] t_window Time window before t_max for which flows will be computed.
* @param[in] dt_reference Reference time step size (coming from ODE simulation).
* @param[in] dt_comparison Time step size. Default is dt_reference. 
*/
void get_flows_from_ode_compartments(mio::osecir::Model<ScalarType>& model_ode,
                                     mio::TimeSeries<ScalarType> compartments, mio::TimeSeries<ScalarType>& flows,
                                     ScalarType t_max, ScalarType t_window, ScalarType dt_reference,
                                     ScalarType dt_comparison = 0.);

/**
* @brief Computes the inital flows that are needed for an IDE simulation given we have the compartments from an ODE 
* simulation for an adequate time window before t0_ide. 
*
* Using model_ode we get simulated results for the compartments stored in secihurd_ode. Then we compute the needed iniital flows 
* for the IDE model model_ide until time t0_ide. Here we assumen, that model_ode and model_ide are matching, i.e. that the parameters
* of model_ideare chosen such that model_ide shopuld reduce to model_ode by choosing exponentially distributed transitions with the
* the corresponding mean stay times etc. The time step size of ODE and IDE simulation can be chosen independently. However, we assume
* that dt_is a multiple of dt_ide. 
*
* @param[in] model_ode ODE model that is used.
* @param[in] model_ide IDE model that is used.
* @param[in] secihurd_ode TimeSeries containing compartments from ODE simulation. 
* @param[in] t0_ide Start time of IDE simulation until which we need to compute flows. 
* @param[in] dt_ode Time step size of ODE simulation.
* @param[in] dt_ide Time step size of IDE simulation. 
*/
void compute_initial_flows_for_ide_from_ode(mio::osecir::Model<ScalarType>& model_ode, mio::isecir::Model& model_ide,
                                            mio::TimeSeries<ScalarType> secihurd_ode, ScalarType t0_ide,
                                            ScalarType dt_ode, ScalarType dt_ide);

} // namespace isecir
} // namespace mio

#endif //IDE_SECIR_INITIALIZE_FROM_ODE_H