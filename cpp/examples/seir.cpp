/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#include "seir/seir.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 1;
    double dt   = 0.001;

    mio::log_info("Simulating SEIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::SeirModel model;

    double total_population = 10000;
    model.populations[{mio::Index<mio::SeirInfType>(mio::SeirInfType::E)}] = 100;
    model.populations[{mio::Index<mio::SeirInfType>(mio::SeirInfType::I)}] = 100;
    model.populations[{mio::Index<mio::SeirInfType>(mio::SeirInfType::R)}] = 100;
    model.populations[{mio::Index<mio::SeirInfType>(mio::SeirInfType::S)}] = total_population - model.populations[{mio::Index<mio::SeirInfType>(mio::SeirInfType::E)}]
                                                                                              - model.populations[{mio::Index<mio::SeirInfType>(mio::SeirInfType::I)}]
                                                                                              - model.populations[{mio::Index<mio::SeirInfType>(mio::SeirInfType::R)}];
    // suscetible now set with every other update
    // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
    model.parameters.set<mio::StageTimeIncubationInv>(1./5.2);
    model.parameters.set<mio::StageTimeInfectiousInv>(1./6);
    model.parameters.set<mio::TransmissionRisk>(0.04);
    model.parameters.get<mio::ContactFrequency>().get_baseline()(0, 0) = 10;

    print_seir_params(model);

    auto seir = simulate(t0, tmax, dt, model);

    printf("\n number total: %f\n",
           seir.get_last_value()[0] + seir.get_last_value()[1] + seir.get_last_value()[2] + seir.get_last_value()[3]);
}
