/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn
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

namespace epi
{
void print_seir_params(const SeirModel& model)
{
    printf("\n SEIR model set.\n Parameters:\n\t Time incubation:\t %.4f \n\t Time infectious:\t %.4f \n\t transmission "
           "risk:\t %.4f \n\t N:\t %d \n\t E0:\t %d \n\t I0:\t %d \n\t R0:\t %d\n",
           1.0 / model.parameters.get<StageTimeIncubationInv>(), 1.0 / model.parameters.get<StageTimeInfectiousInv>(),
           model.parameters.get<TransmissionRisk>(), (int)model.populations.get_total(),
           (int)model.populations[{Index<SeirInfType>(SeirInfType::E)}], (int)model.populations[{Index<SeirInfType>(SeirInfType::I)}],
           (int)model.populations[{Index<SeirInfType>(SeirInfType::R)}]);
}

} // namespace epi
