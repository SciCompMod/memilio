#include "epidemiology/secir/seir.h"

namespace epi
{
void print_seir_params(const SeirModel& model)
{
    printf("\n SEIR model set.\n Parameters:\n\t Time incubation:\t %.4f \n\t Time infectious:\t %.4f \n\t transmission "
           "risk:\t %.4f \n\t N:\t %d \n\t E0:\t %d \n\t I0:\t %d \n\t R0:\t %d\n",
           1.0 / model.parameters.get<StageTimeIncubationInv>(), 1.0 / model.parameters.get<StageTimeInfectiousInv>(),
           model.parameters.get<TransmissionRisk>(), (int)model.populations.get_total(),
           (int)model.populations[{SeirInfType::E}], (int)model.populations[{SeirInfType::I}],
           (int)model.populations[{SeirInfType::R}]);
}

} // namespace epi
