#include "epidemiology/secir/seir.h"

namespace epi
{
void print_seir_params(const SeirModel& model)
{
    printf("\n SEIR model set.\n Parameters:\n\t Time incubation:\t %.4f \n\t Time infectious:\t %.4f \n\t contact "
           "rate:\t %.4f \n\t N:\t %d \n\t E0:\t %d \n\t I0:\t %d \n\t R0:\t %d\n",
           1.0 / model.parameters.times.get_incubation_inv(), 1.0 / model.parameters.times.get_infectious_inv(),
           model.parameters.times.get_cont_freq(), (int)model.populations.get_total(),
           (int)model.populations.get(SeirInfType::E), (int)model.populations.get(SeirInfType::I),
           (int)model.populations.get(SeirInfType::R));
}
} // namespace epi
