#ifndef SEIR_H
#define SEIR_H

#include "epidemiology/model/compartmentalmodel.h"
#include "epidemiology/model/populations.h"
#include "epidemiology/secir/seir_params.h"

namespace epi
{

// Create template specialization for the simple SEIR model

enum class InfectionType
{
    S,
    E,
    I,
    R,
    Count = 4
};

using SeirModel = CompartmentalModel<Populations<InfectionType>, SeirParams>;

/**
 * prints given parameters
 * @param[in] params the SeirParams parameter object
 */
void print_seir_params(const SeirModel& model)
{
    printf("\n SEIR model set.\n Parameters:\n\t Time incubation:\t %.4f \n\t Time infectious:\t %.4f \n\t contact "
           "rate:\t %.4f \n\t N:\t %d \n\t E0:\t %d \n\t I0:\t %d \n\t R0:\t %d\n",
           1.0 / model.parameters.times.get_incubation_inv(), 1.0 / model.parameters.times.get_infectious_inv(),
           model.parameters.times.get_cont_freq(), (int)model.populations.get_total(),
           (int)model.populations.get(InfectionType::E), (int)model.populations.get(InfectionType::I),
           (int)model.populations.get(InfectionType::R));
}

//TODO: Once the flows can be defined at compile time, we might not need
// a factory function anymore. Its main purpose is to instantiate the
// intercompartmental flows
SeirModel create_seir_model()
{
    using Pa = SeirParams;
    using Po = Populations<InfectionType>;
    SeirModel model;

    //S to E
    model.add_flow(std::make_tuple(InfectionType::S), std::make_tuple(InfectionType::E),
                   [](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                       ScalarType cont_freq_eff = p.times.get_cont_freq() * p.dampings.get_factor(t);

                       //TODO: We should probably write a static Po::get_total_from function
                       ScalarType divN = 1.0 / (Po::get_from(y, InfectionType::S) + Po::get_from(y, InfectionType::E) +
                                                Po::get_from(y, InfectionType::I) + Po::get_from(y, InfectionType::R));
                       return cont_freq_eff * Po::get_from(y, InfectionType::S) * Po::get_from(y, InfectionType::I) *
                              divN;
                   });

    //E to I
    model.add_flow(std::make_tuple(InfectionType::E), std::make_tuple(InfectionType::I),
                   [](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                       return p.times.get_incubation_inv() * Po::get_from(y, InfectionType::E);
                   });

    //I to R
    model.add_flow(std::make_tuple(InfectionType::I), std::make_tuple(InfectionType::R),
                   [](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                       return p.times.get_infectious_inv() * Po::get_from(y, InfectionType::I);
                   });
}

} // namespace epi

#endif
