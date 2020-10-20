#ifndef SEIR_H
#define SEIR_H

#include "epidemiology/model/compartmentalmodel.h"
#include "epidemiology/model/populations.h"
#include "epidemiology/secir/seir_params.h"

namespace epi
{

// Create template specialization for the simple SEIR model

enum class SeirInfType
{
    S,
    E,
    I,
    R,
    Count = 4
};

class SeirModel : public CompartmentalModel<Populations<SeirInfType>, SeirParams>
{
    using Pa = SeirParams;
    using Po = Populations<SeirInfType>;

public:
    SeirModel()
    {
        //S to E
        this->add_flow(std::make_tuple(SeirInfType::S), std::make_tuple(SeirInfType::E),
                       [](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           ScalarType cont_freq_eff = p.times.get_cont_freq() * p.dampings.get_factor(t);

                           //TODO: We should probably write a static Po::get_total_from function
                           ScalarType divN = 1.0 / (Po::get_from(y, SeirInfType::S) + Po::get_from(y, SeirInfType::E) +
                                                    Po::get_from(y, SeirInfType::I) + Po::get_from(y, SeirInfType::R));
                           return cont_freq_eff * Po::get_from(y, SeirInfType::S) * Po::get_from(y, SeirInfType::I) *
                                  divN;
                       });

        //E to I
        this->add_flow(std::make_tuple(SeirInfType::E), std::make_tuple(SeirInfType::I),
                       [](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           return p.times.get_incubation_inv() * Po::get_from(y, SeirInfType::E);
                       });

        //I to R
        this->add_flow(std::make_tuple(SeirInfType::I), std::make_tuple(SeirInfType::R),
                       [](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           return p.times.get_infectious_inv() * Po::get_from(y, SeirInfType::I);
                       });
    }
};

/**
 * prints given parameters
 * @param[in] params the SeirParams parameter object
 */
void print_seir_params(const SeirModel& model);

} // namespace epi

#endif
