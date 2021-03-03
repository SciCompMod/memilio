#ifndef SEIR_H
#define SEIR_H

#include "epidemiology/model/compartmentalmodel.h"
#include "epidemiology/model/populations.h"
#include "epidemiology/utils/parameter_set.h"
#include "epidemiology/secir/contact_matrix.h"

namespace epi
{

/****************************************
 * Define SEIR categories using an enum *
 ****************************************/

enum class SeirInfType
{
    S,
    E,
    I,
    R,
    Count = 4
};

using SeirPopulations = Populations<SeirInfType>;

/**************************
 * define some parameters *
 **************************/


struct TransmissionRisk
{
    using Type=double;
    static constexpr Type get_default() {
        return 1.0;
    }
};
struct StageTimeIncubationInv
{
    using Type=double;
    static constexpr Type get_default() {
        return 1.0 / 5.2;
    }
};
struct StageTimeInfectiousInv
{
    using Type=double;
    static constexpr Type get_default() {
        return 1.0 / 6.0;
    }
};
struct ContactFrequency
{
    using Type=ContactMatrix;
    static Type get_default() {
        return ContactMatrix{1};
    }
};

using SeirParameters = ParameterSet<TransmissionRisk,
                                    StageTimeIncubationInv,
                                    StageTimeInfectiousInv,
                                    ContactFrequency>;


/********************
 * define the model *
 ********************/

class SeirModel : public CompartmentalModel<SeirPopulations, SeirParameters>
{
    using Po = SeirPopulations;
    using Pa = SeirParameters;

public:
    SeirModel()
    {
#if !USE_DERIV_FUNC
        //S to E
        this->add_flow(std::make_tuple(SeirInfType::S), std::make_tuple(SeirInfType::E),
                       [](Pa const& p, Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           ScalarType cont_freq_eff = params.contact_frequency.get_matrix_at(t)(0, 0);

                           //TODO: We should probably write a static Po::get_total_from function
                           ScalarType divN = 1.0 / (Po::get_from(y, SeirInfType::S) + Po::get_from(y, SeirInfType::E) +
                                                    Po::get_from(y, SeirInfType::I) + Po::get_from(y, SeirInfType::R));
                           return cont_freq_eff * Po::get_from(y, SeirInfType::S) * Po::get_from(pop, SeirInfType::I) *
                                  divN;
                       });

        //E to I
        this->add_flow(std::make_tuple(SeirInfType::E), std::make_tuple(SeirInfType::I),
                       [](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                           return p.times.get_incubation_inv() * Po::get_from(y, SeirInfType::E);
                       });

        //I to R
        this->add_flow(std::make_tuple(SeirInfType::I), std::make_tuple(SeirInfType::R),
                       [](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                           return p.times.get_infectious_inv() * Po::get_from(y, SeirInfType::I);
                       });
#endif
    }

#if USE_DERIV_FUNC
    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop,
                         Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override
    {
        auto& params = this->parameters;
        double S2E_coeff = params.get<ContactFrequency>().get_matrix_at(t)(0, 0) * params.get<TransmissionRisk>()
                / populations.get_total();

        dydt[(size_t)SeirInfType::S] = -S2E_coeff * y[(size_t)SeirInfType::S] * pop[(size_t)SeirInfType::I];
        dydt[(size_t)SeirInfType::E] = S2E_coeff * y[(size_t)SeirInfType::S] * pop[(size_t)SeirInfType::I] -
                                    params.get<StageTimeIncubationInv>() * y[(size_t)SeirInfType::E];
        dydt[(size_t)SeirInfType::I] = params.get<StageTimeIncubationInv>() * y[(size_t)SeirInfType::E] -
                                    params.get<StageTimeInfectiousInv>() * y[(size_t)SeirInfType::I];
        dydt[(size_t)SeirInfType::R] = params.get<StageTimeInfectiousInv>() * y[(size_t)SeirInfType::I];
    }

#endif // USE_DERIV_FUNC
};

/**
 * prints given parameters
 * @param[in] params the SeirParams parameter object
 */
void print_seir_params(const SeirModel& model);

} // namespace epi

#endif
