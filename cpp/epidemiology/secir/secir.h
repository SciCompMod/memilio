#ifndef SECIR_H
#define SECIR_H

#include "epidemiology/model/compartmentalmodel.h"
#include "epidemiology/model/populations.h"
#include "epidemiology/secir/secir_params.h"

namespace epi
{

// Create template specializations for the age resolved
// SECIHURD model

enum class InfectionType
{
    S,
    E,
    C,
    I,
    H,
    U,
    R,
    D,
    Count
};

enum class AgeGroup5
{
    Group1,
    Group2,
    Group3,
    Group4,
    Group5,
    Count = 5
};

enum class AgeGroup3
{
    Group1,
    Group2,
    Group3,
    Count = 1
};

//TODO: Once the secir_params have been replaced with a
// more abstract compile time map, this is not needed anymore
enum class AgeGroup1
{
    All,
    Count = 1
};

// some convenience typedefs
using SecirModel1 = CompartmentalModel<Populations<AgeGroup1, InfectionType>, SecirParams<(size_t)AgeGroup1::Count>>;
using SecirModel3 = CompartmentalModel<Populations<AgeGroup3, InfectionType>, SecirParams<(size_t)AgeGroup3::Count>>;
using SecirModel5 = CompartmentalModel<Populations<AgeGroup5, InfectionType>, SecirParams<(size_t)AgeGroup5::Count>>;

//TODO: Once the flows can be defined at compile time, we might not need
// a factory function anymore. Its main purpose is to instantiate the
// intercompartmental flows
template <class AgeGroup>
CompartmentalModel<Populations<AgeGroup, InfectionType>, SecirParams<(size_t)AgeGroup::Count>> create_secir_model()
{
    using Pa = SecirParams<(size_t)AgeGroup::Count>;
    using Po = Populations<AgeGroup, InfectionType>;
    CompartmentalModel<Po, Pa> model;

    // alpha  // percentage of asymptomatic cases
    // beta // risk of infection from the infected symptomatic patients
    // rho   // hospitalized per infectious
    // theta // icu per hospitalized
    // delta  // deaths per ICUs
    // 0: S,      1: E,     2: C,     3: I,     4: H,     5: U,     6: R,     7: D

    size_t n_agegroups = model.parameters.get_num_groups();
    for (size_t i = 0; i < n_agegroups; i++) {
        for (size_t j = 0; j < n_agegroups; j++) {

            // Si to Ei individually for each age group j
            model.add_flow(std::make_tuple((AgeGroup)i, InfectionType::S),
                           std::make_tuple((AgeGroup)i, InfectionType::E),
                           [i, j](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                               ScalarType cont_freq_eff =
                                   p.get_contact_patterns().get_cont_freq_mat().get_cont_freq(i, j) *
                                   p.get_contact_patterns().get_cont_freq_mat().get_dampings(i, j).get_factor(t);
                               ScalarType Nj = Po::get_from(y, (AgeGroup)j, InfectionType::S) +
                                               Po::get_from(y, (AgeGroup)j, InfectionType::E) +
                                               Po::get_from(y, (AgeGroup)j, InfectionType::C) +
                                               Po::get_from(y, (AgeGroup)j, InfectionType::I) +
                                               Po::get_from(y, (AgeGroup)j, InfectionType::H) +
                                               Po::get_from(y, (AgeGroup)j, InfectionType::U) +
                                               Po::get_from(y, (AgeGroup)j, InfectionType::R); // without died people
                               ScalarType divNj = 1.0 / Nj; // precompute 1.0/Nj
                               ScalarType Si    = Po::get_from(y, (AgeGroup)i, InfectionType::S);
                               ScalarType Cj    = Po::get_from(y, (AgeGroup)j, InfectionType::C);
                               ScalarType Ij    = Po::get_from(y, (AgeGroup)i, InfectionType::I);
                               return Si * cont_freq_eff * divNj * p.probabilities[i].get_infection_from_contact() *
                                      (p.probabilities[j].get_carrier_infectability() * Cj +
                                       p.probabilities[j].get_risk_from_symptomatic() * Ij);
                           });
        }

        // Ei to Ci
        model.add_flow(std::make_tuple((AgeGroup)i, InfectionType::E), std::make_tuple((AgeGroup)i, InfectionType::C),
                       [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           return Po::get_from(y, (AgeGroup)i, InfectionType::E) /
                                  (2 * p.times[i].get_serialinterval() - p.times[i].get_incubation());
                       });

        // Ci to Ii
        model.add_flow(std::make_tuple((AgeGroup)i, InfectionType::C), std::make_tuple((AgeGroup)i, InfectionType::I),
                       [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           double dummy_R3 = 0.5 / (p.times[i].get_incubation() - p.times[i].get_serialinterval());
                           double alpha    = p.probabilities[i].get_asymp_per_infectious();
                           return ((1 - alpha) * dummy_R3) * Po::get_from(y, (AgeGroup)i, InfectionType::C);
                       });

        // Ci to Ri
        model.add_flow(std::make_tuple((AgeGroup)i, InfectionType::C), std::make_tuple((AgeGroup)i, InfectionType::R),
                       [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           double alpha = p.probabilities[i].get_asymp_per_infectious();
                           return (alpha / p.times[i].get_infectious_asymp()) *
                                  Po::get_from(y, (AgeGroup)i, InfectionType::C);
                       });

        // Ii to Ri
        model.add_flow(std::make_tuple((AgeGroup)i, InfectionType::I), std::make_tuple((AgeGroup)i, InfectionType::R),
                       [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           return Po::get_from(y, (AgeGroup)i, InfectionType::I) *
                                  (1 - p.probabilities[i].get_hospitalized_per_infectious()) /
                                  p.times[i].get_infectious_mild();
                       });

        // Ii to Hi
        model.add_flow(std::make_tuple((AgeGroup)i, InfectionType::I), std::make_tuple((AgeGroup)i, InfectionType::H),
                       [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           return Po::get_from(y, (AgeGroup)i, InfectionType::I) *
                                  p.probabilities[i].get_hospitalized_per_infectious() /
                                  p.times[i].get_home_to_hospitalized();
                       });

        // Hi to Ui
        model.add_flow(std::make_tuple((AgeGroup)i, InfectionType::H), std::make_tuple((AgeGroup)i, InfectionType::U),
                       [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           return Po::get_from(y, (AgeGroup)i, InfectionType::H) *
                                  p.probabilities[i].get_icu_per_hospitalized() / p.times[i].get_hospitalized_to_icu();
                       });

        // Hi to Ri
        model.add_flow(std::make_tuple((AgeGroup)i, InfectionType::H), std::make_tuple((AgeGroup)i, InfectionType::R),
                       [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           return Po::get_from(y, (AgeGroup)i, InfectionType::H) *
                                  (1 - p.probabilities[i].get_icu_per_hospitalized()) /
                                  p.times[i].get_hospitalized_to_home();
                       });

        // Ui to Ri
        model.add_flow(std::make_tuple((AgeGroup)i, InfectionType::U), std::make_tuple((AgeGroup)i, InfectionType::R),
                       [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           return Po::get_from(y, (AgeGroup)i, InfectionType::U) *
                                  (1 - p.probabilities[i].get_dead_per_icu()) / p.times[i].get_icu_to_home();
                       });

        // Ui to Di
        model.add_flow(std::make_tuple((AgeGroup)i, InfectionType::U), std::make_tuple((AgeGroup)i, InfectionType::D),
                       [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                           return Po::get_from(y, (AgeGroup)i, InfectionType::U) *
                                  p.probabilities[i].get_dead_per_icu() / p.times[i].get_icu_to_dead();
                       });
    }

    return model;
}

} // namespace epi

#endif
