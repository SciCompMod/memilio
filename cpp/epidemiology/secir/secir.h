#ifndef SECIR_H
#define SECIR_H

#include "epidemiology/model/compartmentalmodel.h"
#include "epidemiology/model/populations.h"
#include "epidemiology/secir/secir_params.h"
#include "epidemiology/math/smoother.h"

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

enum class AgeGroup2
{
    Group0,
    Group1,
    Count = 2
};

enum class AgeGroup3
{
    Group0,
    Group1,
    Group2,
    Count = 3
};

enum class AgeGroup8
{
    Group0,
    Group1,
    Group2,
    Group3,
    Group4,
    Group5,
    Group6,
    Group7,
    Count = 8
};

//TODO: Once the secir_params have been replaced with a
// more abstract compile time map, this is not needed anymore
enum class AgeGroup1
{
    Group0,
    Count = 1
};

template <class AG>
class SecirModel : public CompartmentalModel<Populations<AG, InfectionType>, SecirParams<(size_t)AG::Count>>
{
    using Pa = SecirParams<(size_t)AG::Count>;
    using Po = Populations<AG, InfectionType>;

public:
    using AgeGroup = AG;

    SecirModel()
    {
        size_t n_agegroups = (size_t)AgeGroup::Count;
        for (size_t i = 0; i < n_agegroups; i++) {
            for (size_t j = 0; j < n_agegroups; j++) {

                // Si to Ei individually for each age group j
                this->add_flow(
                    std::make_tuple((AgeGroup)i, InfectionType::S), std::make_tuple((AgeGroup)i, InfectionType::E),
                    [i, j](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double t) {
                        // effective contact rate by contact rate between groups i and j and damping j
                        ScalarType season_val =
                            (1 + p.get_seasonality() * sin(3.141592653589793 *
                                                           (std::fmod((p.get_start_day() + t), 365.0) / 182.5 + 0.5)));
                        ScalarType cont_freq_eff =
                            season_val * p.get_contact_patterns().get_cont_freq_mat().get_cont_freq((int)i, (int)j) *
                            p.get_contact_patterns().get_cont_freq_mat().get_dampings((int)i, (int)j).get_factor(t);
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
                        ScalarType Ij    = Po::get_from(y, (AgeGroup)j, InfectionType::I);
                        return Si * cont_freq_eff * divNj * p.probabilities[i].get_infection_from_contact() *
                               (p.probabilities[j].get_carrier_infectability() * Cj +
                                p.probabilities[j].get_risk_from_symptomatic() * Ij);
                    });
            }

            // Ei to Ci
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionType::E),
                           std::make_tuple((AgeGroup)i, InfectionType::C),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionType::E) /
                                      (2 * p.times[i].get_serialinterval() - p.times[i].get_incubation());
                           });

            // Ci to Ii
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionType::C),
                           std::make_tuple((AgeGroup)i, InfectionType::I),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               double dummy_R3 = 0.5 / (p.times[i].get_incubation() - p.times[i].get_serialinterval());
                               double alpha    = p.probabilities[i].get_asymp_per_infectious();
                               return ((1 - alpha) * dummy_R3) * Po::get_from(y, (AgeGroup)i, InfectionType::C);
                           });

            // Ci to Ri
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionType::C), std::make_tuple((AgeGroup)i, InfectionType::R),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                    double alpha = p.probabilities[i].get_asymp_per_infectious();
                    return (alpha / p.times[i].get_infectious_asymp()) * Po::get_from(y, (AgeGroup)i, InfectionType::C);
                });

            // Ii to Ri
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionType::I),
                           std::make_tuple((AgeGroup)i, InfectionType::R),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionType::I) *
                                      (1 - p.probabilities[i].get_hospitalized_per_infectious()) /
                                      p.times[i].get_infectious_mild();
                           });

            // Ii to Hi
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionType::I), std::make_tuple((AgeGroup)i, InfectionType::H),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                    return Po::get_from(y, (AgeGroup)i, InfectionType::I) *
                           p.probabilities[i].get_hospitalized_per_infectious() / p.times[i].get_home_to_hospitalized();
                });

            // Hi to Ui
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionType::H),
                           std::make_tuple((AgeGroup)i, InfectionType::U),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               ScalarType icu_occupancy = 0;
                               for (size_t j = 0; j < (size_t)AgeGroup::Count; ++j) {
                                   icu_occupancy += Po::get_from(y, (AgeGroup)j, InfectionType::U);
                               }

                               ScalarType prob_hosp2icu =
                                   smoother_cosine(icu_occupancy, 0.90 * p.get_icu_capacity(), p.get_icu_capacity(),
                                                   p.probabilities[i].get_icu_per_hospitalized(), 0);

                               return Po::get_from(y, (AgeGroup)i, InfectionType::H) * prob_hosp2icu /
                                      p.times[i].get_hospitalized_to_icu();
                           });

            // Hi to Di
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionType::H), std::make_tuple((AgeGroup)i, InfectionType::D),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                    ScalarType icu_occupancy = 0;
                    for (size_t j = 0; j < (size_t)AgeGroup::Count; ++j) {
                        icu_occupancy += Po::get_from(y, (AgeGroup)j, InfectionType::U);
                    }

                    ScalarType prob_hosp2icu =
                        smoother_cosine(icu_occupancy, 0.90 * p.get_icu_capacity(), p.get_icu_capacity(),
                                        p.probabilities[i].get_icu_per_hospitalized(), 0);
                    ScalarType prob_hosp2dead = p.probabilities[i].get_icu_per_hospitalized() - prob_hosp2icu;

                    return Po::get_from(y, (AgeGroup)i, InfectionType::H) * prob_hosp2dead /
                           p.times[i].get_hospitalized_to_icu();
                });

            // Hi to Ri
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionType::H), std::make_tuple((AgeGroup)i, InfectionType::R),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                    return Po::get_from(y, (AgeGroup)i, InfectionType::H) *
                           (1 - p.probabilities[i].get_icu_per_hospitalized()) / p.times[i].get_hospitalized_to_home();
                });

            // Ui to Ri
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionType::U),
                           std::make_tuple((AgeGroup)i, InfectionType::R),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionType::U) *
                                      (1 - p.probabilities[i].get_dead_per_icu()) / p.times[i].get_icu_to_home();
                           });

            // Ui to Di
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionType::U),
                           std::make_tuple((AgeGroup)i, InfectionType::D),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionType::U) *
                                      p.probabilities[i].get_dead_per_icu() / p.times[i].get_icu_to_dead();
                           });
        }
    }
};

} // namespace epi

#endif
