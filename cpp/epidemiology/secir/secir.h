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

enum class AgeGroup6
{
    Group0,
    Group1,
    Group2,
    Group3,
    Group4,
    Group5,
    Count = 6
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

#if USE_DEPRECATED_SECIR_DERIV_FUNC

void get_derivatives(Eigen::Ref<const Eigen::VectorXd> y, double t,
                           Eigen::Ref<Eigen::VectorXd> dydt) const override
{
    // alpha  // percentage of asymptomatic cases
    // beta // risk of infection from the infected symptomatic patients
    // rho   // hospitalized per infectious
    // theta // icu per hospitalized
    // delta  // deaths per ICUs
    // 0: S,      1: E,     2: C,     3: I,     4: H,     5: U,     6: R,     7: D
    auto& params = this->parameters;
    size_t n_agegroups = params.get_num_groups();

    ContactFrequencyMatrix const& cont_freq_matrix = params.get_contact_patterns();

    for (size_t i = 0; i < n_agegroups; i++) {

        size_t Si = Po::get_flat_index((AgeGroup)i, InfectionType::S);
        size_t Ei = Po::get_flat_index((AgeGroup)i, InfectionType::E);
        size_t Ci = Po::get_flat_index((AgeGroup)i, InfectionType::C);
        size_t Ii = Po::get_flat_index((AgeGroup)i, InfectionType::I);
        size_t Hi = Po::get_flat_index((AgeGroup)i, InfectionType::H);
        size_t Ui = Po::get_flat_index((AgeGroup)i, InfectionType::U);
        size_t Ri = Po::get_flat_index((AgeGroup)i, InfectionType::R);
        size_t Di = Po::get_flat_index((AgeGroup)i, InfectionType::D);

        dydt[Si] = 0;
        dydt[Ei] = 0;

        double icu_occupancy = 0;

        for (size_t j = 0; j < n_agegroups; j++) {
            size_t Sj = Po::get_flat_index((AgeGroup)j, InfectionType::S);
            size_t Ej = Po::get_flat_index((AgeGroup)j, InfectionType::E);
            size_t Cj = Po::get_flat_index((AgeGroup)j, InfectionType::C);
            size_t Ij = Po::get_flat_index((AgeGroup)j, InfectionType::I);
            size_t Hj = Po::get_flat_index((AgeGroup)j, InfectionType::H);
            size_t Uj = Po::get_flat_index((AgeGroup)j, InfectionType::U);
            size_t Rj = Po::get_flat_index((AgeGroup)j, InfectionType::R);

            // effective contact rate by contact rate between groups i and j and damping j
            double season_val =
                (1 + params.get_seasonality() *
                         sin(3.141592653589793 * (std::fmod((params.get_start_day() + t), 365.0) / 182.5 + 0.5)));
            double cont_freq_eff = // get effective contact rate between i and j
                season_val * cont_freq_matrix.get_cont_freq(static_cast<int>(i), static_cast<int>(j)) *
                cont_freq_matrix.get_dampings(static_cast<int>(i), static_cast<int>(j)).get_factor(t);
            double Nj      = y[Sj] + y[Ej] + y[Cj] + y[Ij] + y[Hj] + y[Uj] + y[Rj]; // without died people
            double divNj   = 1.0 / Nj; // precompute 1.0/Nj
            double dummy_S = y[Si] * cont_freq_eff * divNj * params.probabilities[i].get_infection_from_contact() *
                             (params.probabilities[j].get_carrier_infectability() * y[Cj] +
                              params.probabilities[j].get_risk_from_symptomatic() * y[Ij]);

            dydt[Si] -= dummy_S; // -R1*(C+beta*I)*S/N0
            dydt[Ei] += dummy_S; // R1*(C+beta*I)*S/N0-R2*E

            icu_occupancy += y[Uj];
        }

        // ICU capacity shortage is close
        double prob_hosp2icu =
            smoother_cosine(icu_occupancy, 0.90 * params.get_icu_capacity(), params.get_icu_capacity(),
                            params.probabilities[i].get_icu_per_hospitalized(), 0);

        double prob_hosp2dead = params.probabilities[i].get_icu_per_hospitalized() - prob_hosp2icu;

        double dummy_R2 =
            1.0 / (2 * params.times[i].get_serialinterval() - params.times[i].get_incubation()); // R2 = 1/(2SI-TINC)
        double dummy_R3 =
            0.5 / (params.times[i].get_incubation() - params.times[i].get_serialinterval()); // R3 = 1/(2(TINC-SI))

        dydt[Ei] -= dummy_R2 * y[Ei]; // only exchange of E and C done here
        dydt[Ci] = dummy_R2 * y[Ei] -
                   ((1 - params.probabilities[i].get_asymp_per_infectious()) * dummy_R3 +
                    params.probabilities[i].get_asymp_per_infectious() / params.times[i].get_infectious_asymp()) *
                       y[Ci];
        dydt[Ii] =
            (1 - params.probabilities[i].get_asymp_per_infectious()) * dummy_R3 * y[Ci] -
            ((1 - params.probabilities[i].get_hospitalized_per_infectious()) / params.times[i].get_infectious_mild() +
             params.probabilities[i].get_hospitalized_per_infectious() / params.times[i].get_home_to_hospitalized()) *
                y[Ii];
        dydt[Hi] =
            params.probabilities[i].get_hospitalized_per_infectious() / params.times[i].get_home_to_hospitalized() *
                y[Ii] -
            ((1 - params.probabilities[i].get_icu_per_hospitalized()) / params.times[i].get_hospitalized_to_home() +
             params.probabilities[i].get_icu_per_hospitalized() / params.times[i].get_hospitalized_to_icu()) *
                y[Hi];
        dydt[Ui] = -((1 - params.probabilities[i].get_dead_per_icu()) / params.times[i].get_icu_to_home() +
                     params.probabilities[i].get_dead_per_icu() / params.times[i].get_icu_to_dead()) *
                   y[Ui];
        // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
        dydt[Ui] += prob_hosp2icu / params.times[i].get_hospitalized_to_icu() * y[Hi];

        dydt[Ri] = params.probabilities[i].get_asymp_per_infectious() / params.times[i].get_infectious_asymp() * y[Ci] +
                   (1 - params.probabilities[i].get_hospitalized_per_infectious()) /
                       params.times[i].get_infectious_mild() * y[Ii] +
                   (1 - params.probabilities[i].get_icu_per_hospitalized()) /
                       params.times[i].get_hospitalized_to_home() * y[Hi] +
                   (1 - params.probabilities[i].get_dead_per_icu()) / params.times[i].get_icu_to_home() * y[Ui];

        dydt[Di] = params.probabilities[i].get_dead_per_icu() / params.times[i].get_icu_to_dead() * y[Ui];
        // add potential, additional deaths due to ICU overflow
        dydt[Di] += prob_hosp2dead / params.times[i].get_hospitalized_to_icu() * y[Hi];
    }
}

#endif // USE_DEPRECATED_SECIR_DERIV_FUNC

};

} // namespace epi

#endif
