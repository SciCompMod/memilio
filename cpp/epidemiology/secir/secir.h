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

enum class InfectionState
{
    Susceptible  = 0,
    Exposed      = 1,
    Carrier      = 2,
    Infected     = 3,
    Hospitalized = 4,
    ICU          = 5,
    Recovered    = 6,
    Dead         = 7,
    Count = 8
};

struct AgeGroup : public Index<AgeGroup> {
    AgeGroup(size_t val) : Index<AgeGroup>(val){}
};

class SecirModel : public CompartmentalModel<Populations<AgeGroup, InfectionState>, SecirParams>
{
    using Pa = SecirParams;
    using Po = Populations<AgeGroup, InfectionState>;

public:

    SecirModel(size_t num_agegroups)
     : CompartmentalModel<Populations<AgeGroup, InfectionState>, SecirParams>(Po({AgeGroup(num_agegroups), InfectionState::Count}),
                                                                             Pa(num_agegroups))
    {
#if !USE_DERIV_FUNC
        size_t n_agegroups = (size_t)AgeGroup::Count;
        for (size_t i = 0; i < n_agegroups; i++) {
            for (size_t j = 0; j < n_agegroups; j++) {

                // Si to Ei individually for each age group j
                this->add_flow(
                    std::make_tuple((AgeGroup)i, InfectionState::S), std::make_tuple((AgeGroup)i, InfectionState::E),
                    [i, j](Pa const& p, Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t) {

                        //symptomatic are less well quarantined when testing and tracing is overwhelmed so they infect more people
                        auto test_and_trace_required = (1 - params.probabilities[i].get_asymp_per_infectious()) * dummy_R3 * Po::get_from(pop, (AgeGroup)i, InfectionState::C);
                        auto risk_from_symptomatic   = smoother_cosine(
                            test_and_trace_required, params.get_test_and_trace_capacity(), params.get_test_and_trace_capacity() * 5,
                            params.probabilities[i].get_risk_from_symptomatic(),
                            params.probabilities[i].get_test_and_trace_max_risk_from_symptomatic());

                        // effective contact rate by contact rate between groups i and j and damping j
                        ScalarType season_val =
                            (1 + p.get_seasonality() * sin(3.141592653589793 *
                                                           (std::fmod((p.get_start_day() + t), 365.0) / 182.5 + 0.5)));
                        ScalarType cont_freq_eff = season_val * contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>(i),
                                                                                                static_cast<Eigen::Index>(j));
                        ScalarType Nj = Po::get_from(pop, (AgeGroup)j, InfectionState::S) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::E) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::C) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::I) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::H) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::U) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::R); // without died people
                        ScalarType divNj = 1.0 / Nj; // precompute 1.0/Nj
                        ScalarType Si    = Po::get_from(y, (AgeGroup)i, InfectionState::S);
                        ScalarType Cj    = Po::get_from(pop, (AgeGroup)j, InfectionState::C);
                        ScalarType Ij    = Po::get_from(pop, (AgeGroup)j, InfectionState::I);
                        return Si * cont_freq_eff * divNj * p.probabilities[i].get_infection_from_contact() *
                               (p.probabilities[j].get_carrier_infectability() * Cj +
                                risk_from_symptomatic * Ij);
                    });
            }

            // Ei to Ci
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::E),
                           std::make_tuple((AgeGroup)i, InfectionState::C),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionState::E) /
                                      (2 * p.times[i].get_serialinterval() - p.times[i].get_incubation());
                           });

            // Ci to Ii
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::C),
                           std::make_tuple((AgeGroup)i, InfectionState::I),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               double dummy_R3 = 0.5 / (p.times[i].get_incubation() - p.times[i].get_serialinterval());
                               double alpha    = p.probabilities[i].get_asymp_per_infectious();
                               return ((1 - alpha) * dummy_R3) * Po::get_from(y, (AgeGroup)i, InfectionState::C);
                           });

            // Ci to Ri
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionState::C), std::make_tuple((AgeGroup)i, InfectionState::R),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                    double alpha = p.probabilities[i].get_asymp_per_infectious();
                    return (alpha / p.times[i].get_infectious_asymp()) * Po::get_from(y, (AgeGroup)i, InfectionState::C);
                });

            // Ii to Ri
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::I),
                           std::make_tuple((AgeGroup)i, InfectionState::R),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionState::I) *
                                      (1 - p.probabilities[i].get_hospitalized_per_infectious()) /
                                      p.times[i].get_infectious_mild();
                           });

            // Ii to Hi
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionState::I), std::make_tuple((AgeGroup)i, InfectionState::H),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                    return Po::get_from(y, (AgeGroup)i, InfectionState::I) *
                           p.probabilities[i].get_hospitalized_per_infectious() / p.times[i].get_home_to_hospitalized();
                });

            // Hi to Ui
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::H),
                           std::make_tuple((AgeGroup)i, InfectionState::U),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               ScalarType icu_occupancy = 0;
                               for (size_t j = 0; j < (size_t)AgeGroup::Count; ++j) {
                                   icu_occupancy += Po::get_from(y, (AgeGroup)j, InfectionState::U);
                               }

                               ScalarType prob_hosp2icu =
                                   smoother_cosine(icu_occupancy, 0.90 * p.get_icu_capacity(), p.get_icu_capacity(),
                                                   p.probabilities[i].get_icu_per_hospitalized(), 0);

                               return Po::get_from(y, (AgeGroup)i, InfectionState::H) * prob_hosp2icu /
                                      p.times[i].get_hospitalized_to_icu();
                           });

            // Hi to Di
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionState::H), std::make_tuple((AgeGroup)i, InfectionState::D),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                    ScalarType icu_occupancy = 0;
                    for (size_t j = 0; j < (size_t)AgeGroup::Count; ++j) {
                        icu_occupancy += Po::get_from(y, (AgeGroup)j, InfectionState::U);
                    }

                    ScalarType prob_hosp2icu =
                        smoother_cosine(icu_occupancy, 0.90 * p.get_icu_capacity(), p.get_icu_capacity(),
                                        p.probabilities[i].get_icu_per_hospitalized(), 0);
                    ScalarType prob_hosp2dead = p.probabilities[i].get_icu_per_hospitalized() - prob_hosp2icu;

                    return Po::get_from(y, (AgeGroup)i, InfectionState::H) * prob_hosp2dead /
                           p.times[i].get_hospitalized_to_icu();
                });

            // Hi to Ri
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionState::H), std::make_tuple((AgeGroup)i, InfectionState::R),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                    return Po::get_from(y, (AgeGroup)i, InfectionState::H) *
                           (1 - p.probabilities[i].get_icu_per_hospitalized()) / p.times[i].get_hospitalized_to_home();
                });

            // Ui to Ri
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::U),
                           std::make_tuple((AgeGroup)i, InfectionState::R),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionState::U) *
                                      (1 - p.probabilities[i].get_dead_per_icu()) / p.times[i].get_icu_to_home();
                           });

            // Ui to Di
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::U),
                           std::make_tuple((AgeGroup)i, InfectionState::D),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionState::U) *
                                      p.probabilities[i].get_dead_per_icu() / p.times[i].get_icu_to_dead();
                           });
        }
#endif
    }

#if USE_DERIV_FUNC

void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop,
                     Eigen::Ref<const Eigen::VectorXd> y, double t,
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

    ContactMatrixGroup const& contact_matrix = params.get_contact_patterns();

    for (size_t i = 0; i < n_agegroups; i++) {

        size_t Si = this->populations.get_flat_index({AgeGroup(i), InfectionState::Susceptible});
        size_t Ei = this->populations.get_flat_index({AgeGroup(i), InfectionState::Exposed});
        size_t Ci = this->populations.get_flat_index({AgeGroup(i), InfectionState::Carrier});
        size_t Ii = this->populations.get_flat_index({AgeGroup(i), InfectionState::Infected});
        size_t Hi = this->populations.get_flat_index({AgeGroup(i), InfectionState::Hospitalized});
        size_t Ui = this->populations.get_flat_index({AgeGroup(i), InfectionState::ICU});
        size_t Ri = this->populations.get_flat_index({AgeGroup(i), InfectionState::Recovered});
        size_t Di = this->populations.get_flat_index({AgeGroup(i), InfectionState::Dead});

        dydt[Si] = 0;
        dydt[Ei] = 0;

        double dummy_R2 =
            1.0 / (2 * params.times[i].get_serialinterval() - params.times[i].get_incubation()); // R2 = 1/(2SI-TINC)
        double dummy_R3 =
            0.5 / (params.times[i].get_incubation() - params.times[i].get_serialinterval()); // R3 = 1/(2(TINC-SI))

        //symptomatic are less well quarantined when testing and tracing is overwhelmed so they infect more people
        auto test_and_trace_required = (1 - params.probabilities[i].get_asymp_per_infectious()) * dummy_R3 * this->populations.get_from(pop, {AgeGroup(i), InfectionState::Carrier});
        auto risk_from_symptomatic   = smoother_cosine(
            test_and_trace_required, params.get_test_and_trace_capacity(), params.get_test_and_trace_capacity() * 5,
            params.probabilities[i].get_risk_from_symptomatic(),
            params.probabilities[i].get_test_and_trace_max_risk_from_symptomatic());

        double icu_occupancy = 0;

        for (size_t j = 0; j < n_agegroups; j++) {
            size_t Sj = this->populations.get_flat_index({AgeGroup(j), InfectionState::Susceptible});
            size_t Ej = this->populations.get_flat_index({AgeGroup(j), InfectionState::Exposed});
            size_t Cj = this->populations.get_flat_index({AgeGroup(j), InfectionState::Carrier});
            size_t Ij = this->populations.get_flat_index({AgeGroup(j), InfectionState::Infected});
            size_t Hj = this->populations.get_flat_index({AgeGroup(j), InfectionState::Hospitalized});
            size_t Uj = this->populations.get_flat_index({AgeGroup(j), InfectionState::ICU});
            size_t Rj = this->populations.get_flat_index({AgeGroup(j), InfectionState::Recovered});

            // effective contact rate by contact rate between groups i and j and damping j
            double season_val =
                (1 + params.get_seasonality() *
                         sin(3.141592653589793 * (std::fmod((params.get_start_day() + t), 365.0) / 182.5 + 0.5)));
            double cont_freq_eff = season_val * contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>(i),
                                                                                static_cast<Eigen::Index>(j));
            double Nj      = pop[Sj] + pop[Ej] + pop[Cj] + pop[Ij] + pop[Hj] + pop[Uj] + pop[Rj]; // without died people
            double divNj   = 1.0 / Nj; // precompute 1.0/Nj
            double dummy_S = y[Si] * cont_freq_eff * divNj * params.probabilities[i].get_infection_from_contact() *
                             (params.probabilities[j].get_carrier_infectability() * pop[Cj] +
                              risk_from_symptomatic * pop[Ij]);

            dydt[Si] -= dummy_S; // -R1*(C+beta*I)*S/N0
            dydt[Ei] += dummy_S; // R1*(C+beta*I)*S/N0-R2*E

            icu_occupancy += y[Uj];
        }

        // ICU capacity shortage is close
        double prob_hosp2icu =
            smoother_cosine(icu_occupancy, 0.90 * params.get_icu_capacity(), params.get_icu_capacity(),
                            params.probabilities[i].get_icu_per_hospitalized(), 0);

        double prob_hosp2dead = params.probabilities[i].get_icu_per_hospitalized() - prob_hosp2icu;


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

#endif // USE_DERIV_FUNC

};

} // namespace epi

#endif
