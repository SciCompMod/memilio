/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Julia Bicker, René Schmieding
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

#ifndef MIO_D_ABM_SINGLE_WELL_H
#define MIO_D_ABM_SINGLE_WELL_H

#include "memilio/config.h"
#include "memilio/math/eigen.h"
#include "d_abm/parameters.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/epidemiology/adoption_rate.h"

/**
 * @brief A Sampler for sampling a position for the singlewell potential F(x,y) = (x^4 + y^4)/2, see SingleWell.
 */
class SWPositionSampler
{
public:
    using Position = Eigen::Vector2d;

    /**
     * @brief Create a sampler.
     * @param[in] bottom_left Coordinates of the bottom left corner of the range in which should be samples.
     * @param[in] top_right Coordinates of the top right corner of the range in which should be samples.
     * @param[in] margin Margin defining the distance that should be kept from the borders defined by bottom_left and top_right when sampling a position.
     */
    SWPositionSampler(const Position& bottom_left, const Position& top_right, ScalarType margin)
    {
        auto assign_range = [&](Position range_x, Position range_y) {
            range_x += Position{margin, -margin};
            range_y += Position{margin, -margin};
            m_range = {range_x, range_y};
        };

        assign_range({bottom_left.x(), top_right.x()}, {bottom_left.y(), top_right.y()});
    }

    /**
     * @brief Sampling a position within [bottom_left.x + margin, bottom_left.y + margin] x [top_right.x - margin, top_right.y - margin]
     */
    Position operator()() const
    {
        return {mio::UniformDistribution<ScalarType>::get_instance()(mio::thread_local_rng(), m_range.first[0],
                                                                     m_range.first[1]),
                mio::UniformDistribution<ScalarType>::get_instance()(mio::thread_local_rng(), m_range.second[0],
                                                                     m_range.second[1])};
    }

private:
    // stores pairs of (x-range, y-range)
    std::pair<Position, Position> m_range;
};

/**
 * @brief Implementation of diffusive ABM, see dabm::Model. 
 * This implementation defines a diffusion process for the potential F(x,y) = (x^4 + y^4)/2.
 * @tparam InfectionState An infection state enum.
 */
template <class InfectionState>
class SingleWell
{
public:
    using Position = Eigen::Vector2d;
    using Status   = InfectionState;

    struct Agent {
        Position position;
        Status status;
    };

    /**
     * @brief Set up a diffusive ABM using the quadwell potential F(x,y) = (x^4 + y^4)/2.
     * @param[in] agents A vector of Agent%s representing the population.
     * @param[in] rates AdoptionRate%s defining InfectionState adoptions, see mio::AdoptionRate.
     * @param[in] contact_radius Contact radius for second-order adoptions.
     * @param[in] sigma Noise term for the diffusion process.
     * @param[in] non_moving_state InfectionStates that are excluded from movement e.g. Dead.
     */
    SingleWell(const std::vector<Agent>& agents, const std::vector<mio::AdoptionRate<ScalarType, Status>>& rates,
               ScalarType contact_radius = 0.4, ScalarType sigma = 0.4, std::vector<Status> non_moving_states = {})
        : populations(agents)
        , m_contact_radius(contact_radius)
        , m_sigma(sigma)
        , m_non_moving_states(non_moving_states)
    {
        for (auto& agent : populations) {
            mio::unused(agent);
            assert(is_in_domain(agent.position));
        }
        for (auto& r : rates) {
            m_adoption_rates.emplace(std::forward_as_tuple(r.region, r.from, r.to), r);
        }
    }

    /**
     * @brief Perform infection state adoption for an Agent.
     * @param[in, out] agent Agent whose infection state is changed.
     * @param[in] new_status Agent's new infection state.
     */
    inline static constexpr void adopt(Agent& agent, const Status& new_status)
    {
        agent.status = new_status;
    }

    /**
     * @brief Calculate adoption rate for an Agent.
     * @param[in] agent Agent for whom the adoption rate is calculated.
     * @param[in] new_status Target infection state of the adoption rate, see mio::AdoptionRate.to.
     * @return Value of agent-dependent AdoptionRate.
     */
    ScalarType adoption_rate(const Agent& agent, const Status& new_status) const
    {
        ScalarType rate = 0;
        // get the correct adoption rate
        const size_t well = 0;
        auto map_itr      = m_adoption_rates.find({well, agent.status, new_status});
        if (map_itr != m_adoption_rates.end()) {
            const auto& adoption_rate = map_itr->second;
            // calculate the current rate, depending on order
            if (adoption_rate.influences.size() == 0) { // first order adoption
                // contact independant rate
                rate = adoption_rate.factor;
            }
            else { // second order adoption
                // accumulate rate per contact with a status in influences
                size_t num_contacts   = 0;
                ScalarType influences = 0;
                for (auto& contact : populations) {
                    // check if contact is indeed a contact
                    if (is_contact(agent, contact)) {
                        num_contacts++;
                        for (size_t i = 0; i < adoption_rate.influences.size(); i++) {
                            if (contact.status == adoption_rate.influences[i].status) {
                                influences += adoption_rate.influences[i].factor;
                            }
                        }
                    }
                }
                // rate = factor * "concentration of contacts with status new_status"
                if (num_contacts > 0) {
                    rate = adoption_rate.factor * (influences / num_contacts);
                }
            }
        }
        // else: no adoption from agent.status to new_status exist
        return rate;
    }

    /**
     * @brief Perform one integration step of the diffusion process for a given Agent using the Euler-Maruyama method.
     * @param[in] dt Step size.
     * @param[in] agent Agent to be moved.
     */
    void move(const ScalarType /*t*/, const ScalarType dt, Agent& agent)
    {
        if (std::find(m_non_moving_states.begin(), m_non_moving_states.end(), agent.status) ==
            m_non_moving_states.end()) {
            Position p = {
                mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(m_rng, 0.0, 1.0),
                mio::DistributionAdapter<std::normal_distribution<ScalarType>>::get_instance()(m_rng, 0.0, 1.0)};

            agent.position = agent.position - dt * grad_U(agent.position) + (m_sigma * std::sqrt(dt)) * p;
        }
        //else{agent has non-moving status}
    }

    /**
     * @brief Calculate the current system state i.e. the populations for each infection state.
     * @return Vector containing the number of agents per infection state.
     */
    Eigen::VectorX<ScalarType> time_point() const
    {
        Eigen::VectorX<ScalarType> val = Eigen::VectorX<ScalarType>::Zero(static_cast<size_t>(Status::Count));
        for (auto& agent : populations) {
            val[static_cast<size_t>(agent.status)] += 1;
        }
        return val;
    }

    std::map<std::tuple<mio::regions::Region, Status, Status>, mio::AdoptionRate<ScalarType, Status>>&
    get_adoption_rates()
    {
        return m_adoption_rates;
    }

    /**
    * Get the RandomNumberGenerator used by this Model for random events.
    * @return The random number generator.
    */
    mio::RandomNumberGenerator& get_rng()
    {
        return m_rng;
    }

    std::vector<Agent> populations; ///< Vector containing all Agent%s in the model.

private:
    /**
     * @brief Calculate the gradient of the potential at a given position.
     * @param[in] p Position at which the gradient is evaluated.
     * @return Value of potential gradient.
     */
    static Position grad_U(const Position p)
    {
        // U is a single well potential
        // U(x0,x1) = (x0^4 + x1^4)/2
        return {2 * p[0] * p[0] * p[0], 2 * p[1] * p[1] * p[1]};
    }

    /**
     * @brief Evaluate whether two agents are within their contact radius.
     * @param[in] agent Agent whose position specifies the contact area.
     * @param[in] contact Potential contact of agent.
     * @return Boolean specifying whether are within each others contact radius.
     */
    bool is_contact(const Agent& agent, const Agent& contact) const
    {
        //      test if contact is in the contact radius                     and test if agent and contact are different objects
        return (agent.position - contact.position).norm() < m_contact_radius && (&agent != &contact);
    }

    /** 
     * @brief Restrict domain to [-2, 2]^2 where "escaping" is impossible.
     * @param[in] p Position to check.
     * @return Boolean specifying whether p is in [-2, 2]^2.
    */
    bool is_in_domain(const Position& p, const ScalarType lower_domain_border = -2,
                      const ScalarType upper_domain_border = 2) const
    {
        // restrict domain to [lower_domain_border, upper_domain_border]^2 where "escaping" is impossible, i.e. it holds x <= grad_U(x) for dt <= 0.1
        return lower_domain_border <= p[0] && p[0] <= upper_domain_border && lower_domain_border <= p[1] &&
               p[1] <= upper_domain_border;
    }

    std::map<std::tuple<mio::regions::Region, Status, Status>, mio::AdoptionRate<ScalarType, Status>>
        m_adoption_rates; ///< Map of AdoptionRates according to their region index and their from -> to infection states.
    ScalarType m_contact_radius; ///< Agents' interaction radius. Within this radius agents are considered as contacts.
    ScalarType m_sigma; ///< Noise term of the diffusion process.
    std::vector<Status> m_non_moving_states; ///< Infection states within which agents do not change their location.
    mio::RandomNumberGenerator m_rng; ///< Model's random number generator.
};

#endif //MIO_D_ABM_SINGLE_WELL_H
