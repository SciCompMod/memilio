/* 
* Copyright (C) 2020-2024 German Aerospace Center (DLR-SC)
*
* Authors: René Schmieding, Julia Bicker
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

#ifndef QUAD_WELL_H
#define QUAD_WELL_H

#include "memilio/math/eigen.h"
#include "d_abm/parameters.h"
#include "memilio/utils/random_number_generator.h"
#include <string>
#include <vector>

using Position = Eigen::Vector2d;

inline size_t well_index(const Position& p)
{
    // 0|1
    // -+-
    // 2|3
    return (p.x() >= 0) + 2 * (p.y() < 0);
}

template <class InfectionState>
class QuadWellModel
{

public:
    using Status = InfectionState;

    struct Agent {
        Position position;
        Status status;
    };

    QuadWellModel(const std::vector<Agent>& agents, const typename mio::dabm::AdoptionRates<Status>::Type& rates,
                  double contact_radius = 0.4, double sigma = 0.4, std::vector<Status> non_moving_states = {})
        : populations(agents)
        , m_contact_radius(contact_radius)
        , m_sigma(sigma)
        , m_non_moving_states(non_moving_states)
        , m_number_transitions(static_cast<size_t>(Status::Count), Eigen::MatrixXd::Zero(4, 4))
    {
        for (auto& agent : populations) {
            mio::unused(agent);
            assert(is_in_domain(agent.position));
        }
        for (auto& r : rates) {
            m_adoption_rates.emplace(std::forward_as_tuple(r.region, r.from, r.to), r);
        }
    }

    inline static constexpr void adopt(Agent& agent, const Status& new_status)
    {
        agent.status = new_status;
    }

    double adoption_rate(const Agent& agent, const Status& new_status) const
    {
        double rate = 0;
        // get the correct adoption rate
        const size_t well = well_index(agent.position);
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
                            if (contact.status == adoption_rate.influences[i]) {
                                influences += adoption_rate.factors[i];
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

    void move(const double /*t*/, const double dt, Agent& agent)
    {
        const auto old_well = well_index(agent.position);
        if (std::find(m_non_moving_states.begin(), m_non_moving_states.end(), agent.status) ==
                m_non_moving_states.end() &&
            std::find(m_non_moving_regions.begin(), m_non_moving_regions.end(), old_well) ==
                m_non_moving_regions.end()) {
            Position p = {mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(m_rng, 0.0, 1.0),
                          mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(m_rng, 0.0, 1.0)};

            agent.position      = agent.position - dt * grad_U(agent.position) + (m_sigma * std::sqrt(dt)) * p;
            const auto new_well = well_index(agent.position);
            if (old_well != new_well) {
                m_number_transitions[static_cast<size_t>(agent.status)](old_well, new_well)++;
            }
        }
        //else{agent has non-moving status or region}
    }

    Eigen::VectorXd time_point() const
    {
        Eigen::VectorXd val = Eigen::VectorXd::Zero(4 * static_cast<size_t>(Status::Count));
        for (auto& agent : populations) {
            // split population into the wells given by grad_U
            auto position =
                static_cast<size_t>(agent.status) + well_index(agent.position) * static_cast<size_t>(Status::Count);
            val[position] += 1;
        }
        return val;
    }

    const std::vector<Eigen::MatrixXd>& number_transitions() const
    {
        return m_number_transitions;
    }

    std::vector<Eigen::MatrixXd>& number_transitions()
    {
        return m_number_transitions;
    }

    std::map<std::tuple<mio::dabm::Region, Status, Status>, mio::dabm::AdoptionRate<Status>>& get_adoption_rates()
    {
        return m_adoption_rates;
    }

    void set_non_moving_regions(std::vector<size_t> non_moving_regions)
    {
        m_non_moving_regions = non_moving_regions;
    }

    mio::RandomNumberGenerator& get_rng()
    {
        return m_rng;
    }

    std::vector<Agent> populations;

private:
    static Position grad_U(const Position x)
    {
        // U is a quad well potential
        // U(x0,x1) = (x0^2 - 1)^2 + (x1^2 - 1)^2
        return {4 * x[0] * (x[0] * x[0] - 1), 4 * x[1] * (x[1] * x[1] - 1)};
    }

    bool is_contact(const Agent& agent, const Agent& contact) const
    {
        //      test if contact is in the contact radius                     and test if agent and contact are different objects
        return (agent.position - contact.position).norm() < m_contact_radius && (&agent != &contact) &&
               well_index(agent.position) == well_index(contact.position);
    }

    bool is_in_domain(const Position& p) const
    {
        // restrict domain to [-2, 2]^2 where "escaping" is impossible, i.e. it holds x <= grad_U(x) for dt <= 0.1
        return -2 <= p[0] && p[0] <= 2 && -2 <= p[1] && p[1] <= 2;
    }

    std::map<std::tuple<mio::dabm::Region, Status, Status>, mio::dabm::AdoptionRate<Status>> m_adoption_rates;
    double m_contact_radius;
    double m_sigma;
    std::vector<Status> m_non_moving_states;
    std::vector<size_t> m_non_moving_regions{};
    std::vector<Eigen::MatrixXd> m_number_transitions;
    mio::RandomNumberGenerator m_rng;
};

#endif
