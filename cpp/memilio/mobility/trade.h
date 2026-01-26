/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Kilian Volmer
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
#ifndef MIO_TRADE_H
#define MIO_TRADE_H

#include "memilio/mobility/graph.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/compartments/feedback_simulation.h"
#include "memilio/geography/regions.h"
#include <limits>
#include <queue>

namespace mio
{
template <typename ScalarType>
struct Farm {
    ScalarType region;
    ScalarType grown_date;
};

template <typename ScalarType>
struct CompareFarms {
    bool operator()(const Farm<ScalarType>& f1, const Farm<ScalarType>& f2)
    {
        return f1.grown_date > f2.grown_date;
    }
};

template <typename ScalarType, class Graph>
class Trade
{
public:
    Trade(Graph& graph)
        : m_farm_queue{}
        , m_graph(graph)
    {
    }

    void add_farm(ScalarType region, ScalarType grown_date)
    {
        m_farm_queue.push(Farm<ScalarType>{region, grown_date});
    }

    const auto next_trade_time()
    {
        if (!m_farm_queue.empty()) {
            return m_farm_queue.top().grown_date;
        }
        else {
            return std::numeric_limits<ScalarType>::max();
        }
    }

    void next_trade()
    {
        m_farm_queue.pop();
    }

    void execute_trades(ScalarType current_time)
    {
        while (!m_farm_queue.empty() && m_farm_queue.top().grown_date <= current_time) {
            auto farm        = m_farm_queue.top();
            auto source      = farm.region;
            auto destination = get_destination();
            m_graph.nodes()[source].property.get_result().get_last_value()[0] -= 1;
            m_graph.nodes()[destination].property.get_result().get_last_value()[0] += 1;
            m_farm_queue.pop();
        }
    }

private:
    ScalarType get_destination()
    {
        auto options = std::vector<ScalarType>{};
        for (auto& node : m_graph.nodes()) {
            if (node.property.get_result().get_last_value()[0] > 0) {
                options.push_back(node.id);
            }
        }
        return options[0];
    }

    std::priority_queue<Farm<ScalarType>, std::vector<Farm<ScalarType>>, CompareFarms<ScalarType>> m_farm_queue;
    Graph& m_graph;
};

} // namespace mio
#endif //MIO_TRADE_H