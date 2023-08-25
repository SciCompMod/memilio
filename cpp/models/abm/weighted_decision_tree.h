/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Sascha Korf
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

#include <vector>
#include <map>
#include <cmath>
#include <memory>
#include <algorithm>

struct TreeNode {
    int value;
};

struct TreeEdge {
    double time_length;
};

class weighted_decision_tree
{
public:
    weighted_decision_tree(std::vector<int> nodes, std::vector<std::tuple<int, int, double>> edges)
    {
        for (auto edge : edges) {
            auto parent_node = std::find(nodes.begin(), nodes.end(), std::get<0>(edge));
            auto child_node  = std::find(nodes.begin(), nodes.end(), std::get<1>(edge));
            build_tree(TreeNode{*parent_node}, TreeEdge{std::get<2>(edge)}, TreeNode{*child_node});
        }
    }

private:
    void build_tree(TreeNode node, TreeEdge edge, TreeNode child)
    {
        m_tree.insert({node, std::make_pair(edge, child)});
    }

    std::map<TreeNode, std::pair<TreeEdge, TreeNode>> m_tree;
};
