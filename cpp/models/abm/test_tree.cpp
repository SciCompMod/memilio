#include "weighted_decision_tree.h"
#include <vector>

int main(int argc, char const* argv[])
{
    auto nodes         = std::vector<int>({1, 4, 10});
    auto edges         = std::vector<std::tuple<int, int, double>>({{1, 4, 0.5}, {4, 10, 0.5}});
    auto decision_tree = weighted_decision_tree(nodes, edges);
    return 0;
}
