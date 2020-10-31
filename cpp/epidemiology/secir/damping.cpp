#include "epidemiology/secir/damping.h"
#include "epidemiology/utils/stl_util.h"
#include "epidemiology/math/smoother.h"

#include <algorithm>
#include <cassert>
#include <stdio.h>
#include <cmath>

namespace epi
{

void Dampings::finalize() const
{
    using std::get;

    if (m_accumulated_dampings_cached.empty()) {
        m_accumulated_dampings_cached.emplace_back(Eigen::MatrixXd::Zero(m_num_groups, m_num_groups),
                                                   SimulationTime(std::numeric_limits<double>::lowest()));

        std::vector<std::tuple<std::reference_wrapper<const Eigen::MatrixXd>, DampingLevel, DampingType>>
            active_by_type;
        std::vector<std::tuple<Eigen::MatrixXd, DampingLevel>> sum_by_level;
        for (auto& damping : m_dampings) {
            update_active_dampings(damping, active_by_type, sum_by_level);
            m_accumulated_dampings_cached.emplace_back(inclusive_exclusive_sum(sum_by_level),
                                                       get<SimulationTime>(damping));
            assert((get<Eigen::MatrixXd>(m_accumulated_dampings_cached.back()).array() <= 1).all() &&
                   "accumulated damping must be below 1.");
        }

        m_accumulated_dampings_cached.emplace_back(get<Eigen::MatrixXd>(m_accumulated_dampings_cached.back()),
                                                   SimulationTime(std::numeric_limits<double>::max()));
    }
}

void Dampings::add_(const Damping& damping)
{
    assert(damping.get_coeffs().rows() == m_num_groups && damping.get_coeffs().cols() == m_num_groups);
    insert_sorted_replace(m_dampings, damping, [](auto& tup1, auto& tup2) {
        return double(std::get<SimulationTime>(tup1)) < double(std::get<SimulationTime>(tup2));
    });
    m_accumulated_dampings_cached.clear();
}

void Dampings::update_active_dampings(
    const Damping& damping,
    std::vector<std::tuple<std::reference_wrapper<const Eigen::MatrixXd>, DampingLevel, DampingType>>& active_by_type,
    std::vector<std::tuple<Eigen::MatrixXd, DampingLevel>>& sum_by_level)
{
    using std::get;

    const int MatrixIdx = 0;

    auto iter_active_same_type = std::find_if(active_by_type.begin(), active_by_type.end(), [&damping](auto& active) {
        return get<DampingLevel>(active) == get<DampingLevel>(damping) &&
               get<DampingType>(active) == get<DampingType>(damping);
    });
    if (iter_active_same_type != active_by_type.end()) {
        //replace active of the same type and level
        auto& active_same_type = *iter_active_same_type;
        auto& sum_same_level   = *std::find_if(sum_by_level.begin(), sum_by_level.end(), [&damping](auto& sum) {
            return get<DampingLevel>(sum) == get<DampingLevel>(damping);
        });
        get<MatrixIdx>(sum_same_level) += get<MatrixIdx>(damping) - get<MatrixIdx>(active_same_type).get();
        get<MatrixIdx>(active_same_type) = get<MatrixIdx>(damping);
    }
    else {
        //add new type
        active_by_type.emplace_back(get<MatrixIdx>(damping), get<DampingLevel>(damping), get<DampingType>(damping));

        auto iter_sum_same_level = std::find_if(sum_by_level.begin(), sum_by_level.end(), [&damping](auto& sum) {
            return get<DampingLevel>(sum) == get<DampingLevel>(damping);
        });
        if (iter_sum_same_level != sum_by_level.end()) {
            //add to existing level
            get<MatrixIdx>(*iter_sum_same_level) += get<MatrixIdx>(damping);
        }
        else {
            //add new level
            sum_by_level.emplace_back(get<MatrixIdx>(damping), get<DampingLevel>(damping));
        }
    }
}

} // namespace epi
