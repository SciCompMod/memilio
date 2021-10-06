/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#include "memilio/epidemiology/uncertain_matrix.h"

namespace mio
{
UncertainContactMatrix::UncertainContactMatrix(size_t num_matrices, Eigen::Index num_groups)
    : UncertainContactMatrix(ContactMatrixGroup(num_matrices, num_groups))
{
}

UncertainContactMatrix::UncertainContactMatrix(const ContactMatrixGroup& cont_freq)
    : m_cont_freq(cont_freq)
    , m_dampings()
    , m_school_holiday_damping(0.0, mio::DampingLevel(0), mio::DampingType(0), mio::SimulationTime(0), {}, Eigen::VectorXd::Zero(cont_freq.get_num_groups()))
    , m_school_holidays()
{
}

UncertainContactMatrix::operator ContactMatrixGroup const &() const
{
    return m_cont_freq;
}

UncertainContactMatrix::operator ContactMatrixGroup&()
{
    return m_cont_freq;
}

UncertainContactMatrix& UncertainContactMatrix::operator=(const ContactMatrixGroup& cont_freq)
{
    m_cont_freq = cont_freq;
    return *this;
}

ContactMatrixGroup& UncertainContactMatrix::get_cont_freq_mat()
{
    return m_cont_freq;
}

ContactMatrixGroup const& UncertainContactMatrix::get_cont_freq_mat() const
{
    return m_cont_freq;
}

ContactMatrixGroup UncertainContactMatrix::draw_sample(bool accum)
{
    draw_sample_dampings();
    return make_matrix(accum);
}

void UncertainContactMatrix::draw_sample_dampings()
{
    for (auto& d : m_dampings) {
        d.draw_sample();
    }
    m_school_holiday_damping.draw_sample();
}

ContactMatrixGroup UncertainContactMatrix::make_matrix(bool accum)
{
    if (!accum) {
        m_cont_freq.clear_dampings();
    }

    auto make_matrix = [](auto&& v) {
        return mio::make_contact_damping_matrix(v);
    };
    mio::apply_dampings(m_cont_freq, m_dampings, make_matrix);

    for (auto h : m_school_holidays) {
        //enable damping at the start of the period
        auto damping = m_school_holiday_damping;
        damping.set_time(h.first);
        mio::apply_dampings(m_cont_freq, make_range(&damping, &damping + 1), make_matrix);

        //disable damping at the end of the period
        damping.get_value() = 0.0;
        damping.set_time(h.second);
        mio::apply_dampings(m_cont_freq, make_range(&damping, &damping + 1), make_matrix);
    }

    return m_cont_freq;
}

} // namespace mio