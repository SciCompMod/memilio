#include "epidemiology/secir/uncertain_matrix.h"

namespace epi
{
UncertainContactMatrix::UncertainContactMatrix(size_t num_matrices, Eigen::Index num_groups)
    : UncertainContactMatrix(ContactMatrixGroup(num_matrices, num_groups))
{
}

UncertainContactMatrix::UncertainContactMatrix(const ContactMatrixGroup& cont_freq)
    : m_cont_freq(cont_freq)
    , m_dampings()
    , m_school_holiday_damping(0.0, epi::DampingLevel(0), epi::DampingType(0), epi::SimulationTime(0), {}, Eigen::VectorXd::Zero(cont_freq.get_num_groups()))
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
        return epi::make_contact_damping_matrix(v);
    };
    epi::apply_dampings(m_cont_freq, m_dampings, make_matrix);

    for (auto h : m_school_holidays) {
        //enable damping at the start of the period
        auto damping = m_school_holiday_damping;
        damping.set_time(h.first);
        epi::apply_dampings(m_cont_freq, make_range(&damping, &damping + 1), make_matrix);

        //disable damping at the end of the period
        damping.get_value() = 0.0;
        damping.set_time(h.second);
        epi::apply_dampings(m_cont_freq, make_range(&damping, &damping + 1), make_matrix);
    }

    return m_cont_freq;
}

} // namespace epi