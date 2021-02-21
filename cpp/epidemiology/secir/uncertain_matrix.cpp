#include "epidemiology/secir/uncertain_matrix.h"
#include "epidemiology/utils/memory.h"
#include "epidemiology/utils/parameter_distributions.h"

#include <memory>

namespace epi
{
UncertainContactMatrix::UncertainContactMatrix(size_t num_matrices, Eigen::Index num_groups)
    : m_cont_freq(num_matrices, num_groups)
    , m_dampings()
{
}

UncertainContactMatrix::UncertainContactMatrix(const ContactMatrixGroup& cont_freq)
    : m_cont_freq(cont_freq)
    , m_dampings()
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
    if (!accum)
    {
        m_cont_freq.clear_dampings();
    }

    for (auto& d : m_dampings)
    {
        d.draw_sample();
    }

    epi::apply_dampings(m_cont_freq, m_dampings, [](auto&& v) {
        return epi::make_contact_damping_sampling_mask(v);
    });

    return m_cont_freq;
}

} // namespace epi