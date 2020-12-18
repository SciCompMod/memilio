#include "epidemiology/secir/uncertain_matrix.h"
#include "epidemiology/utils/memory.h"
#include "epidemiology/utils/parameter_distributions.h"

#include <memory>

namespace epi
{
UncertainContactMatrix::UncertainContactMatrix(Eigen::Index num_groups, size_t num_matrices)
    : m_cont_freq(num_groups, num_matrices)
{
}

UncertainContactMatrix::UncertainContactMatrix(const ContactMatrixGroup& cont_freq)
    : m_cont_freq(cont_freq)
{
}

UncertainContactMatrix::UncertainContactMatrix(UncertainContactMatrix const& other)
    : m_cont_freq(other.m_cont_freq)
{
    if (other.m_damp_diag_base) {
        m_damp_diag_base.reset(other.m_damp_diag_base->clone());
    }
    if (other.m_damp_diag_rel) {
        m_damp_diag_rel.reset(other.m_damp_diag_rel->clone());
    }
    if (other.m_damp_offdiag_rel) {
        m_damp_offdiag_rel.reset(other.m_damp_offdiag_rel->clone());
    }
    if (other.m_damp_nb) {
        m_damp_nb.reset(other.m_damp_nb->clone());
    }
    if (other.m_damp_days) {
        m_damp_days.reset(other.m_damp_days->clone());
    }
}

UncertainContactMatrix& UncertainContactMatrix::operator=(const UncertainContactMatrix& other)
{
    UncertainContactMatrix tmp(other);
    std::swap(m_cont_freq, tmp.m_cont_freq);
    std::swap(m_damp_nb, tmp.m_damp_nb);
    std::swap(m_damp_days, tmp.m_damp_days);
    std::swap(m_damp_diag_base, tmp.m_damp_diag_base);
    std::swap(m_damp_diag_rel, tmp.m_damp_diag_rel);
    std::swap(m_damp_offdiag_rel, tmp.m_damp_offdiag_rel);

    return *this;
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

void UncertainContactMatrix::set_distribution_damp_nb(const ParameterDistribution& dist)
{
    m_damp_nb.reset(dist.clone());
}

void UncertainContactMatrix::set_distribution_damp_days(const ParameterDistribution& dist)
{
    m_damp_days.reset(dist.clone());
}

void UncertainContactMatrix::set_distribution_damp_diag_base(const ParameterDistribution& dist)
{
    m_damp_diag_base.reset(dist.clone());
}

void UncertainContactMatrix::set_distribution_damp_diag_rel(const ParameterDistribution& dist)
{
    m_damp_diag_rel.reset(dist.clone());
}

void UncertainContactMatrix::set_distribution_damp_offdiag_rel(const ParameterDistribution& dist)
{
    m_damp_offdiag_rel.reset(dist.clone());
}

observer_ptr<ParameterDistribution> UncertainContactMatrix::get_distribution_damp_nb()
{
    return m_damp_nb.get();
}

observer_ptr<const ParameterDistribution> UncertainContactMatrix::get_distribution_damp_nb() const
{
    return m_damp_nb.get();
}

observer_ptr<ParameterDistribution> UncertainContactMatrix::get_distribution_damp_days()
{
    return m_damp_days.get();
}

observer_ptr<const ParameterDistribution> UncertainContactMatrix::get_distribution_damp_days() const
{
    return m_damp_days.get();
}

observer_ptr<ParameterDistribution> UncertainContactMatrix::get_distribution_damp_diag_base()
{
    return m_damp_diag_base.get();
}

observer_ptr<const ParameterDistribution> UncertainContactMatrix::get_distribution_damp_diag_base() const
{
    return m_damp_diag_base.get();
}

observer_ptr<ParameterDistribution> UncertainContactMatrix::get_distribution_damp_diag_rel()
{
    return m_damp_diag_rel.get();
}

observer_ptr<const ParameterDistribution> UncertainContactMatrix::get_distribution_damp_diag_rel() const
{
    return m_damp_diag_rel.get();
}

observer_ptr<ParameterDistribution> UncertainContactMatrix::get_distribution_damp_offdiag_rel()
{
    return m_damp_offdiag_rel.get();
}

observer_ptr<const ParameterDistribution> UncertainContactMatrix::get_distribution_damp_offdiag_rel() const
{
    return m_damp_offdiag_rel.get();
}

ContactMatrixGroup UncertainContactMatrix::draw_sample(bool accum)
{
    if (m_damp_nb && m_damp_days && m_damp_diag_base && m_damp_diag_rel && m_damp_offdiag_rel) {
        if (!accum) {
            //matrices with same baseline/minimum but no dampings
            auto cont_freq =
                ContactMatrixGroup{m_cont_freq.get_num_groups(), m_cont_freq.get_num_matrices()};
            for (size_t i = 0; i < m_cont_freq.get_num_matrices(); ++i) {
                cont_freq[i].get_baseline() = m_cont_freq[i].get_baseline();
                cont_freq[i].get_minimum()  = m_cont_freq[i].get_minimum();
            }
            m_cont_freq = cont_freq;
        }

        int nb_dampings = (int)(m_damp_nb->get_sample() + 0.5);
        for (int i = 0; i < nb_dampings; i++) {

            double day            = m_damp_days->get_sample();
            double damp_diag_base = m_damp_diag_base->get_sample();

            // diagonal entries
            Eigen::MatrixXd damping(m_cont_freq.get_num_groups(), m_cont_freq.get_num_groups());
            for (Eigen::Index j = 0; j < m_cont_freq.get_num_groups(); j++) {
                damping(j, j) = damp_diag_base * m_damp_diag_rel->get_sample();
                for (Eigen::Index k = 0; k < m_cont_freq.get_num_groups(); k++) {
                    if (j != k) {
                        damping(j, k) = damping(j, j) * m_damp_offdiag_rel->get_sample();
                    }
                }
            }

            m_cont_freq.add_damping(damping, SimulationTime(day));
        }
    }
    else {
        epi::log_warning("UncertainContactMatrix distributions not set, no sampling conducted.");
    }

    return m_cont_freq;
}

} // namespace epi