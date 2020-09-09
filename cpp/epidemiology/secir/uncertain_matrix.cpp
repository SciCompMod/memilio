#include "epidemiology/secir/uncertain_matrix.h"
#include "epidemiology/utils/memory.h"
#include "epidemiology/utils/parameter_distributions.h"
#include "epidemiology/secir/damping.h"

#include <memory>

namespace epi
{

ContactFrequencyMatrix::ContactFrequencyMatrix()
    : m_cont_freq{{1.0}}
    , m_dampings{{Dampings{}}}
{
}

ContactFrequencyMatrix::ContactFrequencyMatrix(size_t const nb_groups)
    : m_cont_freq{nb_groups, std::vector<double>(nb_groups, 0)}
    , m_dampings{nb_groups, std::vector<Dampings>(nb_groups, Dampings{})}
{
}

int ContactFrequencyMatrix::get_size() const
{
    return static_cast<int>(m_cont_freq.size());
}

void ContactFrequencyMatrix::set_cont_freq(double const cont_freq, int const self_group, int const contact_group)
{
    if (self_group <= contact_group) {
        m_cont_freq[self_group][contact_group] = cont_freq;
    }
    else {
        m_cont_freq[contact_group][self_group] = cont_freq;
    }
}

double ContactFrequencyMatrix::get_cont_freq(int self_group, int contact_group) const
{
    // prevent erroneous nonsymmetry
    return self_group <= contact_group ? m_cont_freq[self_group][contact_group]
                                       : m_cont_freq[contact_group][self_group];
}

void ContactFrequencyMatrix::set_dampings(Dampings const& damping, int self_group, int contact_group)
{
    if (self_group <= contact_group) {
        m_dampings[self_group][contact_group] = damping;
    }
    else {
        m_dampings[contact_group][self_group] = damping;
    }
}

const Dampings& ContactFrequencyMatrix::get_dampings(int self_group, int contact_group) const
{
    // prevent erroneous nonsymmetry
    return self_group <= contact_group ? m_dampings[self_group][contact_group] : m_dampings[contact_group][self_group];
}

void ContactFrequencyMatrix::add_damping(Damping const& damping, int self_group, int contact_group)
{
    if (self_group <= contact_group) {
        m_dampings[self_group][contact_group].add(damping);
    }
    else {
        m_dampings[contact_group][self_group].add(damping);
    }
}

void ContactFrequencyMatrix::clear_dampings()
{
    size_t nb_groups = m_dampings.size();
    m_dampings       = {nb_groups, std::vector<Dampings>(nb_groups, Dampings{})};
}

UncertainContactMatrix::UncertainContactMatrix()
    : m_cont_freq(ContactFrequencyMatrix{1})
{
}

UncertainContactMatrix::UncertainContactMatrix(ContactFrequencyMatrix cont_freq)
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

UncertainContactMatrix::operator ContactFrequencyMatrix const &() const
{
    return m_cont_freq;
}

UncertainContactMatrix::operator ContactFrequencyMatrix&()
{
    return m_cont_freq;
}

UncertainContactMatrix& UncertainContactMatrix::operator=(ContactFrequencyMatrix cont_freq)
{
    m_cont_freq = cont_freq;
    return *this;
}

ContactFrequencyMatrix& UncertainContactMatrix::get_cont_freq_mat()
{
    return m_cont_freq;
}

ContactFrequencyMatrix const& UncertainContactMatrix::get_cont_freq_mat() const
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

ContactFrequencyMatrix UncertainContactMatrix::draw_sample(bool accum)
{
    if (!accum) {
        m_cont_freq.clear_dampings();
    }

    if (m_damp_nb && m_damp_days && m_damp_diag_base && m_damp_diag_rel && m_damp_offdiag_rel) {

        int nb_dampings = (int)(m_damp_nb->get_sample() + 0.5);
        for (int i = 0; i < nb_dampings; i++) {

            double day            = m_damp_days->get_sample();
            double damp_diag_base = m_damp_diag_base->get_sample();

            // diagonal entries
            std::vector<double> damp_diag_val(m_cont_freq.get_size(), 0);
            for (int j = 0; j < m_cont_freq.get_size(); j++) {
                damp_diag_val[j] = damp_diag_base * m_damp_diag_rel->get_sample();
                m_cont_freq.add_damping(Damping(day, damp_diag_val[j]), j, j);
            }

            // offdiagonal entries
            for (int j = 0; j < m_cont_freq.get_size(); j++) {

                for (int k = j + 1; k < m_cont_freq.get_size(); k++) {
                    double damp_offdiag_val = 0.5 * damp_diag_val[j] * m_damp_offdiag_rel->get_sample() +
                                              0.5 * damp_diag_val[k] * m_damp_offdiag_rel->get_sample();
                    m_cont_freq.add_damping(Damping(day, damp_offdiag_val), j, k);
                }
            }
        }
    }
    else {
        epi::log_warning("UncertainContactMatrix distributions not set, no sampling conducted.");
    }

    return m_cont_freq;
}

} // namespace epi