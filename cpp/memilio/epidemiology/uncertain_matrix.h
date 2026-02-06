/* 
* Copyright (C) 2020-2026 MEmilio
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
#ifndef MIO_EPI_ODE_UNCERTAINMATRIX_H
#define MIO_EPI_ODE_UNCERTAINMATRIX_H

#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/damping_sampling.h"

#include <vector>

namespace mio
{

/**
 * @brief The UncertainContactMatrix class consists of a
 *        ContactMatrix with fixed baseline and uncertain Dampings. 
 * 
 * The UncertainContactMatrix class represents a matrix-style model parameter 
 * that can take a ContactMatrix value but that is subjected to a uncertainty,
 * based on contact pattern changes realized by zero or more dampings with uncertain coefficients
 * that are sampled to modify the contacts at some points in time.
 * @see UncertainValue
 */
template <typename FP>
class UncertainContactMatrix
{
public:
    UncertainContactMatrix(size_t num_matrices = 1, Eigen::Index num_groups = 1)
        : UncertainContactMatrix<FP>(ContactMatrixGroup<FP>(num_matrices, num_groups))
    {
    }

    UncertainContactMatrix(const ContactMatrixGroup<FP>& cont_freq)
        : m_cont_freq(cont_freq)
        , m_dampings()
        , m_school_holiday_damping(mio::UncertainValue<FP>(0.0), mio::DampingLevel(0), mio::DampingType(0),
                                   mio::SimulationTime<FP>(0.0), {},
                                   Eigen::VectorX<FP>::Zero(cont_freq.get_num_groups()))
        , m_school_holidays()
    {
    }

    /**
     * @brief Conversion to const ContactMatrix reference by returning the 
     *        ContactMatrix contained in UncertainContactMatrix
     */
    operator ContactMatrixGroup<FP> const&() const
    {
        return m_cont_freq;
    }

    /**
     * @brief Conversion to ContactMatrix reference by returning the 
     *        ContactMatrix contained in UncertainContactMatrix
     */
    operator ContactMatrixGroup<FP>&()
    {
        return m_cont_freq;
    }

    /**
     * @brief Set an UncertainContactMatrix from a ContactMatrix, 
     *        all distributions remain unchanged.
     */
    UncertainContactMatrix<FP>& operator=(const ContactMatrixGroup<FP>& cont_freq)
    {
        m_cont_freq = cont_freq;
        return *this;
    }

    /**
     * @brief Returns the ContactMatrix reference 
     *        of the UncertainContactMatrix object
     */
    ContactMatrixGroup<FP>& get_cont_freq_mat()
    {
        return m_cont_freq;
    }

    /**
     * @brief Returns the const ContactMatrix reference 
     *        of the UncertainContactMatrix object
     */
    ContactMatrixGroup<FP> const& get_cont_freq_mat() const
    {
        return m_cont_freq;
    }

    /**
     * @brief Get a list of uncertain Dampings that are sampled and added to the contact matrix.
     * @return list of damping samplings.
     * @{
     */
    const std::vector<DampingSampling<FP>>& get_dampings() const
    {
        return m_dampings;
    }
    std::vector<DampingSampling<FP>>& get_dampings()
    {
        return m_dampings;
    }
    /**@}*/

    /**
     * Damping that is active during school holiday periods.
     * time is ignored and taken from holidays instead.
     * @{
     */
    const DampingSampling<FP>& get_school_holiday_damping() const
    {
        return m_school_holiday_damping;
    }
    DampingSampling<FP>& get_school_holiday_damping()
    {
        return m_school_holiday_damping;
    }
    /**@}*/

    /**
     * list of school holiday periods.
     * one period is a pair of start and end dates.
     * @{
     */
    std::vector<std::pair<SimulationTime<FP>, SimulationTime<FP>>>& get_school_holidays()
    {
        return m_school_holidays;
    }
    const std::vector<std::pair<SimulationTime<FP>, SimulationTime<FP>>>& get_school_holidays() const
    {
        return m_school_holidays;
    }
    /**@}*/

    /**
     * @brief Samples dampings and adds them to the contact matrix.
     * @param accum accumulating current and newly sampled dampings if true;
     *              default: false; removing all previously set dampings
     */
    ContactMatrixGroup<FP> draw_sample(bool accum = false)
    {
        draw_sample_dampings();
        return make_matrix(accum);
    }

    /**
     * draw sample of all dampings.
     */
    void draw_sample_dampings()
    {
        for (auto& d : m_dampings) {
            d.draw_sample();
        }
        m_school_holiday_damping.draw_sample();
    }

    /**
     * create the contact matrix using the sampled dampings.
     * @param accum accumulating current and newly dampings if true;
     *              default: false; removing all previously set dampings
     */
    ContactMatrixGroup<FP> make_matrix(bool accum = false)
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

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("UncertainContactMatrix");
        obj.add_element("ContactMatrix", m_cont_freq);
        if (!(io.flags() & IOF_OmitDistributions)) {
            obj.add_element("SchoolHolidayDamping", m_school_holiday_damping);
            obj.add_list("SchoolHolidays", m_school_holidays.begin(), m_school_holidays.end());
            obj.add_list("Dampings", m_dampings.begin(), m_dampings.end());
        }
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<UncertainContactMatrix<FP>> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("UncertainContactMatrix");
        if (!(io.flags() & IOF_OmitDistributions)) {
            auto c = obj.expect_element("ContactMatrix", Tag<ContactMatrixGroup<FP>>{});
            auto e = obj.expect_element("SchoolHolidayDamping", Tag<DampingSampling<FP>>{});
            auto f = obj.expect_list("SchoolHolidays", Tag<std::pair<SimulationTime<FP>, SimulationTime<FP>>>{});
            auto d = obj.expect_list("Dampings", Tag<DampingSampling<FP>>{});
            return apply(
                io,
                [](auto&& c_, auto&& d_, auto&& e_, auto&& f_) {
                    auto m                         = UncertainContactMatrix<FP>{c_};
                    m.get_dampings()               = d_;
                    m.get_school_holiday_damping() = e_;
                    m.get_school_holidays()        = f_;
                    return m;
                },
                c, d, e, f);
        }
        else {
            auto c = obj.expect_element("ContactMatrix", Tag<ContactMatrixGroup<FP>>{});
            return apply(
                io,
                [](auto&& c_) {
                    return UncertainContactMatrix<FP>{c_};
                },
                c);
        }
    }

private:
    ContactMatrixGroup<FP> m_cont_freq;
    std::vector<DampingSampling<FP>> m_dampings;
    DampingSampling<FP> m_school_holiday_damping;
    std::vector<std::pair<SimulationTime<FP>, SimulationTime<FP>>> m_school_holidays;
};

} // namespace mio

#endif // MIO_EPI_ODE_UNCERTAINMATRIX_H
