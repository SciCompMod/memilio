/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

#include "memilio/utils/date.h"
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
class UncertainContactMatrix
{
public:
    UncertainContactMatrix(size_t num_matrices = 1, Eigen::Index num_groups = 1);

    UncertainContactMatrix(const ContactMatrixGroup& cont_freq);

    /**
     * @brief Conversion to const ContactMatrix reference by returning the 
     *        ContactMatrix contained in UncertainContactMatrix
     */
    operator ContactMatrixGroup const &() const;

    /**
     * @brief Conversion to ContactMatrix reference by returning the 
     *        ContactMatrix contained in UncertainContactMatrix
     */
    operator ContactMatrixGroup&();

    /**
     * @brief Set an UncertainContactMatrix from a ContactMatrix, 
     *        all distributions remain unchanged.
     */
    UncertainContactMatrix& operator=(const ContactMatrixGroup& cont_freq);

    /**
     * @brief Returns the ContactMatrix reference 
     *        of the UncertainContactMatrix object
     */
    ContactMatrixGroup& get_cont_freq_mat();

    /**
     * @brief Returns the const ContactMatrix reference 
     *        of the UncertainContactMatrix object
     */
    ContactMatrixGroup const& get_cont_freq_mat() const;

    /**
     * @brief Get a list of uncertain Dampings that are sampled and added to the contact matrix.
     * @return list of damping samplings.
     * @{
     */
    const std::vector<DampingSampling>& get_dampings() const
    {
        return m_dampings;
    }
    std::vector<DampingSampling>& get_dampings()
    {
        return m_dampings;
    }
    /**@}*/

    /**
     * Damping that is active during school holiday periods.
     * time is ignored and taken from holidays instead.
     * @{
     */
    const DampingSampling& get_school_holiday_damping() const
    {
        return m_school_holiday_damping;
    }
    DampingSampling& get_school_holiday_damping()
    {
        return m_school_holiday_damping;
    }
    /**@}*/

    /**
     * list of school holiday periods.
     * one period is a pair of start and end dates.
     * @{
     */
    std::vector<std::pair<SimulationTime, SimulationTime>>& get_school_holidays()
    {
        return m_school_holidays;
    }
    const std::vector<std::pair<SimulationTime, SimulationTime>>& get_school_holidays() const
    {
        return m_school_holidays;
    }
    /**@}*/

    /**
     * @brief Samples dampings and adds them to the contact matrix.
     * @param accum accumulating current and newly sampled dampings if true;
     *              default: false; removing all previously set dampings
     */
    ContactMatrixGroup draw_sample(bool accum = false);

    /**
     * draw sample of all dampings.
     */
    void draw_sample_dampings();

    /**
     * create the contact matrix using the sampled dampings.
     * @param accum accumulating current and newly dampings if true;
     *              default: false; removing all previously set dampings
     */
    ContactMatrixGroup make_matrix(bool accum = false);

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
    static IOResult<UncertainContactMatrix> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("UncertainContactMatrix");
        if (!(io.flags() & IOF_OmitDistributions)) {
            auto c = obj.expect_element("ContactMatrix", Tag<ContactMatrixGroup>{});
            auto e = obj.expect_element("SchoolHolidayDamping", Tag<DampingSampling>{});
            auto f = obj.expect_list("SchoolHolidays", Tag<std::pair<SimulationTime, SimulationTime>>{});
            auto d = obj.expect_list("Dampings", Tag<DampingSampling>{});
            return apply(
                io,
                [](auto&& c_, auto&& d_, auto&& e_, auto&& f_) {
                    auto m           = UncertainContactMatrix{c_};
                    m.get_dampings() = d_;
                    m.get_school_holiday_damping() = e_;
                    m.get_school_holidays() = f_;
                    return m;
                },
                c, d, e, f);
        }
        else {
            auto c = obj.expect_element("ContactMatrix", Tag<ContactMatrixGroup>{});
            return apply(
                io,
                [](auto&& c_) {
                    return UncertainContactMatrix{c_};
                },
                c);
        }
    }

private:
    ContactMatrixGroup m_cont_freq;
    std::vector<DampingSampling> m_dampings;
    DampingSampling m_school_holiday_damping;
    std::vector<std::pair<SimulationTime, SimulationTime>> m_school_holidays;
};

} // namespace mio

#endif // MIO_EPI_ODE_UNCERTAINMATRIX_H
