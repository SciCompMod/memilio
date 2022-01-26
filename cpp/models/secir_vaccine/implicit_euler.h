/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Alexander Ruettgers, Martin J. Kuehn, Daniel Abele
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
#ifndef EPI_SECIR_IMPLICIT_EULER_H
#define EPI_SECIR_IMPLICIT_EULER_H

#include "memilio/math/euler.h"
#include "secir/secir.h"

namespace mio
{
namespace vaccinated
{

    /**
 * @brief Implicit Euler integration (not generalized, adapted to SECIHURD-model)
 */
    class ImplicitEulerIntegratorCore : public IntegratorCore
    {
    public:
        /**
     * @brief Setting up the implicit Euler integrator
     * @param params Paramters of the SECIR/SECIHURD model
     */
        ImplicitEulerIntegratorCore(SecirModel const& params);

        /**
     * @brief Fixed step width of the time implicit Euler time integration scheme
     *
     * @param[in] yt value of y at t, y(t)
     * @param[in,out] t current time step h=dt
     * @param[in,out] dt current time step h=dt
     * @param[out] ytp1 approximated value y(t+1)
     */
        bool step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
                  Eigen::Ref<Eigen::VectorXd> ytp1) const override;

        SecirModel const& get_secir_params() const
        {
            return m_model;
        }

        /**
     *  @param tol the required absolute tolerance for the comparison with the Fehlberg approximation (actually not really required but used in SecirSimulation constructor)
     */
        void set_abs_tolerance(double tol)
        {
            m_abs_tol = tol;
        }

    private:
        const SecirModel& m_model;
        double m_abs_tol = 1e-4;
    };

} // namespace vaccinated
} // namespace mio

#endif
