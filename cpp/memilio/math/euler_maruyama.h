/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Martin J. Kuehn, Daniel Abele
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
#ifndef MIO_MATH_EULER_MARUYAMA_H
#define MIO_MATH_EULER_MARUYAMA_H

#include "memilio/config.h"
#include "memilio/math/integrator.h"
#include "memilio/utils/random_number_generator.h"

namespace mio
{

/**
 * @brief Simple explicit euler integration y(t+1) = y(t) + h*f(t,y) for ODE y'(t) = f(t,y)
 * @tparam FP A floating point type, e.g., ScalarType.
 */
template <typename FP>
class EulerMaruyamaIntegratorCore : public IntegratorCore<FP, 2>
{
public:
    EulerMaruyamaIntegratorCore()
        : IntegratorCore<FP, 2>(FP{}, FP{})
    {
    }

    /**
     * @brief Fixed step width of the integration
     *
     * @param[in] yt value of y at t, y(t)
     * @param[in,out] t current time step h=dt
     * @param[in,out] dt current time step h=dt
     * @param[out] ytp1 approximated value y(t+1)
     */
    bool step(const DerivFunction<FP> (&fs)[2], Eigen::Ref<const Eigen::VectorX<FP>> yt, FP& t, FP& dt,
              Eigen::Ref<Eigen::VectorX<FP>> ytp1) const override
    {
        using std::sqrt;

        auto noise = Eigen::VectorX<FP>(yt.size());
        // we are misusing the next step y as temporary space to store the derivative
        fs[0](yt, t, ytp1);
        fs[1](yt, t, noise);

        const auto white_noise = Eigen::VectorX<FP>::NullaryExpr(noise.size(), 1, [this]() {
            return sample_standart_normal_distribution();
        });

        noise = white_noise.array() * noise.array();
        ytp1  = yt + dt * ytp1 + sqrt(dt) * noise;
        rescale(ytp1);
        t += dt;
        return true;
    }

    void rescale(Eigen::Ref<Eigen::VectorX<FP>> x) const
    {
        FP positive{0.0}, negative{0.0};
        for (auto& v : x) {
            if (v >= FP{0.0}) {
                positive += v;
            }
            else {
                negative -= v;
                v = FP{0.0};
            }
        }
        if (negative > Limits<FP>::zero_tolerance()) {
            x = x * (1 + negative / positive);
        }
    }

    FP sample_standart_normal_distribution() const
    {
        return FP{mio::DistributionAdapter<std::normal_distribution<double>>::get_instance()(m_rng, 0.0, 1.0)};
    }

    RandomNumberGenerator& get_rng()
    {
        return m_rng;
    }

    const RandomNumberGenerator& get_rng() const
    {
        return m_rng;
    }

private:
    mutable RandomNumberGenerator m_rng;
};

} // namespace mio

#endif // MIO_MATH_EULER_MARUYAMA_H
