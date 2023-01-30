#ifndef FLOW_CALC_H
#define FLOW_CALC_H

#include "Eigen/src/Core/util/Meta.h"
#include "memilio/config.h"
#include "memilio/math/eigen.h"
#include "memilio/math/integrator.h"
#include "memilio/utils/time_series.h"
#include <memory>

namespace mio
{

class flow_calculator : public mio::IntegratorCore
{
public:
    template <class Integrator, class... Args>
    explicit flow_calculator(std::shared_ptr<TimeSeries<ScalarType>> flows, Args... args)
        : m_dt()
        , m_deriv_integrator(std::make_shared<Integrator>(args...))
        , m_flow_integrator(std::make_shared<Integrator>(args...))
        , m_flows(flows)
    {
    }

    inline bool step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, ScalarType& t, ScalarType& dt,
                     Eigen::Ref<Eigen::VectorXd> ytp1) const override final
    {
        m_dt                       = dt;
        Eigen::VectorXd yt_old     = yt;
        auto t_old                 = t;
        Eigen::VectorXd ytp1_dummy = ytp1;
        const auto rtval           = m_deriv_integrator->step(f, yt, t, dt, ytp1);
        // clean up TimeSeries
        auto flow_itr = m_flows->begin();
        for (Eigen::Index i = m_flows->get_num_time_points() - 1; i > 1; i--) {
            if (m_flows->get_time(i) < m_flows->get_time(i - 1)) {
                flow_itr + i - 1;
                break;
            }
        }
        // calculate flow
        ScalarType last_accepted_dt = t - t_old;
        m_flow_integrator->step(
            [&](Eigen::Ref<const Eigen::VectorXd> /*y*/, double /*t*/, Eigen::Ref<Eigen::VectorXd> dydt) {
                dydt = *flow_itr;
                flow_itr++;
            },
            yt_old, t_old, last_accepted_dt, ytp1_dummy);
        // store ytp1_dummy

        // delete m_flows
        return rtval;
    }

    void set_integrator(std::shared_ptr<IntegratorCore> integrator)
    {
        m_deriv_integrator = std::move(integrator);
    }

    IntegratorCore& get_integrator()
    {
        return *m_deriv_integrator;
    }

    const IntegratorCore& get_integrator() const
    {
        return *m_deriv_integrator;
    }

    ScalarType get_dt() const
    {
        return m_dt;
    }

private:
    mutable ScalarType m_dt;
    std::shared_ptr<IntegratorCore> m_deriv_integrator;
    std::shared_ptr<IntegratorCore> m_flow_integrator;
    std::shared_ptr<TimeSeries<ScalarType>> m_flows;
};

} // namespace mio

#endif