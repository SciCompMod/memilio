#ifndef FLOW_CALC_H
#define FLOW_CALC_H

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
    template <class IntegratorCore>
    flow_calculator(std::shared_ptr<TimeSeries<ScalarType>>& flows, IntegratorCore&& integrator)
        : m_dt()
        , m_deriv_integrator(std::make_shared<IntegratorCore>(integrator))
        , m_flow_integrator(std::make_shared<IntegratorCore>(integrator))
        , m_flows(flows)
    {
        static_cast<IntegratorCore&>(*m_flow_integrator).set_abs_tolerance(1);
        static_cast<IntegratorCore&>(*m_flow_integrator).set_rel_tolerance(1);
    }

    inline bool step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, ScalarType& t, ScalarType& dt,
                     Eigen::Ref<Eigen::VectorXd> ytp1) const override final
    {
        // m_dt                   = dt;
        // Eigen::VectorXd yt_old = yt;
        // auto t_old             = t;
        // Eigen::Index flow_itr  = 0;
        // Eigen::VectorXd flow_integrated(m_flows->get_num_elements());
        // flow_integrated.setZero();
        // const Eigen::VectorXd flow_zero = flow_integrated;

        // const auto rtval                = m_deriv_integrator->step(f, yt, t, dt, ytp1);

        // g() {
        //     temp = get_flows()
        //     return {get_deriv_from_flows(temp), temp}
        // }

        // m_deriv_integrator->step(g, {yt, flows}, t, dt, {ytp1, flows+1});

        // clean up TimeSeries
        /* for (Eigen::Index i = m_flows->get_num_time_points() - 1; i > 1; i--) {
            if (m_flows->get_time(i) < m_flows->get_time(i - 1)) {
                flow_itr + i - 1;
                break;
            }
        } */
        // // calculate flow
        // ScalarType last_accepted_dt = t - t_old;
        // m_flow_integrator->step(
        //     [&](Eigen::Ref<const Eigen::VectorXd> /*y*/, double /*t*/, Eigen::Ref<Eigen::VectorXd> dydt) {
        //         dydt = m_flows->get_value(flow_itr);
        //         flow_itr++;
        //     },
        //     flow_zero, t_old, last_accepted_dt, flow_integrated);
        // // store ytp1_dummy

        // // delete m_flows
        // return rtval;
        return false;
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