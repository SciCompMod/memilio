#include "optimization_settings.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

SecirOptimization::SecirOptimization(const OptimizationModel& optimization_model, size_t num_control_intervals,
                                     size_t pc_resolution, bool random_start, IntegratorType integrator_type,
                                     size_t integrator_resolution, ADType ad_eval_f, ADType ad_eval_jac,
                                     std::vector<DynamicNPIControlParameter> control_parameters,
                                     std::vector<Constraint> path_constraints,
                                     std::vector<Constraint> terminal_constraints)
    : m_optimization_model(optimization_model)
    , m_t0(0.0)
    , m_tmax(m_optimization_model.simulation_days())
    , m_num_control_intervals(num_control_intervals)
    , m_pc_resolution(pc_resolution)
    , m_random_start(random_start)
    , m_integrator_type(integrator_type)
    , m_integrator_resolution(integrator_resolution)
    , m_ad_eval_f(ad_eval_f)
    , m_ad_eval_jac(ad_eval_jac)
    , m_control_parameters(control_parameters)
    , m_path_constraints(path_constraints)
    , m_terminal_constraints(terminal_constraints)
{
    // Number of intervals to check for path constraints (each day).
    m_num_intervals = m_num_control_intervals * m_pc_resolution;
    // Time step size used in the integrator.
    m_dt = (m_tmax - m_t0) / (m_num_intervals * m_integrator_resolution);
    // Number of parameter and constraints.
    m_num_control_parameters = 0;
    for (const auto& cp : control_parameters) {
        m_num_control_parameters += cp.dampings().size();
    }
    m_num_path_constraints     = path_constraints.size();
    m_num_terminal_constraints = terminal_constraints.size();
}

const OptimizationModel& SecirOptimization::optimization_model() const
{
    return this->m_optimization_model;
}

double SecirOptimization::t0() const
{
    return m_t0;
}

double SecirOptimization::tmax() const
{
    return m_tmax;
}

size_t SecirOptimization::num_control_intervals() const
{
    return this->m_num_control_intervals;
}

size_t SecirOptimization::pc_resolution() const
{
    return this->m_pc_resolution;
}

bool SecirOptimization::random_start() const
{
    return this->m_random_start;
}

IntegratorType SecirOptimization::integrator_type() const
{
    return this->m_integrator_type;
}

size_t SecirOptimization::integrator_resolution() const
{
    return this->m_integrator_resolution;
}

ADType SecirOptimization::ad_eval_f() const
{
    return this->m_ad_eval_f;
}

ADType SecirOptimization::ad_eval_jac() const
{
    return this->m_ad_eval_jac;
}

std::vector<DynamicNPIControlParameter>& SecirOptimization::control_parameters()
{
    return this->m_control_parameters;
}

const std::vector<DynamicNPIControlParameter>& SecirOptimization::control_parameters() const
{
    return this->m_control_parameters;
}

std::vector<Constraint>& SecirOptimization::path_constraints()
{
    return this->m_path_constraints;
}

const std::vector<Constraint>& SecirOptimization::path_constraints() const
{
    return this->m_path_constraints;
}

std::vector<Constraint>& SecirOptimization::terminal_constraints()
{
    return this->m_terminal_constraints;
}

const std::vector<Constraint>& SecirOptimization::terminal_constraints() const
{
    return this->m_terminal_constraints;
}

size_t SecirOptimization::num_intervals() const
{
    return this->m_num_intervals;
}

double SecirOptimization::dt() const
{
    return this->m_dt;
}

size_t SecirOptimization::num_control_parameters() const
{
    return this->m_num_control_parameters;
}

size_t SecirOptimization::num_path_constraints() const
{
    return this->m_num_path_constraints;
}

size_t SecirOptimization::num_terminal_constraints() const
{
    return this->m_num_terminal_constraints;
}
