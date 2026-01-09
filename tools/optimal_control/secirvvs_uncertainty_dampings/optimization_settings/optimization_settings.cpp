#include "optimization_settings.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

SecirvvsOptimization::SecirvvsOptimization(const OptimizationModel& optimization_model, size_t num_control_intervals,
                                           size_t pc_resolution, bool random_start, IntegratorType integrator_type,
                                           size_t integrator_resolution, ADType ad_eval_f, ADType ad_eval_jac,
                                           std::vector<DampingControlParameter> control_parameters,
                                           std::vector<Constraint> path_constraints,
                                           std::vector<Constraint> terminal_constraints, ControlActivation activation,
                                           size_t num_simulation_runs, unsigned int base_seed)
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
    , m_activation_function(ControlActivationFunction(activation))
    , m_num_simulation_runs(num_simulation_runs)
    , m_base_seed(base_seed)
{
    // Number of intervals to check for path constraints (each day).
    m_num_intervals = m_num_control_intervals * m_pc_resolution;
    // Time step size used in the integrator.
    m_dt = (m_tmax - m_t0) / (m_num_intervals * m_integrator_resolution);
    // Number of parameter and constraints.
    m_num_control_parameters   = control_parameters.size();
    m_num_path_constraints     = path_constraints.size();
    m_num_terminal_constraints = terminal_constraints.size();
}

const OptimizationModel& SecirvvsOptimization::optimization_model() const
{
    return this->m_optimization_model;
}

double SecirvvsOptimization::t0() const
{
    return m_t0;
}

double SecirvvsOptimization::tmax() const
{
    return m_tmax;
}

size_t SecirvvsOptimization::num_control_intervals() const
{
    return this->m_num_control_intervals;
}

size_t SecirvvsOptimization::pc_resolution() const
{
    return this->m_pc_resolution;
}

bool SecirvvsOptimization::random_start() const
{
    return this->m_random_start;
}

IntegratorType SecirvvsOptimization::integrator_type() const
{
    return this->m_integrator_type;
}

size_t SecirvvsOptimization::integrator_resolution() const
{
    return this->m_integrator_resolution;
}

ADType SecirvvsOptimization::ad_eval_f() const
{
    return this->m_ad_eval_f;
}

ADType SecirvvsOptimization::ad_eval_jac() const
{
    return this->m_ad_eval_jac;
}

std::vector<DampingControlParameter>& SecirvvsOptimization::control_parameters()
{
    return this->m_control_parameters;
}

const std::vector<DampingControlParameter>& SecirvvsOptimization::control_parameters() const
{
    return this->m_control_parameters;
}

std::vector<Constraint>& SecirvvsOptimization::path_constraints()
{
    return this->m_path_constraints;
}

const std::vector<Constraint>& SecirvvsOptimization::path_constraints() const
{
    return this->m_path_constraints;
}

std::vector<Constraint>& SecirvvsOptimization::terminal_constraints()
{
    return this->m_terminal_constraints;
}

const std::vector<Constraint>& SecirvvsOptimization::terminal_constraints() const
{
    return this->m_terminal_constraints;
}

ControlActivationFunction SecirvvsOptimization::activation_function() const
{
    return this->m_activation_function;
}

size_t SecirvvsOptimization::num_simulation_runs() const
{
    return this->m_num_simulation_runs;
}

unsigned int SecirvvsOptimization::base_seed() const
{
    return this->m_base_seed;
}

size_t SecirvvsOptimization::num_intervals() const
{
    return this->m_num_intervals;
}

double SecirvvsOptimization::dt() const
{
    return this->m_dt;
}

size_t SecirvvsOptimization::num_control_parameters() const
{
    return this->m_num_control_parameters;
}

size_t SecirvvsOptimization::num_path_constraints() const
{
    return this->m_num_path_constraints;
}

size_t SecirvvsOptimization::num_terminal_constraints() const
{
    return this->m_num_terminal_constraints;
}
