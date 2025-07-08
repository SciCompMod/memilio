// clang-format off
#include "secirvvs_optimization.h"

// ----------- //
// Constructor //
// ----------- //
SecirvvsOptimization::SecirvvsOptimization(
    const mio::Graph<mio::SimulationNode<double, mio::Simulation<double, mio::osecirvvs::Model<double>>>,
                     mio::MobilityEdge<double>>& graph,
    double t0, double tmax)
    : m_graph(graph)
    , m_t0(t0)
    , m_tmax(tmax)
{}

void SecirvvsOptimization::configure()
{
    m_num_graph_nodes = graph().nodes().size();
    m_num_graph_edges = graph().edges().size();

    m_num_intervals = numControlIntervals() * pcResolution();
    m_dt            = (tmax() - t0()) / (numIntervals() * integratorResolution());

    // Vector of control policies for various parameters: { name, {lower, upper}, effectiveness, cost_penalty }
    m_control_policy = {
        {"SchoolClosure",            {0.0, 1.0}, 1.00, 1.0},
        {"HomeOffice",               {0.0, 1.0}, 0.25, 1.0},
        {"PhysicalDistancingSchool", {0.0, 1.0}, 0.25, 1.0},
        {"PhysicalDistancingWork",   {0.0, 1.0}, 0.25, 1.0},
        {"PhysicalDistancingOther",  {0.0, 1.0}, 0.35, 1.0}
    };

    m_num_controls = m_control_policy.size();
    m_num_path_constraints =  m_path_constraints.size();
    m_num_terminal_constraints = m_terminal_constraints.size();

    validateConstraints();
}


const mio::Graph<mio::SimulationNode<double, mio::Simulation<double, mio::osecirvvs::Model<double>>>,
                 mio::MobilityEdge<double>>&
SecirvvsOptimization::graph() const
{
    return m_graph;
}
double SecirvvsOptimization::t0() const
{
    return m_t0;
}
double SecirvvsOptimization::tmax() const
{
    return m_tmax;
}

// ----------- //
// Constraints //
// ----------- //
void SecirvvsOptimization::addPathConstraint(const std::string& name, const std::pair<double, double>& bounds)
{
    m_path_constraints.emplace_back(name, bounds, std::nullopt);
}
void SecirvvsOptimization::addPathConstraint(const std::string& name, const std::pair<double, double>& bounds,
                                             size_t node_index)
{
    m_path_constraints.emplace_back(name, bounds, node_index);
}

void SecirvvsOptimization::addTerminalConstraint(const std::string& name, const std::pair<double, double>& bounds)
{
    m_terminal_constraints.emplace_back(name, bounds, std::nullopt);
}
void SecirvvsOptimization::addTerminalConstraint(const std::string& name, const std::pair<double, double>& bounds,
                                                 size_t node_index)
{
    m_terminal_constraints.emplace_back(name, bounds, node_index);
}

size_t SecirvvsOptimization::numControlIntervals() const
{
    return m_num_control_intervals;
}
void SecirvvsOptimization::setNumControlIntervals(size_t num_control_intervals)
{
    m_num_control_intervals = num_control_intervals;
}

size_t SecirvvsOptimization::pcResolution() const
{
    return m_pc_resolution;
}
void SecirvvsOptimization::setPcResolution(size_t pc_resolution)
{
    m_pc_resolution = pc_resolution;
}

bool SecirvvsOptimization::randomStart() const
{
    return m_random_start;
}
void SecirvvsOptimization::setRandomStart(bool random_start)
{
    m_random_start = random_start;
}

IntegratorType SecirvvsOptimization::integratorType() const
{
    return m_integrator_type;
}

void SecirvvsOptimization::setIntegratorType(IntegratorType integrator_type)
{
    m_integrator_type = integrator_type;
}

size_t SecirvvsOptimization::integratorResolution() const
{
    return m_integrator_resolution;
}
void SecirvvsOptimization::setIntegratorResolution(size_t integrator_resolution)
{
    m_integrator_resolution = integrator_resolution;
}

ADType SecirvvsOptimization::adEvalF() const
{
    return m_ad_eval_f;
}
void SecirvvsOptimization::setAdEvalF(ADType ad_eval_f)
{
    m_ad_eval_f = ad_eval_f;
}

ADType SecirvvsOptimization::adEvalJac() const
{
    return m_ad_eval_jac;
}
void SecirvvsOptimization::setAdEvalJac(ADType ad_eval_jac)
{
    m_ad_eval_jac = ad_eval_jac;
}

// --- IPOPT Options ---
int SecirvvsOptimization::ipopt_verbose() const
{
    return m_ipopt_verbose;
}
void SecirvvsOptimization::setIpoptVerbose(int ipopt_verbose)
{
    m_ipopt_verbose = ipopt_verbose;
}

int SecirvvsOptimization::ipopt_printFrequency() const
{
    return m_ipopt_print_frequency;
}
void SecirvvsOptimization::setIpoptPrintFrequency(int ipopt_print_frequency)
{
    m_ipopt_print_frequency = ipopt_print_frequency;
}

bool SecirvvsOptimization::ipopt_printTimings() const
{
    return m_ipopt_print_timings;
}
void SecirvvsOptimization::setIpoptPrintTimings(bool ipopt_print_timings)
{
    m_ipopt_print_timings = ipopt_print_timings;
}

int SecirvvsOptimization::ipopt_maxIter() const
{
    return m_ipopt_max_iter;
}
void SecirvvsOptimization::setIpoptMaxIter(int ipopt_max_iter)
{
    m_ipopt_max_iter = ipopt_max_iter;
}

double SecirvvsOptimization::ipopt_tolerance() const
{
    return m_ipopt_tolerance;
}
void SecirvvsOptimization::setIpoptTolerance(double ipopt_tolerance)
{
    m_ipopt_tolerance = ipopt_tolerance;
}

size_t SecirvvsOptimization::numGraphNodes() const{
    return m_num_graph_nodes;
}
size_t SecirvvsOptimization::numGraphEdges() const{
    return m_num_graph_edges;
}

size_t SecirvvsOptimization::numIntervals() const{
    return m_num_intervals;
}
double SecirvvsOptimization::dt() const{
    return m_dt;
}

size_t SecirvvsOptimization::numControls() const{
    return m_num_controls;
}
size_t SecirvvsOptimization::numPathConstraints() const{
    return m_num_path_constraints;
}
size_t SecirvvsOptimization::numTerminalConstraints() const{
    return m_num_terminal_constraints;
}

std::vector<ConstraintPolicy<double>>& SecirvvsOptimization::pathConstraints() {
    return m_path_constraints;
}
std::vector<ConstraintPolicy<double>>& SecirvvsOptimization::terminalConstraints() {
    return m_terminal_constraints;
}