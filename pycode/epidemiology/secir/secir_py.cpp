#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

#include <epidemiology/secir/secir.h>
#include <epidemiology/secir/damping.h>
#include <epidemiology/utils/time_series.h>
#include <epidemiology/secir/parameter_studies.h>
#include <epidemiology/secir/analyze_result.h>

#include <Eigen/Core>
#include <vector>

namespace py = pybind11;

namespace
{

std::vector<epi::TimeSeries<double>> filter_graph_results(
    const std::vector<epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge>>& graph_results)
{
    std::vector<epi::TimeSeries<double>> results;
    results.reserve(graph_results.size());
    std::transform(graph_results.begin(), graph_results.end(), std::back_inserter(results), [](auto&& graph) {
        return graph.nodes()[0].model.get_result();
    });
    return results;
}

} // namespace

PYBIND11_MODULE(_secir, m)
{
    py::enum_<epi::SecirCategory>(m, "SecirCategory")
        .value("InfectionType", epi::SecirCategory::InfectionType)
        .value("AgeGroup", epi::SecirCategory::AgeGroup)
        .value("CategoryCount", epi::SecirCategory::CategoryCount)
        .export_values();

    py::enum_<epi::SecirCompartments>(m, "SecirCompartments")
        .value("S", epi::SecirCompartments::S)
        .value("E", epi::SecirCompartments::E)
        .value("C", epi::SecirCompartments::C)
        .value("I", epi::SecirCompartments::I)
        .value("H", epi::SecirCompartments::H)
        .value("U", epi::SecirCompartments::U)
        .value("R", epi::SecirCompartments::R)
        .value("D", epi::SecirCompartments::D)
        .value("SecirCount", epi::SecirCompartments::SecirCount)
        .export_values();

    py::class_<epi::Damping>(m, "Damping")
        .def(py::init<double, double>(), py::arg("day"), py::arg("factor"))
        .def_readwrite("day", &epi::Damping::day)
        .def_readwrite("factor", &epi::Damping::factor);

    py::class_<epi::Dampings>(m, "Dampings")
        .def(py::init<>())
        .def("add", &epi::Dampings::add)
        .def("get_factor", &epi::Dampings::get_factor);

    py::class_<epi::TimeSeries<double>>(m, "TimeSeries")
        .def(py::init<Eigen::Index>(), py::arg("num_elements"))
        .def("get_num_time_points", &epi::TimeSeries<double>::get_num_time_points)
        .def("get_num_elements", &epi::TimeSeries<double>::get_num_elements)
        .def("get_time", py::overload_cast<Eigen::Index>(&epi::TimeSeries<double>::get_time), py::arg("index"))
        .def("get_last_time", py::overload_cast<>(&epi::TimeSeries<double>::get_last_time))
        .def("get_value", py::overload_cast<Eigen::Index>(&epi::TimeSeries<double>::get_value), py::arg("index"))
        .def("get_last_value", py::overload_cast<>(&epi::TimeSeries<double>::get_last_value))
        .def(
            "__getitem__",
            [](epi::TimeSeries<double>& self, Eigen::Index i) {
                if (i >= 0 && i < self.get_num_time_points()) {
                    return self[i];
                }
                else {
                    throw pybind11::index_error("Index out of range."); //needs to throw exception for iterable
                }
            },
            py::is_operator(), py::arg("index"))
        .def(
            "__setitem__",
            [](epi::TimeSeries<double>& self, Eigen::Index i, Eigen::Ref<const epi::TimeSeries<double>::Vector> expr) {
                if (i >= 0 && i < self.get_num_time_points()) {
                    self[i] = expr;
                }
                else {
                    throw pybind11::index_error("Index out of range."); //needs to throw exception for iterable
                }
            },
            py::is_operator(), py::arg("index"), py::arg("v"))
        .def("add_time_point",
             [](epi::TimeSeries<double>& self) {
                 return self.add_time_point();
             })
        .def("add_time_point",
             [](epi::TimeSeries<double>& self, double t) {
                 return self.add_time_point(t);
             })
        .def("add_time_point",
             [](epi::TimeSeries<double>& self, double t, Eigen::Ref<const epi::TimeSeries<double>::Vector> expr) {
                 return self.add_time_point(t, expr);
             })
        .def("as_ndarray", [](epi::TimeSeries<double>& self) {
            auto m = Eigen::Map<epi::TimeSeries<double>::Matrix>(self.data(), self.get_num_rows(),
                                                                 self.get_num_time_points());
            return Eigen::Ref<epi::TimeSeries<double>::Matrix>(m);
        });

    py::class_<epi::ParameterDistribution>(m, "ParameterDistribution")
        .def_property("lower_bound", &epi::ParameterDistribution::get_lower_bound,
                      &epi::ParameterDistribution::set_lower_bound)
        .def_property("upper_bound", &epi::ParameterDistribution::get_upper_bound,
                      &epi::ParameterDistribution::set_upper_bound)
        .def("add_predefined_sample", &epi::ParameterDistribution::add_predefined_sample)
        .def("remove_predefined_samples", &epi::ParameterDistribution::remove_predefined_samples)
        .def("get_sample", &epi::ParameterDistribution::get_sample);

    py::class_<epi::ParameterDistributionNormal, epi::ParameterDistribution>(m, "ParameterDistributionNormal")
        .def(py::init<double, double, double, double>(), py::arg("lb"), py::arg("ub"), py::arg("mean"),
             py::arg("std_dev"))
        .def(py::init<double, double, double>(), py::arg("lb"), py::arg("ub"), py::arg("mean"))
        .def_property("mean", &epi::ParameterDistributionNormal::get_mean, &epi::ParameterDistributionNormal::set_mean)
        .def_property("standard_dev", &epi::ParameterDistributionNormal::get_standard_dev,
                      &epi::ParameterDistributionNormal::set_standard_dev);

    py::class_<epi::ParameterDistributionUniform, epi::ParameterDistribution>(m, "ParameterDistributionUniform")
        .def(py::init<>())
        .def(py::init<double, double>(), py::arg("lb"), py::arg("ub"));

    py::class_<epi::UncertainValue>(m, "UncertainValue")
        .def(py::init<double>(), py::arg("value") = 0.0)
        .def_property(
            "value",
            [](epi::UncertainValue& self) {
                return double(self);
            },
            [](epi::UncertainValue& self, double v) {
                self = v;
            })
        .def("set_distribution", //a property would be nicer but getter and setter use a different type
             &epi::UncertainValue::set_distribution)
        .def(
            "get_distribution",
            [](const epi::UncertainValue& self) {
                return self.get_distribution().get();
            },
            py::return_value_policy::reference_internal)
        .def(
            "get_distribution",
            [](epi::UncertainValue& self) {
                return self.get_distribution().get();
            },
            py::return_value_policy::reference_internal)
        .def("draw_sample", &epi::UncertainValue::draw_sample);

    py::class_<epi::SecirParams::StageTimes>(m, "StageTimes")
        .def(py::init<>())
        .def("set_incubation", py::overload_cast<double>(&epi::SecirParams::StageTimes::set_incubation))
        .def("set_incubation",
             py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams::StageTimes::set_incubation))
        .def("set_infectious_mild", py::overload_cast<double>(&epi::SecirParams::StageTimes::set_infectious_mild))
        .def("set_infectious_mild",
             py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams::StageTimes::set_infectious_mild))
        .def("set_serialinterval", py::overload_cast<double>(&epi::SecirParams::StageTimes::set_serialinterval))
        .def("set_serialinterval",
             py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams::StageTimes::set_serialinterval))
        .def("set_hospitalized_to_home",
             py::overload_cast<double>(&epi::SecirParams::StageTimes::set_hospitalized_to_home))
        .def("set_hospitalized_to_home", py::overload_cast<const epi::ParameterDistribution&>(
                                             &epi::SecirParams::StageTimes::set_hospitalized_to_home))
        .def("set_home_to_hospitalized",
             py::overload_cast<double>(&epi::SecirParams::StageTimes::set_home_to_hospitalized))
        .def("set_home_to_hospitalized", py::overload_cast<const epi::ParameterDistribution&>(
                                             &epi::SecirParams::StageTimes::set_home_to_hospitalized))
        .def("set_hospitalized_to_icu",
             py::overload_cast<double>(&epi::SecirParams::StageTimes::set_hospitalized_to_icu))
        .def("set_hospitalized_to_icu", py::overload_cast<const epi::ParameterDistribution&>(
                                            &epi::SecirParams::StageTimes::set_hospitalized_to_icu))
        .def("set_icu_to_home", py::overload_cast<double>(&epi::SecirParams::StageTimes::set_icu_to_home))
        .def("set_icu_to_home",
             py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams::StageTimes::set_icu_to_home))
        .def("set_infectious_asymp", py::overload_cast<double>(&epi::SecirParams::StageTimes::set_infectious_asymp))
        .def("set_infectious_asymp",
             py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams::StageTimes::set_infectious_asymp))
        .def("set_icu_to_death", py::overload_cast<double>(&epi::SecirParams::StageTimes::set_icu_to_death))
        .def("set_icu_to_death",
             py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams::StageTimes::set_icu_to_death))

        .def("get_incubation", py::overload_cast<>(&epi::SecirParams::StageTimes::get_incubation),
             py::return_value_policy::reference_internal)
        .def("get_incubation", py::overload_cast<>(&epi::SecirParams::StageTimes::get_incubation, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_infectious_mild", py::overload_cast<>(&epi::SecirParams::StageTimes::get_infectious_mild),
             py::return_value_policy::reference_internal)
        .def("get_infectious_mild", py::overload_cast<>(&epi::SecirParams::StageTimes::get_infectious_mild, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_serialinterval", py::overload_cast<>(&epi::SecirParams::StageTimes::get_serialinterval),
             py::return_value_policy::reference_internal)
        .def("get_serialinterval", py::overload_cast<>(&epi::SecirParams::StageTimes::get_serialinterval, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_hospitalized_to_home", py::overload_cast<>(&epi::SecirParams::StageTimes::get_hospitalized_to_home),
             py::return_value_policy::reference_internal)
        .def("get_hospitalized_to_home",
             py::overload_cast<>(&epi::SecirParams::StageTimes::get_hospitalized_to_home, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_home_to_hospitalized", py::overload_cast<>(&epi::SecirParams::StageTimes::get_home_to_hospitalized),
             py::return_value_policy::reference_internal)
        .def("get_home_to_hospitalized",
             py::overload_cast<>(&epi::SecirParams::StageTimes::get_home_to_hospitalized, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_hospitalized_to_icu", py::overload_cast<>(&epi::SecirParams::StageTimes::get_hospitalized_to_icu),
             py::return_value_policy::reference_internal)
        .def("get_hospitalized_to_icu",
             py::overload_cast<>(&epi::SecirParams::StageTimes::get_hospitalized_to_icu, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_icu_to_home", py::overload_cast<>(&epi::SecirParams::StageTimes::get_icu_to_home),
             py::return_value_policy::reference_internal)
        .def("get_icu_to_home", py::overload_cast<>(&epi::SecirParams::StageTimes::get_icu_to_home, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_infectious_asymp", py::overload_cast<>(&epi::SecirParams::StageTimes::get_infectious_asymp),
             py::return_value_policy::reference_internal)
        .def("get_infectious_asymp",
             py::overload_cast<>(&epi::SecirParams::StageTimes::get_infectious_asymp, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_icu_to_dead", py::overload_cast<>(&epi::SecirParams::StageTimes::get_icu_to_dead),
             py::return_value_policy::reference_internal)
        .def("get_icu_to_dead", py::overload_cast<>(&epi::SecirParams::StageTimes::get_icu_to_dead, py::const_),
             py::return_value_policy::reference_internal);

    py::class_<epi::Populations>(m, "Populations")
        .def(py::init<std::vector<size_t>&>())
        .def("get_num_compartments", &epi::Populations::get_num_compartments)
        .def("get_category_sizes", &epi::Populations::get_category_sizes)
        .def("get_compartments", &epi::Populations::get_compartments)
        .def("get", py::overload_cast<std::vector<size_t> const&>(&epi::Populations::get),
             py::return_value_policy::reference_internal)
        .def("get", py::overload_cast<std::vector<size_t> const&>(&epi::Populations::get, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_group_total", &epi::Populations::get_group_total)
        .def("get_total", &epi::Populations::get_total)
        .def("set", py::overload_cast<std::vector<size_t> const&, double>(&epi::Populations::set))
        .def("set",
             py::overload_cast<std::vector<size_t> const&, const epi::ParameterDistribution&>(&epi::Populations::set))
        .def("set_group_total", &epi::Populations::set_group_total)
        .def("set_total", &epi::Populations::set_total)
        .def("set_difference_from_total", &epi::Populations::set_difference_from_total)
        .def("set_difference_from_group_total", &epi::Populations::set_difference_from_group_total)
        .def("get_flat_index", &epi::Populations::get_flat_index);

    py::class_<epi::SecirParams::Probabilities>(m, "Probabilities")
        .def(py::init<>())
        .def("set_infection_from_contact",
             py::overload_cast<double>(&epi::SecirParams::Probabilities::set_infection_from_contact))
        .def("set_infection_from_contact", py::overload_cast<const epi::ParameterDistribution&>(
                                               &epi::SecirParams::Probabilities::set_infection_from_contact))
        .def("set_carrier_infectability",
             py::overload_cast<double>(&epi::SecirParams::Probabilities::set_carrier_infectability))
        .def("set_carrier_infectability", py::overload_cast<const epi::ParameterDistribution&>(
                                              &epi::SecirParams::Probabilities::set_carrier_infectability))
        .def("set_asymp_per_infectious",
             py::overload_cast<double>(&epi::SecirParams::Probabilities::set_asymp_per_infectious))
        .def("set_asymp_per_infectious", py::overload_cast<const epi::ParameterDistribution&>(
                                             &epi::SecirParams::Probabilities::set_asymp_per_infectious))
        .def("set_risk_from_symptomatic",
             py::overload_cast<double>(&epi::SecirParams::Probabilities::set_risk_from_symptomatic))
        .def("set_risk_from_symptomatic", py::overload_cast<const epi::ParameterDistribution&>(
                                              &epi::SecirParams::Probabilities::set_risk_from_symptomatic))
        .def("set_hospitalized_per_infectious",
             py::overload_cast<double>(&epi::SecirParams::Probabilities::set_hospitalized_per_infectious))
        .def("set_hospitalized_per_infectious", py::overload_cast<const epi::ParameterDistribution&>(
                                                    &epi::SecirParams::Probabilities::set_hospitalized_per_infectious))
        .def("set_icu_per_hospitalized",
             py::overload_cast<double>(&epi::SecirParams::Probabilities::set_icu_per_hospitalized))
        .def("set_icu_per_hospitalized", py::overload_cast<const epi::ParameterDistribution&>(
                                             &epi::SecirParams::Probabilities::set_icu_per_hospitalized))
        .def("set_dead_per_icu", py::overload_cast<double>(&epi::SecirParams::Probabilities::set_dead_per_icu))
        .def("set_dead_per_icu",
             py::overload_cast<const epi::ParameterDistribution&>(&epi::SecirParams::Probabilities::set_dead_per_icu))

        .def("get_infection_from_contact",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_infection_from_contact),
             py::return_value_policy::reference_internal)
        .def("get_infection_from_contact",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_infection_from_contact, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_carrier_infectability",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_carrier_infectability),
             py::return_value_policy::reference_internal)
        .def("get_carrier_infectability",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_carrier_infectability, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_asymp_per_infectious",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_asymp_per_infectious),
             py::return_value_policy::reference_internal)
        .def("get_asymp_per_infectious",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_asymp_per_infectious, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_risk_from_symptomatic",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_risk_from_symptomatic),
             py::return_value_policy::reference_internal)
        .def("get_risk_from_symptomatic",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_risk_from_symptomatic, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_hospitalized_per_infectious",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_hospitalized_per_infectious),
             py::return_value_policy::reference_internal)
        .def("get_hospitalized_per_infectious",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_hospitalized_per_infectious, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_icu_per_hospitalized",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_icu_per_hospitalized),
             py::return_value_policy::reference_internal)
        .def("get_icu_per_hospitalized",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_icu_per_hospitalized, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_dead_per_icu", py::overload_cast<>(&epi::SecirParams::Probabilities::get_dead_per_icu),
             py::return_value_policy::reference_internal)
        .def("get_dead_per_icu", py::overload_cast<>(&epi::SecirParams::Probabilities::get_dead_per_icu, py::const_),
             py::return_value_policy::reference_internal);

    py::class_<epi::ContactFrequencyMatrix>(m, "ContactFrequencyMatrix")
        .def(py::init<>())
        .def(py::init<size_t>())
        .def("set_cont_freq", &epi::ContactFrequencyMatrix::set_cont_freq)
        .def("get_cont_freq", &epi::ContactFrequencyMatrix::get_cont_freq)
        .def("set_dampings", &epi::ContactFrequencyMatrix::set_dampings)
        .def("get_dampings", &epi::ContactFrequencyMatrix::get_dampings, py::return_value_policy::reference_internal)
        .def("add_damping", &epi::ContactFrequencyMatrix::add_damping);

    py::class_<epi::UncertainContactMatrix>(m, "UncertainContactMatrix")
        .def(py::init<>())
        .def(py::init<epi::ContactFrequencyMatrix>())
        .def("get_cont_freq_mat", py::overload_cast<>(&epi::UncertainContactMatrix::get_cont_freq_mat),
             py::return_value_policy::reference_internal);

    py::class_<epi::SecirParams>(m, "SecirParams")
        .def(py::init<size_t>(), py::arg("num_groups") = 1)
        .def_readwrite("times", &epi::SecirParams::times)
        .def_readwrite("populations", &epi::SecirParams::populations)
        .def_readwrite("probabilities", &epi::SecirParams::probabilities)
        .def("check_constraints", &epi::SecirParams::check_constraints)
        .def("apply_constraints", &epi::SecirParams::apply_constraints)
        .def("get_contact_patterns", py::overload_cast<>(&epi::SecirParams::get_contact_patterns),
             py::return_value_policy::reference_internal)
        .def("get_contact_patterns", py::overload_cast<>(&epi::SecirParams::get_contact_patterns, py::const_),
             py::return_value_policy::reference_internal);

    m.def("simulate", py::overload_cast<double, double, double, const epi::SecirParams&>(&epi::simulate),
          "Simulates the SECIR model from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"),
          py::arg("params"));

    py::class_<epi::SecirSimulation>(m, "SecirSimulation")
        .def(py::init<const epi::SecirParams&, double, double>(), py::arg("params"), py::arg("t0") = 0,
             py::arg("dt") = 0.1)
        .def_property_readonly("result", py::overload_cast<>(&epi::SecirSimulation::get_result, py::const_),
                               py::return_value_policy::reference_internal)
        .def("advance", &epi::SecirSimulation::advance, py::arg("tmax"));

    py::class_<epi::MigrationEdge>(m, "MigrationParams")
        .def(py::init<const Eigen::VectorXd&>(), py::arg("coeffs"))
        .def_property(
            "coefficients", [](const epi::MigrationEdge& self) -> auto& { return self.coefficients; },
            [](epi::MigrationEdge& self, const Eigen::VectorXd& v) {
                self.coefficients = v;
            },
            py::return_value_policy::reference_internal);

    py::class_<epi::Edge<epi::MigrationEdge>>(m, "MigrationEdge")
        .def_property_readonly("start_node_idx",
                               [](const epi::Edge<epi::MigrationEdge>& self) {
                                   return self.start_node_idx;
                               })
        .def_property_readonly("end_node_idx",
                               [](const epi::Edge<epi::MigrationEdge>& self) {
                                   return self.end_node_idx;
                               })
        .def_property_readonly(
            "property", [](const epi::Edge<epi::MigrationEdge>& self) -> auto& { return self.property; },
            py::return_value_policy::reference_internal);

    using SecirParamsGraph = epi::Graph<epi::SecirParams, epi::MigrationEdge>;
    py::class_<SecirParamsGraph>(m, "SecirParamsGraph")
        .def(py::init<>())
        .def("add_node", &SecirParamsGraph::add_node<const epi::SecirParams&>,
             py::return_value_policy::reference_internal)
        .def("add_edge", &SecirParamsGraph::add_edge<const epi::MigrationEdge&>,
             py::return_value_policy::reference_internal)
        .def("add_edge", &SecirParamsGraph::add_edge<const Eigen::VectorXd&>,
             py::return_value_policy::reference_internal)
        .def_property_readonly("num_nodes",
                               [](const SecirParamsGraph& self) {
                                   return self.nodes().size();
                               })
        .def(
            "get_node", [](const SecirParamsGraph& self, size_t node_idx) -> auto& { return self.nodes()[node_idx]; },
            py::return_value_policy::reference_internal)
        .def_property_readonly("num_edges",
                               [](const SecirParamsGraph& self) {
                                   return self.edges().size();
                               })
        .def(
            "get_edge", [](const SecirParamsGraph& self, size_t edge_idx) -> auto& { return self.edges()[edge_idx]; },
            py::return_value_policy::reference_internal)
        .def("get_num_out_edges",
                               [](const SecirParamsGraph& self, size_t node_idx) {
                                   return self.out_edges(node_idx).size();
                               })
        .def(
            "get_out_edge", [](const SecirParamsGraph& self, size_t node_idx, size_t edge_idx) -> auto& {
                return self.out_edges(node_idx)[edge_idx];
            },
            py::return_value_policy::reference_internal);

    using MigrationGraph = epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge>;
    py::class_<MigrationGraph>(m, "MigrationGraph")
        .def(py::init<>())
        .def(
            "add_node", [](MigrationGraph & self, const epi::SecirParams& p, double t0, double dt) -> auto& {
                return self.add_node(p, t0, dt).model;
            },
            py::arg("params"), py::arg("t0") = 0.0, py::arg("dt") = 0.1, py::return_value_policy::reference_internal)
        .def("add_edge", &MigrationGraph::add_edge<const epi::MigrationEdge&>,
             py::return_value_policy::reference_internal)
        .def("add_edge", &MigrationGraph::add_edge<const Eigen::VectorXd&>, py::return_value_policy::reference_internal)
        .def_property_readonly("num_nodes",
                               [](const MigrationGraph& self) {
                                   return self.nodes().size();
                               })
        .def(
            "get_node",
            [](const MigrationGraph& self, size_t node_idx) -> auto& { return self.nodes()[node_idx].model; },
            py::return_value_policy::reference_internal)
        .def_property_readonly("num_edges",
                               [](const MigrationGraph& self) {
                                   return self.edges().size();
                               })
        .def(
            "get_edge", [](const MigrationGraph& self, size_t edge_idx) -> auto& { return self.edges()[edge_idx]; },
            py::return_value_policy::reference_internal)
        .def("get_num_out_edges",
                               [](const MigrationGraph& self, size_t node_idx) {
                                   return self.out_edges(node_idx).size();
                               })
        .def(
            "get_out_edge", [](const MigrationGraph& self, size_t node_idx, size_t edge_idx) -> auto& {
                return self.out_edges(node_idx)[edge_idx];
            },
            py::return_value_policy::reference_internal);

    py::class_<epi::GraphSimulation<MigrationGraph>>(m, "MigrationSimulation")
        .def(py::init([](const MigrationGraph& graph, double t0, double dt) {
            return std::make_unique<epi::GraphSimulation<MigrationGraph>>(epi::make_migration_sim(t0, dt, graph));
        }), py::arg("graph"), py::arg("t0") = 0.0, py::arg("dt") = 1.0)
        .def_property_readonly("graph",
                               py::overload_cast<>(&epi::GraphSimulation<MigrationGraph>::get_graph, py::const_),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("t", &epi::GraphSimulation<MigrationGraph>::get_t)
        .def("advance", &epi::GraphSimulation<MigrationGraph>::advance, py::arg("tmax"));

    py::class_<epi::ParameterStudy>(m, "ParameterStudy")
        .def(py::init<const epi::SecirParams&, double, double, size_t>(), py::arg("params"), py::arg("t0"),
             py::arg("tmax"), py::arg("num_runs"))
        .def(py::init<const epi::SecirParams&, double, double, double, size_t>(), py::arg("params"), py::arg("t0"),
             py::arg("tmax"), py::arg("dev_rel"), py::arg("num_runs"))
        .def(py::init<const SecirParamsGraph&, double, double, double, size_t>(), py::arg("params_graph"),
             py::arg("t0"), py::arg("tmax"), py::arg("graph_dt"), py::arg("num_runs"))
        .def_property("num_runs", &epi::ParameterStudy::get_num_runs, &epi::ParameterStudy::set_num_runs)
        .def_property("tmax", &epi::ParameterStudy::get_tmax, &epi::ParameterStudy::set_tmax)
        .def_property("t0", &epi::ParameterStudy::get_t0, &epi::ParameterStudy::set_t0)
        .def_property_readonly("secir_params", py::overload_cast<>(&epi::ParameterStudy::get_secir_params),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("secir_params", py::overload_cast<>(&epi::ParameterStudy::get_secir_params, py::const_),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("secir_params_graph", py::overload_cast<>(&epi::ParameterStudy::get_secir_params_graph),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("secir_params_graph",
                               py::overload_cast<>(&epi::ParameterStudy::get_secir_params_graph, py::const_),
                               py::return_value_policy::reference_internal)
        .def("run", &epi::ParameterStudy::run)
        .def("run",
             [](epi::ParameterStudy& self) { //default argument doesn't seem to work with functions
                 return self.run();
             })
        .def(
            "run_single",
            [](epi::ParameterStudy& self, epi::HandleSimulationResultFunction handle_result_func) {
                return filter_graph_results(self.run(handle_result_func));
            },
            py::arg("handle_result_func"))
        .def("run_single", [](epi::ParameterStudy& self) {
            return filter_graph_results(self.run());
        });

    m.def("set_params_distributions_normal", &epi::set_params_distributions_normal, py::arg("params"), py::arg("t0"),
          py::arg("tmax"), py::arg("dev_rel"));

    m.def("draw_sample", &epi::draw_sample, py::arg("params"));

    m.def("interpolate_simulation_result", py::overload_cast<const epi::TimeSeries<double>&>(&epi::interpolate_simulation_result));
    m.def("interpolate_simulation_result", py::overload_cast<const MigrationGraph&>(&epi::interpolate_simulation_result));

    m.def("interpolate_ensemble_results", &epi::interpolate_ensemble_results<MigrationGraph>);
    m.def("interpolate_ensemble_results", &epi::interpolate_ensemble_results<epi::TimeSeries<double>>);

    m.def("ensemble_mean", &epi::ensemble_mean);
    m.def("ensemble_percentile", &epi::ensemble_percentile);

    m.attr("__version__") = "dev";
}
