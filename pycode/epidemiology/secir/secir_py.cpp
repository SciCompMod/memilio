#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <epidemiology/secir/secir.h>
#include <epidemiology/secir/damping.h>
#include <vector>

namespace py = pybind11;

class SecirResult
{
public:
    void resize(size_t size)
    {
        t.resize(size, 0.);
        nb_sus.resize(size, 0.);
        nb_exp.resize(size, 0.);
        nb_car.resize(size, 0.);
        nb_inf.resize(size, 0.);
        nb_hosp.resize(size, 0.);
        nb_icu.resize(size, 0.);
        nb_rec.resize(size, 0.);
        nb_dead.resize(size, 0.);
    }

    std::vector<double> t;
    std::vector<double> nb_sus;
    std::vector<double> nb_exp;
    std::vector<double> nb_car;
    std::vector<double> nb_inf;
    std::vector<double> nb_hosp;
    std::vector<double> nb_icu;
    std::vector<double> nb_rec;
    std::vector<double> nb_dead;
};

std::vector<SecirResult> simulate_secir(double t0, double tmax, double dt, epi::SecirParams const& params)
{
    std::vector<Eigen::VectorXd> seir(0);
    auto times = simulate(t0, tmax, dt, params, seir);

    SecirResult result;
    printf("seir size is: %zd, %zd", seir.size(), seir[0].size());
    if (seir.size() < 1 || seir[0].size() % epi::SecirCompartments::SecirCount != 0) {
        throw std::runtime_error("Invalid result from secir simulation");
    }

    size_t n_data = seir.size();
    result.resize(n_data);

    int nb_groups = seir[0].size() / epi::SecirCompartments::SecirCount;

    for (size_t irow = 0; irow < seir.size(); ++irow) {

        result.t[irow]       = times[irow];
        result.nb_sus[irow]  = 0;
        result.nb_exp[irow]  = 0;
        result.nb_car[irow]  = 0;
        result.nb_inf[irow]  = 0;
        result.nb_hosp[irow] = 0;
        result.nb_icu[irow]  = 0;
        result.nb_rec[irow]  = 0;
        result.nb_dead[irow] = 0;
        for (size_t groups = 0; groups < nb_groups; ++groups) {

            result.nb_sus[irow] += seir[irow][0 + groups * epi::SecirCompartments::SecirCount];
            result.nb_exp[irow] += seir[irow][1 + groups * epi::SecirCompartments::SecirCount];
            result.nb_car[irow] += seir[irow][2 + groups * epi::SecirCompartments::SecirCount];
            result.nb_inf[irow] += seir[irow][3 + groups * epi::SecirCompartments::SecirCount];
            result.nb_hosp[irow] += seir[irow][4 + groups * epi::SecirCompartments::SecirCount];
            result.nb_icu[irow] += seir[irow][5 + groups * epi::SecirCompartments::SecirCount];
            result.nb_rec[irow] += seir[irow][6 + groups * epi::SecirCompartments::SecirCount];
            result.nb_dead[irow] += seir[irow][7 + groups * epi::SecirCompartments::SecirCount];
        }
    }

    std::vector<SecirResult> GroupResult;

    for (int i = 0; i < nb_groups; ++i) {
        SecirResult temp_result;
        size_t n_data = seir.size();
        temp_result.resize(n_data);
        for (size_t irow = 0; irow < seir.size(); ++irow) {
            temp_result.t[irow]       = times[irow];
            temp_result.nb_sus[irow]  = seir[irow][0 + i * epi::SecirCompartments::SecirCount];
            temp_result.nb_exp[irow]  = seir[irow][1 + i * epi::SecirCompartments::SecirCount];
            temp_result.nb_car[irow]  = seir[irow][2 + i * epi::SecirCompartments::SecirCount];
            temp_result.nb_inf[irow]  = seir[irow][3 + i * epi::SecirCompartments::SecirCount];
            temp_result.nb_hosp[irow] = seir[irow][4 + i * epi::SecirCompartments::SecirCount];
            temp_result.nb_icu[irow]  = seir[irow][5 + i * epi::SecirCompartments::SecirCount];
            temp_result.nb_rec[irow]  = seir[irow][6 + i * epi::SecirCompartments::SecirCount];
            temp_result.nb_dead[irow] = seir[irow][7 + i * epi::SecirCompartments::SecirCount];
        }
        GroupResult.push_back(temp_result);
    }

    return GroupResult;
}

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

    py::class_<SecirResult>(m, "SecirResult")
        .def_readonly("t", &SecirResult::t)
        .def_readonly("nb_sus", &SecirResult::nb_sus)
        .def_readonly("nb_exp", &SecirResult::nb_exp)
        .def_readonly("nb_car", &SecirResult::nb_car)
        .def_readonly("nb_inf", &SecirResult::nb_inf)
        .def_readonly("nb_hosp", &SecirResult::nb_hosp)
        .def_readonly("nb_icu", &SecirResult::nb_icu)
        .def_readonly("nb_rec", &SecirResult::nb_rec)
        .def_readonly("nb_dead", &SecirResult::nb_dead);

    py::class_<epi::SecirParams::StageTimes>(m, "StageTimes")
        .def(py::init<>())
        .def("set_incubation", py::overload_cast<double>(&epi::SecirParams::StageTimes::set_incubation))
        // .def("set_incubation", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::StageTimes::set_incubation)) // invalid template argument?? To do..
        .def("set_infectious_mild", py::overload_cast<double>(&epi::SecirParams::StageTimes::set_infectious_mild))
        // .def("set_infectious_mild", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::StageTimes::set_infectious_mild))
        .def("set_serialinterval", py::overload_cast<double>(&epi::SecirParams::StageTimes::set_serialinterval))
        // .def("set_serialinterval", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::StageTimes::set_serialinterval))
        .def("set_hospitalized_to_home",
             py::overload_cast<double>(&epi::SecirParams::StageTimes::set_hospitalized_to_home))
        // .def("set_hospitalized_to_home", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::StageTimes::set_hospitalized_to_home))
        .def("set_home_to_hospitalized",
             py::overload_cast<double>(&epi::SecirParams::StageTimes::set_home_to_hospitalized))
        // .def("set_home_to_hospitalized", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::StageTimes::set_home_to_hospitalized))
        .def("set_hospitalized_to_icu",
             py::overload_cast<double>(&epi::SecirParams::StageTimes::set_hospitalized_to_icu))
        // .def("set_hospitalized_to_icu", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::StageTimes::set_hospitalized_to_icu))
        .def("set_icu_to_home", py::overload_cast<double>(&epi::SecirParams::StageTimes::set_icu_to_home))
        // .def("set_icu_to_home", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::StageTimes::set_icu_to_home))
        .def("set_infectious_asymp", py::overload_cast<double>(&epi::SecirParams::StageTimes::set_infectious_asymp))
        // .def("set_infectious_asymp", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::StageTimes::set_infectious_asymp))
        .def("set_icu_to_death", py::overload_cast<double>(&epi::SecirParams::StageTimes::set_icu_to_death))
        // .def("set_icu_to_death", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::StageTimes::set_icu_to_death))

        .def("get_incubation", py::overload_cast<>(&epi::SecirParams::StageTimes::get_incubation))
        .def("get_incubation", py::overload_cast<>(&epi::SecirParams::StageTimes::get_incubation, py::const_))
        .def("get_infectious_mild", py::overload_cast<>(&epi::SecirParams::StageTimes::get_infectious_mild))
        .def("get_infectious_mild", py::overload_cast<>(&epi::SecirParams::StageTimes::get_infectious_mild, py::const_))
        .def("get_serialinterval", py::overload_cast<>(&epi::SecirParams::StageTimes::get_serialinterval))
        .def("get_serialinterval", py::overload_cast<>(&epi::SecirParams::StageTimes::get_serialinterval, py::const_))
        .def("get_hospitalized_to_home", py::overload_cast<>(&epi::SecirParams::StageTimes::get_hospitalized_to_home))
        .def("get_hospitalized_to_home",
             py::overload_cast<>(&epi::SecirParams::StageTimes::get_hospitalized_to_home, py::const_))
        .def("get_home_to_hospitalized", py::overload_cast<>(&epi::SecirParams::StageTimes::get_home_to_hospitalized))
        .def("get_home_to_hospitalized",
             py::overload_cast<>(&epi::SecirParams::StageTimes::get_home_to_hospitalized, py::const_))
        .def("get_hospitalized_to_icu", py::overload_cast<>(&epi::SecirParams::StageTimes::get_hospitalized_to_icu))
        .def("get_hospitalized_to_icu",
             py::overload_cast<>(&epi::SecirParams::StageTimes::get_hospitalized_to_icu, py::const_))
        .def("get_icu_to_home", py::overload_cast<>(&epi::SecirParams::StageTimes::get_icu_to_home))
        .def("get_icu_to_home", py::overload_cast<>(&epi::SecirParams::StageTimes::get_icu_to_home, py::const_))
        .def("get_infectious_asymp", py::overload_cast<>(&epi::SecirParams::StageTimes::get_infectious_asymp))
        .def("get_infectious_asymp",
             py::overload_cast<>(&epi::SecirParams::StageTimes::get_infectious_asymp, py::const_))
        .def("get_icu_to_dead", py::overload_cast<>(&epi::SecirParams::StageTimes::get_icu_to_dead))
        .def("get_icu_to_dead", py::overload_cast<>(&epi::SecirParams::StageTimes::get_icu_to_dead, py::const_));

    py::class_<epi::Populations>(m, "Populations")
        .def(py::init<std::vector<size_t>&>())
        .def("get_num_compartments", &epi::Populations::get_num_compartments)
        .def("get_category_sizes", &epi::Populations::get_category_sizes)
        .def("get_compartments", &epi::Populations::get_compartments)
        .def("get", py::overload_cast<std::vector<size_t> const&>(&epi::Populations::get))
        .def("get", py::overload_cast<std::vector<size_t> const&>(&epi::Populations::get, py::const_))
        .def("get_group_total", &epi::Populations::get_group_total)
        .def("get_total", &epi::Populations::get_total)
        .def("set", py::overload_cast<std::vector<size_t> const&, double>(&epi::Populations::set))
        // .def("set", py::overload_cast<std::vector<size_t> const&, &epi::ParameterDistribution>(&epi::Populations::set))
        .def("set_group_total", &epi::Populations::set_group_total)
        .def("set_total", &epi::Populations::set_total)
        .def("set_difference_from_total", &epi::Populations::set_difference_from_total)
        .def("set_difference_from_group_total", &epi::Populations::set_difference_from_group_total)
        .def("get_flat_index", &epi::Populations::get_flat_index);

    py::class_<epi::SecirParams::Probabilities>(m, "Probabilities")
        .def(py::init<>())
        .def("set_infection_from_contact",
             py::overload_cast<double>(&epi::SecirParams::Probabilities::set_infection_from_contact))
        // .def("set_infection_from_contact", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::Probabilities::set_infection_from_contact))
        .def("set_asymp_per_infectious",
             py::overload_cast<double>(&epi::SecirParams::Probabilities::set_asymp_per_infectious))
        // .def("set_asymp_per_infectious", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::Probabilities::set_asymp_per_infectious))
        .def("set_risk_from_symptomatic",
             py::overload_cast<double>(&epi::SecirParams::Probabilities::set_risk_from_symptomatic))
        // .def("set_risk_from_symptomatic", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::Probabilities::set_risk_from_symptomatic))
        .def("set_hospitalized_per_infectious",
             py::overload_cast<double>(&epi::SecirParams::Probabilities::set_hospitalized_per_infectious))
        // .def("set_hospitalized_per_infectious", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::Probabilities::set_hospitalized_per_infectious))
        .def("set_icu_per_hospitalized",
             py::overload_cast<double>(&epi::SecirParams::Probabilities::set_icu_per_hospitalized))
        // .def("set_icu_per_hospitalized", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::Probabilities::set_icu_per_hospitalized))
        .def("set_dead_per_icu", py::overload_cast<double>(&epi::SecirParams::Probabilities::set_dead_per_icu))
        // .def("set_dead_per_icu", py::overload_cast<&epi::ParameterDistribution>(&epi::SecirParams::Probabilities::set_dead_per_icu))

        .def("get_infection_from_contact",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_infection_from_contact))
        .def("get_infection_from_contact",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_infection_from_contact, py::const_))
        .def("get_asymp_per_infectious",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_asymp_per_infectious))
        .def("get_asymp_per_infectious",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_asymp_per_infectious, py::const_))
        .def("get_risk_from_symptomatic",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_risk_from_symptomatic))
        .def("get_risk_from_symptomatic",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_risk_from_symptomatic, py::const_))
        .def("get_hospitalized_per_infectious",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_hospitalized_per_infectious))
        .def("get_hospitalized_per_infectious",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_hospitalized_per_infectious, py::const_))
        .def("get_icu_per_hospitalized",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_icu_per_hospitalized))
        .def("get_icu_per_hospitalized",
             py::overload_cast<>(&epi::SecirParams::Probabilities::get_icu_per_hospitalized, py::const_))
        .def("get_dead_per_icu", py::overload_cast<>(&epi::SecirParams::Probabilities::get_dead_per_icu))
        .def("get_dead_per_icu", py::overload_cast<>(&epi::SecirParams::Probabilities::get_dead_per_icu, py::const_));

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
        .def("get_cont_freq_mat", py::overload_cast<>(&epi::UncertainContactMatrix::get_cont_freq_mat));

    py::class_<epi::SecirParams>(m, "SecirParams")
        .def(py::init<size_t>())
        .def_readwrite("times", &epi::SecirParams::times)
        .def_readwrite("populations", &epi::SecirParams::populations)
        .def_readwrite("probabilities", &epi::SecirParams::probabilities)
        .def("get_contact_patterns", py::overload_cast<>(&epi::SecirParams::get_contact_patterns))
        .def("get_contact_patterns", py::overload_cast<>(&epi::SecirParams::get_contact_patterns, py::const_));

    m.def("simulate", &simulate_secir, "Simulates the SECIR model from t0 to tmax.", py::arg("t0"), py::arg("tmax"),
          py::arg("dt"), py::arg("params"));

    m.attr("__version__") = "dev";
}
