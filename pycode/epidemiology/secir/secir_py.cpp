#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <epidemiology/secir.h>
#include <epidemiology/damping.h>
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

SecirResult simulate_secir(double t0, double tmax, double dt, const epi::ContactFrequencyMatrix& cont_freq_matrix,
                           std::vector<epi::SecirParams> const& params)
{
    std::vector<Eigen::VectorXd> seir(0);
    auto times = simulate(t0, tmax, dt, cont_freq_matrix, params, seir);

    SecirResult result;
    if (seir.size() < 1 || seir[0].size() != 8) {
        throw std::runtime_error("Invalid result from secir simulation");
    }

    size_t n_data = seir.size();
    result.resize(n_data);

    for (size_t irow = 0; irow < seir.size(); ++irow) {

        result.t[irow]       = times[irow];
        result.nb_sus[irow]  = seir[irow][0];
        result.nb_exp[irow]  = seir[irow][1];
        result.nb_car[irow]  = seir[irow][2];
        result.nb_inf[irow]  = seir[irow][3];
        result.nb_hosp[irow] = seir[irow][4];
        result.nb_icu[irow]  = seir[irow][5];
        result.nb_rec[irow]  = seir[irow][6];
        result.nb_dead[irow] = seir[irow][7];
    }

    return result;
}

PYBIND11_MODULE(_secir, m)
{
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
        .def("set_incubation", &epi::SecirParams::StageTimes::set_incubation)
        .def("set_infectious_mild", &epi::SecirParams::StageTimes::set_infectious_mild)
        .def("set_serialinterval", &epi::SecirParams::StageTimes::set_serialinterval)
        .def("set_hospitalized_to_home", &epi::SecirParams::StageTimes::set_hospitalized_to_home)
        .def("set_home_to_hospitalized", &epi::SecirParams::StageTimes::set_home_to_hospitalized)
        .def("set_hospitalized_to_icu", &epi::SecirParams::StageTimes::set_hospitalized_to_icu)
        .def("set_icu_to_home", &epi::SecirParams::StageTimes::set_icu_to_home)
        .def("set_infectious_asymp", &epi::SecirParams::StageTimes::set_infectious_asymp)
        .def("set_icu_to_death", &epi::SecirParams::StageTimes::set_icu_to_death)

        .def("get_incubation_inv", &epi::SecirParams::StageTimes::get_incubation_inv)
        .def("get_infectious_mild_inv", &epi::SecirParams::StageTimes::get_infectious_mild_inv)
        .def("get_serialinterval_inv", &epi::SecirParams::StageTimes::get_serialinterval_inv)
        .def("get_hospitalized_to_home_inv", &epi::SecirParams::StageTimes::get_hospitalized_to_home_inv)
        .def("get_home_to_hospitalized_inv", &epi::SecirParams::StageTimes::get_home_to_hospitalized_inv)
        .def("get_hospitalized_to_icu_inv", &epi::SecirParams::StageTimes::get_hospitalized_to_icu_inv)
        .def("get_icu_to_home_inv", &epi::SecirParams::StageTimes::get_icu_to_home_inv)
        .def("get_infectious_asymp_inv", &epi::SecirParams::StageTimes::get_infectious_asymp_inv)
        .def("get_icu_to_dead_inv", &epi::SecirParams::StageTimes::get_icu_to_dead_inv);

    py::class_<epi::SecirParams::Populations>(m, "Populations")
        .def(py::init<>())
        .def("set_total_t0", &epi::SecirParams::Populations::set_total_t0)
        .def("set_exposed_t0", &epi::SecirParams::Populations::set_exposed_t0)
        .def("set_carrier_t0", &epi::SecirParams::Populations::set_carrier_t0)
        .def("set_infectious_t0", &epi::SecirParams::Populations::set_infectious_t0)
        .def("set_hospital_t0", &epi::SecirParams::Populations::set_hospital_t0)
        .def("set_icu_t0", &epi::SecirParams::Populations::set_icu_t0)
        .def("set_recovered_t0", &epi::SecirParams::Populations::set_recovered_t0)
        .def("set_dead_t0", &epi::SecirParams::Populations::set_dead_t0)

        .def("get_total_t0", &epi::SecirParams::Populations::get_total_t0)
        .def("get_exposed_t0", &epi::SecirParams::Populations::get_exposed_t0)
        .def("get_carrier_t0", &epi::SecirParams::Populations::get_carrier_t0)
        .def("get_infectious_t0", &epi::SecirParams::Populations::get_infectious_t0)
        .def("get_hospitalized_t0", &epi::SecirParams::Populations::get_hospitalized_t0)
        .def("get_icu_t0", &epi::SecirParams::Populations::get_icu_t0)
        .def("get_recovered_t0", &epi::SecirParams::Populations::get_recovered_t0)
        .def("get_dead_t0", &epi::SecirParams::Populations::get_dead_t0);

    py::class_<epi::SecirParams::Probabilities>(m, "Probabilities")
        .def(py::init<>())
        .def("set_asymp_per_infectious", &epi::SecirParams::Probabilities::set_asymp_per_infectious)
        .def("set_risk_from_symptomatic", &epi::SecirParams::Probabilities::set_risk_from_symptomatic)
        .def("set_hospitalized_per_infectious", &epi::SecirParams::Probabilities::set_hospitalized_per_infectious)
        .def("set_icu_per_hospitalized", &epi::SecirParams::Probabilities::set_icu_per_hospitalized)
        .def("set_dead_per_icu", &epi::SecirParams::Probabilities::set_dead_per_icu)

        .def("get_asymp_per_infectious", &epi::SecirParams::Probabilities::get_asymp_per_infectious)
        .def("get_risk_from_symptomatic", &epi::SecirParams::Probabilities::get_risk_from_symptomatic)
        .def("get_hospitalized_per_infectious", &epi::SecirParams::Probabilities::get_hospitalized_per_infectious)
        .def("get_icu_per_hospitalized", &epi::SecirParams::Probabilities::get_icu_per_hospitalized)
        .def("get_dead_per_icu", &epi::SecirParams::Probabilities::get_dead_per_icu);

    py::class_<epi::SecirParams>(m, "SecirParams")
        .def(py::init<>())
        .def_readwrite("times", &epi::SecirParams::times)
        .def_readwrite("populations", &epi::SecirParams::populations)
        .def_readwrite("probabilities", &epi::SecirParams::probabilities);

    py::class_<epi::ContactFrequencyMatrix>(m, "ContactFrequencyMatrix")
        .def(py::init<>())
        .def("set_cont_freq", &epi::ContactFrequencyMatrix::set_cont_freq)
        .def("get_cont_freq", &epi::ContactFrequencyMatrix::get_cont_freq)
        .def("set_dampings", &epi::ContactFrequencyMatrix::set_dampings)
        .def("get_dampings", &epi::ContactFrequencyMatrix::get_dampings)
        .def("add_damping", &epi::ContactFrequencyMatrix::add_damping);

    m.def("print_secir_params", &epi::print_secir_params);
    m.def("simulate", &simulate_secir, "Simulates the SECIR model from t0 to tmax.", py::arg("t0"), py::arg("tmax"),
          py::arg("dt"), py::arg("cont_freq_matrix"), py::arg("params"));

    m.attr("__version__") = "dev";
}
