#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <epidemiology/seirParam.h>
#include <epidemiology/secir.h>
#include <vector>

int add(int i, int j) {
    return i + j;
}

namespace py = pybind11;


class SecirResult
{
public:
    void resize(size_t size) {
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

SecirResult simulate_secir(double t0, double tmax, double dt, const seirParam<double>& params)
{
    std::vector<std::vector<double> > seir(0);
    std::vector<double> times = simulate(t0, tmax, dt, params, seir);

    SecirResult result;
    if (seir.size() < 1 || seir[0].size() != 8) {
        throw std::runtime_error("Invalid result from secir simulation");
    }

    size_t n_data = seir.size();
    result.resize(n_data);

    for (size_t irow = 0; irow < seir.size(); ++irow) {

        result.t[irow] = times[irow];
        result.nb_sus[irow] = seir[irow][0]; 
        result.nb_exp[irow] = seir[irow][1]; 
        result.nb_car[irow] = seir[irow][2]; 
        result.nb_inf[irow] = seir[irow][3]; 
        result.nb_hosp[irow] = seir[irow][4]; 
        result.nb_icu[irow] = seir[irow][5]; 
        result.nb_rec[irow] = seir[irow][6]; 
        result.nb_dead[irow] = seir[irow][7]; 
    }

    return result;
}

PYBIND11_MODULE(secir, m) {
    py::class_<damping<double>>(m, "Damping")
        .def(py::init<double, double>(), py::arg("day"), py::arg("factor"))
        .def_readwrite("day", &damping<double>::day)
        .def_readwrite("factor", &damping<double>::factor);

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

    // T tinc, T tinfmild, T tserint, T thosp2home, T thome2hosp, T thosp2icu, T ticu2home, T tinfasy, T ticu2death,
    // T cont_freq_in, T alpha_in, T beta_in, T delta_in, T rho_in, T theta_in, 
    //  T nb_total_t0_in, T nb_exp_t0_in, T nb_car_t0_in, T nb_inf_t0_in, T nb_hosp_t0_in, T nb_icu_t0_in, T nb_rec_t0_in, T nb_dead_t0_in) 
    py::class_<seirParam<double>>(m, "SeirParam")
        .def(py::init<double, double, double, double, double, double, double, double, double,
                      double, double, double, double, double, double,
                      double, double, double, double, double, double, double, double>(),
                      py::arg("tinc"), py::arg("tinfmild"), py::arg("tserint"), py::arg("thosp2home"), py::arg("thome2hosp"),
                      py::arg("thosp2icu"), py::arg("ticu2home"), py::arg("tinfasy"), py::arg("ticu2death"),
                      py::arg("cont_freq_in"), py::arg("alpha_in"), py::arg("beta_in"), py::arg("delta_in"), py::arg("rho_in"), py::arg("theta_in"),
                      py::arg("nb_total_t0_in"), py::arg("nb_exp_t0_in"), py::arg("nb_car_t0_in"), py::arg("nb_inf_t0_in"),
                      py::arg("nb_hosp_t0_in"), py::arg("nb_icu_t0_in"), py::arg("nb_rec_t0_in"), py::arg("nb_dead_t0_in"))
        .def("add_damping", &seirParam<double>::add_damping);

    m.def("print_seir_params", &printSeirParams<double>);
    m.def("simulate", &simulate_secir, "Simulates the SECIR model from t0 to tmax.",
        py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("params"));

    m.attr("__version__") = "dev";  
}
