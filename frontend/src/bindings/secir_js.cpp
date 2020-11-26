#include <emscripten/bind.h>

#include <epidemiology/secir/secir.h>
#include <vector>

#ifdef DEBUG
#include <sanitizer/lsan_interface.h>
#endif

namespace js = emscripten;

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

SecirResult simulate_secir(double t0, double tmax, double dt, epi::SecirParams const& params)
{
    auto result_timeseries = simulate(t0, tmax, dt, params);

    if (result_timeseries.get_num_time_points() < 1 ||
        result_timeseries.get_num_elements() != epi::SecirCompartments::SecirCount) {
        throw std::runtime_error("Invalid result from secir simulation");
    }

    SecirResult result;
    result.resize(result_timeseries.get_num_time_points());

    for (size_t irow = 0; irow < result_timeseries.get_num_time_points(); ++irow) {
        result.t[irow]       = result_timeseries.get_times()[irow];
        result.nb_sus[irow]  = result_timeseries[irow][0];
        result.nb_exp[irow]  = result_timeseries[irow][1];
        result.nb_car[irow]  = result_timeseries[irow][2];
        result.nb_inf[irow]  = result_timeseries[irow][3];
        result.nb_hosp[irow] = result_timeseries[irow][4];
        result.nb_icu[irow]  = result_timeseries[irow][5];
        result.nb_rec[irow]  = result_timeseries[irow][6];
        result.nb_dead[irow] = result_timeseries[irow][7];
    }

    return result;
}

EMSCRIPTEN_BINDINGS(secirjs)
{
    js::class_<epi::Damping>("Damping")
        .constructor<epi::Damping(Eigen::Index, double, int, int, double)>(
            [](auto n, auto d, auto lv, auto ty, auto t) {
                return epi::Damping(n, d, epi::DampingLevel(lv), epi::DampingType(ty), epi::SimulationTime(t));
            })
        .constructor<epi::Damping(Eigen::Index, double, double)>(
            [](auto n, auto d, auto t) {
                return epi::Damping(n, d, epi::SimulationTime(t));
            })
        .property<double>(
            "time",
            [](auto& self) {
                return double(self.get_time());
            },
            [](auto& self, auto v) {
                self.get_time() = epi::SimulationTime(v);
            })
        .property<int>(
            "type",
            [](auto& self) {
                return int(self.get_type());
            },
            [](auto& self, auto v) {
                self.get_type() = epi::DampingType(v);
            })
        .property<int>(
            "level",
            [](auto& self) {
                return int(self.get_level());
            },
            [](auto& self, auto v) {
                self.get_level() = epi::DampingLevel(v);
            })
        .function("get_coeff",
                  std::function<double(epi::Damping&, Eigen::Index, Eigen::Index)>([](auto& self, auto i, auto j) {
                      return self.get_coeffs()(i, j);
                  }))
        .function("set_coeff", std::function<void(epi::Damping&, Eigen::Index, Eigen::Index, double)>(
                                   [](auto& self, auto i, auto j, double d) {
                                       self.get_coeffs()(i, j) = d;
                                   }));

    js::class_<SecirResult>("SecirResult")
        .property("t", &SecirResult::t)
        .property("nb_sus", &SecirResult::nb_sus)
        .property("nb_exp", &SecirResult::nb_exp)
        .property("nb_car", &SecirResult::nb_car)
        .property("nb_inf", &SecirResult::nb_inf)
        .property("nb_hosp", &SecirResult::nb_hosp)
        .property("nb_icu", &SecirResult::nb_icu)
        .property("nb_rec", &SecirResult::nb_rec)
        .property("nb_dead", &SecirResult::nb_dead);

    js::class_<epi::SecirParams::StageTimes>("StageTimes")
        .constructor<>()
        .function("set_incubation", js::select_overload<void(double)>(&epi::SecirParams::StageTimes::set_incubation))
        .function("set_infectious_mild",
                  js::select_overload<void(double)>(&epi::SecirParams::StageTimes::set_infectious_mild))
        .function("set_serialinterval",
                  js::select_overload<void(double)>(&epi::SecirParams::StageTimes::set_serialinterval))
        .function("set_hospitalized_to_home",
                  js::select_overload<void(double)>(&epi::SecirParams::StageTimes::set_hospitalized_to_home))
        .function("set_home_to_hospitalized",
                  js::select_overload<void(double)>(&epi::SecirParams::StageTimes::set_home_to_hospitalized))
        .function("set_hospitalized_to_icu",
                  js::select_overload<void(double)>(&epi::SecirParams::StageTimes::set_hospitalized_to_icu))
        .function("set_icu_to_home", js::select_overload<void(double)>(&epi::SecirParams::StageTimes::set_icu_to_home))
        .function("set_infectious_asymp",
                  js::select_overload<void(double)>(&epi::SecirParams::StageTimes::set_infectious_asymp))
        .function("set_icu_to_death",
                  js::select_overload<void(double)>(&epi::SecirParams::StageTimes::set_icu_to_death))

        // at the moment, no const getters available in JS
        .function("get_incubation",
                  js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::StageTimes::get_incubation))
        .function("get_infectious_mild",
                  js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::StageTimes::get_infectious_mild))
        .function("get_serialinterval",
                  js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::StageTimes::get_serialinterval))
        .function("get_hospitalized_to_home",
                  js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::StageTimes::get_hospitalized_to_home))
        .function("get_home_to_hospitalized",
                  js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::StageTimes::get_home_to_hospitalized))
        .function("get_hospitalized_to_icu",
                  js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::StageTimes::get_hospitalized_to_icu))
        .function("get_icu_to_home",
                  js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::StageTimes::get_icu_to_home))
        .function("get_infectious_asymp",
                  js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::StageTimes::get_infectious_asymp))
        .function("get_icu_to_dead",
                  js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::StageTimes::get_icu_to_dead));

    js::class_<epi::Populations>("Populations")
        .constructor<std::vector<size_t>&>()
        .function("get_num_compartments", &epi::Populations::get_num_compartments)
        .function("get_category_sizes", &epi::Populations::get_category_sizes)
        .function("get_compartments", &epi::Populations::get_compartments)
        .function("get", js::select_overload<epi::UncertainValue&(std::vector<size_t> const&)>(&epi::Populations::get))
        .function("get_group_total", &epi::Populations::get_group_total)
        .function("get_total", &epi::Populations::get_total)
        .function("set", js::select_overload<void(std::vector<size_t> const&, double)>(&epi::Populations::set))
        // .function("set", js::select_overload<void(std::vector<size_t> const&, &epi::ParameterDistribution)>(&epi::Populations::set))
        .function("set_group_total", &epi::Populations::set_group_total)
        .function("set_difference_from_total", &epi::Populations::set_difference_from_total)
        .function("set_difference_from_group_total", &epi::Populations::set_difference_from_group_total)
        .function("set_total", &epi::Populations::set_total)
        .function("get_flat_index", &epi::Populations::get_flat_index);

    js::class_<epi::SecirParams::Probabilities>("Probabilities")
        .constructor<>()
        .function("set_infection_from_contact",
                  js::select_overload<void(double)>(&epi::SecirParams::Probabilities::set_infection_from_contact))
        .function("set_carrier_infectability",
                  js::select_overload<void(double)>(&epi::SecirParams::Probabilities::set_carrier_infectability))
        .function("set_asymp_per_infectious",
                  js::select_overload<void(double)>(&epi::SecirParams::Probabilities::set_asymp_per_infectious))
        .function("set_risk_from_symptomatic",
                  js::select_overload<void(double)>(&epi::SecirParams::Probabilities::set_risk_from_symptomatic))
        .function("set_test_and_trace_max_risk_from_symptomatic",
                  js::select_overload<void(double)>(
                      &epi::SecirParams::Probabilities::set_test_and_trace_max_risk_from_symptomatic))
        .function("set_hospitalized_per_infectious",
                  js::select_overload<void(double)>(&epi::SecirParams::Probabilities::set_hospitalized_per_infectious))
        .function("set_icu_per_hospitalized",
                  js::select_overload<void(double)>(&epi::SecirParams::Probabilities::set_icu_per_hospitalized))
        .function("set_dead_per_icu",
                  js::select_overload<void(double)>(&epi::SecirParams::Probabilities::set_dead_per_icu))

        // at the moment, no const getters available in JS
        .function("get_infection_from_contact", js::select_overload<epi::UncertainValue&()>(
                                                    &epi::SecirParams::Probabilities::get_infection_from_contact))
        .function("get_carrier_infectability", js::select_overload<epi::UncertainValue&()>(
                                                   &epi::SecirParams::Probabilities::get_carrier_infectability))
        .function("get_asymp_per_infectious", js::select_overload<epi::UncertainValue&()>(
                                                  &epi::SecirParams::Probabilities::get_asymp_per_infectious))
        .function("get_risk_from_symptomatic", js::select_overload<epi::UncertainValue&()>(
                                                   &epi::SecirParams::Probabilities::get_risk_from_symptomatic))
        .function("get_test_and_trace_max_risk_from_symptomatic",
                  js::select_overload<epi::UncertainValue&()>(
                      &epi::SecirParams::Probabilities::get_test_and_trace_max_risk_from_symptomatic))
        .function("get_hospitalized_per_infectious",
                  js::select_overload<epi::UncertainValue&()>(
                      &epi::SecirParams::Probabilities::get_hospitalized_per_infectious))
        .function("get_icu_per_hospitalized", js::select_overload<epi::UncertainValue&()>(
                                                  &epi::SecirParams::Probabilities::get_icu_per_hospitalized))
        .function("get_dead_per_icu",
                  js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::Probabilities::get_dead_per_icu));

    js::class_<epi::ContactMatrix>("ContactMatrix")
        .constructor<Eigen::Index>()
        .function("add_damping", std::function<void(epi::ContactMatrix&, const epi::Damping&)>([](auto& self, auto& d) {
                      self.add_damping(d);
                  }))
        .function("get_num_dampings", std::function<size_t(epi::ContactMatrix&)>([](auto& self) {
                      return self.get_dampings().size();
                  }))
        .function("get_damping", std::function<epi::Damping&(epi::ContactMatrix&, size_t)>(
                                     [](auto& self, auto i) -> auto& { return self.get_dampings()[i]; }))
        .function("get_factor", std::function<double(epi::Dampings&, double, Eigen::Index, Eigen::Index)>(
                                    [](auto& self, auto t, auto i, auto j) {
                                        return self.get_matrix_at(t)(i, j);
                                    }))
        .function("get_baseline",
                  std::function<double(epi::ContactMatrix&, Eigen::Index, Eigen::Index)>([](auto& self, auto i, auto j) {
                      return self.get_baseline()(i, j);
                  }))
        .function("set_baseline",
                  std::function<void(epi::ContactMatrix&, Eigen::Index, Eigen::Index, double)>([](auto& self, auto i, auto j, auto d) {
                      self.get_baseline()(i, j) = d;
                  }))
        .function("get_minimum",
                  std::function<double(epi::ContactMatrix&, Eigen::Index, Eigen::Index)>([](auto& self, auto i, auto j) {
                      return self.get_minimum()(i, j);
                  }))
        .function("set_minimum",
                  std::function<void(epi::ContactMatrix&, Eigen::Index, Eigen::Index, double)>([](auto& self, auto i, auto j, auto d) {
                      self.get_minimum()(i, j) = d;
                  }));

    js::class_<epi::ContactMatrixGroup>("ContactMatrixGroup")
        .constructor<Eigen::Index, size_t>()
        .function("get_num_matrices", std::function<size_t(epi::ContactMatrixGroup&)>([](auto& self) {
                      return self.get_num_matrices();
                  }))
        .function("get_matrix", std::function<epi::ContactMatrix&(epi::ContactMatrixGroup&, size_t)>(
                                    [](auto& self, auto i) -> auto& { return self[i]; }))
        .function("get_factor", std::function<double(epi::ContactMatrixGroup&, double, Eigen::Index, Eigen::Index)>(
                                    [](auto& self, auto t, auto i, auto j) {
                                        return self.get_matrix_at(t)(i, j);
                                    }));

    js::class_<epi::UncertainContactMatrix>("UncertainContactMatrix")
        .constructor<const epi::ContactMatrixGroup&>()
        .function("get_cont_freq_mat",
                  js::select_overload<epi::ContactMatrixGroup&()>(&epi::UncertainContactMatrix::get_cont_freq_mat));

    js::class_<epi::SecirParams>("SecirParams")
        .constructor<const epi::ContactMatrixGroup&>()
        .property("times", &epi::SecirParams::times)
        .property("populations", &epi::SecirParams::populations)
        .property("probabilities", &epi::SecirParams::probabilities)
        .function("set_icu_capacity", js::select_overload<void(double)>(&epi::SecirParams::set_icu_capacity))
        .function("get_icu_capacity", js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::get_icu_capacity))
        .function("set_test_and_trace_capacity",
                  js::select_overload<void(double)>(&epi::SecirParams::set_test_and_trace_capacity))
        .function("get_test_and_trace_capacity",
                  js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::get_test_and_trace_capacity))
        .function("set_start_day", &epi::SecirParams::set_start_day)
        .function("get_start_day", &epi::SecirParams::get_start_day)
        .function("set_seasonality", js::select_overload<void(double)>(&epi::SecirParams::set_seasonality))
        .function("get_seasonality", js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::get_seasonality))
        .function("get_contact_patterns",
                  js::select_overload<epi::UncertainContactMatrix&()>(&epi::SecirParams::get_contact_patterns))
        .function("set_contact_patterns", &epi::SecirParams::set_contact_patterns);

    js::function("simulate", &simulate_secir);

    js::enum_<epi::SecirCompartments>("SecirCompartments")
        .value("E", epi::SecirCompartments::E)
        .value("S", epi::SecirCompartments::S)
        .value("C", epi::SecirCompartments::C)
        .value("I", epi::SecirCompartments::I)
        .value("H", epi::SecirCompartments::H)
        .value("U", epi::SecirCompartments::U)
        .value("R", epi::SecirCompartments::R)
        .value("D", epi::SecirCompartments::D)
        .value("SecirCount", epi::SecirCompartments::SecirCount);

    js::register_vector<double>("VectorDouble");
    js::register_vector<size_t>("VectorSizeT");
    js::register_vector<epi::SecirParams>("VectorSecirParams");
    js::register_vector<epi::SecirParams::Probabilities>("VectorSecirParamsProbabilities");
    js::register_vector<epi::SecirParams::StageTimes>("VectorSecirParamsStageTimes");

#ifdef DEBUG
    /**
     * Checks if all memory allocated by the WebAssembly module is freed. If not it will emit a warning with relevant
     * information about a potential leak.
     */
    js::function("doLeakCheck", &__lsan_do_recoverable_leak_check);
#endif
}
