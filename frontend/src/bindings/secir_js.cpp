#include <emscripten/bind.h>

#include <epidemiology/secir/secir.h>
#include <epidemiology/secir/damping.h>
#include <vector>

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
    std::vector<Eigen::VectorXd> seir(0);
    auto times = simulate(t0, tmax, dt, params, seir);

    SecirResult result;
    if (seir.size() < 1 || seir[0].size() != epi::SecirCompartments::SecirCount) {
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

EMSCRIPTEN_BINDINGS(secirjs)
{
    js::class_<epi::Damping>("Damping")
        .constructor<double, double>()
        .property("day", &epi::Damping::day)
        .property("factor", &epi::Damping::factor);

    js::class_<epi::Dampings>("Dampings")
        .constructor<>()
        .function("add", &epi::Dampings::add)
        .function("get_factor", &epi::Dampings::get_factor);

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
        .function("set_asymp_per_infectious",
                  js::select_overload<void(double)>(&epi::SecirParams::Probabilities::set_asymp_per_infectious))
        .function("set_risk_from_symptomatic",
                  js::select_overload<void(double)>(&epi::SecirParams::Probabilities::set_risk_from_symptomatic))
        .function("set_hospitalized_per_infectious",
                  js::select_overload<void(double)>(&epi::SecirParams::Probabilities::set_hospitalized_per_infectious))
        .function("set_icu_per_hospitalized",
                  js::select_overload<void(double)>(&epi::SecirParams::Probabilities::set_icu_per_hospitalized))
        .function("set_dead_per_icu",
                  js::select_overload<void(double)>(&epi::SecirParams::Probabilities::set_dead_per_icu))

        // at the moment, no const getters available in JS
        .function("get_infection_from_contact", js::select_overload<epi::UncertainValue&()>(
                                                    &epi::SecirParams::Probabilities::get_infection_from_contact))
        .function("get_asymp_per_infectious", js::select_overload<epi::UncertainValue&()>(
                                                  &epi::SecirParams::Probabilities::get_asymp_per_infectious))
        .function("get_risk_from_symptomatic", js::select_overload<epi::UncertainValue&()>(
                                                   &epi::SecirParams::Probabilities::get_risk_from_symptomatic))
        .function("get_hospitalized_per_infectious",
                  js::select_overload<epi::UncertainValue&()>(
                      &epi::SecirParams::Probabilities::get_hospitalized_per_infectious))
        .function("get_icu_per_hospitalized", js::select_overload<epi::UncertainValue&()>(
                                                  &epi::SecirParams::Probabilities::get_icu_per_hospitalized))
        .function("get_dead_per_icu",
                  js::select_overload<epi::UncertainValue&()>(&epi::SecirParams::Probabilities::get_dead_per_icu));

    js::class_<epi::ContactFrequencyMatrix>("ContactFrequencyMatrix")
        .constructor<>()
        .function("set_cont_freq", &epi::ContactFrequencyMatrix::set_cont_freq)
        .function("get_cont_freq", &epi::ContactFrequencyMatrix::get_cont_freq)
        .function("set_dampings", &epi::ContactFrequencyMatrix::set_dampings)
        .function("get_dampings", &epi::ContactFrequencyMatrix::get_dampings)
        .function("add_damping", &epi::ContactFrequencyMatrix::add_damping);

    js::class_<epi::UncertainContactMatrix>("UncertainContactMatrix")
        .constructor<epi::ContactFrequencyMatrix>()
        .function("get_cont_freq_mat",
                  js::select_overload<epi::ContactFrequencyMatrix&()>(&epi::UncertainContactMatrix::get_cont_freq_mat));

    js::class_<epi::SecirParams>("SecirParams")
        .constructor<epi::ContactFrequencyMatrix&>()
        .property("times", &epi::SecirParams::times)
        .property("populations", &epi::SecirParams::populations)
        .property("probabilities", &epi::SecirParams::probabilities)
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
}
