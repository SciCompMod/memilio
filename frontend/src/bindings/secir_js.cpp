#include <emscripten/bind.h>

#include <epidemiology/secir/secir.h>
#include <epidemiology/model/populations.h>
#include <epidemiology/model/simulation.h>
#include <epidemiology/secir/damping.h>
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

template <class Model>
SecirResult simulate_secir(double t0, double tmax, double dt, Model const& model)
{
    auto result_timeseries = simulate(t0, tmax, dt, model);

    if (result_timeseries.get_num_time_points() < 1 ||
        result_timeseries.get_num_elements() != (size_t)epi::InfectionType::Count) {
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

// the following functions help bind class template realizations
//https://stackoverflow.com/questions/64552878/how-can-i-automatically-bind-templated-member-functions-of-variadic-class-templa
template <typename T>
std::string pretty_name()
{
    std::ostringstream o;
    o << typeid(T).name();
    return o.str();
}
template <> std::string pretty_name<epi::AgeGroup1>(){ return "AgeGroup"; }
template <> std::string pretty_name<epi::AgeGroup2>(){ return "AgeGroup"; }
template <> std::string pretty_name<epi::AgeGroup3>(){ return "AgeGroup"; }
template <> std::string pretty_name<epi::InfectionType>(){ return "InfectionType"; }

template <class C>
void bind_populations_members_for_all_cats(js::class_<C>&){}

template <class C, class T, class... Ts>
void bind_populations_members_for_all_cats(js::class_<C>& c)
{
    std::string tname = pretty_name<T>();
    c.function(("set_difference_from_group_total_" + tname).c_str(), &C::template set_difference_from_group_total<T>)
     .function(("set_group_total_" + tname).c_str(), &C::template set_group_total<T>)
     .function(("get_group_total_" + tname).c_str(), &C::template get_group_total<T>);

    // recursively bind the member for each type
    bind_populations_members_for_all_cats<C, Ts...>(c);
}

/*
 * @brief bind Populations class template for any choice of categories
 */
template<class... Cats>
void bind_Populations(std::string const& name)
{
    using Class = epi::Populations<Cats...>;
    js::class_<Class> c(name.c_str());
    c.template constructor<>()
        .class_function("get_num_compartments", &Class::get_num_compartments)
        .function("get_compartments", &Class::get_compartments)
        .function("get", js::select_overload<typename Class::Type&(Cats...)>(&Class::get))
        .function("get_total", &Class::get_total)
        .function("set", js::select_overload<void(ScalarType, Cats...)>(&Class::set))
        .function("set_difference_from_total", &Class::set_difference_from_total)
        .function("set_total", &Class::set_total)
        .function("get_flat_index", &Class::get_flat_index);

        //get_group_total, set_group_total and set_difference_from_group_total
        bind_populations_members_for_all_cats<Class, Cats...>(c);
}

template<int N>
void bind_StageTimes(std::string const& name)
{
    using Class = typename epi::SecirParams<N>::StageTimes;
    js::class_<Class>(name.c_str())
        .template constructor<>()
        .function("set_incubation", js::select_overload<void(double)>(&Class::set_incubation))
        .function("set_infectious_mild",
                  js::select_overload<void(double)>(&Class::set_infectious_mild))
        .function("set_serialinterval",
                  js::select_overload<void(double)>(&Class::set_serialinterval))
        .function("set_hospitalized_to_home",
                  js::select_overload<void(double)>(&Class::set_hospitalized_to_home))
        .function("set_home_to_hospitalized",
                  js::select_overload<void(double)>(&Class::set_home_to_hospitalized))
        .function("set_hospitalized_to_icu",
                  js::select_overload<void(double)>(&Class::set_hospitalized_to_icu))
        .function("set_icu_to_home", js::select_overload<void(double)>(&Class::set_icu_to_home))
        .function("set_infectious_asymp",
                  js::select_overload<void(double)>(&Class::set_infectious_asymp))
        .function("set_icu_to_death",
                  js::select_overload<void(double)>(&Class::set_icu_to_death))

        // at the moment, no const getters available in JS
        .function("get_incubation",
                  js::select_overload<epi::UncertainValue&()>(&Class::get_incubation))
        .function("get_infectious_mild",
                  js::select_overload<epi::UncertainValue&()>(&Class::get_infectious_mild))
        .function("get_serialinterval",
                  js::select_overload<epi::UncertainValue&()>(&Class::get_serialinterval))
        .function("get_hospitalized_to_home",
                  js::select_overload<epi::UncertainValue&()>(&Class::get_hospitalized_to_home))
        .function("get_home_to_hospitalized",
                  js::select_overload<epi::UncertainValue&()>(&Class::get_home_to_hospitalized))
        .function("get_hospitalized_to_icu",
                  js::select_overload<epi::UncertainValue&()>(&Class::get_hospitalized_to_icu))
        .function("get_icu_to_home",
                  js::select_overload<epi::UncertainValue&()>(&Class::get_icu_to_home))
        .function("get_infectious_asymp",
                  js::select_overload<epi::UncertainValue&()>(&Class::get_infectious_asymp))
        .function("get_icu_to_dead",
                  js::select_overload<epi::UncertainValue&()>(&Class::get_icu_to_dead));
}

template<int N>
void bind_Probabilities(std::string const& name)
{
    using Class = typename epi::SecirParams<N>::Probabilities;
    js::class_<Class>(name.c_str())
        .template constructor<>()
        .function("set_infection_from_contact",
                  js::select_overload<void(double)>(&Class::set_infection_from_contact))
        .function("set_carrier_infectability",
                  js::select_overload<void(double)>(&Class::set_carrier_infectability))
        .function("set_asymp_per_infectious",
                  js::select_overload<void(double)>(&Class::set_asymp_per_infectious))
        .function("set_risk_from_symptomatic",
                  js::select_overload<void(double)>(&Class::set_risk_from_symptomatic))
        .function("set_hospitalized_per_infectious",
                  js::select_overload<void(double)>(&Class::set_hospitalized_per_infectious))
        .function("set_icu_per_hospitalized",
                  js::select_overload<void(double)>(&Class::set_icu_per_hospitalized))
        .function("set_dead_per_icu",
                  js::select_overload<void(double)>(&Class::set_dead_per_icu))

        // at the moment, no const getters available in JS
        .function("get_infection_from_contact", js::select_overload<epi::UncertainValue&()>(
                                                    &Class::get_infection_from_contact))
        .function("get_carrier_infectability", js::select_overload<epi::UncertainValue&()>(
                                                   &Class::get_carrier_infectability))
        .function("get_asymp_per_infectious", js::select_overload<epi::UncertainValue&()>(
                                                  &Class::get_asymp_per_infectious))
        .function("get_risk_from_symptomatic", js::select_overload<epi::UncertainValue&()>(
                                                   &Class::get_risk_from_symptomatic))
        .function("get_hospitalized_per_infectious",
                  js::select_overload<epi::UncertainValue&()>(
                      &Class::get_hospitalized_per_infectious))
        .function("get_icu_per_hospitalized", js::select_overload<epi::UncertainValue&()>(
                                                  &Class::get_icu_per_hospitalized))
        .function("get_dead_per_icu",
                  js::select_overload<epi::UncertainValue&()>(&Class::get_dead_per_icu));
}

template <int N>
void bind_SecirParams(std::string const& name)
{
    using Class = typename epi::SecirParams<N>;
    js::class_<Class>(name.c_str())
        .template constructor<epi::ContactFrequencyMatrix&>()
        .property("times", &Class::times)
        .property("probabilities", &Class::probabilities)
        .function("set_icu_capacity", js::select_overload<void(double)>(&Class::set_icu_capacity))
        .function("get_icu_capacity", js::select_overload<epi::UncertainValue&()>(&Class::get_icu_capacity))
        .function("set_start_day", &Class::set_start_day)
        .function("get_start_day", &Class::get_start_day)
        .function("set_seasonality", js::select_overload<void(double)>(&Class::set_seasonality))
        .function("get_seasonality", js::select_overload<epi::UncertainValue&()>(&Class::get_seasonality))
        .function("get_contact_patterns",
                  js::select_overload<epi::UncertainContactMatrix&()>(&Class::get_contact_patterns))
        .function("set_contact_patterns", &Class::set_contact_patterns);

}

template<class Populations, class Parameters>
void bind_CompartmentalModel(std::string const& name)
{
    using Class = epi::CompartmentalModel<Populations, Parameters>;
    js::class_<Class>(name.c_str())
        .template constructor<>()
        .function("apply_constraints", &Class::apply_constraints)
        .function("check_constraints", &Class::check_constraints)
        .function("get_initial_values", &Class::get_initial_values)
        .property("populations", &Class::populations)
        .property("parameters", &Class::parameters);
}

template<class AgeGroup>
void bind_SecirModel(std::string const& name)
{

    using Pa = epi::SecirParams<(size_t)AgeGroup::Count>;
    using Po = epi::Populations<AgeGroup, epi::InfectionType>;
    using Model = epi::CompartmentalModel<Po, Pa>;
    js::class_<epi::SecirModel<AgeGroup>, js::base<Model>>(name.c_str())
            .template constructor<>();
}

template <class AgeGroup>
auto bind_secir_ageres()
{
    size_t constexpr N = (size_t)AgeGroup::Count;

    js::enum_<AgeGroup> agegroup_enum(("AgeGroup" + std::to_string(N)).c_str());
    for (size_t i=0; i < (size_t)AgeGroup::Count; ++i) {
        agegroup_enum.value(("Group" + std::to_string(i)).c_str(), (AgeGroup)i);
    }
    agegroup_enum.value("Count", AgeGroup::Count);

    bind_Populations<AgeGroup, epi::InfectionType>("Populations" + std::to_string(N));
    bind_SecirParams<N>("SecirParams" + std::to_string(N));

    //TODO: Currently, StageTimes and Probabilies are subclasses of the class template SecirParams.
    // If this were not the case, we would not really need a StageTimes class per
    // number of AgeGroups, but since we are planning to remove the SecirParams class
    // this is a valid workaround for now.
    bind_StageTimes<N>("StageTimes" + std::to_string(N));
    bind_Probabilities<N>("Probabilities" + std::to_string(N));

    using Populations = epi::Populations<AgeGroup, epi::InfectionType>;
    bind_CompartmentalModel<Populations, epi::SecirParams<N>>("SecirModelBase" + std::to_string(N));  //<- no flows
    bind_SecirModel<AgeGroup>("SecirModel" + std::to_string(N));  //<-- flows defined in derived constructor at rt

    using SecirModel = epi::SecirModel<AgeGroup>;
    js::function("simulate", &simulate_secir<SecirModel>);

    js::register_vector<SecirModel>(("VectorSecirModel" + std::to_string(N)).c_str());
    js::register_vector<typename epi::SecirParams<N>>(("VectorSecirParams" + std::to_string(N)).c_str());
    js::register_vector<typename epi::SecirParams<N>::Probabilities>(("VectorSecirParamsProbabilities" + std::to_string(N)).c_str());
    js::register_vector<typename epi::SecirParams<N>::StageTimes>(("VectorSecirParamsStageTimes" + std::to_string(N)).c_str());
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

    js::enum_<epi::InfectionType>("InfectionType")
        .value("E", epi::InfectionType::E)
        .value("S", epi::InfectionType::S)
        .value("C", epi::InfectionType::C)
        .value("I", epi::InfectionType::I)
        .value("H", epi::InfectionType::H)
        .value("U", epi::InfectionType::U)
        .value("R", epi::InfectionType::R)
        .value("D", epi::InfectionType::D)
        .value("Count", epi::InfectionType::Count);

    bind_secir_ageres<epi::AgeGroup1>();
    bind_secir_ageres<epi::AgeGroup2>();
    bind_secir_ageres<epi::AgeGroup3>();

    js::register_vector<double>("VectorDouble");
    js::register_vector<size_t>("VectorSizeT");


#ifdef DEBUG
    /**
     * Checks if all memory allocated by the WebAssembly module is freed. If not it will emit a warning with relevant
     * information about a potential leak.
     */
    js::function("doLeakCheck", &__lsan_do_recoverable_leak_check);
#endif
}
