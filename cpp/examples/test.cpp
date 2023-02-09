#include "models/ode_seir/model.h"
#include <cassert>
#include <cstddef>
#include <tuple>
#include <type_traits>

mio::oseir::ParametersBase params;
mio::oseir::ParametersBase& parameters = params;
mio::oseir::Model::Populations populations{mio::oseir::InfectionState::Count};

void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                     Eigen::Ref<Eigen::VectorXd> dydt)
{
    using namespace mio::oseir;
    double coeffStoE = params.get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                       params.get<TransmissionProbabilityOnContact>() / populations.get_total();

    double SToE = coeffStoE * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected];
    double EToI = (1.0 / params.get<TimeExposed>()) * y[(size_t)InfectionState::Exposed];
    double IToR = (1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected];

    dydt[(size_t)InfectionState::Susceptible] = -SToE;
    dydt[(size_t)InfectionState::Exposed]     = SToE - EToI;
    dydt[(size_t)InfectionState::Infected]    = EToI - IToR;
    dydt[(size_t)InfectionState::Recovered]   = IToR;
}

void get_flows(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
               Eigen::Ref<Eigen::VectorXd> flows)
{
    using namespace mio::oseir;
    double coeffStoE = params.get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                       params.get<TransmissionProbabilityOnContact>() / populations.get_total();

    flows[0] = coeffStoE * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected];
    flows[1] = (1.0 / params.get<TimeExposed>()) * y[(size_t)InfectionState::Exposed];
    flows[2] = (1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected];
}

template <class Status, Status Source, Status Target>
struct Flow {
    const static Status source = Source;
    const static Status target = Target;
};

template <class... Flows>
class FlowChart
{
public:
    template <size_t i>
    constexpr auto get()
    {
        return std::get<i>(flows);
    }

    template <class Flow>
    constexpr size_t get()
    {
        return get_impl(typelist<Flow>(), typelist<>(), typelist<Flows...>());
    }

private:
    template <class... Types>
    struct typelist {
    };

    template <class Flow, class... Indices>
    constexpr size_t get_impl(typelist<Flow> f, typelist<Indices...>, typelist<>)
    {
        static_assert(sizeof...(Indices) != sizeof...(Flows), "ERROR IN get<Flow>() - could not find flow");
        static_assert(sizeof...(Indices) == sizeof...(Flows), "ERROR IN get<Flow>() - could not find flow");
        return get_impl_(f, typelist<Indices..., int>());
    }

    template <class Flow, class... Indices, class F, class... Fs>
    inline constexpr std::enable_if_t<!std::is_same<Flow, F>::value, size_t>
    get_impl(typelist<Flow> f, typelist<Indices...>, typelist<F, Fs...>)
    {
        return get_impl_(f, typelist<Indices..., int>(), typelist<Fs...>());
    }

    template <class Flow, class... Indices, class F, class... Fs>
    inline constexpr std::enable_if_t<std::is_same<Flow, F>::value, size_t>
    get_impl(typelist<Flow> f, typelist<Indices...>, typelist<F, Fs...>)
    {
        return sizeof...(Indices);
    }

    std::tuple<Flows...> flows;
};

void get_rhs(Eigen::Ref<const Eigen::VectorXd> flows, Eigen::Ref<Eigen::VectorXd> rhs)
{
    get_rhs_impl(flows, rhs);
}

int main()
{
    using I = mio::oseir::InfectionState;
    FlowChart<Flow<I, I::Susceptible, I::Exposed>, Flow<I, I::Exposed, I::Infected>, Flow<I, I::Infected, I::Recovered>>
        f;
    std::cout << (size_t)f.get<0>().target << " ";
    std::cout << f.get<Flow<I, I::Susceptible, I::Exposed>>() << "\n";
    ;
    double total_population                            = 10000;
    populations[mio::oseir::InfectionState::Exposed]   = 2000;
    populations[mio::oseir::InfectionState::Infected]  = 2000;
    populations[mio::oseir::InfectionState::Recovered] = 2000;
    populations[mio::oseir::InfectionState::Susceptible] =
        total_population - populations[mio::oseir::InfectionState::Exposed] -
        populations[mio::oseir::InfectionState::Infected] - populations[mio::oseir::InfectionState::Recovered];
    // suscetible now set with every other update
    // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
    parameters.set<mio::oseir::TimeExposed>(5.2);
    parameters.set<mio::oseir::TimeInfected>(6);
    parameters.set<mio::oseir::TransmissionProbabilityOnContact>(0.04);
    parameters.get<mio::oseir::ContactPatterns>().get_baseline()(0, 0) = 10;

    Eigen::VectorXd dxdt = populations.get_compartments();
    get_derivatives(populations.get_compartments(), populations.get_compartments(), 10, dxdt);
}