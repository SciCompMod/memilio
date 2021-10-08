#ifndef GEN_SECIHURD_EXAMPLE_H_
#define GEN_SECIHURD_EXAMPLE_H_

#include "epidemiology/math/smoother.h"
#include "epidemiology/model/compartmentalmodel.h"
#include "epidemiology/model/populations.h"
#include "epidemiology/model/simulation.h"
#include "epidemiology/utils/parameter_set.h"
#include "epidemiology/utils/eigen_util.h"
#include "epidemiology/utils/index.h"

namespace gen {

enum Compartments {
    Susceptible,
    Exposed,
    Carrier,
    Infected,
    Hospitalized,
    ICU,
    Recovered,
    Dead,
    Count
};

using Populations = epi::Populations<Compartments>;

struct StartDay {
    using Type = double;
    static Type get_default() {
        return 0;
    }
};
struct Seasonality {
    using Type = double;
    static Type get_default() {
        return 0;
    }
};
struct IncubationTime {
    using Type = double;
    static Type get_default() {
        return 5.2;
    }
};
struct InfectiousTimeMild {
    using Type = double;
    static Type get_default() {
        return 6;
    }
};
struct SerialInterval {
    using Type = double;
    static Type get_default() {
        return 4.2;
    }
};
struct HospitalizedToHomeTime {
    using Type = double;
    static Type get_default() {
        return 1;
    }
};
struct HomeToHospitalizedTime {
    using Type = double;
    static Type get_default() {
        return 5;
    }
};
struct HospitalizedToICUTime {
    using Type = double;
    static Type get_default() {
        return 2;
    }
};
struct ICUToHomeTime {
    using Type = double;
    static Type get_default() {
        return 8;
    }
};
struct ICUToDeathTime {
    using Type = double;
    static Type get_default() {
        return 5;
    }
};
struct InfectionProbabilityFromContact {
    using Type = double;
    static Type get_default() {
        return 0.05;
    }
};
struct RelativeCarrierInfectability {
    using Type = double;
    static Type get_default() {
        return 1;
    }
};
struct AsymptoticCasesPerInfectious {
    using Type = double;
    static Type get_default() {
        return 0.09;
    }
};
struct RiskOfInfectionFromSympomatic {
    using Type = double;
    static Type get_default() {
        return 0.25;
    }
};
struct HospitalizedCasesPerInfectious {
    using Type = double;
    static Type get_default() {
        return 0.2;
    }
};
struct ICUCasesPerHospitalized {
    using Type = double;
    static Type get_default() {
        return 0.25;
    }
};
struct DeathsPerHospitalized {
    using Type = double;
    static Type get_default() {
        return 0.3;
    }
};
struct TestAndTraceCapacity {
    using Type = double;
    static Type get_default() {
        return double(std::numeric_limits<double>::max());
    }
};
struct ICUCapacity {
    using Type = double;
    static Type get_default() {
        return 0;
    }
};
struct InfectiousTimeAsymptomatic {
    using Type = double;
    static Type get_default() {
        return 5.0;
    }
};

using Parameters = epi::ParameterSet<StartDay, Seasonality, IncubationTime, InfectiousTimeMild, SerialInterval, HospitalizedToHomeTime, HomeToHospitalizedTime, HospitalizedToICUTime, ICUToHomeTime, ICUToDeathTime, InfectionProbabilityFromContact, RelativeCarrierInfectability, AsymptoticCasesPerInfectious, RiskOfInfectionFromSympomatic, HospitalizedCasesPerInfectious, ICUCasesPerHospitalized, DeathsPerHospitalized, TestAndTraceCapacity, ICUCapacity, InfectiousTimeAsymptomatic>;

class secihurd_example : public epi::CompartmentalModel<Populations, Parameters> {
public:
    secihurd_example() : epi::CompartmentalModel<Populations, Parameters>(Populations(epi::Index<Compartments>(Compartments::Count)), Parameters()) {};

    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop,
                         Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override
    {
        const Parameters& par = this->parameters;
                auto icu_occupancy = 0;
        auto test_and_trace_required = 0;
        auto dummy_R3 = 0.5/(par.get<IncubationTime>() - par.get<SerialInterval>());
        test_and_trace_required = dummy_R3*pop[Compartments::Carrier]*(1 - par.get<AsymptoticCasesPerInfectious>());
        icu_occupancy = pop[Compartments::ICU];
        dydt[Compartments::Susceptible] = 0;
        dydt[Compartments::Exposed] = 0;
        auto dummy_R2 = 1.0/(-par.get<IncubationTime>() + 2*par.get<SerialInterval>());
        dummy_R3 = 0.5/(par.get<IncubationTime>() - par.get<SerialInterval>());
        auto risk_from_symptomatic = epi::smoother_cosine(test_and_trace_required, par.get<TestAndTraceCapacity>(), 5*par.get<TestAndTraceCapacity>(), par.get<RiskOfInfectionFromSympomatic>(), 0);
        auto season_val = par.get<Seasonality>()*sin(0.017214206321039961*std::fmod(par.get<StartDay>() + t, 365) + 1.5707963267948966) + 1;
        auto cont_freq_eff = season_val;
        auto N = pop[Compartments::Carrier] + pop[Compartments::Exposed] + pop[Compartments::Hospitalized] + pop[Compartments::ICU] + pop[Compartments::Infected] + pop[Compartments::Recovered] + pop[Compartments::Susceptible];
        auto divN = 1.0/N;
        auto dummy_S = cont_freq_eff*divN*par.get<InfectionProbabilityFromContact>()*pop[Compartments::Susceptible]*(par.get<RelativeCarrierInfectability>()*pop[Compartments::Carrier] + pop[Compartments::Infected]*risk_from_symptomatic);
        dydt[Compartments::Susceptible] -= dummy_S;
        dydt[Compartments::Exposed] += dummy_S;
        auto prob_hosp2icu = epi::smoother_cosine(icu_occupancy, 0.9*par.get<ICUCapacity>(), par.get<ICUCapacity>(), par.get<ICUCasesPerHospitalized>(), 0);
        auto prob_hosp2dead = par.get<ICUCasesPerHospitalized>() - prob_hosp2icu;
        dydt[Compartments::Exposed] -= dummy_R2*pop[Compartments::Exposed];
        dydt[Compartments::Carrier] = dummy_R2*pop[Compartments::Exposed] - pop[Compartments::Carrier]*(dummy_R3*(1 - par.get<AsymptoticCasesPerInfectious>()) + par.get<AsymptoticCasesPerInfectious>()/par.get<InfectiousTimeAsymptomatic>());
        dydt[Compartments::Infected] = dummy_R3*pop[Compartments::Carrier]*(1 - par.get<AsymptoticCasesPerInfectious>()) - pop[Compartments::Infected]*((1 - par.get<HospitalizedCasesPerInfectious>())/par.get<InfectiousTimeMild>() + par.get<HospitalizedCasesPerInfectious>()/par.get<HomeToHospitalizedTime>());
        dydt[Compartments::Hospitalized] = -pop[Compartments::Hospitalized]*(par.get<ICUCasesPerHospitalized>()/par.get<HospitalizedToICUTime>() + (1 - par.get<ICUCasesPerHospitalized>())/par.get<HospitalizedToHomeTime>()) + par.get<HospitalizedCasesPerInfectious>()*pop[Compartments::Infected]/par.get<HomeToHospitalizedTime>();
        dydt[Compartments::ICU] = pop[Compartments::ICU]*(-par.get<DeathsPerHospitalized>()/par.get<ICUToDeathTime>() - (1 - par.get<DeathsPerHospitalized>())/par.get<ICUToHomeTime>());
        dydt[Compartments::ICU] += pop[Compartments::Hospitalized]*prob_hosp2icu/par.get<HospitalizedToICUTime>();
        dydt[Compartments::Recovered] = par.get<AsymptoticCasesPerInfectious>()*pop[Compartments::Carrier]/par.get<InfectiousTimeAsymptomatic>() + pop[Compartments::Infected]*(1 - par.get<HospitalizedCasesPerInfectious>())/par.get<InfectiousTimeMild>() + pop[Compartments::ICU]*(1 - par.get<DeathsPerHospitalized>())/par.get<ICUToHomeTime>() + pop[Compartments::Hospitalized]*(1 - par.get<ICUCasesPerHospitalized>())/par.get<HospitalizedToHomeTime>();
        dydt[Compartments::Dead] = par.get<DeathsPerHospitalized>()*pop[Compartments::ICU]/par.get<ICUToDeathTime>();
        dydt[Compartments::Dead] += pop[Compartments::Hospitalized]*prob_hosp2dead/par.get<HospitalizedToICUTime>();
    };
};

} // end namespace gen

#endif