#ifndef GEN_SIR_EXAMPLE_H_
#define GEN_SIR_EXAMPLE_H_

#include "epidemiology/math/smoother.h"
#include "epidemiology/model/compartmentalmodel.h"
#include "epidemiology/model/populations.h"
#include "epidemiology/model/simulation.h"
#include "epidemiology/utils/parameter_set.h"
#include "epidemiology/utils/eigen_util.h"
#include "epidemiology/utils/index.h"

namespace gen {

enum Compartments {
    S,
    I,
    R,
    Count
};

using Populations = epi::Populations<Compartments>;

struct ContactsPerDay {
    using Type = double;
    static Type get_default() {
        return 15;
    }
};
struct RecoveryRate {
    using Type = double;
    static Type get_default() {
        return 0.5;
    }
};

using Parameters = epi::ParameterSet<ContactsPerDay, RecoveryRate>;

class sir_example : public epi::CompartmentalModel<Populations, Parameters> {
public:
    sir_example() : epi::CompartmentalModel<Populations, Parameters>(Populations(epi::Index<Compartments>(Compartments::Count)), Parameters()) {};

    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop,
                         Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override
    {
        const Parameters& par = this->parameters;
                auto N = pop[Compartments::I] + pop[Compartments::R] + pop[Compartments::S];
        dydt[Compartments::S] = -par.get<ContactsPerDay>()*pop[Compartments::I]*pop[Compartments::S]/N;
        dydt[Compartments::I] = -par.get<RecoveryRate>()*pop[Compartments::I] + par.get<ContactsPerDay>()*pop[Compartments::I]*pop[Compartments::S]/N;
        dydt[Compartments::R] = par.get<RecoveryRate>()*pop[Compartments::I];
    };
};

} // end namespace gen

#endif