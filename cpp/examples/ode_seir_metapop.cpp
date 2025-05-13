#include "memilio/compartments/simulation.h"
#include "memilio/utils/custom_index_array.h"
#include "models/ode_seir_metapop/infection_state.h"
#include "models/ode_seir_metapop/model.h"
#include "models/ode_seir_metapop/parameters.h"
#include "models/ode_seir_metapop/regions.h"

int main()
{
    const ScalarType t0   = 0.;
    const ScalarType tmax = 10;
    ScalarType dt         = 0.1;

    mio::oseirmetapop::Model<ScalarType> model(3, 1);

    for (size_t i = 0; i < (size_t)model.parameters.get_num_regions(); i++) {
        model.populations[{mio::oseirmetapop::Region(i), mio::AgeGroup(0),
                           mio::oseirmetapop::InfectionState::Susceptible}] = 10000;
    }

    model.populations[{mio::oseirmetapop::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Exposed}] +=
        100;
    model.populations[{mio::oseirmetapop::Region(0), mio::AgeGroup(0),
                       mio::oseirmetapop::InfectionState::Susceptible}] -= 100;

    Eigen::MatrixXd mobility_data_commuter(3, 3);
    mobility_data_commuter << 0.4, 0.3, 0.3, 0.2, 0.7, 0.1, 0.4, 0.1, 0.5;

    model.set_commuting_strengths(mobility_data_commuter);

    model.parameters.template get<mio::oseirmetapop::ContactPatterns<>>()
        .get_cont_freq_mat()[0]
        .get_baseline()
        .setConstant(2.7);

    model.parameters.set<mio::oseirmetapop::TimeExposed<>>(3.335);
    model.parameters.set<mio::oseirmetapop::TimeInfected<>>(8.097612257);
    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<>>(0.07333);

    auto result = simulate(t0, tmax, dt, model);
    return 0;
}
