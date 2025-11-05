#include "memilio/compartments/simulation.h"
#include "memilio/utils/custom_index_array.h"
#include "models/ode_seir_metapop/model.h"
#include "models/ode_seir_metapop/parameters.h"
#include "memilio/geography/regions.h"

int main()
{
    const ScalarType t0   = 0.;
    const ScalarType tmax = 10;
    ScalarType dt         = 0.1;

    mio::oseirmetapop::Model<ScalarType> model(3, 1);

    for (size_t i = 0; i < (size_t)model.parameters.get_num_regions(); i++) {
        model.populations[{mio::regions::Region(i), mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] =
            10000;
    }

    model.populations[{mio::regions::Region(0), mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] += 100;
    model.populations[{mio::regions::Region(0), mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] -=
        100;

    Eigen::MatrixXd mobility_data_commuter(3, 3);
    mobility_data_commuter << 0.4, 0.3, 0.3, 0.2, 0.7, 0.1, 0.4, 0.1, 0.5;

    model.set_commuting_strengths(mobility_data_commuter);

    mio::oseir::Parameters<ScalarType> local_params(1);

    local_params.template get<mio::oseir::ContactPatterns<>>()
        .get_cont_freq_mat()[0]
        .get_baseline()
        .setConstant(2.7);

    local_params.set<mio::oseir::TimeExposed<>>(3.335);
    local_params.set<mio::oseir::TimeInfected<>>(8.097612257);
    local_params.set<mio::oseir::TransmissionProbabilityOnContact<>>(0.07333);

    model.set_local_parameters(local_params);

    auto result = simulate(t0, tmax, dt, model);

    result.print_table();
    return 0;
}
